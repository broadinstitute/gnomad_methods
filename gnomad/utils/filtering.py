# noqa: D100

import functools
import logging
import operator
from typing import Callable, Dict, List, Optional, Tuple, Union

import hail as hl

import gnomad.utils.annotations as annotate_utils
from gnomad.resources.resource_utils import DataException
from gnomad.utils.reference_genome import get_reference_genome

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Filter genotypes to adj criteria."""
    if "adj" not in list(mt.entry):
        mt = annotate_utils.annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)
    return mt.drop(mt.adj)


def filter_by_frequency(
    t: Union[hl.MatrixTable, hl.Table],
    direction: str,
    frequency: float = None,
    allele_count: int = None,
    population: str = None,
    subpop: str = None,
    downsampling: int = None,
    keep: bool = True,
    adj: bool = True,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter MatrixTable or Table with gnomAD-format frequency data (assumed bi-allelic/split).

    gnomAD frequency data format expectation is: Array[Struct(Array[AC], Array[AF], AN, homozygote_count, meta)].

    At least one of frequency or allele_count is required.

    Subpop can be specified without a population if desired.

    :param t: Input MatrixTable or Table
    :param direction: One of "above", "below", and "equal" (how to apply the filter)
    :param frequency: Frequency to filter by (one of frequency or allele_count is required)
    :param allele_count: Allele count to filter by (one of frequency or allele_count is required)
    :param population: Population in which to filter frequency
    :param subpop: Sub-population in which to filter frequency
    :param downsampling: Downsampling in which to filter frequency
    :param keep: Whether to keep rows passing this frequency (passed to filter_rows)
    :param adj: Whether to use adj frequency
    :return: Filtered MatrixTable or Table
    """
    if frequency is None and allele_count is None:
        raise ValueError("At least one of frequency or allele_count must be specified")
    if direction not in ("above", "below", "equal"):
        raise ValueError('direction needs to be one of "above", "below", or "equal"')
    group = "adj" if adj else "raw"
    criteria = [lambda f: f.meta.get("group", "") == group]
    if frequency is not None:
        if direction == "above":
            criteria.append(lambda f: f.AF[1] > frequency)
        elif direction == "below":
            criteria.append(lambda f: f.AF[1] < frequency)
        else:
            criteria.append(lambda f: f.AF[1] == frequency)
    if allele_count is not None:
        if direction == "above":
            criteria.append(lambda f: f.AC[1] > allele_count)
        elif direction == "below":
            criteria.append(lambda f: f.AC[1] < allele_count)
        else:
            criteria.append(lambda f: f.AC[1] == allele_count)
    size = 1
    if population:
        criteria.append(lambda f: f.meta.get("pop", "") == population)
        size += 1
    if subpop:
        criteria.append(lambda f: f.meta.get("subpop", "") == subpop)
        size += 1
        # If one supplies a subpop but not a population, this will ensure this
        # gets it right
        if not population:
            size += 1
    if downsampling:
        criteria.append(lambda f: f.meta.get("downsampling", "") == str(downsampling))
        size += 1
        if not population:
            size += 1
            criteria.append(lambda f: f.meta.get("pop", "") == "global")
        if subpop:
            raise ValueError("No downsampling data for subpopulations implemented")
    criteria.append(lambda f: f.meta.size() == size)

    filt = lambda x: combine_functions(criteria, x)
    criteria = hl.any(filt, t.freq)
    return (
        t.filter_rows(criteria, keep=keep)
        if isinstance(t, hl.MatrixTable)
        else t.filter(criteria, keep=keep)
    )


def combine_functions(
    func_list: List[Callable[[bool], bool]],
    x: hl.expr.StructExpression,
    operator_func: Callable[[bool, bool], bool] = operator.iand,
) -> bool:
    """
    Combine a list of boolean functions to an Expression using the specified operator.

    .. note::

        The `operator_func` is applied cumulatively from left to right of the `func_list`.

    :param func_list: A list of boolean functions that can be applied to `x`.
    :param x: Expression to be passed to each function in `func_list`.
    :param operator_func: Operator function to combine the functions in `func_list`. Default is `operator.iand`.
    :return: A boolean from the combined operations.
    """
    cond = func_list[0](x)
    for c in func_list[1:]:
        cond = operator_func(cond, c(x))
    return cond


def low_conf_regions_expr(
    locus_expr: hl.expr.LocusExpression,
    filter_lcr: bool = True,
    filter_decoy: bool = True,
    filter_segdup: bool = True,
    filter_exome_low_coverage_regions: bool = False,
    filter_telomeres_and_centromeres: bool = False,
    high_conf_regions: Optional[List[str]] = None,
) -> hl.expr.BooleanExpression:
    """
    Create an expression to filter low confidence regions.

    :param locus_expr: Locus expression to use for filtering.
    :param filter_lcr: Whether to filter LCR regions
    :param filter_decoy: Whether to filter decoy regions
    :param filter_segdup: Whether to filter Segdup regions
    :param filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
    :param filter_telomeres_and_centromeres: Whether to filter telomeres and centromeres
    :param high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: Bool expression of whether loci are not low confidence (TRUE) or low confidence (FALSE)
    """
    build = get_reference_genome(locus_expr).name
    if build == "GRCh37":
        import gnomad.resources.grch37.reference_data as resources
    elif build == "GRCh38":
        import gnomad.resources.grch38.reference_data as resources
    else:
        raise ValueError(f"Unsupported reference genome build: {build}")

    criteria = []
    if filter_lcr:
        lcr = resources.lcr_intervals.ht()
        criteria.append(hl.is_missing(lcr[locus_expr]))

    if filter_decoy:
        decoy = resources.decoy_intervals.ht()
        criteria.append(hl.is_missing(decoy[locus_expr]))

    if filter_segdup:
        segdup = resources.seg_dup_intervals.ht()
        criteria.append(hl.is_missing(segdup[locus_expr]))

    if filter_exome_low_coverage_regions:
        high_cov = resources.high_coverage_intervals.ht()
        criteria.append(hl.is_missing(high_cov[locus_expr]))

    if filter_telomeres_and_centromeres:
        if build != "GRCh38":
            raise DataException(
                "The telomeres_and_centromeres resource only exists for GRCh38"
            )

        telomeres_and_centromeres = resources.telomeres_and_centromeres.ht()
        criteria.append(hl.is_missing(telomeres_and_centromeres[locus_expr]))

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            criteria.append(hl.is_defined(region[locus_expr]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        return filter_criteria
    else:
        raise ValueError("No low confidence regions requested for filtering.")


def filter_low_conf_regions(
    t: Union[hl.MatrixTable, hl.Table],
    **kwargs,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter low-confidence regions.

    :param t: MatrixTable or Table to filter.
    :param kwargs: Keyword arguments to pass to `low_conf_regions_expr`.
    :return: MatrixTable or Table with low confidence regions removed.
    """
    filter_criteria = low_conf_regions_expr(t.locus, **kwargs)
    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(filter_criteria)
    else:
        t = t.filter(filter_criteria)

    return t


def filter_to_autosomes(
    t: Union[hl.MatrixTable, hl.Table],
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter the Table or MatrixTable to autosomes only.

    This assumes that the input contains a field named `locus` of type Locus

    :param t: Input MT/HT
    :return:  MT/HT autosomes
    """
    reference = get_reference_genome(t.locus)
    autosomes = hl.parse_locus_interval(
        f"{reference.contigs[0]}-{reference.contigs[21]}", reference_genome=reference
    )
    return hl.filter_intervals(t, [autosomes])


def add_filters_expr(
    filters: Dict[str, hl.expr.BooleanExpression],
    current_filters: hl.expr.SetExpression = None,
) -> hl.expr.SetExpression:
    """
    Create an expression to create or add filters.

    For each entry in the `filters` dictionary, if the value evaluates to `True`,
    then the key is added as a filter name.

    Current filters are kept if provided using `current_filters`

    :param filters: The filters and their expressions
    :param current_filters: The set of current filters
    :return: An expression that can be used to annotate the filters
    """
    if current_filters is None:
        current_filters = hl.empty_set(hl.tstr)

    return hl.fold(
        lambda x, y: x.union(y),
        current_filters,
        [
            hl.if_else(filter_condition, hl.set([filter_name]), hl.empty_set(hl.tstr))
            for filter_name, filter_condition in filters.items()
        ],
    )


def subset_samples_and_variants(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    sample_path: str,
    header: bool = True,
    table_key: str = "s",
    sparse: bool = False,
    gt_expr: str = "GT",
    remove_dead_alleles: bool = False,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Subset the MatrixTable or VariantDataset to the provided list of samples and their variants.

    :param mtds: Input MatrixTable or VariantDataset
    :param sample_path: Path to a file with list of samples
    :param header: Whether file with samples has a header. Default is True
    :param table_key: Key to sample Table. Default is "s"
    :param sparse: Whether the MatrixTable is sparse. Default is False
    :param gt_expr: Name of field in MatrixTable containing genotype expression. Default is "GT"
    :param remove_dead_alleles: Remove alleles observed in no samples. This option is currently only relevant when `mtds` is a VariantDataset. Default is False
    :return: MatrixTable or VariantDataset subsetted to specified samples and their variants
    """
    sample_ht = hl.import_table(sample_path, no_header=not header, key=table_key)
    sample_count = sample_ht.count()
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        if remove_dead_alleles:
            raise ValueError(
                "Removal of alleles observed in no samples is currently only"
                " implemented when the input dataset is a VariantDataset."
            )
        mt = mtds
    missing_ht = sample_ht.anti_join(mt.cols())
    missing_ht_count = missing_ht.count()
    full_count = mt.count_cols()

    if missing_ht_count != 0:
        missing_samples = missing_ht.s.collect()
        raise DataException(
            f"Only {sample_count - missing_ht_count} out of"
            f" {sample_count} subsetting-table IDs matched IDs in the"
            f" {'VariantDataset' if is_vds else 'MatrixTable'}.\nIDs that aren't in the"
            f" MT: {missing_samples}\n"
        )

    if is_vds:
        mtds = hl.vds.filter_samples(
            mtds, sample_ht, keep=True, remove_dead_alleles=remove_dead_alleles
        )
        n_cols = mtds.variant_data.count_cols()
    else:
        mtds = mtds.semi_join_cols(sample_ht)
        if sparse:
            mtds = mtds.filter_rows(
                hl.agg.any(mtds[gt_expr].is_non_ref() | hl.is_defined(mtds.END))
            )
        else:
            mtds = mtds.filter_rows(hl.agg.any(mtds[gt_expr].is_non_ref()))
        n_cols = mtds.count_cols()

    logger.info(
        "Finished subsetting samples. Kept %d out of %d samples in %s",
        n_cols,
        full_count,
        "VariantDataset" if is_vds else "MatrixTable",
    )
    return mtds


def filter_to_clinvar_pathogenic(
    t: Union[hl.MatrixTable, hl.Table],
    clnrevstat_field: str = "CLNREVSTAT",
    clnsig_field: str = "CLNSIG",
    clnsigconf_field: str = "CLNSIGCONF",
    remove_no_assertion: bool = True,
    remove_conflicting: bool = True,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Return a MatrixTable or Table that filters the clinvar data to pathogenic and likely pathogenic variants.

    Example use:

    .. code-block:: python

        from gnomad.resources.grch38.reference_data import clinvar
        clinvar_ht = clinvar.ht()
        clinvar_ht = filter_to_clinvar_pathogenic(clinvar_ht)

    :param: t: Input dataset that contains clinvar data, could either be a MatrixTable or Table.
    :param clnrevstat_field: The field string for the expression that contains the review status of the clinical significance of clinvar variants.
    :param clnsig_field: The field string for the expression that contains the clinical signifcance of the clinvar variant.
    :param clnsigconf_field: The field string for the expression that contains the conflicting clinical significance values for the variant. For variants with no conflicting significance, this field should be undefined.
    :param remove_no_assertion: Flag for removing entries in which the clnrevstat (clinical significance) has no assertions (zero stars).
    :param remove_conflicting: Flag for removing entries with conflicting clinical interpretations.
    :return: Filtered MatrixTable or Table
    """
    logger.info(
        "Found %d variants before filtering",
        t.count_rows() if isinstance(t, hl.MatrixTable) else t.count(),
    )
    path_expr = (
        t.info[clnsig_field]
        .map(lambda x: x.lower())
        .map(lambda x: x.contains("pathogenic"))
        .any(lambda x: x)
    )

    if remove_no_assertion:
        logger.info("Variants without assertions will be removed.")
        no_star_assertions = hl.literal(
            {
                "no_assertion_provided",
                "no_assertion_criteria_provided",
                "no_interpretation_for_the_individual_variant",
            }
        )
        path_expr = path_expr & (
            hl.set(t.info[clnrevstat_field]).intersection(no_star_assertions).length()
            == 0
        )

    if remove_conflicting:
        logger.info(
            "Variants with conflicting clinical interpretations will be removed."
        )
        path_expr = path_expr & hl.is_missing(t.info[clnsigconf_field])

    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(path_expr)
    else:
        t = t.filter(path_expr)

    logger.info(
        "Found %d variants after filtering to clinvar pathogenic variants.",
        t.count_rows() if isinstance(t, hl.MatrixTable) else t.count(),
    )
    return t


def filter_to_gencode_cds(
    t: Union[hl.MatrixTable, hl.Table], gencode_ht: Optional[hl.Table] = None
) -> hl.Table:
    """
    Filter a Table/MatrixTable to only Gencode CDS regions in protein coding transcripts.

    Example use:

    .. code-block:: python

        from gnomad.resources.grch37.reference_data import gencode
        gencode_ht = gencode.ht()
        gencode_ht = filter_gencode_to_cds(gencode_ht)

    .. note::

        If no Gencode Table is provided, the default version of the Gencode Table
        resource for the genome build of the input Table/MatrixTable will be used.

    .. warning::

        This Gencode CDS interval filter does not take into account the
        transcript_id, it filters to any locus that is found in a CDS interval for
        any protein coding transcript. Therefore, if downstream analyses require
        filtering to CDS intervals by transcript, an additional step must be taken.
        For example, when filtering VEP transcript consequences, there may be cases
        where a variant is retained with this filter, but is considered outside the
        CDS intervals of the transcript per the VEP predicted consequence of the
        variant.

    :param t: Input Table/MatrixTable to filter.
    :param gencode_ht: Gencode Table to use for filtering the input Table/MatrixTable
        to CDS regions. Default is None, which will use the default version of the
        Gencode Table resource.
    :return: Table/MatrixTable filtered to loci in Gencode CDS intervals.
    """
    if gencode_ht is None:
        build = get_reference_genome(t.locus).name
        if build == "GRCh37":
            from gnomad.resources.grch37.reference_data import gencode
        elif build == "GRCh38":
            from gnomad.resources.grch38.reference_data import gencode
        else:
            raise ValueError(f"Unsupported reference genome build: {build}")

        logger.info(
            "No Gencode Table was supplied, using Gencode version %s",
            gencode.default_version,
        )
        gencode_ht = gencode.ht()

    gencode_ht = gencode_ht.filter(
        (gencode_ht.feature == "CDS") & (gencode_ht.transcript_type == "protein_coding")
    )
    logger.warning(
        "This Gencode CDS interval filter does not filter by transcript! Please see the"
        " documentation for more details to confirm it's being used as intended."
    )
    filter_expr = hl.is_defined(gencode_ht[t.locus])

    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(filter_expr)
    else:
        t = t.filter(filter_expr)

    return t


def remove_fields_from_constant(
    constant: List[str], fields_to_remove: List[str]
) -> List[str]:
    """
    Remove fields from a list and display any field(s) missing from the original list.

    :param constant: List of fields
    :param fields_to_remove: List of fields to remove from `constant`
    """
    for field in fields_to_remove:
        if field in constant:
            constant.remove(field)
        else:
            logger.info("%s missing from %s", field, constant)

    return constant


def filter_x_nonpar(
    t: Union[hl.Table, hl.MatrixTable],
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter to loci that are in non-PAR regions on chromosome X.

    :param t: Input Table or MatrixTable.
    :return: Filtered Table or MatrixTable.
    """
    rg = t.locus.dtype.reference_genome
    t = hl.filter_intervals(
        t,
        [
            hl.parse_locus_interval(contig, reference_genome=rg.name)
            for contig in rg.x_contigs
        ],
    )
    non_par_expr = t.locus.in_x_nonpar()

    return (
        t.filter(non_par_expr)
        if isinstance(t, hl.Table)
        else t.filter_rows(non_par_expr)
    )


def filter_y_nonpar(
    t: Union[hl.Table, hl.MatrixTable],
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter to loci that are in non-PAR regions on chromosome Y.

    :param t: Input Table or MatrixTable.
    :return: Filtered Table or MatrixTable.
    """
    rg = t.locus.dtype.reference_genome
    t = hl.filter_intervals(
        t,
        [
            hl.parse_locus_interval(contig, reference_genome=rg.name)
            for contig in rg.y_contigs
        ],
    )
    non_par_expr = t.locus.in_y_nonpar()

    return (
        t.filter(non_par_expr)
        if isinstance(t, hl.Table)
        else t.filter_rows(non_par_expr)
    )


def filter_by_numeric_expr_range(
    t: Union[hl.MatrixTable, hl.Table],
    filter_expr: hl.NumericExpression,
    filter_range: tuple,
    keep_between: bool = True,
    inclusive: bool = True,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter rows in the Table/MatrixTable based on the range of a numeric expression.

    :param t: Input Table/MatrixTable.
    :param filter_expr: NumericExpression to apply `filter_range` to.
    :param filter_range: Range of values to apply to `filter_expr`.
    :param keep_between: Whether to keep the values between `filter_range` instead of keeping values outside `filter_range`. Default is True.
    :param inclusive: Whether or not to include the `filter_range` values themselves. Default is True.
    :return: Table/MatrixTable filtered to rows with specified criteria.
    """
    if inclusive and keep_between or not inclusive and not keep_between:
        criteria = (filter_expr >= filter_range[0]) & (filter_expr <= filter_range[1])
    else:
        criteria = (filter_expr > filter_range[0]) & (filter_expr < filter_range[1])

    if isinstance(t, hl.MatrixTable):
        return t.filter_rows(criteria, keep=keep_between)
    else:
        return t.filter(criteria, keep=keep_between)


def filter_for_mu(
    ht: hl.Table, gerp_lower_cutoff: float = -3.9885, gerp_upper_cutoff: float = 2.6607
) -> hl.Table:
    """
    Filter to non-coding annotations and remove GERP outliers.

    .. note::

        Values for `gerp_lower_cutoff` and `gerp_upper_cutoff` default to -3.9885 and
        2.6607, respectively. These values were precalculated on the GRCh37 context
        table and define the 5th and 95th percentiles.

    :param ht: Input Table.
    :param gerp_lower_cutoff: Minimum GERP score for variant to be included. Default is -3.9885.
    :param gerp_upper_cutoff: Maximum GERP score for variant to be included. Default is 2.6607.
    :return: Table filtered to intron or intergenic variants with GERP outliers removed.
    """
    ht = filter_by_numeric_expr_range(
        ht,
        filter_expr=ht.gerp,
        filter_range=(gerp_lower_cutoff, gerp_upper_cutoff),
        keep_between=True,
        inclusive=False,
    )
    ht = ht.filter(
        (ht.vep.most_severe_consequence == "intron_variant")
        | (ht.vep.most_severe_consequence == "intergenic_variant")
    )

    return ht


def split_vds_by_strata(
    vds: hl.vds.VariantDataset, strata_expr: hl.expr.Expression
) -> Dict[str, hl.vds.VariantDataset]:
    """
    Split a VDS into multiple VDSs based on `strata_expr`.

    :param vds: Input VDS.
    :param strata_expr: Expression on VDS variant_data MT to split on.
    :return: Dictionary where strata value is key and VDS is value.
    """
    vmt = vds.variant_data
    s_by_strata = vmt.aggregate_cols(
        hl.agg.group_by(strata_expr, hl.agg.collect_as_set(vmt.s))
    )

    return {
        strata: hl.vds.filter_samples(vds, list(s)) for strata, s in s_by_strata.items()
    }


def filter_arrays_by_meta(
    meta_expr: hl.expr.ArrayExpression,
    meta_indexed_exprs: Union[
        Dict[str, hl.expr.ArrayExpression], hl.expr.ArrayExpression
    ],
    items_to_filter: Union[Dict[str, List[str]], List[str]],
    keep: bool = True,
    combine_operator: str = "and",
    exact_match: bool = False,
) -> Tuple[
    hl.expr.ArrayExpression,
    Union[Dict[str, hl.expr.ArrayExpression], hl.expr.ArrayExpression],
]:
    """
    Filter both metadata array expression and meta data indexed expression by `items_to_filter`.

    The `items_to_filter` can be used to filter in the following ways based on
    `meta_expr` items:
    - By a list of keys, e.g. ["sex", "downsampling"].
    - By specific key: value pairs, e.g. to filter where 'pop' is 'han' or 'papuan'
    {"pop": ["han", "papuan"]}, or where 'pop' is 'afr' and/or 'sex' is 'XX'
    {"pop": ["afr"], "sex": ["XX"]}.

    The items can be kept or removed from `meta_indexed_expr` and `meta_expr` based on
    the value of `keep`. For example if `meta_indexed_exprs` is {'freq': ht.freq,
    'freq_meta_sample_count': ht.index_globals().freq_meta_sample_count} and `meta_expr`
    is ht.freq_meta then if `keep` is True, the items specified by `items_to_filter`
    such as  'pop' = 'han' will be kept and all other items will be removed from the
    ht.freq, ht.freq_meta_sample_count, and ht.freq_meta. `meta_indexed_exprs` can also
    be a single array expression such as ht.freq.

    The filtering can also be applied such that all criteria must be met
    (`combine_operator` = "and") by the `meta_expr` item in order to be filtered,
    or at least one of the specified criteria must be met (`combine_operator` = "or")
    by the `meta_expr` item in order to be filtered.

    The `exact_match` parameter can be used to apply the `keep` parameter to only items
    specified in the `items_to_filter` parameter. For example, by default, if `keep` is
    True, `combine_operator` is "and", and `items_to_filter` is ["sex", "downsampling"],
    then all items in `meta_expr` with both "sex" and "downsampling" as keys will be
    kept. However, if `exact_match` is True, then the items
    in `meta_expr` will only be kept if "sex" and "downsampling" are the only keys in
    the meta dict.

    :param meta_expr: Metadata expression that contains the values of the elements in
        `meta_indexed_expr`. The most often used expression is `freq_meta` to index into
        a 'freq' array.
    :param meta_indexed_exprs: Either a Dictionary where the keys are the expression name
        and the values are the expressions indexed by the `meta_expr` such as a 'freq'
        array or just a single expression indexed by the `meta_expr`.
    :param items_to_filter: Items to filter by, either a list or a dictionary.
    :param keep: Whether to keep or remove the items specified by `items_to_filter`.
    :param combine_operator: Whether to use "and" or "or" to combine the items
        specified by `items_to_filter`.
    :param exact_match: Whether to apply the `keep` parameter to only the items
        specified in the `items_to_filter` parameter or to all items in `meta_expr`.
        See the example above for more details. Default is False.
    :return: A Tuple of the filtered metadata expression and a dictionary of metadata
        indexed expressions when meta_indexed_expr is a Dictionary or a single filtered
        array expression when meta_indexed_expr is a single array expression.
    """
    meta_expr = meta_expr.collect(_localize=False)[0]

    if isinstance(meta_indexed_exprs, hl.expr.ArrayExpression):
        meta_indexed_exprs = {"_tmp": meta_indexed_exprs}

    if combine_operator == "and":
        operator_func = hl.all
    elif combine_operator == "or":
        operator_func = hl.any
    else:
        raise ValueError(
            "combine_operator must be one of 'and' or 'or', but found"
            f" {combine_operator}!"
        )

    if isinstance(items_to_filter, list):
        items_to_filter_set = hl.set(items_to_filter)
        items_to_filter = [[k] for k in items_to_filter]
        if exact_match:
            filter_func = lambda m, k: (
                hl.len(hl.set(m.keys()).difference(items_to_filter_set)) == 0
            ) & m.contains(k)
        else:
            filter_func = lambda m, k: m.contains(k)
    elif isinstance(items_to_filter, dict):
        items_to_filter = [
            [(k, v) for v in values] for k, values in items_to_filter.items()
        ]
        items_to_filter_set = hl.set(hl.flatten(items_to_filter))
        if exact_match:
            filter_func = lambda m, k: (
                (hl.len(hl.set(m.items()).difference(items_to_filter_set)) == 0)
                & (m.get(k[0], "") == k[1])
            )
        else:
            filter_func = lambda m, k: (m.get(k[0], "") == k[1])
    else:
        raise TypeError("items_to_filter must be a list or a dictionary!")

    meta_expr = hl.enumerate(meta_expr).filter(
        lambda m: hl.bind(
            lambda x: hl.if_else(keep, x, ~x),
            operator_func(
                [hl.any([filter_func(m[1], v) for v in k]) for k in items_to_filter]
            ),
        ),
    )

    meta_indexed_exprs = {
        k: meta_expr.map(lambda x: v[x[0]]) for k, v in meta_indexed_exprs.items()
    }
    meta_expr = meta_expr.map(lambda x: x[1])

    if "_tmp" in meta_indexed_exprs:
        meta_indexed_exprs = meta_indexed_exprs["_tmp"]

    return meta_expr, meta_indexed_exprs
