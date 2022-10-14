# noqa: D100

import functools
import logging
import operator
from typing import Callable, Dict, List, Optional, Union

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.reference_genome import get_reference_genome

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Filter genotypes to adj criteria."""
    if "adj" not in list(mt.entry):
        mt = annotate_adj(mt)
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
            raise Exception("No downsampling data for subpopulations implemented")
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


def filter_low_conf_regions(
    mt: Union[hl.MatrixTable, hl.Table],
    filter_lcr: bool = True,
    filter_decoy: bool = True,
    filter_segdup: bool = True,
    filter_exome_low_coverage_regions: bool = False,
    filter_telomeres_and_centromeres: bool = False,
    high_conf_regions: Optional[List[str]] = None,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter low-confidence regions.

    :param mt: MatrixTable or Table to filter
    :param filter_lcr: Whether to filter LCR regions
    :param filter_decoy: Whether to filter decoy regions
    :param filter_segdup: Whether to filter Segdup regions
    :param filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
    :param filter_telomeres_and_centromeres: Whether to filter telomeres and centromeres
    :param high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MatrixTable or Table with low confidence regions removed
    """
    build = get_reference_genome(mt.locus).name
    if build == "GRCh37":
        import gnomad.resources.grch37.reference_data as resources
    elif build == "GRCh38":
        import gnomad.resources.grch38.reference_data as resources

    criteria = []
    if filter_lcr:
        lcr = resources.lcr_intervals.ht()
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_decoy:
        decoy = resources.decoy_intervals.ht()
        criteria.append(hl.is_missing(decoy[mt.locus]))

    if filter_segdup:
        segdup = resources.seg_dup_intervals.ht()
        criteria.append(hl.is_missing(segdup[mt.locus]))

    if filter_exome_low_coverage_regions:
        high_cov = resources.high_coverage_intervals.ht()
        criteria.append(hl.is_missing(high_cov[mt.locus]))

    if filter_telomeres_and_centromeres:
        if build != "GRCh38":
            raise DataException(
                "The telomeres_and_centromeres resource only exists for GRCh38"
            )

        telomeres_and_centromeres = resources.telomeres_and_centromeres.ht()
        criteria.append(hl.is_missing(telomeres_and_centromeres[mt.locus]))

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            criteria.append(hl.is_defined(region[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


def filter_to_autosomes(
    t: Union[hl.MatrixTable, hl.Table]
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
            hl.cond(filter_condition, hl.set([filter_name]), hl.empty_set(hl.tstr))
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
    t: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter to loci that are in non-PAR regions on chromosome X.

    :param t: Input Table or MatrixTable.
    :return: Filtered Table or MatrixTable.
    """
    rg = t.locus.dtype.reference_genome
    t = hl.filter_intervals(
        t, [hl.parse_locus_interval(contig) for contig in rg.x_contigs]
    )
    non_par_expr = t.locus.in_x_nonpar()

    return (
        t.filter(non_par_expr)
        if isinstance(t, hl.Table)
        else t.filter_rows(non_par_expr)
    )


def filter_y_nonpar(
    t: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter to loci that are in non-PAR regions on chromosome Y.

    :param t: Input Table or MatrixTable.
    :return: Filtered Table or MatrixTable.
    """
    rg = t.locus.dtype.reference_genome
    t = hl.filter_intervals(
        t, [hl.parse_locus_interval(contig) for contig in rg.y_contigs]
    )
    non_par_expr = t.locus.in_y_nonpar()

    return (
        t.filter(non_par_expr)
        if isinstance(t, hl.Table)
        else t.filter_rows(non_par_expr)
    )
