import functools
import logging
import operator
from typing import Dict, List, Optional, Union

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.reference_genome import get_reference_genome

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter genotypes to adj criteria
    """
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
    Filter MatrixTable or Table with gnomAD-format frequency data (assumed bi-allelic/split)
    (i.e. Array[Struct(Array[AC], Array[AF], AN, homozygote_count, meta)])

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
        # If one supplies a subpop but not a population, this will ensure this gets it right
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

    def combine_functions(func_list, x):
        cond = func_list[0](x)
        for c in func_list[1:]:
            cond &= c(x)
        return cond

    filt = lambda x: combine_functions(criteria, x)
    criteria = hl.any(filt, t.freq)
    return (
        t.filter_rows(criteria, keep=keep)
        if isinstance(t, hl.MatrixTable)
        else t.filter(criteria, keep=keep)
    )


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
    Filters low-confidence regions

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
    Filters the Table or MatrixTable to autosomes only.
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
    Creates an expression to create or add filters.
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
    mt: hl.MatrixTable,
    sample_path: str,
    header: bool = True,
    table_key: str = "s",
    sparse: bool = False,
    gt_expr: str = "GT",
) -> hl.MatrixTable:
    """
    Subsets the MatrixTable to the provided list of samples and their variants

    :param mt: Input MatrixTable
    :param sample_path: Path to a file with list of samples
    :param header: Whether file with samples has a header. Default is True
    :param table_key: Key to sample Table. Default is "s"
    :param sparse: Whether the MatrixTable is sparse. Default is False
    :param gt_expr: Name of field in MatrixTable containing genotype expression. Default is "GT"
    :return: MatrixTable subsetted to specified samples and their variants
    """
    sample_ht = hl.import_table(sample_path, no_header=not header, key=table_key)
    sample_count = sample_ht.count()
    missing_ht = sample_ht.anti_join(mt.cols())
    missing_ht_count = missing_ht.count()
    full_count = mt.count_cols()

    if missing_ht_count != 0:
        missing_samples = missing_ht.s.collect()
        raise DataException(
            f"Only {sample_count - missing_ht_count} out of {sample_count} "
            "subsetting-table IDs matched IDs in the MT.\n"
            f"IDs that aren't in the MT: {missing_samples}\n"
        )

    mt = mt.semi_join_cols(sample_ht)
    if sparse:
        mt = mt.filter_rows(
            hl.agg.any(mt[gt_expr].is_non_ref() | hl.is_defined(mt.END))
        )
    else:
        mt = mt.filter_rows(hl.agg.any(mt[gt_expr].is_non_ref()))

    logger.info(
        f"Finished subsetting samples. Kept {mt.count_cols()} "
        f"out of {full_count} samples in MT"
    )
    return mt


def filter_to_clinvar_pathogenic(
    t: Union[hl.MatrixTable, hl.Table],
    clnrevstat_field: str = "CLNREVSTAT",
    clnsig_field: str = "CLNSIG",
    clnsigconf_field: str = "CLNSIGCONF",
    remove_no_assertion: bool = True,
    remove_conflicting: bool = True,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Returns a MatrixTable or Table that filters the clinvar data to pathogenic and likely pathogenic variants.

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
        f"Found {t.count_rows() if isinstance(t, hl.MatrixTable) else t.count()} variants before filtering"
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
            f"Variants with conflicting clinical interpretations will be removed."
        )
        path_expr = path_expr & hl.is_missing(t.info[clnsigconf_field])

    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(path_expr)
    else:
        t = t.filter(path_expr)

    logger.info(
        f"Found {t.count_rows() if isinstance(t, hl.MatrixTable) else t.count()} variants after filtering to clinvar pathogenic variants."
    )
    return t
