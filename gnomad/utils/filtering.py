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
    high_conf_regions: Optional[List[str]] = None,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters low-confidence regions

    :param mt: MatrixTable or Table to filter
    :param filter_lcr: Whether to filter LCR regions
    :param filter_decoy: Whether to filter decoy regions
    :param filter_segdup: Whether to filter Segdup regions
    :param filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
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
