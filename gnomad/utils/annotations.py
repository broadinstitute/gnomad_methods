# noqa: D100

import itertools
import logging
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import hail as hl

from gnomad.utils.gen_stats import to_phred

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

ANNOTATIONS_HISTS = {
    "FS": (0, 50, 50),  # NOTE: in 2.0.2 release this was on (0,20)
    "InbreedingCoeff": (-0.25, 0.25, 50),
    "MQ": (0, 80, 40),
    "RAW_MQ": (2, 13, 33),
    "MQRankSum": (-15, 15, 60),
    "QD": (0, 40, 40),
    "ReadPosRankSum": (-15, 15, 60),
    "SOR": (0, 10, 50),
    "BaseQRankSum": (-15, 15, 60),
    "ClippingRankSum": (-5, 5, 40),
    "DP": (1, 9, 32),  # NOTE: in 2.0.2 release this was on (0,8)
    "VQSLOD": (-30, 30, 60),  # NOTE: in 2.0.2 release this was on (-20,20)
    "AS_VQSLOD": (-30, 30, 60),
    "rf_tp_probability": (0, 1, 50),
    "pab_max": (0, 1, 50),
}


def pop_max_expr(
    freq: hl.expr.ArrayExpression,
    freq_meta: hl.expr.ArrayExpression,
    pops_to_exclude: Optional[Set[str]] = None,
) -> hl.expr.StructExpression:
    """

    Create an expression containing the frequency information about the population that has the highest AF in `freq_meta`.

    Populations specified in `pops_to_exclude` are excluded and only frequencies from adj populations are considered.

    This resulting struct contains the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32
        - pop: str

    :param freq: ArrayExpression of Structs with fields ['AC', 'AF', 'AN', 'homozygote_count']
    :param freq_meta: ArrayExpression of meta dictionaries corresponding to freq (as returned by annotate_freq)
    :param pops_to_exclude: Set of populations to skip for popmax calcluation

    :return: Popmax struct
    """
    _pops_to_exclude = (
        hl.literal(pops_to_exclude)
        if pops_to_exclude is not None
        else hl.empty_set(hl.tstr)
    )

    # pylint: disable=invalid-unary-operand-type
    popmax_freq_indices = hl.range(0, hl.len(freq_meta)).filter(
        lambda i: (hl.set(freq_meta[i].keys()) == {"group", "pop"})
        & (freq_meta[i]["group"] == "adj")
        & (~_pops_to_exclude.contains(freq_meta[i]["pop"]))
    )
    freq_filtered = popmax_freq_indices.map(
        lambda i: freq[i].annotate(pop=freq_meta[i]["pop"])
    ).filter(lambda f: f.AC > 0)

    sorted_freqs = hl.sorted(freq_filtered, key=lambda x: x.AF, reverse=True)
    return hl.or_missing(hl.len(sorted_freqs) > 0, sorted_freqs[0])


def project_max_expr(
    project_expr: hl.expr.StringExpression,
    gt_expr: hl.expr.CallExpression,
    alleles_expr: hl.expr.ArrayExpression,
    n_projects: int = 5,
) -> hl.expr.ArrayExpression:
    """
    Create an expression that computes allele frequency information by project for the `n_projects` with the largest AF at this row.

    Will return an array with one element per non-reference allele.

    Each of these elements is itself an array of structs with the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32
        - project: str

    .. note::

        Only projects with AF > 0 are returned.
        In case of ties, the project ordering is not guaranteed, and at most `n_projects` are returned.

    :param project_expr: column expression containing the project
    :param gt_expr: entry expression containing the genotype
    :param alleles_expr: row expression containing the alleles
    :param n_projects: Maximum number of projects to return for each row
    :return: projectmax expression
    """
    n_alleles = hl.len(alleles_expr)

    # compute call stats by  project
    project_cs = hl.array(
        hl.agg.group_by(project_expr, hl.agg.call_stats(gt_expr, alleles_expr))
    )

    return hl.or_missing(
        n_alleles > 1,  # Exclude monomorphic sites
        hl.range(1, n_alleles).map(
            lambda ai: hl.sorted(
                project_cs.filter(
                    # filter to projects with AF > 0
                    lambda x: x[1].AF[ai]
                    > 0
                ),
                # order the callstats computed by AF in decreasing order
                lambda x: -x[1].AF[ai],
                # take the n_projects projects with largest AF
            )[:n_projects].map(
                # add the project in the callstats struct
                lambda x: x[1].annotate(
                    AC=x[1].AC[ai],
                    AF=x[1].AF[ai],
                    AN=x[1].AN,
                    homozygote_count=x[1].homozygote_count[ai],
                    project=x[0],
                )
            )
        ),
    )


def faf_expr(
    freq: hl.expr.ArrayExpression,
    freq_meta: hl.expr.ArrayExpression,
    locus: hl.expr.LocusExpression,
    pops_to_exclude: Optional[Set[str]] = None,
    faf_thresholds: List[float] = [0.95, 0.99],
) -> Tuple[hl.expr.ArrayExpression, List[Dict[str, str]]]:
    """
    Calculate the filtering allele frequency (FAF) for each threshold specified in `faf_thresholds`.

    See http://cardiodb.org/allelefrequencyapp/ for more information.

    The FAF is computed for each of the following population stratification if found in `freq_meta`:

        - All samples, with adj criteria
        - For each population, with adj criteria
        - For all sex/population on the non-PAR regions of sex chromosomes (will be missing on autosomes and PAR regions of sex chromosomes)

    Each of the FAF entry is a struct with one entry per threshold specified in `faf_thresholds` of type float64.

    This returns a tuple with two expressions:

        1. An array of FAF expressions as described above
        2. An array of dict containing the metadata for each of the array elements, in the same format as that produced by `annotate_freq`.

    :param freq: ArrayExpression of call stats structs (typically generated by hl.agg.call_stats)
    :param freq_meta: ArrayExpression of meta dictionaries corresponding to freq (typically generated using annotate_freq)
    :param locus: locus
    :param pops_to_exclude: Set of populations to exclude from faf calculation (typically bottlenecked or consanguineous populations)
    :param faf_thresholds: List of FAF thresholds to compute
    :return: (FAF expression, FAF metadata)
    """
    _pops_to_exclude = (
        hl.literal(pops_to_exclude)
        if pops_to_exclude is not None
        else hl.empty_set(hl.tstr)
    )

    # pylint: disable=invalid-unary-operand-type
    faf_freq_indices = hl.range(0, hl.len(freq_meta)).filter(
        lambda i: (freq_meta[i].get("group") == "adj")
        & (
            (freq_meta[i].size() == 1)
            | (
                (hl.set(freq_meta[i].keys()) == {"pop", "group"})
                & (~_pops_to_exclude.contains(freq_meta[i]["pop"]))
            )
        )
    )
    sex_faf_freq_indices = hl.range(0, hl.len(freq_meta)).filter(
        lambda i: (freq_meta[i].get("group") == "adj")
        & (freq_meta[i].contains("sex"))
        & (
            (freq_meta[i].size() == 2)
            | (
                (hl.set(freq_meta[i].keys()) == {"pop", "group", "sex"})
                & (~_pops_to_exclude.contains(freq_meta[i]["pop"]))
            )
        )
    )

    faf_expr = faf_freq_indices.map(
        lambda i: hl.struct(
            **{
                f"faf{str(threshold)[2:]}": hl.experimental.filtering_allele_frequency(
                    freq[i].AC, freq[i].AN, threshold
                )
                for threshold in faf_thresholds
            }
        )
    )

    faf_expr = faf_expr.extend(
        sex_faf_freq_indices.map(
            lambda i: hl.or_missing(
                ~locus.in_autosome_or_par(),
                hl.struct(
                    **{
                        f"faf{str(threshold)[2:]}": (
                            hl.experimental.filtering_allele_frequency(
                                freq[i].AC, freq[i].AN, threshold
                            )
                        )
                        for threshold in faf_thresholds
                    }
                ),
            )
        )
    )

    faf_meta = faf_freq_indices.extend(sex_faf_freq_indices).map(lambda i: freq_meta[i])
    return faf_expr, hl.eval(faf_meta)


def qual_hist_expr(
    gt_expr: Optional[hl.expr.CallExpression] = None,
    gq_expr: Optional[hl.expr.NumericExpression] = None,
    dp_expr: Optional[hl.expr.NumericExpression] = None,
    ad_expr: Optional[hl.expr.ArrayNumericExpression] = None,
    adj_expr: Optional[hl.expr.BooleanExpression] = None,
    ab_expr: Optional[hl.expr.NumericExpression] = None,
) -> hl.expr.StructExpression:
    """
    Return a struct expression with genotype quality histograms based on the arguments given (dp, gq, ad, ab).

    .. note::

        - If `gt_expr` is provided, will return histograms for non-reference samples only as well as all samples.
        - `gt_expr` is required for the allele-balance histogram, as it is only computed on het samples.
        - If `ab_expr` is provided, the allele-balance histogram is computed using this expression instead of the ad_expr.
        - If `adj_expr` is provided, additional histograms are computed using only adj samples.

    :param gt_expr: Entry expression containing genotype
    :param gq_expr: Entry expression containing genotype quality
    :param dp_expr: Entry expression containing depth
    :param ad_expr: Entry expression containing allelic depth (bi-allelic here)
    :param adj_expr: Entry expression containing adj (high quality) genotype status
    :param ab_expr: Entry expression containing allele balance (bi-allelic here)
    :return: Genotype quality histograms expression
    """
    qual_hists = {}
    if gq_expr is not None:
        qual_hists["gq_hist"] = hl.agg.hist(gq_expr, 0, 100, 20)
    if dp_expr is not None:
        qual_hists["dp_hist"] = hl.agg.hist(dp_expr, 0, 100, 20)

    if gt_expr is not None:
        qual_hists = {
            **{
                f"{qual_hist_name}_all": qual_hist_expr
                for qual_hist_name, qual_hist_expr in qual_hists.items()
            },
            **{
                f"{qual_hist_name}_alt": hl.agg.filter(
                    gt_expr.is_non_ref(), qual_hist_expr
                )
                for qual_hist_name, qual_hist_expr in qual_hists.items()
            },
        }
        ab_hist_msg = "Using the %s to compute allele balance histogram..."
        if ab_expr is not None:
            logger.info(ab_hist_msg, "ab_expr")
            qual_hists["ab_hist_alt"] = hl.agg.filter(
                gt_expr.is_het(), hl.agg.hist(ab_expr, 0, 1, 20)
            )
        elif ad_expr is not None:
            logger.info(ab_hist_msg, "ad_expr")
            qual_hists["ab_hist_alt"] = hl.agg.filter(
                gt_expr.is_het(), hl.agg.hist(ad_expr[1] / hl.sum(ad_expr), 0, 1, 20)
            )

    else:
        qual_hists = {
            f"{qual_hist_name}_all": qual_hist_expr
            for qual_hist_name, qual_hist_expr in qual_hists.items()
        }

    if adj_expr is not None:
        qual_hists.update(
            {
                f"{qual_hist_name}_adj": hl.agg.filter(adj_expr, qual_hist_expr)
                for qual_hist_name, qual_hist_expr in qual_hists.items()
            }
        )

    return hl.struct(**qual_hists)


def age_hists_expr(
    adj_expr: hl.expr.BooleanExpression,
    gt_expr: hl.expr.CallExpression,
    age_expr: hl.expr.NumericExpression,
    lowest_boundary: int = 30,
    highest_boundary: int = 80,
    n_bins: int = 10,
) -> hl.expr.StructExpression:
    """
    Return a StructExpression with the age histograms for hets and homs.

    :param adj_expr: Entry expression containing whether a genotype is high quality (adj) or not
    :param gt_expr: Entry expression containing the genotype
    :param age_expr: Col expression containing the sample's age
    :param lowest_boundary: Lowest bin boundary (any younger sample will be binned in n_smaller)
    :param highest_boundary: Highest bin boundary (any older sample will be binned in n_larger)
    :param n_bins: Total number of bins
    :return: A struct with `age_hist_het` and `age_hist_hom`
    """
    return hl.struct(
        age_hist_het=hl.agg.filter(
            adj_expr & gt_expr.is_het(),
            hl.agg.hist(age_expr, lowest_boundary, highest_boundary, n_bins),
        ),
        age_hist_hom=hl.agg.filter(
            adj_expr & gt_expr.is_hom_var(),
            hl.agg.hist(age_expr, lowest_boundary, highest_boundary, n_bins),
        ),
    )


def annotate_freq(
    mt: hl.MatrixTable,
    sex_expr: Optional[hl.expr.StringExpression] = None,
    pop_expr: Optional[hl.expr.StringExpression] = None,
    subpop_expr: Optional[hl.expr.StringExpression] = None,
    additional_strata_expr: Optional[
        Union[
            List[Dict[str, hl.expr.StringExpression]],
            Dict[str, hl.expr.StringExpression],
        ]
    ] = None,
    downsamplings: Optional[List[int]] = None,
    ds_pop_counts: Optional[Dict[str, int]] = None,
    downsampling_expr: Optional[hl.expr.StructExpression] = None,
    entry_agg_funcs: Optional[Dict[str, Tuple[Callable, Callable]]] = None,
    annotate_mt: bool = True,
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Annotate `mt` with stratified allele frequencies.

    The output Matrix table will include:
        - row annotation `freq` containing the stratified allele frequencies
        - global annotation `freq_meta` with metadata
        - global annotation `freq_sample_count` with sample count information

    .. note::

        Currently this only supports bi-allelic sites.

        The input `mt` needs to have the following entry fields:
          - GT: a CallExpression containing the genotype
          - adj: a BooleanExpression containing whether the genotype is of high quality
            or not.

        All expressions arguments need to be expression on the input `mt`.

    .. rubric:: `freq` row annotation

    The `freq` row annotation is an Array of Structs, with each Struct containing the
    following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data, and the
    metadata about these annotations is stored in the globals.

    .. rubric:: Global `freq_meta` metadata annotation

    The global annotation `freq_meta` is added to the input `mt`. It is a list of dict.
    Each element of the list contains metadata on a frequency stratification and the
    index in the list corresponds to the index of that frequency stratification in the
    `freq` row annotation.

    .. rubric:: Global `freq_sample_count` annotation

    The global annotation `freq_sample_count` is added to the input `mt`. This is a
    sample count per sample grouping defined in the `freq_meta` global annotation.

    .. rubric:: The `downsamplings` parameter

    If the `downsamplings` parameter is used, frequencies will be computed for all
    samples and by population (if `pop_expr` is specified) by downsampling the number
    of samples without replacement to each of the numbers specified in the
    `downsamplings` array, provided that there are enough samples in the dataset. In
    addition, if `pop_expr` is specified, a downsampling to each of the exact number
    of samples present in each population is added. Note that samples are randomly
    sampled only once, meaning that the lower downsamplings are subsets of the higher
    ones.

    .. rubric:: The `additional_strata_expr` parameter

    If the `additional_strata_expr` parameter is used, frequencies will be computed for
    each of the strata dictionaries across all values. For example, if
    `additional_strata_expr` is set to `[{'platform': mt.platform},
    {'platform':mt.platform, 'pop': mt.pop}, {'age_bin': mt.age_bin}]`, then
    frequencies will be computed for each of the values of `mt.platform`, each of the
    combined values of `mt.platform` and `mt.pop`, and each of the values of
    `mt.age_bin`.

    :param mt: Input MatrixTable
    :param sex_expr: When specified, frequencies are stratified by sex. If `pop_expr`
        is also specified, then a pop/sex stratifiction is added.
    :param pop_expr: When specified, frequencies are stratified by population. If
        `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param subpop_expr: When specified, frequencies are stratified by sub-continental
        population. Note that `pop_expr` is required as well when using this option.
    :param additional_strata_expr: When specified, frequencies are stratified by the
        given additional strata. This can e.g. be used to stratify by platform,
        platform-pop, platform-pop-sex.
    :param downsamplings: When specified, frequencies are computed by downsampling the
        data to the number of samples given in the list. Note that if `pop_expr` is
        specified, downsamplings by population is also computed.
    :param ds_pop_counts: When specified, frequencies are computed by downsampling the
        data to the number of samples per pop in the dict. The key is the population
        and the value is the number of samples.
    :param downsampling_expr: When specified, frequencies are computed using the
        downsampling indices in the provided StructExpression. Note that if `pop_idx`
        is specified within the struct, downsamplings by population is also computed.
    :param entry_agg_funcs: When specified, additional annotations are added to the
        output Table/MatrixTable. The keys of the dict are the names of the annotations
        and the values are tuples of functions. The first function is used to transform
        the `mt` entries in some way, and the second function is used to aggregate the
        output from the first function.
    :param annotate_mt: Whether to return the full MatrixTable with annotations added
        instead of only a Table with `freq` and other annotations. Default is True.
    :return: MatrixTable or Table with `freq` annotation.
    """
    errors = []
    if subpop_expr is not None and pop_expr is None:
        errors.append("annotate_freq requires pop_expr when using subpop_expr")

    if downsampling_expr is not None:
        if downsamplings is None:
            errors.append(
                "annotate_freq requires `downsamplings` when using `downsampling_expr`"
            )
        if downsampling_expr.get("global_idx") is None:
            errors.append(
                "annotate_freq requires `downsampling_expr` with key 'global_idx'"
            )
        if downsampling_expr.get("pop_idx") is None:
            if pop_expr is not None:
                errors.append(
                    "annotate_freq requires `downsampling_expr` with key 'pop_idx' when"
                    " using `pop_expr`"
                )
        else:
            if pop_expr is None or ds_pop_counts is None:
                errors.append(
                    "annotate_freq requires `pop_expr` and `ds_pop_counts` when using"
                    " `downsampling_expr` with pop_idx"
                )

    if errors:
        raise ValueError("The following errors were found: \n" + "\n".join(errors))

    # Build list of strata expressions based on supplied parameters.
    strata_expr = []
    if pop_expr is not None:
        strata_expr.append({"pop": pop_expr})
    if sex_expr is not None:
        strata_expr.append({"sex": sex_expr})
        if pop_expr is not None:
            strata_expr.append({"pop": pop_expr, "sex": sex_expr})
    if subpop_expr is not None:
        strata_expr.append({"pop": pop_expr, "subpop": subpop_expr})

    # Add downsampling to strata expressions, include pop in the strata if supplied.
    if downsampling_expr is not None:
        downsampling_strata = {"downsampling": downsampling_expr}
        if pop_expr is not None:
            downsampling_strata["pop"] = pop_expr
        strata_expr.append(downsampling_strata)

    # Add additional strata expressions.
    if additional_strata_expr is not None:
        if isinstance(additional_strata_expr, dict):
            additional_strata_expr = [additional_strata_expr]
        strata_expr.extend(additional_strata_expr)

    # Annotate mt with all annotations used in the strata expression list before any
    # modifications are made to mt.
    mt = mt.annotate_cols(**{k: v for d in strata_expr for k, v in d.items()})

    # Get downsampling_expr if it is None, but downsamplings is supplied.
    if downsamplings is not None and downsampling_expr is None:
        mt = annotate_downsamplings(
            mt, downsamplings, pop_expr=None if pop_expr is None else mt.pop
        )
        downsamplings = hl.eval(mt.downsamplings)
        ds_pop_counts = hl.eval(mt.ds_pop_counts)
        downsampling_strata = {"downsampling": mt.downsampling}
        if pop_expr is not None:
            downsampling_strata["pop"] = pop_expr
        strata_expr.append(downsampling_strata)

    strata_expr = [{k: mt[k] for k in d} for d in strata_expr]

    ht = compute_freq_by_strata(
        mt,
        strata_expr,
        downsamplings=downsamplings,
        ds_pop_counts=ds_pop_counts,
        entry_agg_funcs=entry_agg_funcs,
    )

    if annotate_mt:
        return mt.annotate_rows(**ht[mt.row_key])
    else:
        return ht


def get_lowqual_expr(
    alleles: hl.expr.ArrayExpression,
    qual_approx_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
    snv_phred_threshold: int = 30,
    snv_phred_het_prior: int = 30,  # 1/1000
    indel_phred_threshold: int = 30,
    indel_phred_het_prior: int = 39,  # 1/8,000
) -> Union[hl.expr.BooleanExpression, hl.expr.ArrayExpression]:
    """
    Compute lowqual threshold expression for either split or unsplit alleles based on QUALapprox or AS_QUALapprox.

    .. note::

        When running This lowqual annotation using QUALapprox, it differs from the GATK LowQual filter.
        This is because GATK computes this annotation at the site level, which uses the least stringent prior for mixed sites.
        When run using AS_QUALapprox, this implementation can thus be more stringent for certain alleles at mixed sites.

    :param alleles: Array of alleles
    :param qual_approx_expr: QUALapprox or AS_QUALapprox
    :param snv_phred_threshold: Phred-scaled SNV "emission" threshold (similar to GATK emission threshold)
    :param snv_phred_het_prior: Phred-scaled SNV heterozygosity prior (30 = 1/1000 bases, GATK default)
    :param indel_phred_threshold: Phred-scaled indel "emission" threshold (similar to GATK emission threshold)
    :param indel_phred_het_prior: Phred-scaled indel heterozygosity prior (30 = 1/1000 bases, GATK default)
    :return: lowqual expression (BooleanExpression if `qual_approx_expr`is Numeric, Array[BooleanExpression] if `qual_approx_expr` is ArrayNumeric)
    """
    min_snv_qual = snv_phred_threshold + snv_phred_het_prior
    min_indel_qual = indel_phred_threshold + indel_phred_het_prior
    min_mixed_qual = max(min_snv_qual, min_indel_qual)

    if isinstance(qual_approx_expr, hl.expr.ArrayNumericExpression):
        return hl.range(1, hl.len(alleles)).map(
            lambda ai: hl.cond(
                hl.is_snp(alleles[0], alleles[ai]),
                qual_approx_expr[ai - 1] < min_snv_qual,
                qual_approx_expr[ai - 1] < min_indel_qual,
            )
        )
    else:
        return (
            hl.case()
            .when(
                hl.range(1, hl.len(alleles)).all(
                    lambda ai: hl.is_snp(alleles[0], alleles[ai])
                ),
                qual_approx_expr < min_snv_qual,
            )
            .when(
                hl.range(1, hl.len(alleles)).all(
                    lambda ai: hl.is_indel(alleles[0], alleles[ai])
                ),
                qual_approx_expr < min_indel_qual,
            )
            .default(qual_approx_expr < min_mixed_qual)
        )


def get_annotations_hists(
    ht: hl.Table,
    annotations_hists: Dict[str, Tuple],
    log10_annotations: List[str] = ["DP"],
) -> Dict[str, hl.expr.StructExpression]:
    """
    Create histograms for variant metrics in ht.info.

    Used when creating site quality distribution json files.

    :param ht: Table with variant metrics
    :param annotations_hists: Dictionary of metrics names and their histogram values (start, end, bins)
    :param log10_annotations: List of metrics to log scale
    :return: Dictionary of merics and their histograms
    :rtype: Dict[str, hl.expr.StructExpression]
    """
    # Check all fields in ht.info and create histograms if they are in
    # annotations_hists dict
    return {
        field: hl.agg.hist(
            hl.log10(ht.info[field]) if field in log10_annotations else ht.info[field],
            start,
            end,
            bins,
        )
        for field, (start, end, bins) in annotations_hists.items()
        if field in ht.row.info
    }


def create_frequency_bins_expr(
    AC: hl.expr.NumericExpression, AF: hl.expr.NumericExpression
) -> hl.expr.StringExpression:
    """
    Create bins for frequencies in preparation for aggregating QUAL by frequency bin.

    Bins:
        - singleton
        - doubleton
        - 0.00005
        - 0.0001
        - 0.0002
        - 0.0005
        - 0.001,
        - 0.002
        - 0.005
        - 0.01
        - 0.02
        - 0.05
        - 0.1
        - 0.2
        - 0.5
        - 1

    NOTE: Frequencies should be frequencies from raw data.
    Used when creating site quality distribution json files.

    :param AC: Field in input that contains the allele count information
    :param AF: Field in input that contains the allele frequency information
    :return: Expression containing bin name
    :rtype: hl.expr.StringExpression
    """
    bin_expr = (
        hl.case()
        .when(AC == 1, "binned_singleton")
        .when(AC == 2, "binned_doubleton")
        .when((AC > 2) & (AF < 0.00005), "binned_0.00005")
        .when((AF >= 0.00005) & (AF < 0.0001), "binned_0.0001")
        .when((AF >= 0.0001) & (AF < 0.0002), "binned_0.0002")
        .when((AF >= 0.0002) & (AF < 0.0005), "binned_0.0005")
        .when((AF >= 0.0005) & (AF < 0.001), "binned_0.001")
        .when((AF >= 0.001) & (AF < 0.002), "binned_0.002")
        .when((AF >= 0.002) & (AF < 0.005), "binned_0.005")
        .when((AF >= 0.005) & (AF < 0.01), "binned_0.01")
        .when((AF >= 0.01) & (AF < 0.02), "binned_0.02")
        .when((AF >= 0.02) & (AF < 0.05), "binned_0.05")
        .when((AF >= 0.05) & (AF < 0.1), "binned_0.1")
        .when((AF >= 0.1) & (AF < 0.2), "binned_0.2")
        .when((AF >= 0.2) & (AF < 0.5), "binned_0.5")
        .when((AF >= 0.5) & (AF <= 1), "binned_1")
        .default(hl.null(hl.tstr))
    )
    return bin_expr


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.

    Defaults correspond to gnomAD values.
    """
    return (
        (gq_expr >= adj_gq)
        & hl.cond(gt_expr.is_haploid(), dp_expr >= haploid_adj_dp, dp_expr >= adj_dp)
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            .default(
                (ad_expr[gt_expr[0]] / dp_expr >= adj_ab)
                & (ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            )
        )
    )


def annotate_adj(
    mt: hl.MatrixTable,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid).

    Defaults correspond to gnomAD values.
    """
    if "GT" not in mt.entry and "LGT" in mt.entry:
        logger.warning("No GT field found, using LGT instead.")
        gt_expr = mt.LGT
    else:
        gt_expr = mt.GT

    if "AD" not in mt.entry and "LAD" in mt.entry:
        logger.warning("No AD field found, using LAD instead.")
        ad_expr = mt.LAD
    else:
        ad_expr = mt.AD

    return mt.annotate_entries(
        adj=get_adj_expr(
            gt_expr, mt.GQ, mt.DP, ad_expr, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    )


def add_variant_type(alt_alleles: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """Get Struct of variant_type and n_alt_alleles from ArrayExpression of Strings."""
    ref = alt_alleles[0]
    alts = alt_alleles[1:]
    non_star_alleles = hl.filter(lambda a: a != "*", alts)
    return hl.struct(
        variant_type=hl.cond(
            hl.all(lambda a: hl.is_snp(ref, a), non_star_alleles),
            hl.cond(hl.len(non_star_alleles) > 1, "multi-snv", "snv"),
            hl.cond(
                hl.all(lambda a: hl.is_indel(ref, a), non_star_alleles),
                hl.cond(hl.len(non_star_alleles) > 1, "multi-indel", "indel"),
                "mixed",
            ),
        ),
        n_alt_alleles=hl.len(non_star_alleles),
    )


def annotate_allele_info(ht: hl.Table) -> hl.Table:
    """
    Return bi-allelic sites Table with an 'allele_info' annotation.

    .. note::

        This function requires that the input `ht` is unsplit and returns a split `ht`.

    'allele_info' is a struct with the following information:
        - variant_type: Variant type (snv, indel, multi-snv, multi-indel, or mixed).
        - n_alt_alleles: Total number of alternate alleles observed at variant locus.
        - has_star: True if the variant contains a star allele.
        - allele_type: Allele type (snv, insertion, deletion, or mixed).
        - was_mixed: True if the variant was mixed (i.e. contained both SNVs and indels).
        - nonsplit_alleles: Array of alleles before splitting.

    :param Table ht: Unsplit input Table.
    :return: Split Table with allele data annotation added,
    """
    ht = ht.annotate(
        allele_info=hl.struct(
            **add_variant_type(ht.alleles),
            has_star=hl.any(lambda a: a == "*", ht.alleles),
        )
    )

    ht = hl.split_multi(ht)

    ref_expr = ht.alleles[0]
    alt_expr = ht.alleles[1]
    allele_type_expr = (
        hl.case()
        .when(hl.is_snp(ref_expr, alt_expr), "snv")
        .when(hl.is_insertion(ref_expr, alt_expr), "ins")
        .when(hl.is_deletion(ref_expr, alt_expr), "del")
        .default("complex")
    )
    ht = ht.transmute(
        allele_info=ht.allele_info.annotate(
            allele_type=allele_type_expr,
            was_mixed=ht.allele_info.variant_type == "mixed",
            nonsplit_alleles=ht.old_alleles,
        )
    )

    return ht


def annotation_type_is_numeric(t: Any) -> bool:
    """
    Given an annotation type, return whether it is a numerical type or not.

    :param t: Type to test
    :return: If the input type is numeric
    """
    return t in (hl.tint32, hl.tint64, hl.tfloat32, hl.tfloat64)


def annotation_type_in_vcf_info(t: Any) -> bool:
    """
    Given an annotation type, returns whether that type can be natively exported to a VCF INFO field.

    .. note::

        Types that aren't natively exportable to VCF will be converted to String on export.

    :param t: Type to test
    :return: If the input type can be exported to VCF
    """
    return (
        annotation_type_is_numeric(t)
        or t in (hl.tstr, hl.tbool)
        or (
            isinstance(t, (hl.tarray, hl.tset))
            and annotation_type_in_vcf_info(t.element_type)
        )
    )


def bi_allelic_site_inbreeding_expr(
    call: Optional[hl.expr.CallExpression] = None,
    callstats_expr: Optional[hl.expr.StructExpression] = None,
) -> hl.expr.Float32Expression:
    """
    Return the site inbreeding coefficient as an expression to be computed on a MatrixTable.

    This is implemented based on the GATK InbreedingCoeff metric:
    https://software.broadinstitute.org/gatk/documentation/article.php?id=8032

    .. note::

        The computation is run based on the counts of alternate alleles and thus should only be run on bi-allelic sites.

    :param call: Expression giving the calls in the MT
    :param callstats_expr: StructExpression containing only alternate allele AC, AN, and homozygote_count as integers. If passed, used to create expression in place of GT calls.
    :return: Site inbreeding coefficient expression
    """
    if call is None and callstats_expr is None:
        raise ValueError("One of `call` or `callstats_expr` must be passed.")

    def inbreeding_coeff(
        gt_counts: hl.expr.DictExpression,
    ) -> hl.expr.Float32Expression:
        n = gt_counts.get(0, 0) + gt_counts.get(1, 0) + gt_counts.get(2, 0)
        p = (2 * gt_counts.get(0, 0) + gt_counts.get(1, 0)) / (2 * n)
        q = (2 * gt_counts.get(2, 0) + gt_counts.get(1, 0)) / (2 * n)
        return 1 - (gt_counts.get(1, 0) / (2 * p * q * n))

    if callstats_expr is not None:
        # Check that AC, AN, and homozygote count are all ints
        if not (
            (
                (callstats_expr.AC.dtype == hl.tint32)
                | (callstats_expr.AC.dtype == hl.tint64)
            )
            & (
                (callstats_expr.AN.dtype == hl.tint32)
                | (callstats_expr.AN.dtype == hl.tint64)
            )
            & (
                (callstats_expr.homozygote_count.dtype == hl.tint32)
                | (callstats_expr.homozygote_count.dtype == hl.tint64)
            )
        ):
            raise ValueError(
                "callstats_expr must be a StructExpression containing fields 'AC',"
                " 'AN', and 'homozygote_count' of types int32 or int64."
            )
        n = callstats_expr.AN / 2
        q = callstats_expr.AC / callstats_expr.AN
        p = 1 - q
        return 1 - (callstats_expr.AC - (2 * callstats_expr.homozygote_count)) / (
            2 * p * q * n
        )
    else:
        return hl.bind(inbreeding_coeff, hl.agg.counter(call.n_alt_alleles()))


def fs_from_sb(
    sb: Union[hl.expr.ArrayNumericExpression, hl.expr.ArrayExpression],
    normalize: bool = True,
    min_cell_count: int = 200,
    min_count: int = 4,
    min_p_value: float = 1e-320,
) -> hl.expr.Int64Expression:
    """
    Compute `FS` (Fisher strand balance) annotation from  the `SB` (strand balance table) field.

    `FS` is the phred-scaled value of the double-sided Fisher exact test on strand balance.

    Using default values will have the same behavior as the GATK implementation, that is:
    - If sum(counts) > 2*`min_cell_count` (default to GATK value of 200), they are normalized
    - If sum(counts) < `min_count` (default to GATK value of 4), returns missing
    - Any p-value < `min_p_value` (default to GATK value of 1e-320) is truncated to that value

    In addition to the default GATK behavior, setting `normalize` to `False` will perform a chi-squared test
    for large counts (> `min_cell_count`) instead of normalizing the cell values.

    .. note::

        This function can either take
        - an array of length four containing the forward and reverse strands' counts of ref and alt alleles: [ref fwd, ref rev, alt fwd, alt rev]
        - a two dimensional array with arrays of length two, containing the counts: [[ref fwd, ref rev], [alt fwd, alt rev]]

    GATK code here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/FisherStrand.java

    :param sb: Count of ref/alt reads on each strand
    :param normalize: Whether to normalize counts is sum(counts) > min_cell_count (normalize=True), or use a chi sq instead of FET (normalize=False)
    :param min_cell_count: Maximum count for performing a FET
    :param min_count: Minimum total count to output FS (otherwise null it output)
    :return: FS value
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(lambda x: hl.flatten(x), sb)

    sb_sum = hl.bind(lambda x: hl.sum(x), sb)

    # Normalize table if counts get too large
    if normalize:
        fs_expr = hl.bind(
            lambda sb, sb_sum: hl.cond(
                sb_sum <= 2 * min_cell_count,
                sb,
                sb.map(lambda x: hl.int(x / (sb_sum / min_cell_count))),
            ),
            sb,
            sb_sum,
        )

        # FET
        fs_expr = to_phred(
            hl.max(
                hl.fisher_exact_test(
                    fs_expr[0], fs_expr[1], fs_expr[2], fs_expr[3]
                ).p_value,
                min_p_value,
            )
        )
    else:
        fs_expr = to_phred(
            hl.max(
                hl.contingency_table_test(
                    sb[0], sb[1], sb[2], sb[3], min_cell_count=min_cell_count
                ).p_value,
                min_p_value,
            )
        )

    # Return null if counts <= `min_count`
    return hl.or_missing(
        sb_sum > min_count, hl.max(0, fs_expr)  # Needed to avoid -0.0 values
    )


def sor_from_sb(
    sb: Union[hl.expr.ArrayNumericExpression, hl.expr.ArrayExpression]
) -> hl.expr.Float64Expression:
    """
    Compute `SOR` (Symmetric Odds Ratio test) annotation from  the `SB` (strand balance table) field.

    .. note::

        This function can either take
        - an array of length four containing the forward and reverse strands' counts of ref and alt alleles: [ref fwd, ref rev, alt fwd, alt rev]
        - a two dimensional array with arrays of length two, containing the counts: [[ref fwd, ref rev], [alt fwd, alt rev]]

    GATK code here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandOddsRatio.java

    :param sb: Count of ref/alt reads on each strand
    :return: SOR value
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(lambda x: hl.flatten(x), sb)

    sb = sb.map(lambda x: hl.float64(x) + 1)

    ref_fw = sb[0]
    ref_rv = sb[1]
    alt_fw = sb[2]
    alt_rv = sb[3]
    symmetrical_ratio = ((ref_fw * alt_rv) / (alt_fw * ref_rv)) + (
        (alt_fw * ref_rv) / (ref_fw * alt_rv)
    )
    ref_ratio = hl.min(ref_rv, ref_fw) / hl.max(ref_rv, ref_fw)
    alt_ratio = hl.min(alt_fw, alt_rv) / hl.max(alt_fw, alt_rv)
    sor = hl.log(symmetrical_ratio) + hl.log(ref_ratio) - hl.log(alt_ratio)

    return sor


def pab_max_expr(
    gt_expr: hl.expr.CallExpression,
    ad_expr: hl.expr.ArrayExpression,
    la_expr: Optional[hl.expr.ArrayExpression] = None,
    n_alleles_expr: Optional[hl.expr.Int32Expression] = None,
) -> hl.expr.ArrayExpression:
    """
    Compute the maximum p-value of the binomial test for the alternate allele balance (PAB) for each allele.

    .. note::

        This function can take a `gt_expr` and `ad_expr` that use local or global
        alleles. If they use local alleles, `la_expr` and `n_alleles_expr` should be
        provided to transform `gt_expr` and `ad_expr` to global alleles.

    :param gt_expr: Genotype call expression.
    :param ad_expr: Allele depth expression.
    :param la_expr: Allele local index expression. When provided `gt_expr` and
        `ad_expr` are transformed from using local alleles to global alleles using
        `la_expr`.
    :param n_alleles_expr: Number of alleles expression. Required when 'la_expr' is
        provided.
    :return: Array expression of maximum p-values.
    """
    if la_expr is not None:
        if n_alleles_expr is None:
            raise ValueError("Must provide `n_alleles_expr` if `la_expr` is provided!")

        ad_expr = hl.vds.local_to_global(
            ad_expr, la_expr, n_alleles_expr, fill_value=0, number="R"
        )
        gt_expr = hl.vds.lgt_to_gt(gt_expr, la_expr)

    expr = hl.agg.array_agg(
        lambda x: hl.agg.filter(
            gt_expr.is_het(),
            hl.agg.max(hl.binom_test(x, hl.sum(ad_expr), 0.5, "two-sided")),
        ),
        ad_expr[1:],  # Skip ref allele
    )

    return expr


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Return a boolean expression selecting bi-allelic sites only, accounting for whether the input MT/HT was split.

    :param t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    """
    return ~t.was_split if "was_split" in t.row else (hl.len(t.alleles) == 2)


def unphase_call_expr(call_expr: hl.expr.CallExpression) -> hl.expr.CallExpression:
    """
    Generate unphased version of a call expression (which can be phased or not).

    :param call_expr: Input call expression
    :return: unphased call expression
    """
    return (
        hl.case()
        .when(call_expr.is_diploid(), hl.call(call_expr[0], call_expr[1], phased=False))
        .when(call_expr.is_haploid(), hl.call(call_expr[0], phased=False))
        .default(hl.null(hl.tcall))
    )


def region_flag_expr(
    t: Union[hl.Table, hl.MatrixTable],
    non_par: bool = True,
    prob_regions: Dict[str, hl.Table] = None,
) -> hl.expr.StructExpression:
    """
    Create a `region_flag` struct that contains flags for problematic regions (i.e., LCR, decoy, segdup, and nonpar regions).

    .. note:: No hg38 resources for decoy or self chain are available yet.

    :param t: Input Table/MatrixTable
    :param non_par: If True, flag loci that occur within pseudoautosomal regions on sex chromosomes
    :param prob_regions: If supplied, flag loci that occur within regions defined in Hail Table(s)
    :return: `region_flag` struct row annotation
    """
    prob_flags_expr = (
        {"non_par": (t.locus.in_x_nonpar() | t.locus.in_y_nonpar())
         } if non_par else {}  # fmt: skip
    )

    if prob_regions is not None:
        prob_flags_expr.update(
            {
                region_name: hl.is_defined(region_table[t.locus])
                for region_name, region_table in prob_regions.items()
            }
        )

    return hl.struct(**prob_flags_expr)


def missing_callstats_expr() -> hl.expr.StructExpression:
    """
    Create a missing callstats struct for insertion into frequency annotation arrays when data is missing.

    :return: Hail Struct with missing values for each callstats element
    """
    return hl.struct(
        AC=hl.missing(hl.tint32),
        AF=hl.missing(hl.tfloat64),
        AN=hl.missing(hl.tint32),
        homozygote_count=hl.missing(hl.tint32),
    )


def set_female_y_metrics_to_na_expr(
    t: Union[hl.Table, hl.MatrixTable]
) -> hl.expr.ArrayExpression:
    """
    Set Y-variant frequency callstats for female-specific metrics to missing structs.

    .. note:: Requires freq, freq_meta, and freq_index_dict annotations to be present in Table or MatrixTable

    :param t: Table or MatrixTable for which to adjust female metrics
    :return: Hail array expression to set female Y-variant metrics to missing values
    """
    female_idx = hl.map(
        lambda x: t.freq_index_dict[x],
        hl.filter(lambda x: x.contains("XX"), t.freq_index_dict.keys()),
    )
    freq_idx_range = hl.range(hl.len(t.freq_meta))

    new_freq_expr = hl.if_else(
        (t.locus.in_y_nonpar() | t.locus.in_y_par()),
        hl.map(
            lambda x: hl.if_else(
                female_idx.contains(x), missing_callstats_expr(), t.freq[x]
            ),
            freq_idx_range,
        ),
        t.freq,
    )

    return new_freq_expr


def hemi_expr(
    locus: hl.expr.LocusExpression,
    sex_expr: hl.expr.StringExpression,
    gt: hl.expr.CallExpression,
    male_str: str = "XY",
) -> hl.expr.BooleanExpression:
    """
    Return whether genotypes are hemizygous.

    Return missing expression if locus is not in chrX/chrY non-PAR regions.

    :param locus: Input locus.
    :param sex_expr: Input StringExpression indicating whether sample is XX or XY.
    :param gt: Input genotype.
    :param xy_str: String indicating whether sample is XY. Default is "XY".
    :return: BooleanExpression indicating whether genotypes are hemizygous.
    """
    return hl.or_missing(
        locus.in_x_nonpar() | locus.in_y_nonpar(),
        # Haploid genotypes have a single integer, so checking if
        # mt.GT[0] is alternate allele
        gt.is_haploid() & (sex_expr == male_str) & (gt[0] == 1),
    )


def merge_freq_arrays(
    farrays: List[hl.expr.ArrayExpression],
    fmeta: List[List[Dict[str, str]]],
    operation: str = "sum",
    set_negatives_to_zero: bool = False,
) -> Tuple[hl.expr.ArrayExpression, List[Dict[str, int]]]:
    """
    Merge a list of frequency arrays based on the supplied `operation`.

    .. warning::
        Arrays must be on the same Table.

    .. note::

        Arrays do not have to contain the same groupings or order of groupings but
        the array indices for a freq array in `farrays` must be the same as its associated
        frequency metadata index in `fmeta` i.e., `farrays = [freq1, freq2]` then `fmeta`
        must equal `[fmeta1, fmeta2]` where fmeta1 contains the metadata information
        for freq1.

        If `operation` is set to "sum", groups in the merged array
        will be the union of groupings found within the arrays' metadata and all arrays
        with be summed by grouping. If `operation` is set to "diff", the merged array
        will contain groups only found in the first array of `fmeta`. Any array containing
        any of these groups will have thier values subtracted from the values of the first array.

    :param farrays: List of frequency arrays to merge. First entry in the list is the primary array to which other arrays will be added or subtracted. All arrays must be on the same Table.
    :param fmeta: List of frequency metadata for arrays being merged.
    :param operation: Merge operation to perform. Options are "sum" and "diff". If "diff" is passed, the first freq array in the list will have the other arrays subtracted from it.
    :param set_negatives_to_zero: If True, set negative array values to 0 for AC, AN, AF, and homozygote_count. If False, raise a ValueError. Default is True.
    :return: Tuple of merged frequency array and its frequency metadata list.
    """
    if len(farrays) < 2:
        raise ValueError("Must provide at least two frequency arrays to merge!")
    if len(farrays) != len(fmeta):
        raise ValueError("Length of farrays and fmeta must be equal!")
    if operation not in ["sum", "diff"]:
        raise ValueError("Operation must be either 'sum' or 'diff'!")

    # Create a list where each entry is a dictionary whose key is an aggregation
    # group and the value is the corresponding index in the freq array.
    fmeta = [hl.dict(hl.enumerate(f).map(lambda x: (x[1], [x[0]]))) for f in fmeta]

    # Merge dictionaries in the list into a single dictionary where key is aggregation
    # group and the value is a list of the group's index in each of the freq arrays, if
    # it exists. For "sum" operation, use keys, aka groups, found in all freq dictionaries.
    # For "diff" operations, only use key_set from the first entry.
    fmeta = hl.fold(
        lambda i, j: hl.dict(
            (
                hl.if_else(operation == "sum", (i.key_set() | j.key_set()), i.key_set())
            ).map(
                lambda k: (
                    k,
                    i.get(k, [hl.missing(hl.tint32)]).extend(
                        j.get(k, [hl.missing(hl.tint32)])
                    ),
                )
            )
        ),
        fmeta[0],
        fmeta[1:],
    )

    # Create a list of tuples from the dictionary, sorted by the list of indices for
    # each aggregation group.
    fmeta = hl.sorted(fmeta.items(), key=lambda f: f[1])

    # Create a list of the aggregation groups, maintaining the sorted order.
    new_freq_meta = fmeta.map(lambda x: x[0])

    # Create array for each aggregation group of arrays containing the group's freq
    # values from each freq array.
    freq_meta_idx = fmeta.map(lambda x: hl.zip(farrays, x[1]).map(lambda i: i[0][i[1]]))

    def _sum_or_diff_fields(
        field_1_expr: str, field_2_expr: str
    ) -> hl.expr.Int32Expression:
        """
        Sum or subtract fields in call statistics struct.

        :param field_1_expr: First field to sum or diff.
        :param field_2_expr: Second field to sum or diff.
        :return: Merged field value.
        """
        return hl.if_else(
            operation == "sum",
            hl.or_else(field_1_expr, 0) + hl.or_else(field_2_expr, 0),
            hl.or_else(field_1_expr, 0) - hl.or_else(field_2_expr, 0),
        )

    # Iterate through the groups and their freq lists to merge callstats.
    callstat_ann = ["AC", "AN", "homozygote_count"]
    new_freq = freq_meta_idx.map(
        lambda x: hl.bind(
            lambda y: y.annotate(AF=hl.if_else(y.AN > 0, y.AC / y.AN, 0)),
            hl.fold(
                lambda i, j: hl.struct(
                    **{ann: _sum_or_diff_fields(i[ann], j[ann]) for ann in callstat_ann}
                ),
                x[0].select("AC", "AN", "homozygote_count"),
                x[1:],
            ),
        )
    )
    # Check and see if any annotation within the merged array is negative. If so,
    # raise an error if set_negatives_to_zero is False or set the value to 0 if
    # set_negatives_to_zero is True.
    if operation == "diff":
        negative_value_error_msg = (
            "Negative values found in merged frequency array. Review data or set"
            " `set_negatives_to_zero` to True to set negative values to 0."
        )
        callstat_ann.append("AF")
        new_freq = new_freq.map(
            lambda x: x.annotate(
                **{
                    ann: (
                        hl.case()
                        .when(set_negatives_to_zero, hl.max(x[ann], 0))
                        .or_error(negative_value_error_msg)
                    )
                    for ann in callstat_ann
                }
            )
        )

    return new_freq, new_freq_meta


def annotate_downsamplings(
    t: Union[hl.MatrixTable, hl.Table],
    downsamplings: List[int],
    pop_expr: Optional[hl.expr.StringExpression] = None,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotate MatrixTable or Table with downsampling groups.

    :param t: Input MatrixTable or Table.
    :param downsamplings: List of downsampling sizes.
    :param pop_expr: Optional expression for population group. When provided, population
        sample sizes are added as values to downsamplings.
    :return: MatrixTable or Table with downsampling annotations.
    """
    if isinstance(t, hl.MatrixTable):
        if pop_expr is not None:
            ht = t.annotate_cols(pop=pop_expr).cols()
        else:
            ht = t.cols()
    else:
        if pop_expr is not None:
            ht = t.annotate(pop=pop_expr)
        else:
            ht = t

    ht = ht.annotate(r=hl.rand_unif(0, 1))
    ht = ht.order_by(ht.r)

    # Add a global index for use in computing frequencies, or other aggregate stats on
    # the downsamplings.
    scan_expr = {"global_idx": hl.scan.count()}

    # If pop_expr is provided, add all pop counts to the downsamplings list.
    if pop_expr is not None:
        pop_counts = ht.aggregate(
            hl.agg.filter(hl.is_defined(ht.pop), hl.agg.counter(ht.pop))
        )
        downsamplings = [x for x in downsamplings if x <= sum(pop_counts.values())]
        downsamplings = sorted(set(downsamplings + list(pop_counts.values())))
        # Add an index by pop for use in computing frequencies, or other aggregate stats
        # on the downsamplings.
        scan_expr["pop_idx"] = hl.scan.counter(ht.pop).get(ht.pop, 0)
    else:
        pop_counts = None
    logger.info("Found %i downsamplings: %s", len(downsamplings), downsamplings)

    ht = ht.annotate(**scan_expr)
    ht = ht.key_by("s").select(*scan_expr)

    if isinstance(t, hl.MatrixTable):
        t = t.annotate_cols(downsampling=ht[t.s])
    else:
        t = t.annotate(downsampling=ht[t.s])

    t = t.annotate_globals(
        downsamplings=downsamplings,
        ds_pop_counts=pop_counts,
    )

    return t


def compute_freq_by_strata(
    mt: hl.MatrixTable,
    strata_expr: List[Dict[str, hl.expr.StringExpression]],
    downsamplings: Optional[List[int]] = None,
    ds_pop_counts: Optional[Dict[str, int]] = None,
    entry_agg_funcs: Optional[Dict[str, Tuple[Callable, Callable]]] = None,
) -> hl.Table:
    """
    Compute allele frequencies by strata and downsamplings, when provided.

    .. note::
        This function is primarily used through annotate_freq but can be used
        independently if desired.

    :param mt: Input MatrixTable.
    :param strata_expr: List of dicts of strata expressions.
    :param downsamplings: Optional list of downsampling groups.
    :param ds_pop_counts: Optional dict of population counts for downsampling groups.
    :param entry_agg_funcs: Optional dict of entry aggregation functions.
    :return: Table or MatrixTable with allele frequencies by strata.
    """
    n_samples = mt.count_cols()

    # Get counters for all strata.
    strata_counts = mt.aggregate_cols(
        hl.struct(
            **{
                k: hl.agg.filter(hl.is_defined(v), hl.agg.counter({k: v}))
                for strata in strata_expr
                for k, v in strata.items()
            }
        )
    )

    # Add all desired strata to sample group filters.
    sample_group_filters = [({}, True)]
    for strata in strata_expr:
        downsampling_expr = strata.get("downsampling")
        strata_values = []
        for s in strata:
            if s == "downsampling":
                v = [("downsampling", d) for d in downsamplings]
            else:
                v = [(s, k[s]) for k in strata_counts.get(s, {})]
                if s == "pop" and downsampling_expr is not None:
                    v.append(("pop", "global"))
            strata_values.append(v)

        # Get all combinations of strata values.
        strata_combinations = itertools.product(*strata_values)
        for combo in strata_combinations:
            combo = dict(combo)
            ds = combo.get("downsampling")
            pop = combo.get("pop")
            # If combo contains downsampling, determine the downsampling index
            # annotation to use.
            downsampling_idx = "global_idx"
            if ds is not None:
                if pop is not None and pop != "global":
                    # Don't include population downsamplings where the downsampling is
                    # larger than the number of samples in the population.
                    if ds > ds_pop_counts[pop]:
                        continue
                    downsampling_idx = "pop_idx"

            # If combo contains downsampling, add downsampling filter expression.
            combo_filter_exprs = []
            for s, v in combo.items():
                if s == "downsampling":
                    combo_filter_exprs.append(downsampling_expr[downsampling_idx] < v)
                else:
                    if s != "pop" or v != "global":
                        combo_filter_exprs.append(strata[s] == v)
            combo = {k: str(v) for k, v in combo.items()}
            sample_group_filters.append((combo, hl.all(combo_filter_exprs)))

    n_groups = len(sample_group_filters)
    logger.info("number of filters: %i", n_groups)

    # Annotate columns with group_membership.
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Get sample count per strata group.
    freq_sample_count = mt.aggregate_cols(
        hl.agg.array_agg(lambda x: hl.agg.count_where(x), mt.group_membership)
    )

    # Create and annotate global expression with meta and sample count information
    freq_meta_expr = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]
    # Add the "raw" group, representing all samples, to the freq_meta_expr list.
    freq_meta_expr.insert(1, {"group": "raw"})
    freq_sample_count.insert(1, freq_sample_count[0])
    mt = mt.annotate_globals(
        freq_meta=freq_meta_expr,
        freq_sample_count=freq_sample_count,
    )

    # Create frequency expression array from the sample groups.
    ht = mt.localize_entries("entries", "cols")
    ht = ht.annotate_globals(
        indices_by_group=hl.range(n_groups).map(
            lambda g_i: hl.range(n_samples).filter(
                lambda s_i: ht.cols[s_i].group_membership[g_i]
            )
        )
    )
    ht = ht.annotate(
        adj_array=ht.entries.map(lambda e: e.adj),
        gt_array=ht.entries.map(lambda e: e.GT),
    )

    def _agg_by_group(ht, agg_func, agg_expr, *args):
        adj_agg_expr = ht.indices_by_group.map(
            lambda s_indices: s_indices.aggregate(
                lambda i: hl.agg.filter(ht.adj_array[i], agg_func(agg_expr[i], *args))
            )
        )
        raw_agg_expr = agg_expr.aggregate(lambda x: agg_func(x, *args))
        # Add the "raw" group, representing all samples, to the adj_agg_expr list.
        return adj_agg_expr[:1].append(raw_agg_expr).extend(adj_agg_expr[1:])

    freq_expr = _agg_by_group(ht, hl.agg.call_stats, ht.gt_array, ht.alleles)

    # Select non-ref allele (assumes bi-allelic).
    ann_expr = {
        "freq": freq_expr.map(
            lambda cs: cs.annotate(
                AC=cs.AC[1],
                # TODO: This is NA in case AC and AN are 0 -- should we set it to 0?
                AF=cs.AF[1],
                homozygote_count=cs.homozygote_count[1],
            )
        )
    }

    # Add annotations for any supplied entry transform and aggregation functions.
    if entry_agg_funcs is not None:
        for ann, f in entry_agg_funcs.items():
            transform_func = f[0]
            agg_func = f[1]
            ann_expr[ann] = _agg_by_group(
                ht,
                agg_func,
                hl.map(lambda e, s: transform_func(e, s), ht.entries, ht.cols),
            )

    ht = ht.select(**ann_expr)

    return ht
