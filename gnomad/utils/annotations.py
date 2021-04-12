import logging
from typing import Any, Dict, List, Optional, Set, Tuple, Union

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
    Creates an expression containing popmax: the frequency information about the population
    that has the highest AF from the populations provided in `freq_meta`,
    excluding those specified in `pops_to_exclude`.
    Only frequencies from adj populations are considered.

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
    _pops_to_exclude = hl.literal(pops_to_exclude)

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
    Creates an expression that computes allele frequency information by project for the `n_projects` with the largest AF at this row.
    This return an array with one element per non-reference allele.

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
                lambda x: -x[1].AF[ai]
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
    Calculates the filtering allele frequency (FAF) for each threshold specified in `faf_thresholds`.
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
        hl.literal(pops_to_exclude) if pops_to_exclude is not None else {}
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
                        f"faf{str(threshold)[2:]}": hl.experimental.filtering_allele_frequency(
                            freq[i].AC, freq[i].AN, threshold
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
) -> hl.expr.StructExpression:
    """
    Returns a struct expression with genotype quality histograms based on the arguments given (dp, gq, ad).

    .. note::

        - If `gt_expr` is provided, will return histograms for non-reference samples only as well as all samples.
        - `gt_expr` is required for the allele-balance histogram, as it is only computed on het samples.
        - If `adj_expr` is provided, additional histograms are computed using only adj samples.

    :param gt_expr: Entry expression containing genotype
    :param gq_expr: Entry expression containing genotype quality
    :param dp_expr: Entry expression containing depth
    :param ad_expr: Entry expression containing allelic depth (bi-allelic here)
    :param adj_expr: Entry expression containing adj (high quality) genotype status
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
        if ad_expr is not None:
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
    Returns a StructExpression with the age histograms for hets and homs.

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
    additional_strata_expr: Optional[Dict[str, hl.expr.StringExpression]] = None,
    downsamplings: Optional[List[int]] = None,
) -> hl.MatrixTable:
    """
    Adds a row annotation `freq` to the input `mt` with stratified allele frequencies, a global annotation `freq_meta`
    with metadata, and a global annotation `freq_sample_count` with sample count information.

    .. note::

        Currently this only supports bi-allelic sites.
        The input `mt` needs to have the following entry fields:
        - GT: a CallExpression containing the genotype
        - adj: a BooleanExpression containing whether the genotype is of high quality or not.
        All expressions arguments need to be expression on the input `mt`.

    .. rubric:: `freq` row annotation

    The `freq` row annotation is an Array of Struct, with each Struct containing the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data, and the metadata about these annotations is
    stored in the globals.

    .. rubric:: Global `freq_meta` metadata annotation

    The global annotation `freq_meta` is added to the input `mt`. It is a list of dict.
    Each element of the list contains metadata on a frequency stratification and the index in the list corresponds
    to the index of that frequency stratification in the `freq` row annotation.

    .. rubric:: Global `freq_sample_count` annotation

    The global annotation `freq_sample_count` is added to the input `mt`. This is a sample count per sample grouping
    defined in the `freq_meta` global annotation.

    .. rubric:: The `downsamplings` parameter

    If the `downsamplings` parameter is used, frequencies will be computed for all samples and by population
    (if `pop_expr` is specified) by downsampling the number of samples without replacement to each of the numbers specified in the
    `downsamplings` array, provided that there are enough samples in the dataset.
    In addition, if `pop_expr` is specified, a downsampling to each of the exact number of samples present in each population is added.
    Note that samples are randomly sampled only once, meaning that the lower downsamplings are subsets of the higher ones.

    :param mt: Input MatrixTable
    :param sex_expr: When specified, frequencies are stratified by sex. If `pop_expr` is also specified, then a pop/sex stratifiction is added.
    :param pop_expr: When specified, frequencies are stratified by population. If `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param subpop_expr: When specified, frequencies are stratified by sub-continental population. Note that `pop_expr` is required as well when using this option.
    :param additional_strata_expr: When specified, frequencies are stratified by the given additional strata found in the dict. This can e.g. be used to stratify by platform.
    :param downsamplings: When specified, frequencies are computed by downsampling the data to the number of samples given in the list. Note that if `pop_expr` is specified, downsamplings by population is also computed.
    :return: MatrixTable with `freq` annotation
    """

    if subpop_expr is not None and pop_expr is None:
        raise NotImplementedError(
            "annotate_freq requires pop_expr when using subpop_expr"
        )

    if additional_strata_expr is None:
        additional_strata_expr = {}

    _freq_meta_expr = hl.struct(**additional_strata_expr)
    if sex_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(sex=sex_expr)
    if pop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(pop=pop_expr)
    if subpop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(subpop=subpop_expr)

    # Annotate cols with provided cuts
    mt = mt.annotate_cols(_freq_meta=_freq_meta_expr)

    # Get counters for sex, pop and subpop if set
    cut_dict = {
        cut: hl.agg.filter(
            hl.is_defined(mt._freq_meta[cut]), hl.agg.counter(mt._freq_meta[cut])
        )
        for cut in mt._freq_meta
        if cut != "subpop"
    }
    if "subpop" in mt._freq_meta:
        cut_dict["subpop"] = hl.agg.filter(
            hl.is_defined(mt._freq_meta.pop) & hl.is_defined(mt._freq_meta.subpop),
            hl.agg.counter(
                hl.struct(subpop=mt._freq_meta.subpop, pop=mt._freq_meta.pop)
            ),
        )

    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    sample_group_filters = []

    # Create downsamplings if needed
    if downsamplings is not None:
        # Add exact pop size downsampling if pops were provided
        if cut_data.get("pop"):
            downsamplings = list(
                set(downsamplings + list(cut_data.get("pop").values()))
            )  # Add the pops values if not in yet
            downsamplings = sorted(
                [x for x in downsamplings if x <= sum(cut_data.get("pop").values())]
            )
        logger.info(f"Found {len(downsamplings)} downsamplings: {downsamplings}")

        # Shuffle the samples, then create a global index for downsampling
        # And a pop-index if pops were provided
        downsampling_ht = mt.cols()
        downsampling_ht = downsampling_ht.annotate(r=hl.rand_unif(0, 1))
        downsampling_ht = downsampling_ht.order_by(downsampling_ht.r)
        scan_expr = {"global_idx": hl.scan.count()}
        if cut_data.get("pop"):
            scan_expr["pop_idx"] = hl.scan.counter(downsampling_ht._freq_meta.pop).get(
                downsampling_ht._freq_meta.pop, 0
            )
        downsampling_ht = downsampling_ht.annotate(**scan_expr)
        downsampling_ht = downsampling_ht.key_by("s").select(*scan_expr)
        mt = mt.annotate_cols(downsampling=downsampling_ht[mt.s])
        mt = mt.annotate_globals(downsamplings=downsamplings)

        # Create downsampled sample groups
        sample_group_filters.extend(
            [
                (
                    {"downsampling": str(ds), "pop": "global"},
                    mt.downsampling.global_idx < ds,
                )
                for ds in downsamplings
            ]
        )
        if cut_data.get("pop"):
            sample_group_filters.extend(
                [
                    (
                        {"downsampling": str(ds), "pop": pop},
                        (mt.downsampling.pop_idx < ds) & (mt._freq_meta.pop == pop),
                    )
                    for ds in downsamplings
                    for pop, pop_count in cut_data.get("pop", {}).items()
                    if ds <= pop_count
                ]
            )

    # Add all desired strata, starting with the full set and ending with downsamplings (if any)
    sample_group_filters = (
        [({}, True)]
        + [({"pop": pop}, mt._freq_meta.pop == pop) for pop in cut_data.get("pop", {})]
        + [({"sex": sex}, mt._freq_meta.sex == sex) for sex in cut_data.get("sex", {})]
        + [
            (
                {"pop": pop, "sex": sex},
                (mt._freq_meta.sex == sex) & (mt._freq_meta.pop == pop),
            )
            for sex in cut_data.get("sex", {})
            for pop in cut_data.get("pop", {})
        ]
        + [
            (
                {"subpop": subpop.subpop, "pop": subpop.pop},
                (mt._freq_meta.pop == subpop.pop)
                & (mt._freq_meta.subpop == subpop.subpop),
            )
            for subpop in cut_data.get("subpop", {})
        ]
        + [
            ({strata: str(s_value)}, mt._freq_meta[strata] == s_value)
            for strata in additional_strata_expr
            for s_value in cut_data.get(strata, {})
        ]
        + sample_group_filters
    )

    freq_sample_count = mt.aggregate_cols(
        [hl.agg.count_where(x[1]) for x in sample_group_filters]
    )

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta and sample count information
    freq_meta_expr = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]
    freq_meta_expr.insert(1, {"group": "raw"})
    freq_sample_count.insert(1, freq_sample_count[0])
    mt = mt.annotate_globals(
        freq_meta=freq_meta_expr, freq_sample_count=freq_sample_count,
    )

    # Create frequency expression array from the sample groups
    # Adding sample_group_filters_range_array to reduce memory usage in this array_agg
    mt = mt.annotate_rows(
        sample_group_filters_range_array=hl.range(len(sample_group_filters))
    )
    freq_expr = hl.agg.array_agg(
        lambda i: hl.agg.filter(
            mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)
        ),
        mt.sample_group_filters_range_array,
    )

    # Insert raw as the second element of the array
    freq_expr = (
        freq_expr[:1]
        .extend([hl.agg.call_stats(mt.GT, mt.alleles)])
        .extend(freq_expr[1:])
    )

    # Select non-ref allele (assumes bi-allelic)
    freq_expr = freq_expr.map(
        lambda cs: cs.annotate(
            AC=cs.AC[1],
            AF=cs.AF[
                1
            ],  # TODO This is NA in case AC and AN are 0 -- should we set it to 0?
            homozygote_count=cs.homozygote_count[1],
        )
    )

    # Return MT with freq row annotation
    return mt.annotate_rows(freq=freq_expr).drop("_freq_meta")


def get_lowqual_expr(
    alleles: hl.expr.ArrayExpression,
    qual_approx_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
    snv_phred_threshold: int = 30,
    snv_phred_het_prior: int = 30,  # 1/1000
    indel_phred_threshold: int = 30,
    indel_phred_het_prior: int = 39,  # 1/8,000
) -> Union[hl.expr.BooleanExpression, hl.expr.ArrayExpression]:
    """
    Computes lowqual threshold expression for either split or unsplit alleles based on QUALapprox or AS_QUALapprox

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
    Creates histograms for variant metrics in ht.info.
    Used when creating site quality distribution json files.

    :param ht: Table with variant metrics
    :param annotations_hists: Dictionary of metrics names and their histogram values (start, end, bins)
    :param log10_annotations: List of metrics to log scale
    :return: Dictionary of merics and their histograms
    :rtype: Dict[str, hl.expr.StructExpression]
    """
    # Check all fields in ht.info and create histograms if they are in annotations_hists dict
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
    Creates bins for frequencies in preparation for aggregating QUAL by frequency bin.

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
    Gets adj genotype annotation.
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
    Annotate genotypes with adj criteria (assumes diploid)
    Defaults correspond to gnomAD values.
    """
    return mt.annotate_entries(
        adj=get_adj_expr(
            mt.GT, mt.GQ, mt.DP, mt.AD, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    )


def add_variant_type(alt_alleles: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """
    Get Struct of variant_type and n_alt_alleles from ArrayExpression of Strings (all alleles)
    """
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


def annotation_type_is_numeric(t: Any) -> bool:
    """
    Given an annotation type, returns whether it is a numerical type or not.

    :param t: Type to test
    :return: If the input type is numeric
    """
    return (
        isinstance(t, hl.tint32)
        or isinstance(t, hl.tint64)
        or isinstance(t, hl.tfloat32)
        or isinstance(t, hl.tfloat64)
    )


def annotation_type_in_vcf_info(t: Any) -> bool:
    """
    Given an annotation type, returns whether that type can be natively exported to a VCF INFO field.
    Note types that aren't natively exportable to VCF will be converted to String on export.

    :param t: Type to test
    :return: If the input type can be exported to VCF
    """
    return (
        annotation_type_is_numeric(t)
        or isinstance(t, hl.tstr)
        or isinstance(t, hl.tarray)
        or isinstance(t, hl.tset)
        or isinstance(t, hl.tbool)
    )


def bi_allelic_site_inbreeding_expr(
    call: hl.expr.CallExpression,
) -> hl.expr.Float32Expression:
    """
    Return the site inbreeding coefficient as an expression to be computed on a MatrixTable.

    This is implemented based on the GATK InbreedingCoeff metric:
    https://software.broadinstitute.org/gatk/documentation/article.php?id=8032

    .. note::

        The computation is run based on the counts of alternate alleles and thus should only be run on bi-allelic sites.

    :param call: Expression giving the calls in the MT
    :return: Site inbreeding coefficient expression
    """

    def inbreeding_coeff(
        gt_counts: hl.expr.DictExpression,
    ) -> hl.expr.Float32Expression:
        n = gt_counts.get(0, 0) + gt_counts.get(1, 0) + gt_counts.get(2, 0)
        p = (2 * gt_counts.get(0, 0) + gt_counts.get(1, 0)) / (2 * n)
        q = (2 * gt_counts.get(2, 0) + gt_counts.get(1, 0)) / (2 * n)
        return 1 - (gt_counts.get(1, 0) / (2 * p * q * n))

    return hl.bind(inbreeding_coeff, hl.agg.counter(call.n_alt_alleles()))


def fs_from_sb(
    sb: Union[hl.expr.ArrayNumericExpression, hl.expr.ArrayExpression],
    normalize: bool = True,
    min_cell_count: int = 200,
    min_count: int = 4,
    min_p_value: float = 1e-320,
) -> hl.expr.Int64Expression:
    """
    Computes `FS` (Fisher strand balance) annotation from  the `SB` (strand balance table) field.
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
    Computes `SOR` (Symmetric Odds Ratio test) annotation from  the `SB` (strand balance table) field.

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


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression selecting bi-allelic sites only,
    accounting for whether the input MT/HT was split.

    :param t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    """
    return ~t.was_split if "was_split" in t.row else (hl.len(t.alleles) == 2)


def unphase_call_expr(call_expr: hl.expr.CallExpression) -> hl.expr.CallExpression:
    """
    Generate unphased version of a call expression (which can be phased or not)

    :param call_expr: Input call expression
    :return: unphased call expression
    """
    return (
        hl.case()
        .when(call_expr.is_diploid(), hl.call(call_expr[0], call_expr[1], phased=False))
        .when(call_expr.is_haploid(), hl.call(call_expr[0], phased=False))
        .default(hl.null(hl.tcall))
    )


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
