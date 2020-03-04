import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import hail as hl

from gnomad.utils.generic import filter_to_autosomes
from gnomad.utils.gnomad_functions import annotate_adj

# TODO: Use import below when relatedness PR goes in
# from gnomad.utils.relatedness import SIBLINGS


logger = logging.getLogger("gnomad.utils")


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
    Adds a row annotation `freq` to the input `mt` with stratified allele frequencies,
    and a global annotation `freq_meta` with metadata.

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

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta information
    freq_meta_expr = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]
    freq_meta_expr.insert(1, {"group": "raw"})
    mt = mt.annotate_globals(freq_meta=freq_meta_expr)

    # Create frequency expression array from the sample groups
    freq_expr = hl.agg.array_agg(
        lambda i: hl.agg.filter(
            mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)
        ),
        hl.range(len(sample_group_filters)),
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

        This lowqual annotation differs from the GATK LowQual filter. GATK uses the SNV prior if the site has any SNV
        so that SNVs are not penalized if there are also indels at those sites. This implementation will use the indel
        threshold for indels and the SNV threshold for the SNVs therefore excluding more indels than the GATK filter.

    :param alleles: Array of alleles
    :param qual_approx_expr: QUALapprox or AS_QUALapprox
    :param snv_phred_threshold: Phred-scaled SNV "emission" threshold (similar to GATK emission threshold)
    :param snv_phred_het_prior: Phred-scaled SNV heterozygosity prior (30 = 1/1000 bases, GATK default)
    :param indel_phred_threshold: Phred-scaled indel "emission" threshold (similar to GATK emission threshold)
    :param indel_phred_het_prior: Phred-scaled indel heterozygosity prior (30 = 1/1000 bases, GATK default)
    :return: lowqual expression (BooleanExpression if `qual_approx_expr`is Numeric, Array[BooleanExpression] if `qual_approx_expr` is ArrayNumeric)
    """

    def low_qual_expr(
        ref: hl.expr.StringExpression,
        alt: hl.expr.StringExpression,
        qual_approx: hl.expr.NumericExpression,
    ) -> hl.expr.BooleanExpression:
        return hl.cond(
            hl.is_snp(ref, alt),
            qual_approx < snv_phred_threshold + snv_phred_het_prior,
            qual_approx < indel_phred_threshold + indel_phred_het_prior,
        )

    if isinstance(qual_approx_expr, hl.expr.ArrayNumericExpression):
        return hl.range(1, hl.len(alleles)).map(
            lambda ai: low_qual_expr(alleles[0], alleles[ai], qual_approx_expr[ai - 1])
        )
    else:
        return low_qual_expr(alleles[0], alleles[1], qual_approx_expr)


def generate_trio_stats_expr(
    trio_mt: hl.MatrixTable,
    transmitted_strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    de_novo_strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    ac_strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    proband_is_female_expr: Optional[hl.expr.BooleanExpression] = None,
) -> hl.expr.StructExpression:
    """
    Generates a row-wise expression containing the following counts:

        - Number of alleles in het parents transmitted to the proband
        - Number of alleles in het parents not transmitted to the proband
        - Number of de novo mutations
        - Parent allele count
        - Proband allele count

    Transmission and de novo mutation metrics and allele counts can be stratified using additional filters.
    `transmitted_strata`, `de_novo_strata`, and `ac_strata` all expect a dictionary of filtering expressions keyed
    by their desired suffix to append for labeling. The default will perform counts using all genotypes and append
    'raw' to the label.

    .. note::

        Expects that `mt` is dense if dealing with a sparse MT `hl.experimental.densify` must be run first.

    :param trio_mt: A trio standard trio MT (with the format as produced by hail.methods.trio_matrix)
    :param transmitted_strata: Strata for the transmission counts
    :param de_novo_strata: Strata for the de novo counts
    :param ac_strata: Strata for the parent and child allele counts
    :param proband_is_female_expr: An optional expression giving the sex the proband. If not given, DNMs are only computed for autosomes.
    :return: An expression with the counts
    """

    # Create map for transmitted, untransmitted and DNM
    hom_ref = 0
    het = 1
    hom_var = 2

    auto_or_par = 2
    hemi_x = 1
    hemi_y = 0

    trans_config_counts = {
        # kid, dad, mom, copy -> t, u
        (hom_ref, het, het, auto_or_par): (0, 2),
        (hom_ref, hom_ref, het, auto_or_par): (0, 1),
        (hom_ref, het, hom_ref, auto_or_par): (0, 1),
        (het, het, het, auto_or_par): (1, 1),
        (het, hom_ref, het, auto_or_par): (1, 0),
        (het, het, hom_ref, auto_or_par): (1, 0),
        (het, hom_var, het, auto_or_par): (0, 1),
        (het, het, hom_var, auto_or_par): (0, 1),
        (hom_var, het, het, auto_or_par): (2, 0),
        (hom_var, het, hom_var, auto_or_par): (1, 0),
        (hom_var, hom_var, het, auto_or_par): (1, 0),
        (hom_ref, hom_ref, het, hemi_x): (0, 1),
        (hom_ref, hom_var, het, hemi_x): (0, 1),
        (hom_var, hom_ref, het, hemi_x): (1, 0),
        (hom_var, hom_var, het, hemi_x): (1, 0),
    }

    trans_count_map = hl.literal(trans_config_counts)

    def _get_copy_state(locus: hl.expr.LocusExpression) -> hl.expr.Int32Expression:
        """
        Helper method to go from LocusExpression to a copy-state int for indexing into the
        trans_count_map.
        """
        return (
            hl.case()
            .when(locus.in_autosome_or_par(), auto_or_par)
            .when(locus.in_x_nonpar(), hemi_x)
            .when(locus.in_y_nonpar(), hemi_y)
            .or_missing()
        )

    def _is_dnm(
        proband_gt: hl.expr.CallExpression,
        father_gt: hl.expr.CallExpression,
        mother_gt: hl.expr.CallExpression,
        locus: hl.expr.LocusExpression,
        proband_is_female: Optional[hl.expr.BooleanExpression],
    ) -> hl.expr.BooleanExpression:
        """
        Helper method to get whether a given genotype combination is a DNM at a given locus with a given proband sex.
        """
        if proband_is_female is None:
            logger.warning(
                "Since no proband sex expression was given to generate_trio_stats_expr, only DNMs in autosomes will be counted."
            )
            return hl.or_missing(
                locus.in_autosome(),
                proband_gt.is_het() & father_gt.is_hom_ref() & mother_gt.is_hom_ref(),
            )
        return hl.cond(
            locus.in_autosome_or_par() | (proband_is_female & locus.in_x_nonpar()),
            proband_gt.is_het() & father_gt.is_hom_ref() & mother_gt.is_hom_ref(),
            hl.or_missing(
                ~proband_is_female, proband_gt.is_hom_var() & father_gt.is_hom_ref()
            ),
        )

    def _ac_an_parent_child_count(
        proband_gt: hl.expr.CallExpression,
        father_gt: hl.expr.CallExpression,
        mother_gt: hl.expr.CallExpression,
    ) -> Dict[str, hl.expr.Int64Expression]:
        """
        Helper method to get AC and AN for parents and children
        """
        ac_parent_expr = hl.agg.sum(
            father_gt.n_alt_alleles() + mother_gt.n_alt_alleles()
        )
        an_parent_expr = hl.agg.sum(
            (hl.is_defined(father_gt) + hl.is_defined(mother_gt)) * 2
        )
        ac_child_expr = hl.agg.sum(proband_gt.n_alt_alleles())
        an_child_expr = hl.agg.sum(hl.is_defined(proband_gt) * 2)

        return {
            f"ac_parents": ac_parent_expr,
            f"an_parents": an_parent_expr,
            f"ac_children": ac_child_expr,
            f"an_children": an_child_expr,
        }

    # Create transmission counters
    trio_stats = hl.struct(
        **{
            f"{name2}_{name}": hl.agg.filter(
                trio_mt.proband_entry.GT.is_non_ref() & expr,
                hl.agg.sum(
                    trans_count_map.get(
                        (
                            trio_mt.proband_entry.GT.n_alt_alleles(),
                            trio_mt.father_entry.GT.n_alt_alleles(),
                            trio_mt.mother_entry.GT.n_alt_alleles(),
                            _get_copy_state(trio_mt.locus),
                        ),
                        default=0,
                    )[i]
                ),
            )
            for name, expr in transmitted_strata.items()
            for i, name2 in enumerate(["n_transmitted", "n_untransmitted"])
        }
    )

    # Create de novo counters
    trio_stats = trio_stats.annotate(
        **{
            f"n_de_novos_{name}": hl.agg.filter(
                _is_dnm(
                    trio_mt.proband_entry.GT,
                    trio_mt.father_entry.GT,
                    trio_mt.mother_entry.GT,
                    trio_mt.locus,
                    proband_is_female_expr,
                )
                & expr,
                hl.agg.count(),
            )
            for name, expr in de_novo_strata.items()
        }
    )

    trio_stats = trio_stats.annotate(
        **{
            f"{name2}_{name}": hl.agg.filter(
                expr,
                _ac_an_parent_child_count(
                    trio_mt.proband_entry.GT,
                    trio_mt.father_entry.GT,
                    trio_mt.mother_entry.GT,
                )[name2],
            )
            for name, expr in ac_strata.items()
            for name2 in ["ac_parents", "an_parents", "ac_children", "an_children"]
        }
    )

    return trio_stats


def filter_mt_to_trios(mt: hl.MatrixTable, fam_ht: hl.Table) -> hl.MatrixTable:
    """
    Filters a MatrixTable to a set of trios in `fam_ht`, filters to autosomes, and annotates with adj.

    :param mt: A Matrix Table to filter to only trios
    :param fam_ht: A Table of trios to filter to, loaded using `hl.import_fam`
    :return: A MT filtered to trios and adj annotated
    """
    # Filter MT to samples present in any of the trios
    fam_ht = fam_ht.annotate(fam_members=[fam_ht.id, fam_ht.pat_id, fam_ht.mat_id])
    fam_ht = fam_ht.explode("fam_members", name="s")
    fam_ht = fam_ht.key_by("s").select().distinct()

    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.col_key]))
    mt = filter_to_autosomes(mt)
    mt = annotate_adj(mt)

    return mt


def default_generate_trio_stats(mt: hl.MatrixTable,) -> hl.Table:
    """
    Default function to run `generate_trio_stats_expr` to get trio stats stratified by raw and adj

    .. note::

        Expects that `mt` is it a trio matrix table that was annotated with adj and if dealing with
        a sparse MT `hl.experimental.densify` must be run first.

    :param mt: A Trio Matrix Table returned from `hl.trio_matrix`. Must be dense
    :return: Table with trio stats
    """
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    logger.info(f"Generating trio stats using {mt.count_cols()} trios.")
    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={"raw": True, "adj": trio_adj},
            de_novo_strata={"raw": True, "adj": trio_adj},
            ac_strata={"raw": True, "adj": trio_adj},
            proband_is_female_expr=mt.is_female,
        )
    ).rows()

    return ht


def generate_sib_stats_expr(
    mt: hl.MatrixTable,
    sib_ht: hl.Table,
    i_col: str = "i",
    j_col: str = "j",
    strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    is_female: Optional[hl.expr.BooleanExpression] = None,
) -> hl.expr.StructExpression:
    """
    Generates a row-wise expression containing the number of alternate alleles in common between sibling pairs.

    The sibling sharing counts can be stratified using additional filters using `stata`.

    .. note::

        This function expects that the `mt` has either been split or filtered to only bi-allelics
        If a sample has multiple sibling pairs, only one pair will be counted

    :param mt: Input matrix table
    :param sib_ht: Table defining sibling pairs with one sample in a col (`i_col`) and the second in another col (`j_col`)
    :param i_col: Column containing the 1st sample of the pair in the relationship table
    :param j_col: Column containing the 2nd sample of the pair in the relationship table
    :param strata: Dict with additional strata to use when computing shared sibling variant counts
    :param is_female: An optional column in mt giving the sample sex. If not given, counts are only computed for autosomes.
    :return: A Table with the sibling shared variant counts
    """

    def get_alt_count(locus, gt, is_female):
        """
        Helper method to calculate alt allele count with sex info if present
        """
        if is_female is None:
            return hl.or_missing(locus.in_autosome(), gt.n_alt_alleles())
        return (
            hl.case()
            .when(locus.in_autosome_or_par(), gt.n_alt_alleles())
            .when(
                ~is_female & (locus.in_x_nonpar() | locus.in_y_nonpar()),
                hl.min(1, gt.n_alt_alleles()),
            )
            .when(is_female & locus.in_y_nonpar(), 0)
            .default(0)
        )

    if is_female is None:
        logger.warning(
            "Since no sex expression was given to generate_sib_stats_expr, only variants in autosomes will be counted."
        )

    # If a sample is in sib_ht more than one time, keep only one of the sibling pairs
    # First filter to only samples found in mt to keep as many pairs as possible
    s_to_keep = mt.aggregate_cols(hl.agg.collect_as_set(mt.s), _localize=False)
    sib_ht = sib_ht.filter(
        s_to_keep.contains(sib_ht[i_col].s) & s_to_keep.contains(sib_ht[j_col].s)
    )
    sib_ht = sib_ht.add_index("sib_idx")
    sib_ht = sib_ht.annotate(sibs=[sib_ht[i_col].s, sib_ht[j_col].s])
    sib_ht = sib_ht.explode("sibs")
    sib_ht = sib_ht.group_by("sibs").aggregate(
        sib_idx=(hl.agg.take(sib_ht.sib_idx, 1, ordering=sib_ht.sib_idx)[0])
    )
    sib_ht = sib_ht.group_by(sib_ht.sib_idx).aggregate(sibs=hl.agg.collect(sib_ht.sibs))
    sib_ht = sib_ht.filter(hl.len(sib_ht.sibs) == 2).persist()

    logger.info(
        f"Generating sibling variant sharing counts using {sib_ht.count()} pairs."
    )
    sib_ht = sib_ht.explode("sibs").key_by("sibs")[mt.s]

    # Create sibling sharing counters
    sib_stats = hl.struct(
        **{
            f"n_sib_shared_variants_{name}": hl.sum(
                hl.agg.filter(
                    expr,
                    hl.agg.group_by(
                        sib_ht.sib_idx,
                        hl.or_missing(
                            hl.agg.sum(hl.is_defined(mt.GT)) == 2,
                            hl.agg.min(get_alt_count(mt.locus, mt.GT, is_female)),
                        ),
                    ),
                ).values()
            )
            for name, expr in strata.items()
        }
    )

    sib_stats = sib_stats.annotate(
        **{
            f"ac_sibs_{name}": hl.agg.filter(
                expr & hl.is_defined(sib_ht.sib_idx), hl.agg.sum(mt.GT.n_alt_alleles())
            )
            for name, expr in strata.items()
        }
    )

    return sib_stats


def default_generate_sib_stats(
    mt: hl.MatrixTable,
    relatedness_ht: hl.Table,
    sex_ht: hl.Table,
    i_col: str = "i",
    j_col: str = "j",
    relationship_col: str = "relationship",
) -> hl.Table:
    """
    This is meant as a default wrapper for `generate_sib_stats_expr`. It returns a hail table with counts of variants
    shared by pairs of siblings in `relatedness_ht`.

    This function takes a hail Table with a row for each pair of individuals i,j in the data that are related (it's OK to have unrelated samples too).
    The `relationship_col` should be a column specifying the relationship between each two samples as defined by
    the constants in `gnomad.utils.relatedness`. This relationship_col will be used to filter to only pairs of
    samples that are annotated as `SIBLINGS`.

    :param mt: Input Matrix table
    :param relatedness_ht: Input relationship table
    :param sex_ht: A Table containing sex information for the samples
    :param i_col: Column containing the 1st sample of the pair in the relationship table
    :param j_col: Column containing the 2nd sample of the pair in the relationship table
    :param relationship_col: Column containing the relationship for the sample pair as defined in this module constants.
    :return: A Table with the sibling shared variant counts
    """
    sex_ht = sex_ht.annotate(
        is_female=hl.case()
        .when(sex_ht.sex_karyotype == "XX", True)
        .when(sex_ht.sex_karyotype == "XY", False)
        .or_missing()
    )

    # TODO: Change to use SIBLINGS constant when relatedness PR goes in
    sib_ht = relatedness_ht.filter(relatedness_ht[relationship_col] == "Siblings")
    s_to_keep = sib_ht.aggregate(
        hl.agg.explode(
            lambda s: hl.agg.collect_as_set(s), [sib_ht[i_col].s, sib_ht[j_col].s]
        ),
        _localize=False,
    )
    mt = mt.filter_cols(s_to_keep.contains(mt.s))
    mt = annotate_adj(mt)

    mt = mt.annotate_cols(is_female=sex_ht[mt.s].is_female)

    sib_stats_ht = mt.select_rows(
        **generate_sib_stats_expr(
            mt,
            sib_ht,
            i_col=i_col,
            j_col=j_col,
            strata={"raw": True, "adj": mt.adj},
            is_female=mt.is_female,
        )
    ).rows()

    return sib_stats_ht
  
 
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
            hl.log10(ht[field]) if field in log10_annotations else ht[field],
            start,
            end,
            bin,
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
    bin_expr=(
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
