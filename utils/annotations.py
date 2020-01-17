from gnomad_hail import *


def pop_max_expr(
        freq: hl.expr.ArrayExpression,
        freq_meta: hl.expr.ArrayExpression,
        pops_to_exclude: Optional[Set[str]] = None
) -> hl.expr.StructExpression:
    """
    Creates an expression containing popmax: the frequency information about the population
    that has the highest AF from the populations provided in `freq_meta`,
    excluding those specified in `pops_to_exclude`.
    Only frequencies from adj populations are considered.
    This resulting struct contains the following fields:
    - AC: int32,
    - AF: float64,
    - AN: int32,
    - homozygote_count: int32,
    - pop: str

    :param ArrayExpression freq: ArrayExpression of Structs with fields ['AC', 'AF', 'AN', 'homozygote_count']
    :param ArrayExpression freq_meta: ArrayExpression of meta dictionaries corresponding to freq (as returned by annotate_freq)
    :param set of str pops_to_exclude: Set of populations to skip for popmax calcluation
    :return: Popmax struct
    :rtype: StructExpression
    """
    _pops_to_exclude = hl.literal(pops_to_exclude)
    popmax_freq_indices = hl.range(0, hl.len(freq_meta)).filter(
        lambda i:
        (hl.set(freq_meta[i].keys()) == {'group', 'pop'}) &
        (freq_meta[i]['group'] == 'adj') &
        (~_pops_to_exclude.contains(freq_meta[i]['pop']))
    )
    freq_filtered = popmax_freq_indices.map(
        lambda i: freq[i].annotate(pop=freq_meta[i]['pop'])
    ).filter(
        lambda f: f.AC > 0
    )

    sorted_freqs = hl.sorted(freq_filtered, key=lambda x: x.AF, reverse=True)
    return hl.or_missing(
        hl.len(sorted_freqs) > 0,
        sorted_freqs[0]
    )


def project_max_expr(
        project_expr: hl.expr.StringExpression,
        gt_expr: hl.expr.CallExpression,
        alleles_expr: hl.expr.ArrayExpression,
        n_projects: int = 5
) -> hl.expr.ArrayExpression:
    """
    Creates an expression that computes allele frequency information by project for the `n_projects` with the largest AF at this row.
    This return an array with one element per non-reference allele.
    Each of these elements is itself an array of structs with the following fields:
    - AC: int32,
    - AF: float64,
    - AN: int32,
    - homozygote_count: int32,
    - project: str

    .. note::

        Only projects with AF > 0 are returned.
        In case of ties, the project ordering is not guaranteed, and at most `n_projects` are returned.

    :param StringExpression project_expr: column expression containing the project
    :param CallExpression gt_expr: entry expression containing the genotype
    :param ArrayExpression alleles_expr: row expression containing the alleles
    :param int n_projects: Maximum number of projects to return for each row
    :return: projectmax expression
    :rtype: ArrayExpression
    """

    n_alleles = hl.len(alleles_expr)

    # compute call stats by  project
    project_cs = hl.array(hl.agg.group_by(project_expr, hl.agg.call_stats(gt_expr, alleles_expr)))

    return hl.or_missing(
        n_alleles > 1, # Exclude monomorphic sites
        hl.range(1, n_alleles).map(
            lambda ai: hl.sorted(
                project_cs.filter(
                    # filter to projects with AF > 0
                    lambda x: x[1].AF[ai] > 0
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
                    project=x[0]
                )
            )
        )
    )


def faf_expr(
        freq: hl.expr.ArrayExpression,
        freq_meta: hl.expr.ArrayExpression,
        locus: hl.expr.LocusExpression,
        pops_to_exclude: Optional[Set[str]] = None,
        faf_thresholds: List[float] = [0.95, 0.99]

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

    :param ArrayExpression freq: ArrayExpression of call stats structs (typically generated by hl.agg.call_stats)
    :param ArrayExpression freq_meta: ArrayExpression of meta dictionaries corresponding to freq (typically generated using annotate_freq)
    :param LocusExpression locus: locus
    :param set of str pops_to_exclude: Set of populations to exclude from faf calculation (typically bottlenecked or consanguineous populations)
    :param list of float faf_thresholds: List of FAF thresholds to compute
    :return: (FAF expression, FAF metadata)
    :rtype: Tuple(ArrayExpression, List[Dict[str, str]])
    """
    _pops_to_exclude = hl.literal(pops_to_exclude) if pops_to_exclude is not None else {}
    faf_freq_indices = hl.range(0, hl.len(freq_meta)).filter(
        lambda i:
        (freq_meta[i].get('group') == 'adj') &
        (
            (freq_meta[i].size() == 1) |
            ((hl.set(freq_meta[i].keys()) == {'pop', 'group'}) & (~_pops_to_exclude.contains(freq_meta[i]['pop'])))
        )
    )
    sex_faf_freq_indices = hl.range(0, hl.len(freq_meta)).filter(
        lambda i:
        (freq_meta[i].get('group') == 'adj') &
        (freq_meta[i].contains('sex')) &
        (
            (freq_meta[i].size() == 2) |
            ((hl.set(freq_meta[i].keys()) == {'pop', 'group', 'sex'}) & (~_pops_to_exclude.contains(freq_meta[i]['pop'])))
        )
    )

    faf_expr = faf_freq_indices.map(
        lambda i: hl.struct(**{
            f'faf{str(threshold)[2:]}': hl.experimental.filtering_allele_frequency(freq[i].AC, freq[i].AN, threshold)
            for threshold in faf_thresholds
        })
    )

    faf_expr = faf_expr.extend(
        sex_faf_freq_indices.map(
            lambda i: hl.or_missing(
                ~locus.in_autosome_or_par(),
                hl.struct(**{
                    f'faf{str(threshold)[2:]}': hl.experimental.filtering_allele_frequency(freq[i].AC, freq[i].AN, threshold)
                    for threshold in faf_thresholds
                })
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
        adj_expr: Optional[hl.expr.BooleanExpression] = None
) -> hl.expr.StructExpression:
    """
    Returns a struct expression with genotype quality histograms based on the arguments given (dp, gq, ad).

    .. note::

        - If `gt_expr` is provided, will return histograms for non-reference samples only as well as all samples.
        - `gt_expr` is required for the allele-balance histogram, as it is only computed on het samples.
        - If `adj_expr` is provided, additional histograms are computed using only adj samples.

    :param CallExpression gt_expr: Entry expression containing genotype
    :param NumericExpression gq_expr: Entry expression containing genotype quality
    :param NumericExpression dp_expr: Entry expression containing depth
    :param ArrayNumericExpression ad_expr: Entry expression containing allelic depth (bi-allelic here)
    :param BooleanExpression adj_expr: Entry expression containing adj (high quality) genotype status
    :return: Genotype quality histograms expression
    :rtype: StructExpression
    """
    qual_hists = {}
    if gq_expr is not None:
        qual_hists['gq_hist'] = hl.agg.hist(gq_expr, 0, 100, 20)
    if dp_expr is not None:
        qual_hists['dp_hist'] = hl.agg.hist(dp_expr, 0, 100, 20)

    if gt_expr is not None:
        qual_hists= {
            **{f'{qual_hist_name}_all': qual_hist_expr for qual_hist_name, qual_hist_expr in qual_hists.items()},
            **{f'{qual_hist_name}_alt': hl.agg.filter(gt_expr.is_non_ref(), qual_hist_expr) for qual_hist_name, qual_hist_expr in qual_hists.items()}
        }
        if ad_expr is not None:
            qual_hists['ab_hist_alt'] = hl.agg.filter(gt_expr.is_het(), hl.agg.hist(ad_expr[1] / hl.sum(ad_expr), 0, 1, 20))

    else:
        qual_hists = {f'{qual_hist_name}_all': qual_hist_expr for qual_hist_name, qual_hist_expr in qual_hists.items()}

    if adj_expr is not None:
        qual_hists.update({
            f'{qual_hist_name}_adj': hl.agg.filter(adj_expr, qual_hist_expr) for qual_hist_name, qual_hist_expr in qual_hists.items()
        })

    return hl.struct(**qual_hists)


def age_hists_expr(
        adj_expr: hl.expr.BooleanExpression,
        gt_expr: hl.expr.CallExpression,
        age_expr: hl.expr.NumericExpression,
        lowest_boundary: int = 30,
        highest_boundary: int = 80,
        n_bins: int = 10
) -> hl.expr.StructExpression:
    """
    Returns a StructExpression with the age histograms for hets and homs.

    :param BooleanExpression adj_expr: Entry expression containing whether a genotype is high quality (adj) or not
    :param CallExpression gt_expr: Entry expression containing the genotype
    :param NumericExpression age_expr: Col expression containing the sample's age
    :param int lowest_boundary: Lowest bin boundary (any younger sample will be binned in n_smaller)
    :param int highest_boundary: Highest bin boundary (any older sample will be binned in n_larger)
    :param int n_bins: Total number of bins
    :return: A struct with `age_hist_het` and `age_hist_hom`
    :rtype: StructExpression.
    """
    return hl.struct(
        age_hist_het=hl.agg.filter(adj_expr & gt_expr.is_het(), hl.agg.hist(age_expr, lowest_boundary, highest_boundary, n_bins)),
        age_hist_hom=hl.agg.filter(adj_expr & gt_expr.is_hom_var(), hl.agg.hist(age_expr, lowest_boundary, highest_boundary, n_bins))
    )


def annotate_freq(
        mt: hl.MatrixTable,
        sex_expr: Optional[hl.expr.StringExpression] = None,
        pop_expr: Optional[hl.expr.StringExpression] = None,
        subpop_expr: Optional[hl.expr.StringExpression] = None,
        additional_strata_expr: Optional[Dict[str, hl.expr.StringExpression]] = None,
        downsamplings: Optional[List[int]] = None
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
    - AC: int32,
    - AF: float64,
    - AN: int32,
    - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data,
    and the metadata about these annotations is stored in the globals.

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

    :param MatrixTable mt: Input MatrixTable
    :param StringExpression sex_expr: When specified, frequencies are stratified by sex. If `pop_expr` is also specified, then a pop/sex stratifiction is added.
    :param StringExpression pop_expr: When specified, frequencies are stratified by population. If `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param StringExpression subpop_expr: When specified, frequencies are stratified by sub-continental population. Note that `pop_expr` is required as well when using this option.
    :param dict of str -> StringExpression additional_strata_expr: When specified, frequencies are stratified by the given additional strata found in the dict. This can e.g. be used to stratify by platform.
    :param list of int downsamplings: When specified, frequencies are computed by downsampling the data to the number of samples given in the list. Note that if `pop_expr` is specified, downsamplings by population is also computed.
    :return: MatrixTable with `freq` annotation
    """

    if subpop_expr is not None and pop_expr is None:
        raise NotImplementedError("annotate_freq requires pop_expr when using subpop_expr")

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
    mt = mt.annotate_cols(
        _freq_meta=_freq_meta_expr
    )

    # Get counters for sex, pop and subpop if set
    cut_dict = {
        cut: hl.agg.filter(hl.is_defined(mt._freq_meta[cut]), hl.agg.counter(mt._freq_meta[cut]))
        for cut in mt._freq_meta if cut != 'subpop'
    }
    if 'subpop' in mt._freq_meta:
        cut_dict['subpop'] = hl.agg.filter(
            hl.is_defined(mt._freq_meta.pop) & hl.is_defined(mt._freq_meta.subpop),
            hl.agg.counter(hl.struct(subpop=mt._freq_meta.subpop, pop=mt._freq_meta.pop))
        )

    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    sample_group_filters = []

    # Create downsamplings if needed
    if downsamplings is not None:
        # Add exact pop size downsampling if pops were provided
        if cut_data.get('pop'):
            downsamplings = list(set(downsamplings + list(cut_data.get('pop').values())))  # Add the pops values if not in yet
            downsamplings = sorted([x for x in downsamplings if x <= sum(cut_data.get('pop').values())])
        logger.info(f'Found {len(downsamplings)} downsamplings: {downsamplings}')

        # Shuffle the samples, then create a global index for downsampling
        # And a pop-index if pops were provided
        downsampling_ht = mt.cols()
        downsampling_ht = downsampling_ht.annotate(r=hl.rand_unif(0, 1))
        downsampling_ht = downsampling_ht.order_by(downsampling_ht.r)
        scan_expr = {'global_idx': hl.scan.count()}
        if cut_data.get('pop'):
            scan_expr['pop_idx'] = hl.scan.counter(downsampling_ht._freq_meta.pop).get(downsampling_ht._freq_meta.pop, 0)
        downsampling_ht = downsampling_ht.annotate(**scan_expr)
        downsampling_ht = downsampling_ht.key_by('s').select(*scan_expr)
        mt = mt.annotate_cols(downsampling=downsampling_ht[mt.s])
        mt = mt.annotate_globals(downsamplings=downsamplings)

        # Create downsampled sample groups
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': 'global'},
             mt.downsampling.global_idx < ds) for ds in downsamplings
        ])
        if cut_data.get('pop'):
            sample_group_filters.extend([
                ({'downsampling': str(ds), 'pop': pop},
                 (mt.downsampling.pop_idx < ds) & (mt._freq_meta.pop == pop))
                for ds in downsamplings for pop, pop_count in cut_data.get('pop', {}).items() if ds <= pop_count
            ])

    # Add all desired strata, starting with the full set and ending with downsamplings (if any)
    sample_group_filters = [({}, True)] + [
        ({'pop': pop}, mt._freq_meta.pop == pop) for pop in cut_data.get('pop', {})
        ] + [
           ({'sex': sex}, mt._freq_meta.sex == sex) for sex in cut_data.get('sex', {})
       ] + [
           ({'pop': pop, 'sex': sex}, (mt._freq_meta.sex == sex) & (mt._freq_meta.pop == pop))
           for sex in cut_data.get('sex', {}) for pop in cut_data.get('pop', {})
       ] + [
           ({'subpop': subpop.subpop, 'pop': subpop.pop}, (mt._freq_meta.pop == subpop.pop) & (mt._freq_meta.subpop == subpop.subpop))
           for subpop in cut_data.get('subpop', {})
       ] + [
           ({strata: str(s_value)}, mt._freq_meta[strata] == s_value)
           for strata in additional_strata_expr for s_value in cut_data.get(strata, {})
       ] + sample_group_filters

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta information
    freq_meta_expr = [dict(**sample_group[0], group='adj') for sample_group in sample_group_filters]
    freq_meta_expr.insert(1, {'group': 'raw'})
    mt = mt.annotate_globals(freq_meta=freq_meta_expr)

    # Create frequency expression array from the sample groups
    freq_expr = hl.agg.array_agg(
        lambda i: hl.agg.filter(mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)),
        hl.range(len(sample_group_filters))
    )

    # Insert raw as the second element of the array
    freq_expr = freq_expr[:1].extend([
        hl.agg.call_stats(mt.GT, mt.alleles)
    ]).extend(
        freq_expr[1:]
    )

    # Select non-ref allele (assumes bi-allelic)
    freq_expr = freq_expr.map(
        lambda cs: cs.annotate(
            AC=cs.AC[1],
            AF=cs.AF[1], #TODO This is NA in case AC and AN are 0 -- should we set it to 0?
            homozygote_count=cs.homozygote_count[1]
        )
    )

    # Return MT with freq row annotation
    return mt.annotate_rows(freq=freq_expr).drop('_freq_meta')
