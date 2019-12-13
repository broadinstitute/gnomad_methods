from gnomad_hail import *


def generate_downsamplings_cumulative(mt: hl.MatrixTable, downsamplings: List[int]) -> Tuple[hl.MatrixTable, List[int]]:
    pop_data = [x[0] for x in get_sample_data(mt, [mt.meta.pop])]
    pops = Counter(pop_data)
    downsamplings = list(set(downsamplings + list(pops.values())))  # Add the pops values if not in yet
    downsamplings = sorted([x for x in downsamplings if x <= sum(pops.values())])
    downsamplings = sorted([x for x in downsamplings if x <= sum(pops.values())])
    ht = mt.cols()
    ht = ht.annotate(r=hl.rand_unif(0, 1))
    ht = ht.order_by(ht.r).add_index('global_idx')

    for i, pop in enumerate(pops):
        pop_ht = ht.filter(ht.meta.pop == pop).add_index('pop_idx')
        if not i:
            global_ht = pop_ht
        else:
            global_ht = global_ht.union(pop_ht)
    global_ht = global_ht.key_by('s')
    logger.info(f'Found {len(downsamplings)} downsamplings: {downsamplings}')
    return mt.annotate_cols(downsampling=global_ht[mt.s]), downsamplings


def popmax_expr(freq: hl.expr.ArrayExpression, freq_meta: hl.expr.ArrayExpression, populations: Set[str]) -> hl.expr.ArrayExpression:
    """
    Calculates popmax (add an additional entry into freq with popmax: pop)

    :param ArrayExpression freq: ArrayExpression of Structs with ['ac', 'an', 'hom']
    :param ArrayExpression freq_meta: ArrayExpression of meta dictionaries corresponding to freq
    :param set of str populations: Set of populations over which to calculate popmax
    :return: Frequency data with annotated popmax
    :rtype: ArrayExpression
    """
    pops_to_use = hl.literal(populations)
    freq = hl.map(lambda x: x[0].annotate(meta=x[1]), hl.zip(freq, freq_meta))
    freq_filtered = hl.filter(
        lambda f: (f.meta.size() == 2) & (f.meta.get('group') == 'adj') &
                  pops_to_use.contains(f.meta.get('pop')) & (f.AC > 0),
        freq
    )
    sorted_freqs = hl.sorted(freq_filtered, key=lambda x: x.AF, reverse=True)
    return hl.or_missing(
        hl.len(sorted_freqs) > 0,
        hl.struct(
            AC=sorted_freqs[0].AC,
            AF=sorted_freqs[0].AF,
            AN=sorted_freqs[0].AN,
            homozygote_count=sorted_freqs[0].homozygote_count,
            pop=sorted_freqs[0].meta['pop']
        )
    )


def project_max_expr(
        project_expr: hl.expr.StringExpression,
        gt_expr: hl.expr.CallExpression,
        alleles_expr: hl.expr.ArrayExpression,
        n_projects: int = 5
) -> hl.expr.ArrayExpression:
    """
    Creates the projectmax annotation, which is an array containing for each non-ref allele
    an array with AC, AN and AF for the `n_projects` with the largest AF at this row.
    Note that only projects with AF > 0 are returned.

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
        n_alleles > 1,
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
                    AC=[x[1].AC[0], x[1].AC[ai]],
                    AF=[x[1].AF[0], x[1].AF[ai]],
                    homozygote_count=[x[1].homozygote_count[0], x[1].homozygote_count[ai]],
                    project=x[0]
                )
            )
        )
    )


def faf_expr(freq: hl.expr.ArrayExpression, freq_meta: hl.expr.ArrayExpression, locus: hl.expr.LocusExpression, populations: Set[str]) -> hl.expr.ArrayExpression:
    """
    Calculates the filtering allele frequency

    :param ArrayExpression freq: ArrayExpression of Structs with ['ac', 'an', 'hom']
    :param ArrayExpression freq_meta: ArrayExpression of meta dictionaries corresponding to freq
    :param LocusExpression locus: LocusExpression
    :param set of str populations: Set of populations over which to calculate popmax
    :return: Frequency data with annotated popmax
    :rtype: ArrayExpression
    """
    pops_to_use = hl.literal(populations)
    freq = hl.map(lambda x: x[0].annotate(meta=x[1]), hl.zip(freq, freq_meta))
    freqs_to_use = hl.filter(lambda f:
                             ((f.meta.size() == 1) & (f.meta.get('group') == 'adj')) |
                             ((f.meta.size() == 2) & (f.meta.get('group') == 'adj') & pops_to_use.contains(f.meta.get('pop'))) |
                             (~locus.in_autosome_or_par() & (
                                     ((f.meta.size() == 2) & (f.meta.get('group') == 'adj') & f.meta.contains('sex')) |
                                     ((f.meta.size() == 3) & (f.meta.get('group') == 'adj') & pops_to_use.contains(f.meta.get('pop')) & f.meta.contains('sex')))),
                             freq)
    return freqs_to_use.map(lambda f: hl.struct(
        meta=f.meta,
        faf95=hl.experimental.filtering_allele_frequency(f.AC, f.AN, 0.95),
        faf99=hl.experimental.filtering_allele_frequency(f.AC, f.AN, 0.99)
    ))


def qual_hist_expr(
        gt_expr: Optional[hl.expr.CallExpression] = None,
        gq_expr: Optional[hl.expr.NumericExpression] = None,
        dp_expr: Optional[hl.expr.NumericExpression] = None,
        ad_expr: Optional[hl.expr.ArrayNumericExpression] = None
) -> hl.expr.StructExpression:
    """
    Returns an expression with genotype quality hsitograms based on the arguments given (dp, gq, ad).
    Notes
    -----
    - If `gt_expr` is provided, will return histograms for non-reference samples only as well as all samples.
    - `gt_expr` is required for the allele-balance histogram, as it is only computed on het samples.

    :param CallExpression gt_expr: Genotype
    :param NumericExpression gq_expr: Genotype quality
    :param NumericExpression dp_expr: Depth
    :param ArrayNumericExpression ad_expr: Allelic Depth (bi-allelic here)
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
            **{f'{qual}_all': qual_hist for qual, qual_hist in qual_hists},
            **{f'{qual}_alt': hl.agg.filter(gt_expr.is_non_ref(), qual_hist) for qual, qual_hist in qual_hists}
        }
        if ad_expr is not None:
            qual_hists['ab_hist_alt'] = hl.agg.filter(gt_expr.is_het(), hl.agg.hist(ad_expr[1] / hl.sum(ad_expr), 0, 1, 20))

    else:
        qual_hists = {f'{qual}_all': qual_hist for qual, qual_hist in qual_hists}

    return hl.struct(**qual_hists)


def generate_frequency_data(
        mt: hl.MatrixTable,
        sex_expr: Optional[hl.expr.StringExpression] = None,
        pop_expr: Optional[hl.expr.StringExpression] = None,
        subpop_expr: Optional[hl.expr.StringExpression] = None,
        platform_expr: Optional[hl.expr.StringExpression] = None,
        project_expr: Optional[hl.expr.StringExpression] = None,
        age_expr: Optional[hl.expr.NumericExpression] = None,
        calculate_faf: bool = True,
        calculate_popmax: bool = True,
        pops_to_remove_for_faf_and_popmax: Optional[List[str]] = None,
        downsamplings: Optional[List[int]] = None
) -> Tuple[hl.Table, hl.Table]:
    """
    Creates a table with allele frequencies by population, sex, subpopulation.
    Additionally, the following can also be computed:
    - age histograms
    - filtering allele frequencies
    - frequencies by platform
    - frequencies by downsampling the data to N samples (incl. by pop)
    - project max

    The input MT needs the following fields:
    - meta.pop
    - meta.sex

    Important note
    --------------
    Currently this only supports bi-allelic sites.

    :param MatrixTable mt: Input MatrixTable
    :param bool calculate_downsampling: Calculate frequencies for downsampled data
    :param bool calculate_by_platform: Calculate frequencies for PCR-free data
    :param bool calculate_age_hists: Calculate age histograms for het and hom_var calls
    :param list of str pops_to_remove_for_faf_and_popmax: Populations to remove for the popmax calculation (typically inbred/bottleneck pops)
    """

    # Annotate cols with provided cuts
    mt = mt.select_cols(
        **{
            name[:-5]: expr for name, expr in locals()
            if expr is not None and name in
               [ 'sex_expr', 'pop_expr', 'subpop_expr', 'platform_expr', 'project_expr', 'age_expr']
        }
    )

    # Create downsamplings if needed
    if downsamplings is not None:
        mt, downsamplings = generate_downsamplings_cumulative(mt, downsamplings)

    # Get counters for sex, pop and subpop if set
    cut_dict = {
        cut: hl.agg.filter(hl.is_defined(mt[cut]), hl.agg.counter(mt[cut]))
        for cut in ['sex', 'pop', 'platform'] if cut in mt.col_value
    }
    if 'subpop' in mt.col_value:
        cut_dict['subpop'] = hl.agg.filter(
            hl.is_defined(mt.pop) & hl.is_defined(mt.subpop),
            hl.agg.counter(hl.struct(subpop=mt.subpop, pop=mt.pop))
        )

    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))

    sample_group_filters = [({}, True)]
    sample_group_filters.extend(
        [
            ({'pop': pop}, mt.pop == pop) for pop in cut_data.get('pop', {})
        ] + [
            ({'sex': sex}, mt.sex == sex) for sex in cut_data.get('sex', {})
        ] + [
            ({'pop': pop, 'sex': sex}, (mt.sex == sex) & (mt.pop == pop))
            for sex in cut_data.get('sex', {}) for pop in cut_data.get('pop', {})
        ] + [
            ({'subpop': subpop.subpop, 'pop': subpop.pop}, mt.subpop == subpop.subpop)
            for subpop in cut_data.get('subpop', {})
        ] + [
            ({'platform': str(platform)}, mt.platform == platform)
            for platform in cut_data.get('platform', {})
        ]
    )

    if downsamplings is not None:
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': 'global'},
             mt.downsampling.global_idx < ds) for ds in downsamplings
        ])
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': pop},
             (mt.downsampling.pop_idx < ds) & (mt.pop == pop))
            for ds in downsamplings for pop, pop_count in cut_data.get('pop', {}).items() if ds <= pop_count
        ])

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta information
    meta_expressions = [dict(**sample_group[0], group='adj') for sample_group in sample_group_filters]
    meta_expressions.insert(1, {'group': 'raw'})
    global_expression = {
        'freq_meta': meta_expressions
    }
    if downsamplings is not None:
        global_expression['downsamplings'] = downsamplings

    mt = mt.annotate_globals(**global_expression)

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
            AF=cs.AF[1],
            homozygote_count=cs.homozygote_count[1]
        )
    )

    # Create row expressions
    row_expression = {'freq': freq_expr}

    if 'age' in mt.col_value:
        row_expression.update({
            'age_hist_het': hl.agg.filter(mt.adj & mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
            'age_hist_hom': hl.agg.filter(mt.adj & mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10))
        })

    if 'project' in mt.col_value:
        # Note that the [0] at the end is because the mt here is bi-allelic
        row_expression['project_max'] = project_max_expr(mt.project, mt.GT, mt.alleles, 5)[0]

    pops = {pop for pop in cut_data.get('pop', []) if pop not in pops_to_remove_for_faf_and_popmax}
    if calculate_popmax:
        row_expression['popmax'] = popmax_expr(freq_expr, mt.freq_meta, populations=pops)

    if calculate_faf:
        row_expression['faf'] = faf_expr(freq_expr, mt.freq_meta, mt.locus, populations=pops)

    mt = mt.select_rows(**row_expression)

    # Select col expressions
    cols = ['group_membership']
    if downsamplings:
        cols.append('downsampling')
    mt = mt.select_cols(*cols)

    return mt.rows(), mt.cols()