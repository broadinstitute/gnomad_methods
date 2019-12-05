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
    return mt.annotate_cols(downsampling=global_ht[mt.s]), downsamplings


def add_faf_expr(freq: hl.expr.ArrayExpression, freq_meta: hl.expr.ArrayExpression, locus: hl.expr.LocusExpression, populations: Set[str]) -> hl.expr.ArrayExpression:
    """
    Calculates popmax (add an additional entry into freq with popmax: pop)

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


def generate_frequency_data(
        mt: hl.MatrixTable,
        calculate_age_hists: bool = True,
        calculate_faf: bool = True,
        calculate_by_platform: bool = False,
        calculate_downsampling: bool = False,
        calculate_project_max: bool = False,
        pops_to_remove_for_popmax: Optional[List[str]] = None
) -> Tuple[hl.Table, hl.Table]:
    """
    :param MatrixTable mt: Input MatrixTable
    :param bool calculate_downsampling: Calculate frequencies for downsampled data
    :param bool calculate_by_platform: Calculate frequencies for PCR-free data
    :param bool calculate_age_hists: Calculate age histograms for het and hom_var calls
    :param list of str pops_to_remove_for_popmax: Populations to remove for the popmax calculation (typically inbred/bottleneck pops)
    """
    if calculate_downsampling:
        mt, downsamplings = generate_downsamplings_cumulative(mt)
        print(f'Got {len(downsamplings)} downsamplings: {downsamplings}')
    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.meta.pop), hl.agg.counter(mt.meta.pop)),
                'sex': hl.agg.filter(hl.is_defined(mt.meta.sex), hl.agg.collect_as_set(mt.meta.sex)),
                'subpop': hl.agg.filter(hl.is_defined(mt.meta.subpop) & hl.is_defined(mt.meta.pop),
                                        hl.agg.collect_as_set(hl.struct(subpop=mt.meta.subpop, pop=mt.meta.pop)))
                }
    if calculate_by_platform:
        cut_dict['platform'] = hl.agg.filter(hl.is_defined(mt.meta.qc_platform),
                                             hl.agg.collect_as_set(mt.meta.qc_platform))
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))

    sample_group_filters = [({}, True)]
    sample_group_filters.extend([
        ({'pop': pop}, mt.meta.pop == pop) for pop in cut_data.pop
    ] + [
        ({'sex': sex}, mt.meta.sex == sex) for sex in cut_data.sex
    ] + [
        ({'pop': pop, 'sex': sex}, (mt.meta.sex == sex) & (mt.meta.pop == pop))
        for sex in cut_data.sex for pop in cut_data.pop
    ] + [
        ({'subpop': subpop.subpop, 'pop': subpop.pop},
         mt.meta.subpop == subpop.subpop)
        for subpop in cut_data.subpop
    ])

    if calculate_by_platform:
        sample_group_filters.extend([
            ({'platform': str(platform)}, mt.meta.qc_platform == platform)
            for platform in cut_data.platform
        ])

    if calculate_downsampling:
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': 'global'},
             mt.downsampling.global_idx < ds) for ds in downsamplings
        ])
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': pop},
             (mt.downsampling.pop_idx < ds) & (mt.meta.pop == pop))
            for ds in downsamplings for pop, pop_count in cut_data.pop.items() if ds <= pop_count
        ])
    mt = mt.select_cols(group_membership=[x[1] for x in sample_group_filters], project_id=mt.meta.project_id, age=mt.meta.age)
    mt = mt.select_rows()

    def get_meta_expressions(sample_group_filters):
        meta_expressions = []
        for i in range(len(sample_group_filters)):
            subgroup_dict = sample_group_filters[i][0]
            subgroup_dict['group'] = 'adj'
            meta_expressions.append(subgroup_dict)
        meta_expressions.insert(1, {'group': 'raw'})
        return meta_expressions

    def get_freq_expressions(mt, n_groups):

        adj_freq_expressions = hl.agg.array_agg(
            lambda i: hl.agg.filter(mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)),
            hl.range(n_groups)
        )

        # Insert raw as the second element of the array
        return adj_freq_expressions[:1].extend([
            hl.agg.call_stats(mt.GT, mt.alleles)
        ]).extend(
            adj_freq_expressions[1:]
        ).map(
            lambda cs: cs.annotate(
                AC=cs.AC[1],
                AF=cs.AF[1],
                homozygote_count=cs.homozygote_count[1]
            )
        )

    frequency_expression = get_freq_expressions(mt, len(sample_group_filters))
    print(f'Calculating {len(sample_group_filters) + 1} aggregators...')
    global_expression = {
        'freq_meta': get_meta_expressions(sample_group_filters)
    }
    mt = mt.annotate_rows(freq=frequency_expression)

    if calculate_age_hists:
        mt = mt.annotate_rows(age_hist_het=hl.agg.filter(mt.adj & mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
                              age_hist_hom=hl.agg.filter(mt.adj & mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10)))
    if calculate_downsampling: global_expression['downsamplings'] = downsamplings
    mt = mt.annotate_globals(**global_expression)
    sample_data = mt.cols()

    pops = set(cut_data.pop.keys())
    [pops.discard(x) for x in pops_to_remove_for_popmax]

    mt = mt.annotate_rows(popmax=add_popmax_expr(mt.freq, mt.freq_meta, populations=pops))

    if calculate_faf:
        mt = mt.annotate_rows(faf=add_faf_expr(mt.freq, mt.freq_meta, mt.locus, populations=pops))

    if calculate_project_max:
        mt = get_projectmax(mt, mt.project_id)

    return mt.rows(), sample_data