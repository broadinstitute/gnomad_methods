from .generic import *


def get_lowqual_expr(
        alleles: hl.expr.ArrayExpression,
        qual_approx_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
        snv_phred_threshold: 30,
        snv_phred_het_prior: 30,  # 1/1000
        indel_phred_threshold: 30,
        indel_phred_het_prior: 39  # 1/8,000
) -> Union[hl.expr.BooleanExpression, hl.expr.ArrayExpression]:
    """
    Computes lowqual threshold expression for either split or unsplit alleles based on
    (AS_)QUALapprox

    :param ArrayExpression alleles: Array of alleles
    :param ArraynumericExpression or NumericExpression qual_approx_expr: QUALapprox or AS_QUALapprox
    :param int snv_phred_threshold: Phred-scaled SNV "emission" threshold (similar to GATK emission threshold)
    :param int snv_phred_het_prior: Phred-scaled SNV heterozygosity prior (30 = 1/1000 bases, GATK default)
    :param int indel_phred_threshold: Phred-scaled indel "emission" threshold (similar to GATK emission threshold)
    :param int indel_phred_het_prior: Phred-scaled indel heterozygosity prior (30 = 1/1000 bases, GATK default)
    :return: lowqual expression (BooleanExpression if `qual_approx_expr`is Numeric, Array[BooleanExpression] if `qual_approx_expr` is ArrayNumeric)
    :rtype: BooleanExpression or ArrayExpression
    """
    def low_qual_expr(ref: hl.expr.StringExpression, alt: hl.expr.StringExpression, qual_approx: hl.expr.NumericExpression) -> BooleanExpression:
        return hl.cond(
            hl.is_snp(ref, alt),
            qual_approx < snv_phred_threshold + snv_phred_het_prior,
            qual_approx < indel_phred_threshold + indel_phred_het_prior
        )
    if isinstance(qual_approx_expr, hl.expr.ArrayNumericExpression):
        return hl.range(1, hl.len(alleles)).map(lambda ai: low_qual_expr(alleles[0], alleles[ai], qual_approx_expr[ai - 1]))
    else:
        return low_qual_expr(alleles[0], alleles[1], qual_approx_expr)


def generate_fam_stats_expr(
        trio_mt: hl.MatrixTable,
        transmitted_strata: Dict[str, Optional[hl.expr.BooleanExpression]] = {'raw': None},
        de_novo_strata: Dict[str, Optional[hl.expr.BooleanExpression]] = {'raw': None},
        proband_is_female_expr: Optional[hl.expr.BooleanExpression] = None
) -> hl.expr.StructExpression:
    """
    Generates a row-wise expression containing the following counts:
    - Number of alleles in het parents transmitted to the proband
    - Number of alleles in het parents not transmitted to the proband
    - Number of de novo mutations

    Both transmission and de novo mutation metrics can be stratified using additional filters.
    If an empty dict is passed as one of the strata arguments, then this metric isn't computed.

    :param MatrixTable trio_mt: A trio standard trio MT (with the format as produced by hail.methods.trio_matrix
    :param dict of str -> BooleanExpression transmitted_strata: Strata for the transmission counts
    :param dict of str -> BooleanExpression de_novo_strata: Strata for the de novo counts
    :param BooleanExpression proband_is_female_expr: An optional expression giving the sex the proband. If not given, DNMs are only computed for autosomes.
    :return: An expression with the counts
    :rtype: StructExpression
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
        (hom_var, hom_var, het, hemi_x): (1, 0)
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

    def _get_composite_filter_expr(
            expr1: hl.expr.BooleanExpression,
            expr2: Optional[hl.expr.BooleanExpression]
    ) -> hl.expr.BooleanExpression:
        """
        Helper method to join two expression with support for None.
        """
        if expr2 is None:
            return expr1
        else:
            return expr1 & expr2

    def _is_dnm(
            proband_gt: hl.expr.CallExpression,
            father_gt: hl.expr.CallExpression,
            mother_gt: hl.expr.CallExpression,
            locus: hl.expr.LocusExpression,
            proband_is_female: Optional[hl.expr.BooleanExpression]
    ) -> hl.expr.BooleanExpression:
        """
        Helper method to get whether a given genotype combination is a DNM at a given locus with a given proband sex.
        """
        if proband_is_female is None:
            logger.warning("Since no proband sex expression was given to generate_fam_stats_expr, only DNMs in autosomes will be counted.")
            return hl.or_missing(
                locus.in_autosome(),
                proband_gt.is_het() & father_gt.is_hom_ref() & mother_gt.is_hom_ref()
            )
        return (
            hl.cond(
                locus.in_autosome_or_par() |
                (proband_is_female & locus.in_x_nonpar()),
                proband_gt.is_het() & father_gt.is_hom_ref() & mother_gt.is_hom_ref(),
                hl.or_missing(
                    ~proband_is_female,
                    proband_gt.is_het() & father_gt.is_hom_ref()
                )
            )
        )

    # Create transmission counters
    fam_stats = hl.struct(
        **{
            name: hl.agg.filter(
                _get_composite_filter_expr(trio_mt.proband_entry.GT.is_non_ref(), expr),
                hl.agg.sum(
                    trans_count_map.get(
                        (
                            trio_mt.proband_entry.GT.n_alt_alleles(),
                            trio_mt.father_entry.GT.n_alt_alleles(),
                            trio_mt.mother_entry.GT.n_alt_alleles(),
                            _get_copy_state(trio_mt.locus)
                        )
                    )
                )
            ) for name, expr in transmitted_strata.items()
        }
    )

    fam_stats = fam_stats.select(
        **{f'n_transmitted_{name}': fam_stats[name][0] for name in fam_stats},
        **{f'n_untransmitted_{name}': fam_stats[name][1] for name in fam_stats}
    )

    # Create de novo counters
    fam_stats = fam_stats.annotate({
        f'n_de_novo_{name}': hl.agg.filter(
            _get_composite_filter_expr(
                _is_dnm(
                    trio_mt.proband_entry.GT,
                    trio_mt.father_entry.GT,
                    trio_mt.mother_entry.GT,
                    trio_mt.locus,
                    proband_is_female_expr
                ), expr
            ),
            hl.agg.count()
        ) for name, expr in de_novo_strata.items()
    })

    return fam_stats


def compute_binned_rank(
        ht: hl.Table,
        score_expr: hl.expr.NumericExpression,
        rank_expr: Dict[str, hl.expr.BooleanExpression] = {'rank': True},
        stratify_snv_indel: bool = True,
        n_bins: int = 100,
        k: int = 1000,
        desc: bool = True
) -> hl.Table:
    """
    Returns a table containing a binned rank for each row.
    The bin is computed by dividing the `score_expr` into `n_bins` bins containing an equal number of elements.
    This is done based on quantiles computed with hl.agg.approx_quantiles.
    If a single value in `score_expr` spans more than one bin, the rows with this value are distributed
    randomly across the bins it spans.

    Notes
    -----
    The `rank_expr` defines which data the rank(s) should be computed on. E.g., to get an SNV rank and an Indel rank,
    the following could be used:
    rank_expr={
       'snv_rank': hl.is_snp(ht.alleles[0], ht.alleles[1]),
       'indels_rank': ~hl.is_snp(ht.alleles[0], ht.alleles[1])
    }

    :param Table ht: Input Table
    :param NumericExpression score_expr: Expression containing the score
    :param dict of str -> BooleanExpression rank_expr: Rank(s) to be computed (see notes)
    :param int n_bins: Number of bins to bin the data into
    :param int k: The `k` parameter of approx_quantiles
    :param bool desc: Whether to bin the score in descending order
    :return: Table with the binned ranks
    :rtype: Table
    """
    import math

    def quantiles_to_bin_boundaries(quantiles: List[int]) -> Dict:
        """
        Merges bins with the same boundaries into a unique bin while keeping track of
        which bins have been merged and the global index of all bins.

        :param quantiles: Original bins boundaries
        :return: (dict of the indices of bins for which multiple bins were collapsed -> number of bins collapsed,
                  Global indices of merged bins,
                  Merged bins boundaries)
        """

        # Pad the quantiles to create boundaries for the first and last bins
        bin_boundaries = [-math.inf] + quantiles + [math.inf]
        merged_bins = defaultdict(int)

        # If every quantile has a unique value, then bin boudaries are unique
        # and can be passed to binary_search as-is
        if len(quantiles) == len(set(quantiles)):
            return dict(
                merged_bins=merged_bins,
                global_bin_indices=list(range(len(bin_boundaries))),
                bin_boundaries=bin_boundaries
            )

        indexed_bins = list(enumerate(bin_boundaries))
        i = 1
        while i < len(indexed_bins):
            if indexed_bins[i - 1][1] == indexed_bins[i][1]:
                merged_bins[i - 1] += 1
                indexed_bins.pop(i)
            else:
                i += 1

        return dict(
            merged_bins=merged_bins,
            global_bin_indices=[x[0] for x in indexed_bins],
            bin_boundaries=[x[1] for x in indexed_bins]
        )

    if stratify_snv_indel:
        # For each bin, add a SNV / indel stratification
        bin_expr = {
            f'{bin_id}_{snv}': (bin_expr & snv_expr)
            for bin_id, bin_expr in bin_expr.items()
            for snv, snv_expr in [
                ('snv', ht.snv),
                ('indel', ~ht.snv)
            ]
        }

    ht = ht.annotate(
        **{f'_filter_{rid}': rexpr for rid, rexpr in bin_expr.items()},
        _score=score_expr
    )

    logger.info(f'Adding rank using approximate_quantiles binned into {n_bins}, using k={k}')
    rank_stats = ht.aggregate(
        hl.struct(
            **{
                rid: hl.agg.filter(
                    ht[f'_filter_{rid}'],
                    hl.struct(
                        n=hl.agg.count(),
                        quantiles=hl.agg.approx_quantiles(ht._score, [x / (n_bins) for x in range(1, n_bins)], k=k)
                    )
                )
                for rid in bin_expr
            }
        )
    )

    # Take care of bins with duplicated boundaries
    rank_stats = rank_stats.annotate(
        **{
            rname: rank_stats[rname].annotate(
                    **quantiles_to_bin_boundaries(rank_stats[rname].quantiles)
            ) for rname in rank_stats
        }
    )

    logger.debug(str(rank_stats))

    ht = ht.annotate_globals(
        rank_stats=hl.literal(
            rank_stats,
            dtype=hl.tstruct(**{
                bin_id: hl.tstruct(
                    n=hl.tint64,
                    quantiles=hl.tarray(hl.tfloat64),
                    bin_boundaries=hl.tarray(hl.tfloat64),
                    global_bin_indices=hl.tarray(hl.tint32),
                    merged_bins=hl.tdict(hl.tint32, hl.tint32)
                ) for bin_id in bin_expr
            })
        )
    )

    # Annotate the rank as the index in the unique boundaries array
    ht = ht.annotate(
        **{
            bin_id: hl.or_missing(
                ht[f'_filter_{bin_id}'],
                hl.binary_search(ht.rank_stats[bin_id].bin_boundaries, ht._score),
            ) for bin_id in bin_expr
        }
    )

    # Convert the rank to global rank by expanding merged bins, that is:
    # If a value falls in a bin that needs expansion, assign it randomly to one of the expanded bins
    # Otherwise, simply modify the rank bin to its global index (with expanded bins that is)
    ht = ht.select(
        **{
            bin_id: hl.cond(
                ht.rank_stats[bin_id].merged_bins.contains(ht[bin_id]),
                ht[bin_id] + hl.int(hl.rand_unif(0, ht.rank_stats[bin_id].merged_bins[ht[bin_id]] + 1)),
                ht.rank_stats[bin_id].global_bin_indices[ht[bin_id]]
            )
            for bin_id in bin_expr
        }
    )

    if desc:
        ht = ht.annotate(
            **{bin_id: n_bins - ht[bin_id] for bin_id in bin_expr}
        )

    # Annotate the HT with the bin
    # Because SNV and indel rows are mutually exclusive, re-combine them into a single bin.
    if stratify_snv_indel:
        ht = ht.annotate(
            **{
                bin_id: hl.cond(
                    ht.snv,
                    bin_ht[f'{bin_id}_snv'],
                    bin_ht[f'{bin_id}_indel']
                )
                for bin_id in bin_expr
            }
        )

    return ht


def create_binned_ht(ht: hl.Table, n_bins: int = 100) -> hl.Table:
    """
    Annotates table with a bin, where variants are binned based on score into `n_bins` equally-sized bins.
    Note that the following fields should be present:
    - score
    - ac - expected that this is the adj filtered allele count
    - ac_raw - expected that this is the raw allele count before adj filtering

    Computes bin numbers stratified by SNV / Indels and with the following sub bins
    - singletons
    - biallelics
    - biallelic singletons
    - adj
    - adj biallelics
    - adj singletons
    - adj biallelic singletons

    :param Table ht: Input table
    :param int n_bins: Number of bins to bin into
    :return table with bin number for each variant
    :rtype: Table
    """

    ht = ht.annotate(
        singleton=ht.ac_raw == 1,
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1])
    )

    ht = ht.filter(
        ht.ac_raw > 0
    ).persist()

    # Desired bins and sub-bins
    # TODO: add to a get default bins and change names to bin from rank and and make other function aggregate
    # TODO: make a parameter for singleton and adj, both boolean criteria
    bin_expr = {
        'bin': True,
        'singleton_bin': ht.singleton,
        'biallelic_bin': ~ht.was_split,
        'biallelic_singleton_bin': ~ht.was_split & ht.singleton,
        'adj_bin': ht.ac > 0,
        'adj_biallelic_bin': ~ht.was_split & (ht.ac > 0),
        'adj_singleton_bin': ht.singleton & (ht.ac > 0),
        'adj_biallelic_singleton_bin': ~ht.was_split & ht.singleton & (ht.ac > 0)
    }

    return compute_binned_rank(ht, ht.score, bin_expr, n_bins)[ht.key]


def compute_aggregate_binned_data(rank_ht: hl.Table, checkpoint_path: Optional[str] = None) -> hl.Table:
    """
    Creates binned data from a rank Table grouped by rank_id (rank, biallelic, etc.), contig, snv, bi_allelic and singleton
    containing the information needed for evaluation plots.

    Requires that `info` be annotate on the Table with a struct that includes QD, FS, and MQ

    Note that the following fields should be present:
    - ac
    - ac_raw

    :param Table rank_ht: Input rank Table
    :param str checkpoint_path: If provided an intermediate checkpoint table is created with all required annotations before shuffling.
    :return Table grouped by rank(s) and with counts of QC metrics
    :rtype Table
    """

    # Load external evaluation data
    clinvar_ht = hl.read_table(clinvar_ht_path)
    info_ht = hl.read_table(get_info_ht_path())
    # TODO: Can stay
    info_ht = info_ht.annotate(
        fail_hard_filters=(info_ht.info.QD < 2) | (info_ht.info.FS > 60) | (info_ht.info.MQ < 30)
    )
    ht_truth_data = hl.read_table(truth_ht_path)
    fam_ht = hl.read_table(fam_stats_ht_path)

    # Annotate rank table with the evaluation data
    rank_ht = rank_ht.annotate(
        indel_length=hl.abs(rank_ht.alleles[0].length() - rank_ht.alleles[1].length())
    )

    # Explode the rank table by rank_id
    rank_ht = rank_ht.annotate(
        rank_bins=hl.array([
            hl.Struct(
                rank_id=rank_name,
                bin=rank_ht[rank_name]
            )
            for rank_name in rank_ht.rank_stats
        ])
    )
    rank_ht = rank_ht.explode(rank_ht.rank_bins)
    rank_ht = rank_ht.transmute(
        rank_id=rank_ht.rank_bins.rank_id,
        bin=rank_ht.rank_bins.bin
    )
    rank_ht = rank_ht.filter(hl.is_defined(rank_ht.bin))

    if checkpoint_path is not None:
        rank_ht.checkpoint(checkpoint_path, overwrite=True)
    else:
        rank_ht = rank_ht.persist()

    # Group by rank_id, bin and additional stratification desired
    # and compute QC metrics per bin
    # TODO: pass a dict for these or pass a function that takes a HT and returns a dict of aggregators with the corresponding a
    #  where we have a default function
    return (
        rank_ht
            .group_by(
            rank_id=rank_ht.rank_id,
            contig=rank_ht.locus.contig,
            snv=hl.is_snp(rank_ht.alleles[0], rank_ht.alleles[1]),
            bi_allelic=~rank_ht.was_split,
            singleton=rank_ht.singleton,
            release_adj=rank_ht.ac > 0,
            bin=rank_ht.bin
        )._set_buffer_size(20000)
            .aggregate(
            min_score=hl.agg.min(rank_ht.score),
            max_score=hl.agg.max(rank_ht.score),
            n=hl.agg.count(),
            n_ins=hl.agg.count_where(hl.is_insertion(rank_ht.alleles[0], rank_ht.alleles[1])),
            n_del=hl.agg.count_where(hl.is_deletion(rank_ht.alleles[0], rank_ht.alleles[1])),
            n_ti=hl.agg.count_where(hl.is_transition(rank_ht.alleles[0], rank_ht.alleles[1])),
            n_tv=hl.agg.count_where(hl.is_transversion(rank_ht.alleles[0], rank_ht.alleles[1])),
            n_1bp_indel=hl.agg.count_where(rank_ht.indel_length == 1),
            n_mod3bp_indel=hl.agg.count_where((rank_ht.indel_length % 3) == 0),
            n_singleton=hl.agg.count_where(rank_ht.singleton),
            fail_hard_filters=hl.agg.count_where(rank_ht.fail_hard_filters),
            n_vqsr_pos_train=hl.agg.count_where(rank_ht.positive_train_site),
            n_vqsr_neg_train=hl.agg.count_where(rank_ht.negative_train_site)
        )
    )


def default_aggregators(ht: hl.Table) -> Dict[str, hl.expr.Aggregator]:
    #    Requires that `info` be annotate on the Table with a struct that includes QD, FS, and MQ
    # Load external evaluation data
    clinvar = hl.read_table(clinvar_ht_path)[ht.key]
    info = hl.read_table(get_info_ht_path())[ht.key]
    # TODO: Can stay
    info = info.annotate(
        fail_hard_filters=(info.info.QD < 2) | (info.info.FS > 60) | (info.info.MQ < 30)
    )
    truth_data = hl.read_table(truth_ht_path)[ht.key]
    fam = hl.read_table(fam_stats_ht_path)[ht.key]

    return dict(
        n_clinvar=hl.agg.count_where(hl.is_defined(clinvar)),
        n_de_novos_hq=hl.agg.sum(fam.n_de_novos_hq),
        n_de_novos_adj=hl.agg.sum(fam.n_de_novos_adj),
        n_de_novo=hl.agg.sum(fam.n_de_novos_raw),
        n_trans_singletons=hl.agg.filter(ht.ac_raw == 2, hl.agg.sum(fam.n_transmitted_raw)),
        n_untrans_singletons=hl.agg.filter((ht.ac_raw < 3) & (fam.unrelated_qc_callstats.AC[1] == 1),
                                           hl.agg.sum(fam.tdt.u)),  # TODO adapt names
        # n_train_trans_singletons=hl.agg.filter((ht.ac_raw == 2) & rank_ht.positive_train_site, hl.agg.sum(fam.n_transmitted_raw)),
        n_omni=hl.agg.count_where(truth_data.omni),
        n_mills=hl.agg.count_where(truth_data.mills),
        n_hapmap=hl.agg.count_where(truth_data.hapmap),
        n_kgp_phase1_hc=hl.agg.count_where(truth_data.kgp_phase1_hc),
        fail_hard_filters=hl.agg.count_where(rank_ht.fail_hard_filters)
    )


def compute_binned_truth_sample_concordance(
        ht: hl.Table,
        binned_rank_score_ht: hl.Table
) -> hl.Table:
    """
    The input HT should contain two row fields:
    * GT: a CallExpression containing the genotype of the evaluation data for the sample
    * truth_GT: a CallExpression containing the genotype of the truth sample

    The table is grouped by global/truth sample rank, bin and variant type and
    contains TP, FP and FN.

    :param Table ht: Input HT
    :param Table binned_rank_score_ht: Table with the binned rank for each variant
    :return: Binned truth sample concordance HT
    :rtype: Table
    """

    # Annotate score and global rank
    indexed_ranked_score_ht = binned_rank_score_ht[ht.key]
    ht = ht.annotate(
    score=indexed_ranked_score_ht.score,
    global_rank=indexed_ranked_score_ht.rank
    )

    # Annotate the truth sample rank
    ht = compute_binned_rank(
        ht,
        score_expr=ht.score,
        rank_expr={'truth_sample_rank': hl.expr.bool(True)},
        n_bins=100
    )

    # Explode the global and truth sample ranks
    ht = ht.annotate(rank=[
        hl.tuple(['global_rank', ht.global_rank]),
        hl.tuple(['truth_sample_rank', ht.truth_sample_rank])
    ])

    ht = ht.explode(ht.rank)
    ht = ht.annotate(
        rank_id=ht.rank[0],
        bin=hl.int(ht.rank[1])
    )

    # Compute TP, FP and FN by rank_id, variant type and bin
    return ht.group_by(
        'rank_id',
        'snv',
        'bin'
    ).aggregate(
        # TP => allele is found in both data sets
        tp=hl.agg.count_where(ht.GT.is_non_ref() & ht.truth_GT.is_non_ref()),
        # FP => allele is found only in test data set
        fp=hl.agg.count_where(ht.GT.is_non_ref() & hl.or_else(ht.truth_GT.is_hom_ref(), True)),
        # FN => allele is found in truth data only
        fn=hl.agg.count_where(ht.GT.is_hom_ref() & hl.or_else(ht.truth_GT.is_non_ref(), True)),
        min_score=hl.agg.min(ht.score),
        max_score=hl.agg.max(ht.score),
        n_alleles=hl.agg.count()
    ).repartition(5)


def create_truth_sample_ht(
        mt: hl.MatrixTable,
        truth_mt: hl.MatrixTable,
        high_confidence_intervals_ht: hl.Table,
        keep_lowqual: bool
) -> hl.Table:
    """
    Computes a table comparing a truth sample in callset vs the truth.

    :param MatrixTable mt: MT of truth sample from callset to be compared to truth
    :param MatrixTable truth_mt: MT of truth sample
    :param Table high_confidence_intervals_ht: High confidence interval HT
    :param bool keep_lowqual: If False lowqual variants are removed
    :return:
    :rtype: Table
    """

    def split_filter_and_flatten_ht(truth_mt: hl.MatrixTable, high_confidence_intervals_ht: hl.Table) -> hl.Table:
        """
        Splits a truth sample MT and filter it to the given high confidence intervals.
        Then "flatten" it as a HT by annotating GT in a row field.

        :param MatrixTable truth_mt: Truth sample MT
        :param Table high_confidence_intervals_ht: High confidence intervals
        :return: Truth sample table with GT as a row annotation
        :rtype: Table
        """
        assert(truth_mt.count_cols() == 1)

        if not 'was_split' in truth_mt.row:
            truth_mt = hl.split_multi_hts(truth_mt)

        truth_mt = truth_mt.filter_rows(
            hl.is_defined(high_confidence_intervals_ht[truth_mt.locus])
        )
        truth_mt = truth_mt.rename({'GT': '_GT'})
        return truth_mt.annotate_rows(GT=hl.agg.take(truth_mt._GT, 1)[0]).rows()

    # Load truth sample MT,
    # restrict it to high confidence intervals
    # and flatten it to a HT by annotating GT in a row annotation
    truth_ht = split_filter_and_flatten_ht(
        truth_mt,
        high_confidence_intervals_ht
    )
    truth_ht = truth_ht.rename({f: f'truth_{f}' for f in truth_ht.row_value})

    #  Similarly load, filter and flatten callset truth sample MT
    ht = split_filter_and_flatten_ht(
        mt,
        high_confidence_intervals_ht
    )


    # Outer join of truth and callset truth and annotate the score and global rank bin
    ht = truth_ht.join(ht, how="outer")
    ht = ht.annotate(
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1])
    )

    return ht
