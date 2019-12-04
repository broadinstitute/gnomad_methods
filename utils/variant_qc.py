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
        proband_is_female_expr: Optional[hl.expr.BooleanExpression] =  None
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
        **{
            f'n_transmitted_{name}': fam_stats[name][0],
            f'n_untransmitted_{name}': fam_stats[name][1]
        } for name in fam_stats
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

    ht = ht.annotate(
        **{f'_filter_{rid}': rexpr for rid, rexpr in rank_expr.items()},
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
                for rid in rank_expr
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
                rank_id: hl.tstruct(
                    n=hl.tint64,
                    quantiles=hl.tarray(hl.tfloat64),
                    bin_boundaries=hl.tarray(hl.tfloat64),
                    global_bin_indices=hl.tarray(hl.tint32),
                    merged_bins=hl.tdict(hl.tint32, hl.tint32)
                ) for rank_id in rank_expr
            })
        )
    )

    # Annotate the rank as the index in the unique boundaries array
    ht = ht.annotate(
        **{
            rank_id: hl.or_missing(
                ht[f'_filter_{rank_id}'],
                hl.binary_search(ht.rank_stats[rank_id].bin_boundaries, ht._score),
            ) for rank_id in rank_expr
        }
    )

    # Convert the rank to global rank by expanding merged bins, that is:
    # If a value falls in a bin that needs expansion, assign it randomly to one of the expanded bins
    # Otherwise, simply modify the rank bin to its global index (with expanded bins that is)
    ht = ht.select(
        **{
            rank_id: hl.cond(
                ht.rank_stats[rank_id].merged_bins.contains(ht[rank_id]),
                ht[rank_id] + hl.int(hl.rand_unif(0, ht.rank_stats[rank_id].merged_bins[ht[rank_id]] + 1)),
                ht.rank_stats[rank_id].global_bin_indices[ht[rank_id]]
            )
            for rank_id in rank_expr
        }
    )

    if desc:
        ht = ht.annotate(
            **{rank_id: n_bins - ht[rank_id] for rank_id in rank_expr}
        )

    return ht