from .generic import *


def get_lowqual_expr(
        alleles: hl.expr.ArrayExpression,
        qual_approx_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
        snp_phred_threshold: 30,
        snp_phred_het_prior: 30,  # 1/1000
        indel_phred_threshold: 30,
        indel_phred_het_prior: 39  # 1/8,000
) -> Union[hl.expr.BooleanExpression, hl.expr.ArrayExpression]:
    """
    Computes lowqual threshold expression for either split or unsplit alleles based on
    (AS_)QUALapprox

    :param ArrayExpression alleles: Expression
    :param qual_approx_expr:
    :return:
    """
    def low_qual_expr(ref: hl.expr.StringExpression, alt: hl.expr.StringExpression, qual_approx: hl.expr.NumericExpression) -> BooleanExpression:
        return hl.cond(
            hl.is_snp(ref, alt),
            qual_approx < snp_phred_threshold + snp_phred_het_prior,
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

    :param MatrixTable trio_mt: A trio stndard trio MT (with the format as produced by hail.methods.trio_matrix
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
                trio_mt.proband_entry.GT.is_het() & trio_mt.father_entry.GT.is_hom_ref() & trio_mt.mother_entry.GT.is_hom_ref()
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
    fam_stats_expr = {
        f'n_transmitted_{name}': hl.agg.filter(
            _get_composite_filter_expr(trio_mt.proband_entry.GT.is_non_ref(), expr),
            hl.agg.sum(
                trans_count_map.get(
                    (
                        trio_mt.proband_entry.GT.n_alt_alleles(),
                        trio_mt.father_entry.GT.n_alt_alleles(),
                        trio_mt.mother_entry.GT.n_alt_alleles(),
                        _get_copy_state(trio_mt.locus)
                    )
                )[0]
            )
        ) for name, expr in transmitted_strata.items()
    }

    fam_stats_expr.update({
        f'n_untransmitted_{name}': hl.agg.filter(
            _get_composite_filter_expr(trio_mt.proband_entry.GT.is_non_ref(), expr),
            hl.agg.sum(
                trans_count_map.get(
                    (
                        trio_mt.proband_entry.GT.n_alt_alleles(),
                        trio_mt.father_entry.GT.n_alt_alleles(),
                        trio_mt.mother_entry.GT.n_alt_alleles(),
                        _get_copy_state(trio_mt.locus)
                    )
                )[1]
            )
        ) for name, expr in transmitted_strata.items()
    })

    # Create de novo counters
    fam_stats_expr.update({
        f'n_de_novo_{name}': hl.agg.filter(
            _get_composite_filter_expr(
                _is_dnm(
                    trio_mt.proband_entry.GT,
                    trio_mt.mother_entry.GT,
                    trio_mt.father_entry.GT,
                    trio_mt.locus,
                    proband_is_female_expr
                ), expr
            ),
            hl.agg.count()
        ) for name, expr in de_novo_strata.items()
    })

    return hl.struct(**fam_stats_expr)