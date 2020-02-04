import itertools
from .generic import *
import gnomad_hail.resources.grch37 as grch37_resources
import gnomad_hail.resources.grch38 as grch38_resources


def get_lowqual_expr(
        alleles: hl.expr.ArrayExpression,
        qual_approx_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
        snv_phred_threshold: int = 30,
        snv_phred_het_prior: int = 30,  # 1/1000
        indel_phred_threshold: int = 30,
        indel_phred_het_prior: int = 39  # 1/8,000
) -> Union[hl.expr.BooleanExpression, hl.expr.ArrayExpression]:
    """
    Computes lowqual threshold expression for either split or unsplit alleles based on QUALapprox or AS_QUALapprox

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
            qual_approx:
            hl.expr.NumericExpression
    ) -> BooleanExpression:
        return hl.cond(
            hl.is_snp(ref, alt),
            qual_approx < snv_phred_threshold + snv_phred_het_prior,
            qual_approx < indel_phred_threshold + indel_phred_het_prior
        )
    if isinstance(qual_approx_expr, hl.expr.ArrayNumericExpression):
        return hl.range(1, hl.len(alleles)).map(
            lambda ai: low_qual_expr(alleles[0], alleles[ai], qual_approx_expr[ai - 1])
        )
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

    :param trio_mt: A trio standard trio MT (with the format as produced by hail.methods.trio_matrix
    :param transmitted_strata: Strata for the transmission counts
    :param de_novo_strata: Strata for the de novo counts
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
        **dict(itertools.chain.from_iterable(
            [
                (f'n_transmitted_{name}', fam_stats[name][0]),
                (f'n_untransmitted_{name}', fam_stats[name][1]),
            ]
            for name in fam_stats
        ))
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


def compute_quantile_bin(
        ht: hl.Table,
        score_expr: hl.expr.NumericExpression,
        bin_expr: Dict[str, hl.expr.BooleanExpression] = {'bin': True},
        compute_snv_indel_separately: bool = True,
        n_bins: int = 100,
        k: int = 1000,
        desc: bool = True
) -> hl.Table:
    """
    Returns a table with a bin for each row based on quantiles of `score_expr`.

    The bin is computed by dividing the `score_expr` into `n_bins` bins containing an equal number of elements.
    This is done based on quantiles computed with hl.agg.approx_quantiles. If a single value in `score_expr` spans more
    than one bin, the rows with this value are distributed randomly across the bins it spans.

    If `compute_snv_indel_separately` is True all items in `bin_expr` will be stratified by snv / indels for the bin
    calculation. Because SNV and indel rows are mutually exclusive, they are re-combined into a single annotation. For
    example if we have the following four variants and scores and `n_bins` of 4:

    ========   =======   ======   ==================================   =================================
    Variant    Type      Score    bin                                  bin
                                  compute_snv_indel_separately=False   compute_snv_indel_separately=True
    ========   =======   ======   ==================================   =================================
    Var1       SNV       0.1      1                                    1
    Var2       SNV       0.2      2                                    2
    Var3       Indel     0.3      3                                    1
    Var4       Indel     0.4      4                                    2
    ========   =======   ======   ==================================   =================================

    .. note::

        The `bin_expr` defines which data the bin(s) should be computed on. E.g., to get a biallelic quantile bin and an
        singleton quantile bin, the following could be used:

        .. code-block:: python

            bin_expr={
                'biallelic_bin': ~ht.was_split,
                'singleton_bin': ht.singleton
            }

    :param ht: Input Table
    :param score_expr: Expression containing the score
    :param bin_expr: Quantile bin(s) to be computed (see notes)
    :param compute_snv_indel_separately: Should all `bin_expr` items be stratified by snv / indels
    :param n_bins: Number of bins to bin the data into
    :param k: The `k` parameter of approx_quantiles
    :param desc: Whether to bin the score in descending order
    :return: Table with the quantile bins
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

    if compute_snv_indel_separately:
        # For each bin, add a SNV / indel stratification
        bin_expr = {
            f'{bin_id}_{snv}': (bin_expr & snv_expr)
            for bin_id, bin_expr in bin_expr.items()
            for snv, snv_expr in [
                ('snv', ht.snv),
                ('indel', ~ht.snv)
            ]
        }

    bin_ht = ht.annotate(
        **{f'_filter_{bin_id}': bin_expr for bin_id, bin_expr in bin_expr.items()},
        _score=score_expr
    )

    logger.info(f'Adding quantile bins using approximate_quantiles binned into {n_bins}, using k={k}')
    bin_stats = bin_ht.aggregate(
        hl.struct(
            **{
                bin_id: hl.agg.filter(
                    bin_ht[f'_filter_{bin_id}'],
                    hl.struct(
                        n=hl.agg.count(),
                        quantiles=hl.agg.approx_quantiles(bin_ht._score, [x / (n_bins) for x in range(1, n_bins)], k=k)
                    )
                )
                for bin_id in bin_expr
            }
        )
    )

    # Take care of bins with duplicated boundaries
    bin_stats = bin_stats.annotate(
        **{
            rname: bin_stats[rname].annotate(
                    **quantiles_to_bin_boundaries(bin_stats[rname].quantiles)
            ) for rname in bin_stats
        }
    )

    bin_ht = bin_ht.annotate_globals(
        bin_stats=hl.literal(
            bin_stats,
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

    # Annotate the bin as the index in the unique boundaries array
    bin_ht = bin_ht.annotate(
        **{
            bin_id: hl.or_missing(
                bin_ht[f'_filter_{bin_id}'],
                hl.binary_search(bin_ht.bin_stats[bin_id].bin_boundaries, bin_ht._score),
            ) for bin_id in bin_expr
        }
    )

    # Convert the bin to global bin by expanding merged bins, that is:
    # If a value falls in a bin that needs expansion, assign it randomly to one of the expanded bins
    # Otherwise, simply modify the bin to its global index (with expanded bins that is)
    bin_ht = bin_ht.select(
        **{
            bin_id: hl.cond(
                bin_ht.bin_stats[bin_id].merged_bins.contains(bin_ht[bin_id]),
                bin_ht[bin_id] + hl.int(hl.rand_unif(0, bin_ht.bin_stats[bin_id].merged_bins[bin_ht[bin_id]] + 1)),
                bin_ht.bin_stats[bin_id].global_bin_indices[bin_ht[bin_id]]
            )
            for bin_id in bin_expr
        }
    )

    if desc:
        bin_ht = bin_ht.annotate(
            **{bin_id: n_bins - bin_ht[bin_id] for bin_id in bin_expr}
        )

    # Because SNV and indel rows are mutually exclusive, re-combine them into a single bin.
    # Update the global bin_stats struct to reflect the change in bin names in the table
    if compute_snv_indel_separately:
        bin_expr_no_snv = {bin_id.rsplit("_", 1)[0] for bin_id in bin_ht.bin_stats}
        bin_ht = bin_ht.annotate_globals(
            bin_stats=hl.struct(
                **{
                    bin_id: hl.struct(
                        **{snv: bin_ht.bin_stats[f"{bin_id}_{snv}"] for snv in ["snv", "indel"]}
                    )
                    for bin_id in bin_expr_no_snv
                }
            )
        )

        bin_ht = bin_ht.transmute(
            **{
                bin_id: hl.cond(
                    ht[bin_ht.key].snv,
                    bin_ht[f'{bin_id}_snv'],
                    bin_ht[f'{bin_id}_indel']
                )
                for bin_id in bin_expr_no_snv
            }
        )

    return bin_ht


def create_binned_ht(
        ht: hl.Table,
        n_bins: int = 100,
        singleton: bool = True,
        biallelic: bool = True,
        adj: bool = True,
        add_substrat: Optional[Dict[str, hl.expr.BooleanExpression]] = None
) -> hl.Table:
    """
    Annotates table with a bin, where variants are binned based on score into `n_bins` equally-sized bins.
    Note that the following fields should be present:
    - score
    - ac - expected that this is the adj filtered allele count
    - ac_raw - expected that this is the raw allele count before adj filtering

    Computes bin numbers stratified by SNV / Indels and with the following optional sub bins
    - singletons
    - biallelics
    - biallelic singletons
    - adj
    - adj biallelics
    - adj singletons
    - adj biallelic singletons

    :param ht: Input table
    :param n_bins: Number of bins to bin into
    :param singleton: Should bins be stratified by singletons
    :param biallelic: Should bins be stratified by bi-alleleic variants
    :param adj: Should bins be stratified by adj filtering
    :param add_substrat: Any additional stratifications for adding bins
    :return: table with bin number for each variant
    """
    def update_bin_expr(
            bin_expr: Dict[str, hl.expr.BooleanExpression],
            new_expr: hl.expr.BooleanExpression,
            new_id: str
    ) -> Dict[str, hl.expr.BooleanExpression]:
        """
        Updates a dictionary of expressions to add another stratification

        :param bin_expr: Dictionary of expressions to add another
        stratification to
        :param new_expr: New Boolean expression to add to `bin_expr`
        :param new_id: Name to add to each current key in `bin_expr` to indicate the new stratification
        :return: Dictionary of `bin_expr` updated with `new_expr` added as an additional stratification to all
        expressions already in `bin_expr`
        """
        bin_expr.update({
            f'{new_id}_{bin_id}': bin_expr & new_expr
            for bin_id, bin_expr in bin_expr.items()
        })
        return bin_expr

    ht = ht.annotate(
        singleton=ht.ac_raw == 1,
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),

    )

    ht = ht.filter(
        ht.ac_raw > 0
    ).persist()

    # Desired bins and sub-bins
    bin_expr = {
        'bin': True
    }

    if singleton:
        bin_expr = update_bin_expr(bin_expr, ht.singleton, "singleton")

    if biallelic:
        bin_expr = update_bin_expr(bin_expr, ~ht.was_split, "biallelic")

    if adj:
        bin_expr = update_bin_expr(bin_expr, (ht.ac > 0), "adj")

    if add_substrat:
        for add_id, add_expr in add_substrat.items():
            bin_expr = update_bin_expr(bin_expr, add_expr, add_id)

    bin_ht = compute_quantile_bin(
        ht,
        score_expr=ht.score,
        bin_expr=bin_expr,
        n_bins=n_bins
    )
    ht = ht.join(bin_ht, how='left')

    return ht


def default_score_bin_agg(
        ht: hl.Table,
        truth_ht: hl.Table,
        fam_stats_ht: hl.Table
) -> Dict[str, hl.expr.Aggregation]:
    """
    Default aggregation function to pass to `compute_aggregate_binned_data` to add aggregations for number of ClinVar
    variants, number of truth variants (omni, mills, hapmap, and kgp_phase1), and family statistics.

    Note that the following fields should be present:

    In ht:
        - ac_raw - expected that this is the raw allele count before adj filtering
    In truth_ht (truth_data annotation):
        - omni
        - mills
        - hapmap
        - kgp_phase1_hc
    In fam_stats_ht:
        - n_de_novos_hq
        - n_de_novos_adj
        - n_de_novos_raw
        - n_transmitted_raw
        - unrelated_qc_callstats
        - tdt

    :param ht: Table that aggregation will be performed on
    :param truth_ht: Path to truth sites HT
    :param fam_stats_ht: Path to family statistics HT
    :return: a dictionary containing aggrecations to perform on ht
    """
    # Load external evaluation data
    build = get_reference_genome(ht.locus).name
    clinvar = (grch37_resources.reference_data.clinvar if build == 'GRCh37' else grch38_resources.reference_data.clinvar).ht()[ht.key]
    truth_data = truth_ht[ht.key].truth_data
    fam = fam_stats_ht[ht.key]

    return dict(
        n_clinvar=hl.agg.count_where(hl.is_defined(clinvar)),
        n_de_novos_hq=hl.agg.sum(fam.n_de_novos_hq),
        n_de_novos_adj=hl.agg.sum(fam.n_de_novos_adj),
        n_de_novo=hl.agg.sum(fam.n_de_novos_raw),
        n_trans_singletons=hl.agg.filter(ht.ac_raw == 2, hl.agg.sum(fam.n_transmitted_raw)),
        n_untrans_singletons=hl.agg.filter((ht.ac_raw < 3) & (fam.unrelated_qc_callstats.AC[1] == 1),
                                           hl.agg.sum(fam.tdt.u)),
        ## n_train_trans_singletons=hl.agg.filter((ht.ac_raw == 2) & rank_ht.positive_train_site, hl.agg.sum(fam.n_transmitted_raw)),
        n_omni=hl.agg.count_where(truth_data.omni),
        n_mills=hl.agg.count_where(truth_data.mills),
        n_hapmap=hl.agg.count_where(truth_data.hapmap),
        n_kgp_phase1_hc=hl.agg.count_where(truth_data.kgp_phase1_hc)
    )


def compute_aggregate_binned_data(
        bin_ht: hl.Table,
        agg_func: Callable = default_score_bin_agg,
        checkpoint_path: Optional[str] = None,
        **kwargs
) -> hl.Table:
    """
    Aggregates a Table that has been annotated with bins based on quantiles (`compute_quantile_bin` or
    `create_binned_ht`). The table will be grouped by bin_id (bin, biallelic, etc.), contig, snv, bi_allelic and
    singleton. Then for each grouping, min/max of `score` will be computed and any other desired aggregations.

    Automatic aggregations that will be done are:
        `min_score` - minimun of score annotation per group
        `max_score` - maiximum of score annotation per group
        `n` - count of variants per group
        `n_ins` - count of insertion per group
        `n_ins` - count of insertion per group
        `n_del` - count of deletions per group
        `n_ti` - count of transitions per group
        `n_tv` - count of trnasversions per group
        `n_1bp_indel` - count of one base pair indels per group
        `n_mod3bp_indel` - count of indels with a length divisible by three per group
        `n_singleton` - count of singletons per group
        `fail_hard_filters` - count of variants per group with QD < 2 | FS > 60 | MQ < 30
        `n_vqsr_pos_train` - count of variants that were a VQSR positive train site per group
        `n_vqsr_neg_train` - count of variants that were a VQSR negative train site per group

    Requires that `bin_ht` be annotated with an `info` struct that includes QD, FS, and MQ in order to add an
    annotation for `fail_hard_filters`

    :param bin_ht: Input Table with a `bin_id` annotation
    :param agg_func: Function that returns a dict of any additional aggregations to perform
    :param checkpoint_path: If provided an intermediate checkpoint table is created with all required annotations before shuffling.
    :return: Table grouped by rank(s) and with counts of QC metrics
    """
    # Annotate binned table with the evaluation data
    bin_ht = bin_ht.annotate(
        indel_length=hl.abs(bin_ht.alleles[0].length() - bin_ht.alleles[1].length()),
        fail_hard_filters=(bin_ht.info.QD < 2) | (bin_ht.info.FS > 60) | (bin_ht.info.MQ < 30)
    )

    # Explode the rank table by bin_id
    bin_ht = bin_ht.annotate(
        quantile_bins=hl.array([
            hl.Struct(
                bin_id=bin_name,
                bin=bin_ht[bin_name]
            )
            for bin_name in bin_ht.bin_stats
        ])
    )
    bin_ht = bin_ht.explode(bin_ht.quantile_bins)
    bin_ht = bin_ht.transmute(
        bin_id=bin_ht.quantile_bins.bin_id,
        bin=bin_ht.quantile_bins.bin
    )
    bin_ht = bin_ht.filter(hl.is_defined(bin_ht.bin))

    if checkpoint_path is not None:
        bin_ht.checkpoint(checkpoint_path, overwrite=True)
    else:
        bin_ht = bin_ht.persist()

    # Group by bin_id, bin and additional stratification desired and compute QC metrics per bin
    return (
        bin_ht.group_by(
            bin_id=bin_ht.bin_id,
            contig=bin_ht.locus.contig,
            snv=hl.is_snp(bin_ht.alleles[0], bin_ht.alleles[1]),
            bi_allelic=~bin_ht.was_split,
            singleton=bin_ht.singleton,
            release_adj=bin_ht.ac > 0,
            bin=bin_ht.bin
        )._set_buffer_size(20000)
            .aggregate(
            min_score=hl.agg.min(bin_ht.score),
            max_score=hl.agg.max(bin_ht.score),
            n=hl.agg.count(),
            n_ins=hl.agg.count_where(hl.is_insertion(bin_ht.alleles[0], bin_ht.alleles[1])),
            n_del=hl.agg.count_where(hl.is_deletion(bin_ht.alleles[0], bin_ht.alleles[1])),
            n_ti=hl.agg.count_where(hl.is_transition(bin_ht.alleles[0], bin_ht.alleles[1])),
            n_tv=hl.agg.count_where(hl.is_transversion(bin_ht.alleles[0], bin_ht.alleles[1])),
            n_1bp_indel=hl.agg.count_where(bin_ht.indel_length == 1),
            n_mod3bp_indel=hl.agg.count_where((bin_ht.indel_length % 3) == 0),
            n_singleton=hl.agg.count_where(bin_ht.singleton),
            fail_hard_filters=hl.agg.count_where(bin_ht.fail_hard_filters),
            n_vqsr_pos_train=hl.agg.count_where(bin_ht.positive_train_site),
            n_vqsr_neg_train=hl.agg.count_where(bin_ht.negative_train_site),
            **agg_func(bin_ht, **kwargs)
        )
    )


def compute_binned_truth_sample_concordance(
        ht: hl.Table,
        binned_score_ht: hl.Table
) -> hl.Table:
    """
    The input HT should contain two row fields:
    * GT: a CallExpression containing the genotype of the evaluation data for the sample
    * truth_GT: a CallExpression containing the genotype of the truth sample

    The table is grouped by global/truth sample bin and variant type and
    contains TP, FP and FN.

    :param ht: Input HT
    :param binned_score_ht: Table with the an annotation for quantile bin for each variant
    :return: Binned truth sample concordance HT
    """

    # Annotate score and global bin
    indexed_binned_score_ht = binned_score_ht[ht.key]
    ht = ht.annotate(
        score=indexed_binned_score_ht.score,
        global_bin=indexed_binned_score_ht.bin
    )

    # Annotate the truth sample quantile bin
    bin_ht = compute_quantile_bin(
        ht,
        score_expr=ht.score,
        bin_expr={'truth_sample_bin': hl.expr.bool(True)},
        n_bins=100
    )
    ht = ht.join(bin_ht, how='left')

    # Explode the global and truth sample bins
    ht = ht.annotate(bin=[
        hl.tuple(['global_bin', ht.global_bin]),
        hl.tuple(['truth_sample_bin', ht.truth_sample_bin])
    ])

    ht = ht.explode(ht.bin)
    ht = ht.annotate(
        bin_id=ht.bin[0],
        bin=hl.int(ht.bin[1])
    )

    # Compute TP, FP and FN by bin_id, variant type and bin
    return ht.group_by(
        'bin_id',
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
        high_confidence_intervals_ht: hl.Table
) -> hl.Table:
    """
    Computes a table comparing a truth sample in callset vs the truth.

    :param mt: MT of truth sample from callset to be compared to truth
    :param truth_mt: MT of truth sample
    :param high_confidence_intervals_ht: High confidence interval HT
    :return: Table containing both the callset truth sample and the truth data
    """

    def split_filter_and_flatten_ht(truth_mt: hl.MatrixTable, high_confidence_intervals_ht: hl.Table) -> hl.Table:
        """
        Splits a truth sample MT and filter it to the given high confidence intervals.
        Then "flatten" it as a HT by annotating GT in a row field.

        :param truth_mt: Truth sample MT
        :param high_confidence_intervals_ht: High confidence intervals
        :return: Truth sample table with GT as a row annotation
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

    # Outer join of truth and callset truth and annotate the score and global bin
    ht = truth_ht.join(ht, how="outer")
    ht = ht.annotate(
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1])
    )

    return ht
