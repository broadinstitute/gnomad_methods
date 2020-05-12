import logging
from typing import Dict, List, Optional, Tuple, Union

import hail as hl
import pyspark.sql

import gnomad.resources.grch37 as grch37_resources
import gnomad.resources.grch38 as grch38_resources
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
)
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.variant_qc.evaluation import compute_quantile_bin
from gnomad.variant_qc.random_forest import (
    get_features_importance,
    test_model,
    train_rf,
)
from gnomad_qc.v2.variant_qc.variantqc import sample_rf_training_examples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def create_binned_ht(
    ht: hl.Table,
    n_bins: int = 100,
    singleton: bool = True,
    biallelic: bool = True,
    adj: bool = True,
    add_substrat: Optional[Dict[str, hl.expr.BooleanExpression]] = None,
) -> hl.Table:
    """
    This is meant as a default wrapper for `compute_quantile_bin`. It annotates table with a bin, where variants are
    binned based on score into `n_bins` equally-sized bins.

    .. note::

        The following fields should be present:
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
        new_id: str,
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
        bin_expr.update(
            {
                f"{new_id}_{bin_id}": bin_expr & new_expr
                for bin_id, bin_expr in bin_expr.items()
            }
        )
        return bin_expr

    ht = ht.annotate(
        singleton=ht.ac_raw == 1, snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
    )

    ht = ht.filter(ht.ac_raw > 0).persist()

    # Desired bins and sub-bins
    bin_expr = {"bin": True}

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
        ht, score_expr=ht.score, bin_expr=bin_expr, n_bins=n_bins
    )
    ht = ht.join(bin_ht, how="left")

    return ht


def score_bin_agg(
    ht: hl.GroupedTable, fam_stats_ht: hl.Table
) -> Dict[str, hl.expr.Aggregation]:
    """
    Default aggregation function to add aggregations for min/max of score, number of ClinVar variants, number of truth
    variants (omni, mills, hapmap, and kgp_phase1), and family statistics.

    .. note::

        This function uses `ht._parent` to get the origin Table from the GroupedTable for the aggregation

    This can easily be combined with the GroupedTable returned by `compute_grouped_binned_ht`

    Example:

    .. code-block:: python

        binned_ht = create_binned_ht(...)
        grouped_binned_ht = compute_grouped_binned_ht(binned_ht)
        agg_ht = grouped_binned_ht.aggregate(score_bin_agg(**grouped_binned_ht, ...))

    .. note::

        The following annotations should be present:

        In ht:
            - score
            - singleton
            - positive_train_site
            - negative_train_site
            - ac_raw - expected that this is the raw allele count before adj filtering
            - info - struct that includes QD, FS, and MQ in order to add an annotation for fail_hard_filters

        In truth_ht:
            - omni
            - mills
            - hapmap
            - kgp_phase1_hc

        In fam_stats_ht:
            - n_de_novos_adj
            - n_de_novos_raw
            - n_transmitted_raw
            - unrelated_qc_callstats
            - tdt

    Automatic aggregations that will be done are:
        - `min_score` - minimun of score annotation per group
        - `max_score` - maiximum of score annotation per group
        - `n` - count of variants per group
        - `n_ins` - count of insertion per group
        - `n_ins` - count of insertion per group
        - `n_del` - count of deletions per group
        - `n_ti` - count of transitions per group
        - `n_tv` - count of trnasversions per group
        - `n_1bp_indel` - count of one base pair indels per group
        - `n_mod3bp_indel` - count of indels with a length divisible by three per group
        - `n_singleton` - count of singletons per group
        - `fail_hard_filters` - count of variants per group with QD < 2 | FS > 60 | MQ < 30
        - `n_vqsr_pos_train` - count of variants that were a VQSR positive train site per group
        - `n_vqsr_neg_train` - count of variants that were a VQSR negative train site per group
        - `n_clinvar` - count of clinvar variants
        - `n_de_novos_adj` - count of adj filtered de dovo variants
        - `n_de_novos` - count of raw unfilterd filtered de dovo variants
        - `n_trans_singletons` - count of transmitted singletons
        - `n_untrans_singletons` - count of untransmitted singletons
        - `n_omni` - count of omni truth variants
        - `n_mills` - count of mills truth variants
        - `n_hapmap` - count of hapmap truth variants
        - `n_kgp_phase1_hc` - count of 1000 genomes phase 1 high confidence truth variants

    :param ht: Table that aggregation will be performed on
    :param fam_stats_ht: Path to family statistics HT
    :return: a dictionary containing aggregations to perform on ht
    """
    # Annotate binned table with the evaluation data
    ht = ht._parent
    indel_length = hl.abs(ht.alleles[0].length() - ht.alleles[1].length())
    # Load external evaluation data
    build = get_reference_genome(ht.locus).name
    clinvar = (
        grch37_resources.reference_data.clinvar
        if build == "GRCh37"
        else grch38_resources.reference_data.clinvar
    ).ht()[ht.key]
    truth_data = (
        grch37_resources.reference_data.get_truth_ht()
        if build == "GRCh37"
        else grch38_resources.reference_data.get_truth_ht()
    )[ht.key]
    fam = fam_stats_ht[ht.key]

    return dict(
        min_score=hl.agg.min(ht.score),
        max_score=hl.agg.max(ht.score),
        n=hl.agg.count(),
        n_ins=hl.agg.count_where(hl.is_insertion(ht.alleles[0], ht.alleles[1])),
        n_del=hl.agg.count_where(hl.is_deletion(ht.alleles[0], ht.alleles[1])),
        n_ti=hl.agg.count_where(hl.is_transition(ht.alleles[0], ht.alleles[1])),
        n_tv=hl.agg.count_where(hl.is_transversion(ht.alleles[0], ht.alleles[1])),
        n_1bp_indel=hl.agg.count_where(indel_length == 1),
        n_mod3bp_indel=hl.agg.count_where((indel_length % 3) == 0),
        n_singleton=hl.agg.count_where(ht.singleton),
        fail_hard_filters=hl.agg.count_where(
            (ht.info.QD < 2) | (ht.info.FS > 60) | (ht.info.MQ < 30)
        ),
        n_pos_train=hl.agg.count_where(ht.positive_train_site),
        n_neg_train=hl.agg.count_where(ht.negative_train_site),
        n_clinvar=hl.agg.count_where(hl.is_defined(clinvar)),
        n_de_novos_adj=hl.agg.sum(fam.n_de_novos_adj),
        n_de_novo=hl.agg.sum(fam.n_de_novos_raw),
        n_trans_singletons=hl.agg.filter(
            ht.ac_raw == 2, hl.agg.sum(fam.n_transmitted_raw)
        ),
        n_untrans_singletons=hl.agg.filter(
            (ht.ac_raw < 3) & (fam.unrelated_qc_callstats.AC[1] == 1),
            hl.agg.sum(fam.tdt.u),
        ),
        # n_train_trans_singletons=hl.agg.filter((ht.ac_raw == 2) & rank_ht.positive_train_site, hl.agg.sum(fam.n_transmitted_raw)),
        n_omni=hl.agg.count_where(truth_data.omni),
        n_mills=hl.agg.count_where(truth_data.mills),
        n_hapmap=hl.agg.count_where(truth_data.hapmap),
        n_kgp_phase1_hc=hl.agg.count_where(truth_data.kgp_phase1_hc),
    )


def generate_trio_stats(mt: hl.MatrixTable,) -> hl.Table:
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


def generate_sib_stats(
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

    sib_ht = relatedness_ht.filter(relatedness_ht[relationship_col] == SIBLINGS)
    s_to_keep = sib_ht.aggregate(
        hl.agg.explode(
            lambda s: hl.agg.collect_as_set(s), [sib_ht[i_col].s, sib_ht[j_col].s]
        ),
        _localize=False,
    )
    mt = mt.filter_cols(s_to_keep.contains(mt.s))
    if "adj" not in mt.entry:
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


def generate_rf_training(
    ht: hl.Table,
    rf_features: List[str],
    tp_expr: hl.expr.BooleanExpression,
    fp_expr: hl.expr.BooleanExpression,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    label_col: str = "rf_label",
    train_col: str = "rf_train",
    test_intervals: Union[str, List[str]] = "chr20",
    filter_expr: Optional[hl.expr.BooleanExpression] = None,
) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    """
    Perform random forest (RF) training using a Table annotated with features and training data.

    :param ht: Table annotated with features for the RF model and the positive and negative training data
    :param rf_features: List of column names to use as features in the RF training
    :param tp_expr: True positive (TP) training expression
    :param fp_expr: False positive (FP) training expression
    :param fp_to_tp Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used
    :param num_trees: Number of trees in the RF model
    :param max_depth: Maxmimum tree depth in the RF model
    :param str label_col: Name of column to store the training label
    :param str train_col: Name of column to store whether the site is used for training or not
    :param test_intervals: The specified interval(s) will be held out for testing and used for evaluation only.
    :param filter_expr: Can be used to filter to specified intervals before training
    :return: Table annotated with TP and FP training sets used in the RF training and the resulting RF model
    """
    ht = ht.annotate(tp=tp_expr, fp=fp_expr, _filter=filter_expr)

    test_intervals_str = (
        []
        if not test_intervals
        else [test_intervals]
        if isinstance(test_intervals, str)
        else test_intervals
    )
    test_intervals_locus = [
        hl.parse_locus_interval(x, reference_genome=get_reference_genome(ht.locus).name)
        for x in test_intervals_str
    ]

    if test_intervals_locus:
        ht = ht.annotate_globals(test_intervals=test_intervals_locus)

    ht = sample_rf_training_examples(
        ht,
        tp_col="tp",
        fp_col="fp",
        fp_to_tp=fp_to_tp,
        label_col=label_col,
        train_col=train_col,
    )
    ht = ht.persist()

    rf_ht = ht

    if filter_expr is not None:
        logger.info("Filtering training set using filter_expr")
        rf_ht = rf_ht.filter(rf_ht._filter)
    rf_ht.drop("_filter")

    summary = rf_ht.group_by("tp", "fp", train_col, label_col).aggregate(
        n=hl.agg.count()
    )
    logger.info("Summary of TP/FP and RF training labels:")
    summary.show(n=20)

    rf_ht = rf_ht.filter(rf_ht[train_col])

    logger.info(
        "Training RF model:\nfeatures: {}\nnum_tree: {}\nmax_depth:{}\nTest intervals: {}".format(
            ",".join(rf_features), num_trees, max_depth, ",".join(test_intervals_str),
        )
    )

    rf_model = train_rf(
        rf_ht,
        rf_features=rf_features,
        label=label_col,
        num_trees=num_trees,
        max_depth=max_depth,
    )

    test_results = None
    if test_intervals:
        logger.info(f"Testing model on intervals {','.join(test_intervals_str)}")
        test_ht = hl.filter_intervals(ht, test_intervals_locus, keep=True)
        test_ht = test_ht.filter(hl.is_defined(test_ht[label_col]))
        test_results = test_model(
            test_ht, rf_model, features=rf_features, label=label_col
        )
    logger.info("Writing RF training HT")
    features_importance = get_features_importance(rf_model)
    ht = ht.annotate_globals(
        features_importance=features_importance,
        features=rf_features,
        test_intervals=test_intervals_str,
        test_results=test_results,
    )

    return ht, rf_model


def generate_final_rf_ht(
    rf_result_ht: hl.Table,
    ac0_filter_expr: hl.expr.BooleanExpression,
    ts_ac_filter_expr: hl.expr.BooleanExpression,
    snp_cutoff: Union[int, float],
    indel_cutoff: Union[int, float],
    determine_cutoff_from_bin: bool = False,
    aggregated_bin_ht: Optional[hl.Table] = None,
    inbreeding_coeff_cutoff: float = -0.3,
) -> hl.Table:
    """
    Prepares finalized RF model given an RF result table from `rf.apply_rf_model` and cutoffs for filtering.

    If `determine_cutoff_from_bin` is True a `aggregated_bin_ht` must be supplied to determine the SNP and indel RF
    probabilities to use as cutoffs from an aggregated quantile bin Table like one created by
    `compute_grouped_binned_ht` in combination with `score_bin_agg`.

    :param rf_result_ht: RF result table from `rf.apply_rf_model` to prepare as the final RF Table
    :param ac0_filter_expr: Expression that indicates if a variant should be filtered as allele count 0 (AC0)
    :param ts_ac_filter_expr: Expression in `rf_result_ht` that indicates if a variant is a transmitted singleton
    :param snp_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for SNP variant QC filter
    :param indel_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for indel variant QC filter
    :param determine_cutoff_from_bin: If True the RF probability will be determined using bin info in `aggregated_bin_ht`
    :param aggregated_bin_ht: File with aggregate counts of variants based on quantile bins
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter to use for variants
    :return: Finalized random forest Table annotated with variant filters
    """
    # Determine SNP and indel RF cutoffs if given bin instead of RF probability
    if determine_cutoff_from_bin:
        snp_cutoff_global = hl.struct(min_score=snp_cutoff)
        indel_cutoff_global = hl.struct(min_score=indel_cutoff)
    else:
        snp_rf_cutoff, indel_rf_cutoff = aggregated_bin_ht.aggregate(
            [
                hl.agg.filter(
                    aggregated_bin_ht.snv & (aggregated_bin_ht.bin == snp_cutoff),
                    hl.agg.min(aggregated_bin_ht.min_score),
                ),
                hl.agg.filter(
                    ~aggregated_bin_ht.snv & (aggregated_bin_ht.bin == indel_cutoff),
                    hl.agg.min(aggregated_bin_ht.min_score),
                ),
            ]
        )
        snp_cutoff_global = hl.struct(bin=snp_cutoff, min_score=snp_rf_cutoff)
        indel_cutoff_global = hl.struct(bin=indel_cutoff, min_score=indel_rf_cutoff)

        logger.info(
            f"Using a SNP RF probability cutoff of {snp_rf_cutoff} and an indel RF probability cutoff of {indel_rf_cutoff}."
        )

    # Add filters to RF HT
    rf_result_ht = rf_result_ht.annotate_globals(
        rf_snv_cutoff=snp_cutoff_global, rf_indel_cutoff=indel_cutoff_global,
    )
    rf_filter_criteria = (
        hl.is_snp(rf_result_ht.alleles[0], rf_result_ht.alleles[1])
        & (rf_result_ht.rf_probability["TP"] < rf_result_ht.rf_snv_cutoff.min_score)
    ) | (
        ~hl.is_snp(rf_result_ht.alleles[0], rf_result_ht.alleles[1])
        & (rf_result_ht.rf_probability["TP"] < rf_result_ht.rf_indel_cutoff.min_score)
    )
    rf_result_ht = rf_result_ht.annotate(
        filters=hl.case()
        .when(rf_filter_criteria, {"RF"})
        .when(~rf_filter_criteria, hl.empty_set(hl.tstr))
        .or_error("Missing RF probability!")
    )

    inbreeding_coeff_filter_criteria = hl.is_defined(rf_result_ht.inbreeding_coeff) & (
        rf_result_ht.inbreeding_coeff < inbreeding_coeff_cutoff
    )
    rf_result_ht = rf_result_ht.annotate(
        filters=hl.cond(
            inbreeding_coeff_filter_criteria,
            rf_result_ht.filters.add("InbreedingCoeff"),
            rf_result_ht.filters,
        )
    )

    rf_result_ht = rf_result_ht.annotate(
        filters=hl.cond(
            ac0_filter_expr, rf_result_ht.filters.add("AC0"), rf_result_ht.filters
        )
    )

    # Fix annotations for release
    annotations_expr = {
        "tp": hl.or_else(rf_result_ht.tp, False),
        "transmitted_singleton": hl.or_missing(
            ts_ac_filter_expr, rf_result_ht.transmitted_singleton
        ),
        "rf_probability": rf_result_ht.rf_probability["TP"],
    }
    if "feature_imputed" in rf_result_ht.row:
        annotations_expr.update(
            {
                x: hl.or_missing(~rf_result_ht.feature_imputed[x], rf_result_ht[x])
                for x in [f for f in rf_result_ht.row.feature_imputed]
            }
        )

    rf_result_ht = rf_result_ht.transmute(**annotations_expr)

    return rf_result_ht
