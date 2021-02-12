import logging
from typing import Dict, List, Optional, Tuple

import hail as hl
import pyspark.sql

import gnomad.resources.grch37 as grch37_resources
import gnomad.resources.grch38 as grch38_resources
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
)
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.variant_qc.evaluation import compute_ranked_bin
from gnomad.variant_qc.random_forest import (
    get_features_importance,
    test_model,
    train_rf,
)
from gnomad.variant_qc.training import sample_training_examples

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

INBREEDING_COEFF_HARD_CUTOFF = -0.3


def create_binned_ht(
    ht: hl.Table,
    n_bins: int = 100,
    singleton: bool = True,
    biallelic: bool = True,
    adj: bool = True,
    add_substrat: Optional[Dict[str, hl.expr.BooleanExpression]] = None,
) -> hl.Table:
    """
    This is meant as a default wrapper for `compute_ranked_bin`. It annotates table with a bin, where variants are
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

    def _update_bin_expr(
        bin_expr: Dict[str, hl.expr.BooleanExpression],
        new_expr: hl.expr.BooleanExpression,
        new_id: str,
    ) -> Dict[str, hl.expr.BooleanExpression]:
        """
        Updates a dictionary of expressions to add another stratification

        :param bin_expr: Dictionary of expressions to add another stratification to
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

    # Desired bins and sub-bins
    bin_expr = {"bin": True}

    if singleton:
        bin_expr = _update_bin_expr(bin_expr, ht.ac_raw == 1, "singleton")

    if biallelic:
        bin_expr = _update_bin_expr(bin_expr, ~ht.was_split, "biallelic")

    if adj:
        bin_expr = _update_bin_expr(bin_expr, (ht.ac > 0), "adj")

    if add_substrat:
        for add_id, add_expr in add_substrat.items():
            bin_expr = _update_bin_expr(bin_expr, add_expr, add_id)

    bin_ht = compute_ranked_bin(
        ht, score_expr=ht.score, bin_expr=bin_expr, n_bins=n_bins
    )

    ht = ht.select_globals()
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
            - ac - expected that this is the allele count after adj filtering
            - ac_qc_samples_unrelated_raw - allele count before adj filtering for unrelated samples passing sample QC
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
            - n_untransmitted_raw

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
        - `n_de_novos_singleton_adj` - count of singleton de novo variants after adj filtration
        - `n_de_novo_singleton` - count of raw unfiltered singleton de novo variants
        - `n_de_novos_adj` - count of adj filtered de novo variants
        - `n_de_novos` - count of raw unfiltered de novo variants
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
        n_de_novos_singleton_adj=hl.agg.filter(
            ht.ac == 1, hl.agg.sum(fam.n_de_novos_adj)
        ),
        n_de_novo_singleton=hl.agg.filter(
            ht.ac_raw == 1, hl.agg.sum(fam.n_de_novos_raw)
        ),
        n_de_novos_adj=hl.agg.sum(fam.n_de_novos_adj),
        n_de_novo=hl.agg.sum(fam.n_de_novos_raw),
        n_trans_singletons=hl.agg.filter(
            ht.ac_raw == 2, hl.agg.sum(fam.n_transmitted_raw)
        ),
        n_untrans_singletons=hl.agg.filter(
            (ht.ac_raw < 3) & (ht.ac_qc_samples_unrelated_raw == 1),
            hl.agg.sum(fam.n_untransmitted_raw),
        ),
        n_train_trans_singletons=hl.agg.filter(
            (ht.ac_raw == 2) & ht.positive_train_site, hl.agg.sum(fam.n_transmitted_raw)
        ),
        n_omni=hl.agg.count_where(truth_data.omni),
        n_mills=hl.agg.count_where(truth_data.mills),
        n_hapmap=hl.agg.count_where(truth_data.hapmap),
        n_kgp_phase1_hc=hl.agg.count_where(truth_data.kgp_phase1_hc),
    )


def generate_trio_stats(
    mt: hl.MatrixTable, autosomes_only: bool = True, bi_allelic_only: bool = True
) -> hl.Table:
    """
    Default function to run `generate_trio_stats_expr` to get trio stats stratified by raw and adj

    .. note::

        Expects that `mt` is it a trio matrix table that was annotated with adj and if dealing with
        a sparse MT `hl.experimental.densify` must be run first.

        By default this pipeline function will filter `mt` to only autosomes and bi-allelic sites.

    :param mt: A Trio Matrix Table returned from `hl.trio_matrix`. Must be dense
    :param autosomes_only: If set, only autosomal intervals are used.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :return: Table with trio stats
    """
    if autosomes_only:
        mt = filter_to_autosomes(mt)
    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

    logger.info(f"Generating trio stats using {mt.count_cols()} trios.")
    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={"raw": True, "adj": trio_adj},
            de_novo_strata={"raw": True, "adj": trio_adj},
            ac_strata={"raw": True, "adj": trio_adj},
        )
    ).rows()

    return ht


def generate_sib_stats(
    mt: hl.MatrixTable,
    relatedness_ht: hl.Table,
    i_col: str = "i",
    j_col: str = "j",
    relationship_col: str = "relationship",
    autosomes_only: bool = True,
    bi_allelic_only: bool = True,
) -> hl.Table:
    """
    This is meant as a default wrapper for `generate_sib_stats_expr`. It returns a hail table with counts of variants
    shared by pairs of siblings in `relatedness_ht`.

    This function takes a hail Table with a row for each pair of individuals i,j in the data that are related (it's OK to have unrelated samples too).
    The `relationship_col` should be a column specifying the relationship between each two samples as defined by
    the constants in `gnomad.utils.relatedness`. This relationship_col will be used to filter to only pairs of
    samples that are annotated as `SIBLINGS`.

    .. note::

        By default this pipeline function will filter `mt` to only autosomes and bi-allelic sites.

    :param mt: Input Matrix table
    :param relatedness_ht: Input relationship table
    :param i_col: Column containing the 1st sample of the pair in the relationship table
    :param j_col: Column containing the 2nd sample of the pair in the relationship table
    :param relationship_col: Column containing the relationship for the sample pair as defined in this module constants.
    :param autosomes_only: If set, only autosomal intervals are used.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :return: A Table with the sibling shared variant counts
    """
    if autosomes_only:
        mt = filter_to_autosomes(mt)
    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

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

    sib_stats_ht = mt.select_rows(
        **generate_sib_stats_expr(
            mt, sib_ht, i_col=i_col, j_col=j_col, strata={"raw": True, "adj": mt.adj},
        )
    ).rows()

    return sib_stats_ht


def train_rf_model(
    ht: hl.Table,
    rf_features: List[str],
    tp_expr: hl.expr.BooleanExpression,
    fp_expr: hl.expr.BooleanExpression,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    test_expr: hl.expr.BooleanExpression = False,
) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    """
    Perform random forest (RF) training using a Table annotated with features and training data.

    .. note::

        This function uses `train_rf` and extends it by:
            - Adding an option to apply the resulting model to test variants which are withheld from training.
            - Uses a false positive (FP) to true positive (TP) ratio to determine what variants to use for RF training.

    The returned Table includes the following annotations:
        - rf_train: indicates if the variant was used for training of the RF model.
        - rf_label: indicates if the variant is a TP or FP.
        - rf_test: indicates if the variant was used in testing of the RF model.
        - features: global annotation of the features used for the RF model.
        - features_importance: global annotation of the importance of each feature in the model.
        - test_results: results from testing the model on variants defined by `test_expr`.

    :param ht: Table annotated with features for the RF model and the positive and negative training data.
    :param rf_features: List of column names to use as features in the RF training.
    :param tp_expr: TP training expression.
    :param fp_expr: FP training expression.
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used.
    :param num_trees: Number of trees in the RF model.
    :param max_depth: Maxmimum tree depth in the RF model.
    :param test_expr: An expression specifying variants to hold out for testing and use for evaluation only.
    :return: Table with TP and FP training sets used in the RF training and the resulting RF model.
    """

    ht = ht.annotate(_tp=tp_expr, _fp=fp_expr, rf_test=test_expr)

    rf_ht = sample_training_examples(
        ht, tp_expr=ht._tp, fp_expr=ht._fp, fp_to_tp=fp_to_tp, test_expr=ht.rf_test
    )
    ht = ht.annotate(rf_train=rf_ht[ht.key].train, rf_label=rf_ht[ht.key].label)

    summary = ht.group_by("_tp", "_fp", "rf_train", "rf_label", "rf_test").aggregate(
        n=hl.agg.count()
    )
    logger.info("Summary of TP/FP and RF training data:")
    summary.show(n=20)

    logger.info(
        "Training RF model:\nfeatures: {}\nnum_tree: {}\nmax_depth:{}".format(
            ",".join(rf_features), num_trees, max_depth
        )
    )

    rf_model = train_rf(
        ht.filter(ht.rf_train),
        features=rf_features,
        label="rf_label",
        num_trees=num_trees,
        max_depth=max_depth,
    )

    test_results = None
    if test_expr is not None:
        logger.info(f"Testing model on specified variants or intervals...")
        test_ht = ht.filter(hl.is_defined(ht.rf_label) & ht.rf_test)
        test_results = test_model(
            test_ht, rf_model, features=rf_features, label="rf_label"
        )

    features_importance = get_features_importance(rf_model)
    ht = ht.select_globals(
        features_importance=features_importance,
        features=rf_features,
        test_results=test_results,
    )

    return ht.select("rf_train", "rf_label", "rf_test"), rf_model
