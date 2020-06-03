import logging
from typing import List, Optional, Tuple

import hail as hl
from pprint import pformat

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def sample_training_examples(
    ht: hl.Table,
    tp_expr: hl.BooleanExpression,
    fp_expr: hl.BooleanExpression,
    fp_to_tp: float = 1.0,
    test_expr: Optional[hl.expr.BooleanExpression] = None,
) -> hl.Table:
    """
    Returns a Table of all positive and negative training examples in `ht` with an annotation for those that should be
    used for training given a true positive (TP) to false positive (FP) ratio.

    The returned Table has the following annotations:
        - train: indicates if the variant should be used for training. The only variants labeled as False will be
          those defined by `test_expr`.
        - label: indicates if a variant is a 'TP' or 'FP' and will also be labeled as such for variants defined by `test_expr`.

    .. note::

        - This function does not support multi-allelic variants.
        - The function will give some stats about the TPs/FPs provided (Ti, Tv, indels).
        - If a site qualifies as both TP and FP, it's excluded from both and therefore not in the returned Table.

    :param ht: Input Table.
    :param tp_expr: Expression for TP examples.
    :param fp_expr: Expression for FP examples.
    :param fp_to_tp: FP to TP ratio. If set to <= 0, all training examples are used.
    :param test_expr: Optional expression to exclude a set of variants from training set. Still contains TP/FP label annotation.
    :return: Table subset with corresponding TP and FP examples with desired FP to TP ratio.
    """

    ht = ht.select(
        _tp=hl.or_else(tp_expr, False),
        _fp=hl.or_else(fp_expr, False),
        _exclude=False if test_expr is None else test_expr,
    )
    ht = ht.filter((ht._tp | ht._fp) & ~(ht._tp & ht._fp)).persist()

    # Get stats about TP / FP sets
    def _get_train_counts(ht: hl.Table) -> Tuple[int, int]:
        """
        Determine the number of TP and FP variants in the input Table and report some stats on Ti, Tv, indels.

        :param ht: Input Table
        :return: Counts of TP and FP variants in the table
        """
        train_stats = hl.struct(n=hl.agg.count())

        if "alleles" in ht.row and ht.row.alleles.dtype == hl.tarray(hl.tstr):
            train_stats = train_stats.annotate(
                ti=hl.agg.count_where(
                    hl.expr.is_transition(ht.alleles[0], ht.alleles[1])
                ),
                tv=hl.agg.count_where(
                    hl.expr.is_transversion(ht.alleles[0], ht.alleles[1])
                ),
                indel=hl.agg.count_where(
                    hl.expr.is_indel(ht.alleles[0], ht.alleles[1])
                ),
            )

        # Sample training examples
        pd_stats = (
            ht.group_by(**{"contig": ht.locus.contig, "tp": ht._tp, "fp": ht._fp})
            .aggregate(**train_stats)
            .to_pandas()
        )

        logger.info(pformat(pd_stats))
        pd_stats = pd_stats.fillna(False)

        # Number of true positive and false positive variants to be sampled for the training set
        n_tp = pd_stats[pd_stats["tp"] & ~pd_stats["fp"]]["n"].sum()
        n_fp = pd_stats[~pd_stats["tp"] & pd_stats["fp"]]["n"].sum()

        return n_tp, n_fp

    n_tp, n_fp = _get_train_counts(ht.filter(~ht._exclude))

    train_expr = True
    if fp_to_tp > 0:
        desired_fp = fp_to_tp * n_tp
        if desired_fp < n_fp:
            prob_tp = 1.0
            prob_fp = desired_fp / n_fp
        else:
            prob_tp = n_fp / desired_fp
            prob_fp = 1.0

        logger.info(
            f"Training examples sampling: tp={prob_tp}*{n_tp}, fp={prob_fp}*{n_fp}"
        )

        if prob_fp < 1.0:
            train_expr = hl.if_else(
                ht._fp & hl.or_else(~ht._tp, True),
                hl.rand_bool(prob_fp),
                hl.if_else(
                    hl.or_else(ht._tp, False), hl.or_else(~ht._fp, True), ht._fp,
                ),
            )
        elif prob_tp < 1.0:
            train_expr = hl.if_else(
                ht._tp & hl.or_else(~ht._fp, True),
                hl.rand_bool(prob_tp),
                hl.if_else(
                    hl.or_else(ht._fp, False), hl.or_else(~ht._tp, True), ht._fp,
                ),
            )
    else:
        logger.info(f"Using all {n_tp} TP and {n_fp} FP training examples.")

    label_expr = (
        hl.case(missing_false=True)
        .when(ht._tp & hl.or_else(~ht._fp, True), "TP")
        .when(ht._fp & hl.or_else(~ht._tp, True), "FP")
        .default(hl.null(hl.tstr))
    )

    return ht.select(train=train_expr & ~ht._exclude, label=label_expr)
