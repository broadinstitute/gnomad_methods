# noqa: D100

import logging
from typing import Dict, Optional

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_ranked_bin(
    ht: hl.Table,
    score_expr: hl.expr.NumericExpression,
    bin_expr: Dict[str, hl.expr.BooleanExpression] = {"bin": True},
    compute_snv_indel_separately: bool = True,
    n_bins: int = 100,
    desc: bool = True,
) -> hl.Table:
    r"""
    Return a table with a bin for each row based on the ranking of `score_expr`.

    The bin is computed by dividing the `score_expr` into `n_bins` bins containing approximately equal numbers of elements.
    This is done by ranking the rows by `score_expr` (and a random number in cases where multiple variants have the same score)
    and then assigning the variant to a bin based on its ranking.

    If `compute_snv_indel_separately` is True all items in `bin_expr` will be stratified by snv / indels for the ranking and
    bin calculation. Because SNV and indel rows are mutually exclusive, they are re-combined into a single annotation. For
    example if we have the following four variants and scores and `n_bins` of 2:

    ========   =======   ======   =================   =================
    Variant    Type      Score    bin - `compute_snv_indel_separately`:
    --------   -------   ------   -------------------------------------
    \          \         \        False               True
    ========   =======   ======   =================   =================
    Var1       SNV       0.1      1                   1
    Var2       SNV       0.2      1                   2
    Var3       Indel     0.3      2                   1
    Var4       Indel     0.4      2                   2
    ========   =======   ======   =================   =================

    .. note::

        The `bin_expr` defines which data the bin(s) should be computed on. E.g., to get biallelic specific binning
        and singleton specific binning, the following could be used:

        .. code-block:: python

            bin_expr={
                'biallelic_bin': ~ht.was_split,
                'singleton_bin': ht.singleton
            }

    :param ht: Input Table
    :param score_expr: Expression containing the score
    :param bin_expr: Specific row grouping(s) to perform ranking and binning on (see note)
    :param compute_snv_indel_separately: Should all `bin_expr` items be stratified by SNVs / indels
    :param n_bins: Number of bins to bin the data into
    :param desc: Whether to bin the score in descending order
    :return: Table with the requested bin annotations
    """
    if compute_snv_indel_separately:
        # For each bin, add a SNV / indel stratification
        bin_expr = {
            f"{bin_id}_{snv}": (bin_expr & snv_expr)
            for bin_id, bin_expr in bin_expr.items()
            for snv, snv_expr in [
                ("snv", hl.is_snp(ht.alleles[0], ht.alleles[1])),
                ("indel", ~hl.is_snp(ht.alleles[0], ht.alleles[1])),
            ]
        }

    bin_ht = ht.select(
        **{f"_filter_{bin_id}": bin_expr for bin_id, bin_expr in bin_expr.items()},
        _score=score_expr,
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
        _rand=hl.rand_unif(0, 1),
    )

    logger.info(
        "Sorting the HT by score_expr followed by a random float between 0 and 1. "
        "Then adding a row index per grouping defined by bin_expr..."
    )
    bin_ht = bin_ht.order_by("_score", "_rand")
    bin_ht = bin_ht.annotate(
        **{
            f"{bin_id}_rank": hl.or_missing(
                bin_ht[f"_filter_{bin_id}"],
                hl.scan.count_where(bin_ht[f"_filter_{bin_id}"]),
            )
            for bin_id in bin_expr
        }
    )
    bin_ht = bin_ht.key_by("locus", "alleles")

    # Annotate globals with variant counts per group defined by bin_expr. This
    # is used to determine bin assignment
    bin_ht = bin_ht.annotate_globals(
        bin_group_variant_counts=bin_ht.aggregate(
            hl.Struct(
                **{
                    bin_id: hl.agg.filter(
                        bin_ht[f"_filter_{bin_id}"],
                        hl.agg.count(),
                    )
                    for bin_id in bin_expr
                }
            )
        )
    )

    logger.info("Binning ranked rows into %d bins...", n_bins)
    bin_ht = bin_ht.select(
        "snv",
        **{
            bin_id: hl.int(
                hl.floor(
                    (
                        n_bins
                        * (
                            bin_ht[f"{bin_id}_rank"]
                            / hl.float64(bin_ht.bin_group_variant_counts[bin_id])
                        )
                    )
                    + 1
                )
            )
            for bin_id in bin_expr
        },
    )

    if desc:
        bin_ht = bin_ht.annotate(
            **{bin_id: n_bins - bin_ht[bin_id] + 1 for bin_id in bin_expr}
        )

    # Because SNV and indel rows are mutually exclusive, re-combine them into a single bin.
    # Update the global bin_group_variant_counts struct to reflect the change
    # in bin names in the table
    if compute_snv_indel_separately:
        bin_expr_no_snv = {
            bin_id.rsplit("_", 1)[0] for bin_id in bin_ht.bin_group_variant_counts
        }
        bin_ht = bin_ht.annotate_globals(
            bin_group_variant_counts=hl.struct(
                **{
                    bin_id: hl.struct(
                        **{
                            snv: bin_ht.bin_group_variant_counts[f"{bin_id}_{snv}"]
                            for snv in ["snv", "indel"]
                        }
                    )
                    for bin_id in bin_expr_no_snv
                }
            )
        )

        bin_ht = bin_ht.transmute(
            **{
                bin_id: hl.if_else(
                    bin_ht.snv,
                    bin_ht[f"{bin_id}_snv"],
                    bin_ht[f"{bin_id}_indel"],
                )
                for bin_id in bin_expr_no_snv
            }
        )

    return bin_ht


def compute_grouped_binned_ht(
    bin_ht: hl.Table,
    checkpoint_path: Optional[str] = None,
) -> hl.GroupedTable:
    """
    Group a Table that has been annotated with bins (`compute_ranked_bin` or `create_binned_ht`).

    The table will be grouped by bin_id (bin, biallelic, etc.), contig, snv, bi_allelic and singleton.

    .. note::

        If performing an aggregation following this grouping (such as `score_bin_agg`) then the aggregation
        function will need to use `ht._parent` to get the origin Table from the GroupedTable for the aggregation

    :param bin_ht: Input Table with a `bin_id` annotation
    :param checkpoint_path: If provided an intermediate checkpoint table is created with all required annotations before shuffling.
    :return: Table grouped by bins(s)
    """
    # Explode the rank table by bin_id
    bin_ht = bin_ht.annotate(
        bin_groups=hl.array(
            [
                hl.Struct(bin_id=bin_name, bin=bin_ht[bin_name])
                for bin_name in bin_ht.bin_group_variant_counts
            ]
        )
    )
    bin_ht = bin_ht.explode(bin_ht.bin_groups)
    bin_ht = bin_ht.transmute(
        bin_id=bin_ht.bin_groups.bin_id, bin=bin_ht.bin_groups.bin
    )
    bin_ht = bin_ht.filter(hl.is_defined(bin_ht.bin))

    if checkpoint_path is not None:
        bin_ht.checkpoint(checkpoint_path, overwrite=True)
    else:
        bin_ht = bin_ht.persist()

    # Group by bin_id, bin and additional stratification desired and compute
    # QC metrics per bin
    return bin_ht.group_by(
        bin_id=bin_ht.bin_id,
        contig=bin_ht.locus.contig,
        snv=hl.is_snp(bin_ht.alleles[0], bin_ht.alleles[1]),
        bi_allelic=~bin_ht.was_split,
        singleton=bin_ht.singleton,
        release_adj=bin_ht.ac > 0,
        bin=bin_ht.bin,
    )._set_buffer_size(20000)


def compute_binned_truth_sample_concordance(
    ht: hl.Table,
    binned_score_ht: hl.Table,
    n_bins: int = 100,
    add_bins: Dict[str, hl.expr.BooleanExpression] = {},
) -> hl.Table:
    """
    Determine the concordance (TP, FP, FN) between a truth sample within the callset and the samples truth data grouped by bins computed using `compute_ranked_bin`.

    .. note::
        The input 'ht` should contain three row fields:
            - score: value to use for binning
            - GT: a CallExpression containing the genotype of the evaluation data for the sample
            - truth_GT: a CallExpression containing the genotype of the truth sample
        The input `binned_score_ht` should contain:
             - score: value used to bin the full callset
             - bin: the full callset bin

    'add_bins` can be used to add additional global and truth sample binning to the final binned truth sample
    concordance HT. The keys in `add_bins` must be present in `binned_score_ht` and the values in `add_bins`
    should be expressions on `ht` that define a subset of variants to bin in the truth sample. An example is if we want
    to look at the global and truth sample binning on only bi-allelic variants. `add_bins` could be set to
    {'biallelic_bin': ht.biallelic}.

    The table is grouped by global/truth sample bin and variant type and contains TP, FP and FN.

    :param ht: Input HT
    :param binned_score_ht: Table with the bin annotation for each variant
    :param n_bins: Number of bins to bin the data into
    :param add_bins: Dictionary of additional global bin columns (key) and the expr to use for binning the truth sample (value)
    :return: Binned truth sample concordance HT
    """
    # Annotate score and global bin
    indexed_binned_score_ht = binned_score_ht[ht.key]
    ht = ht.annotate(
        **{f"global_{bin_id}": indexed_binned_score_ht[bin_id] for bin_id in add_bins},
        **{f"_{bin_id}": bin_expr for bin_id, bin_expr in add_bins.items()},
        score=indexed_binned_score_ht.score,
        global_bin=indexed_binned_score_ht.bin,
    )

    # Annotate the truth sample bin
    bin_ht = compute_ranked_bin(
        ht,
        score_expr=ht.score,
        bin_expr={
            "truth_sample_bin": hl.expr.bool(True),
            **{f"truth_sample_{bin_id}": ht[f"_{bin_id}"] for bin_id in add_bins},
        },
        n_bins=n_bins,
    )
    ht = ht.join(bin_ht, how="left")

    bin_list = [
        hl.tuple(["global_bin", ht.global_bin]),
        hl.tuple(["truth_sample_bin", ht.truth_sample_bin]),
    ]
    bin_list.extend(
        [hl.tuple([f"global_{bin_id}", ht[f"global_{bin_id}"]]) for bin_id in add_bins]
    )
    bin_list.extend(
        [
            hl.tuple([f"truth_sample_{bin_id}", ht[f"truth_sample_{bin_id}"]])
            for bin_id in add_bins
        ]
    )

    # Explode the global and truth sample bins
    ht = ht.annotate(bin=bin_list)

    ht = ht.explode(ht.bin)
    ht = ht.annotate(bin_id=ht.bin[0], bin=hl.int(ht.bin[1]))

    # Compute TP, FP and FN by bin_id, variant type and bin
    return (
        ht.group_by("bin_id", "snv", "bin")
        .aggregate(
            # TP => allele is found in both data sets
            tp=hl.agg.count_where(ht.GT.is_non_ref() & ht.truth_GT.is_non_ref()),
            # FP => allele is found only in test data set
            fp=hl.agg.count_where(
                ht.GT.is_non_ref() & hl.or_else(ht.truth_GT.is_hom_ref(), True)
            ),
            # FN => allele is found in truth data only
            fn=hl.agg.count_where(
                hl.or_else(ht.GT.is_hom_ref(), True) & ht.truth_GT.is_non_ref()
            ),
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n_alleles=hl.agg.count(),
        )
        .repartition(5)
    )


def create_truth_sample_ht(
    mt: hl.MatrixTable, truth_mt: hl.MatrixTable, high_confidence_intervals_ht: hl.Table
) -> hl.Table:
    """
    Compute a table comparing a truth sample in callset vs the truth.

    :param mt: MT of truth sample from callset to be compared to truth
    :param truth_mt: MT of truth sample
    :param high_confidence_intervals_ht: High confidence interval HT
    :return: Table containing both the callset truth sample and the truth data
    """

    def split_filter_and_flatten_ht(
        truth_mt: hl.MatrixTable, high_confidence_intervals_ht: hl.Table
    ) -> hl.Table:
        """
        Split a truth sample MT, filter it to the given high confidence intervals, and then "flatten" it as a HT by annotating GT in a row field.

        :param truth_mt: Truth sample MT
        :param high_confidence_intervals_ht: High confidence intervals
        :return: Truth sample table with GT as a row annotation
        """
        assert truth_mt.count_cols() == 1

        if not "was_split" in truth_mt.row:
            truth_mt = hl.split_multi_hts(truth_mt)

        truth_mt = truth_mt.filter_rows(
            hl.is_defined(high_confidence_intervals_ht[truth_mt.locus])
        )
        rename_entries = {"GT": "_GT"}
        if "adj" in truth_mt.entry:
            rename_entries.update({"adj": "_adj"})

        truth_mt = truth_mt.rename(rename_entries)
        return truth_mt.annotate_rows(
            **{x: hl.agg.take(truth_mt[f"_{x}"], 1)[0] for x in rename_entries}
        ).rows()

    # Load truth sample MT,
    # restrict it to high confidence intervals
    # and flatten it to a HT by annotating GT in a row annotation
    truth_ht = split_filter_and_flatten_ht(truth_mt, high_confidence_intervals_ht)
    truth_ht = truth_ht.rename({f: f"truth_{f}" for f in truth_ht.row_value})

    #  Similarly load, filter and flatten callset truth sample MT
    ht = split_filter_and_flatten_ht(mt, high_confidence_intervals_ht)

    # Outer join of truth and callset truth and annotate the score and global bin
    ht = truth_ht.join(ht, how="outer")
    ht = ht.annotate(snv=hl.is_snp(ht.alleles[0], ht.alleles[1]))

    return ht


def add_rank(
    ht: hl.Table,
    score_expr: hl.expr.NumericExpression,
    subrank_expr: Optional[Dict[str, hl.expr.BooleanExpression]] = None,
) -> hl.Table:
    """
    Add rank based on the `score_expr`. Rank is added for snvs and indels separately.

    If one or more `subrank_expr` are provided, then subrank is added based on all sites for which the boolean expression is true.

    In addition, variant counts (snv, indel separately) is added as a global (`rank_variant_counts`).

    :param ht: input Hail Table containing variants (with QC annotations) to be ranked
    :param score_expr: the Table annotation by which ranking should be scored
    :param subrank_expr: Any subranking to be added in the form name_of_subrank: subrank_filtering_expr
    :return: Table with rankings added
    """
    key = ht.key
    if subrank_expr is None:
        subrank_expr = {}

    temp_expr = {"_score": score_expr}
    temp_expr.update({f"_{name}": expr for name, expr in subrank_expr.items()})
    rank_ht = ht.select(**temp_expr, is_snv=hl.is_snp(ht.alleles[0], ht.alleles[1]))

    rank_ht = rank_ht.key_by("_score").persist()
    scan_expr = {
        "rank": hl.if_else(
            rank_ht.is_snv,
            hl.scan.count_where(rank_ht.is_snv),
            hl.scan.count_where(~rank_ht.is_snv),
        )
    }
    scan_expr.update(
        {
            name: hl.or_missing(
                rank_ht[f"_{name}"],
                hl.if_else(
                    rank_ht.is_snv,
                    hl.scan.count_where(rank_ht.is_snv & rank_ht[f"_{name}"]),
                    hl.scan.count_where(~rank_ht.is_snv & rank_ht[f"_{name}"]),
                ),
            )
            for name in subrank_expr
        }
    )
    rank_ht = rank_ht.annotate(**scan_expr)

    rank_ht = rank_ht.key_by(*key).persist()
    rank_ht = rank_ht.select(*scan_expr.keys())

    ht = ht.annotate(**rank_ht[key])
    return ht
