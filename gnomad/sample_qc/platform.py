# noqa: D100

import logging
from typing import List, Optional, Tuple

import hail as hl
import numpy as np

from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_autosomes

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_callrate_mt(
    mt: hl.MatrixTable,
    intervals_ht: hl.Table,
    bi_allelic_only: bool = True,
    autosomes_only: bool = True,
    match: bool = True,
) -> hl.MatrixTable:
    """
    Compute a sample/interval MT with each entry containing the call rate for that sample/interval.

    This can be used as input for imputing exome sequencing platforms.

    .. note::

        The input interval HT should have a key of type Interval.
        The resulting table will have a key of the same type as the `intervals_ht` table and
        contain an `interval_info` field containing all non-key fields of the `intervals_ht`.

    :param mt: Input MT
    :param intervals_ht: Table containing the intervals. This table has to be keyed by locus.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :param autosomes_only: If set, only autosomal intervals are used.
    :param matches: If set, returns all intervals in intervals_ht that overlap the locus in the input MT.
    :return: Callrate MT
    """
    logger.info("Computing call rate MatrixTable")

    if len(intervals_ht.key) != 1 or not isinstance(
        intervals_ht.key[0], hl.expr.IntervalExpression
    ):
        logger.warning(
            "Call rate matrix computation expects `intervals_ht` with a key of type"
            " Interval. Found: %s",
            intervals_ht.key,
        )

    if autosomes_only:
        callrate_mt = filter_to_autosomes(mt)

    if bi_allelic_only:
        callrate_mt = callrate_mt.filter_rows(bi_allelic_expr(callrate_mt))

    intervals_ht = intervals_ht.annotate(_interval_key=intervals_ht.key)
    callrate_mt = callrate_mt.annotate_rows(
        _interval_key=intervals_ht.index(
            callrate_mt.locus, all_matches=match
        )._interval_key
    )

    if match:
        callrate_mt = callrate_mt.explode_rows("_interval_key")

    callrate_mt = callrate_mt.filter_rows(
        hl.is_defined(callrate_mt._interval_key.interval)
    )
    callrate_mt = callrate_mt.select_entries(
        GT=hl.or_missing(hl.is_defined(callrate_mt.GT), hl.struct())
    )
    callrate_mt = callrate_mt.group_rows_by(**callrate_mt._interval_key).aggregate(
        callrate=hl.agg.fraction(hl.is_defined(callrate_mt.GT))
    )
    intervals_ht = intervals_ht.drop("_interval_key")
    callrate_mt = callrate_mt.annotate_rows(
        interval_info=hl.struct(**intervals_ht[callrate_mt.row_key])
    )
    return callrate_mt


def run_platform_pca(
    callrate_mt: hl.MatrixTable,
    binarization_threshold: Optional[float] = 0.25,
    n_pcs: int = 10,
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Run PCA on a sample/interval MT with each entry containing the call rate.

    When `binzarization_threshold` is set, the callrate is transformed to a 0/1 value based on the threshold.
    E.g. with the default threshold of 0.25, all entries with a callrate < 0.25 are considered as 0s, others as 1s.

    :param callrate_mt: Input callrate MT
    :param binarization_threshold: binzarization_threshold. None is no threshold desired
    :param n_pcs: Number of PCs to compute
    :return: eigenvalues, scores_ht, loadings_ht
    """
    logger.info("Running platform PCA")

    if binarization_threshold is not None:
        callrate_mt = callrate_mt.annotate_entries(
            callrate=hl.int(callrate_mt.callrate > binarization_threshold)
        )
    # Center until Hail's PCA does it for you
    callrate_mt = callrate_mt.annotate_rows(
        mean_callrate=hl.agg.mean(callrate_mt.callrate)
    )
    callrate_mt = callrate_mt.annotate_entries(
        callrate=callrate_mt.callrate - callrate_mt.mean_callrate
    )
    eigenvalues, scores, loadings = hl.pca(
        callrate_mt.callrate,
        compute_loadings=True,
        k=n_pcs,
    )  # TODO:  Evaluate whether computing loadings is a good / worthy thing
    logger.info("Platform PCA eigenvalues: %s", eigenvalues)

    return eigenvalues, scores, loadings


def assign_platform_from_pcs(
    platform_pca_scores_ht: hl.Table,
    pc_scores_ann: str = "scores",
    hdbscan_min_cluster_size: Optional[int] = None,
    hdbscan_min_samples: int = None,
) -> hl.Table:
    """
    Assign platforms using HBDSCAN on the results of call rate PCA.

    :param platform_pca_scores_ht: Input table with the PCA score for each sample
    :param pc_scores_ann: Field containing the scores
    :param hdbscan_min_cluster_size: HDBSCAN `min_cluster_size` parameter. If not specified the smallest of 500 and 0.1*n_samples will be used.
    :param hdbscan_min_samples: HDBSCAN `min_samples` parameter
    :return: A Table with a `qc_platform` annotation containing the platform based on HDBSCAN clustering
    """
    import hdbscan

    logger.info("Assigning platforms based on platform PCA clustering")

    # Read and format data for clustering
    data = platform_pca_scores_ht.to_pandas()
    callrate_data = np.matrix(data[pc_scores_ann].tolist())
    logger.info("Assigning platforms to %d samples.", len(callrate_data))

    # Cluster data
    if hdbscan_min_cluster_size is None:
        hdbscan_min_cluster_size = min(500, 0.1 * data.shape[0])
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples
    )
    cluster_labels = clusterer.fit_predict(callrate_data)
    n_clusters = len(set(cluster_labels)) - (
        -1 in cluster_labels
    )  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info("Found %d unique platforms during platform imputation.", n_clusters)

    data["qc_platform"] = cluster_labels
    ht = hl.Table.from_pandas(data, key=[*platform_pca_scores_ht.key])
    ht = ht.annotate(qc_platform="platform_" + hl.str(ht.qc_platform))
    return ht
