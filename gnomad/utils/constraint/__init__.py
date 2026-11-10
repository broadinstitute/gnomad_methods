"""Script containing generic constraint functions that may be used in the constraint pipeline."""

import logging

from gnomad.utils.constraint.constants import (
    CLASSIC_LOF_ANNOTATIONS,
    COVERAGE_CUTOFF,
    DEFAULT_FIELDS_TO_SUM,
    DEFAULT_GENCODE_ANNOTATIONS,
)
from gnomad.utils.constraint.context import (
    annotate_mutation_type,
    annotate_with_mu,
    assemble_constraint_context_ht,
    calculate_gerp_cutoffs,
    collapse_strand,
    get_mu_annotation_expr,
    transform_methylation_level,
    trimer_from_heptamer,
)
from gnomad.utils.constraint.counts import (
    _build_sum_agg_struct,
    _resolve_row_annotation_expr,
    aggregate_constraint_metrics_expr,
    count_observed_and_possible_by_group,
    counts_agg_expr,
    variant_observed_and_possible_expr,
    variant_observed_expr,
    weighted_sum_agg_expr,
)
from gnomad.utils.constraint.grouping import (
    annotate_exploded_vep_for_constraint_groupings,
    build_constraint_consequence_groups,
    explode_downsamplings_oe,
    get_all_gen_anc_lengths,
    get_constraint_grouping_expr,
)
from gnomad.utils.constraint.metrics import (
    _oe_ci_discretized_poisson,
    _oe_ci_gamma,
    add_gencode_transcript_annotations,
    calculate_raw_z_score,
    calculate_raw_z_score_sd,
    compute_pli,
    get_constraint_flags,
    oe_confidence_interval,
)
from gnomad.utils.constraint.models import (
    apply_models,
    apply_plateau_models,
    build_coverage_model,
    build_models,
    build_plateau_models,
    calibration_model_group_expr,
    coverage_correction_expr,
)
from gnomad.utils.constraint.ranking import (
    annotate_bins_by_threshold,
    compute_percentile_thresholds,
    rank_and_assign_bins,
    rank_array_element_metrics,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)

__all__ = [
    "COVERAGE_CUTOFF",
    "CLASSIC_LOF_ANNOTATIONS",
    "DEFAULT_FIELDS_TO_SUM",
    "DEFAULT_GENCODE_ANNOTATIONS",
    "get_mu_annotation_expr",
    "annotate_with_mu",
    "annotate_mutation_type",
    "trimer_from_heptamer",
    "collapse_strand",
    "transform_methylation_level",
    "assemble_constraint_context_ht",
    "calculate_gerp_cutoffs",
    "variant_observed_expr",
    "variant_observed_and_possible_expr",
    "weighted_sum_agg_expr",
    "counts_agg_expr",
    "count_observed_and_possible_by_group",
    "aggregate_constraint_metrics_expr",
    "build_constraint_consequence_groups",
    "explode_downsamplings_oe",
    "get_all_gen_anc_lengths",
    "get_constraint_grouping_expr",
    "annotate_exploded_vep_for_constraint_groupings",
    "calibration_model_group_expr",
    "build_models",
    "build_plateau_models",
    "build_coverage_model",
    "apply_plateau_models",
    "coverage_correction_expr",
    "apply_models",
    "compute_pli",
    "oe_confidence_interval",
    "calculate_raw_z_score",
    "get_constraint_flags",
    "calculate_raw_z_score_sd",
    "add_gencode_transcript_annotations",
    "rank_and_assign_bins",
    "compute_percentile_thresholds",
    "annotate_bins_by_threshold",
    "rank_array_element_metrics",
]
