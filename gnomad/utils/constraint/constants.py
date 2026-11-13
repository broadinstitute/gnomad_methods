"""Constants shared by the constraint functions."""

COVERAGE_CUTOFF = 30
"""
Minimum median exome coverage differentiating high coverage sites from low coverage sites.

Low coverage sites require an extra calibration when computing the proportion of expected variation.
"""


CLASSIC_LOF_ANNOTATIONS = (
    "stop_gained",
    "splice_donor_variant",
    "splice_acceptor_variant",
)
"""Classic loss-of-function VEP annotations."""


DEFAULT_FIELDS_TO_SUM = (
    "mu_snp",
    "mu",
    "observed_variants",
    "possible_variants",
    "predicted_proportion_observed",
    "coverage_correction",
    "expected_variants",
)
"""Default fields summed when aggregating expected variants by constraint group."""


DEFAULT_GENCODE_ANNOTATIONS = (
    "transcript_id_version",
    "gene_id_version",
    "level",
    "transcript_type",
    "start_position",
    "end_position",
)
"""Default GENCODE annotations added to a Table by transcript id."""
