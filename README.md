# Hail utilities for gnomAD

[![PyPI](https://img.shields.io/pypi/v/gnomad)](https://pypi.org/project/gnomad/)

This repo contains a number of [Hail](https://hail.is/) utility functions and scripts for the [gnomAD project](http://gnomad.broadinstitute.org) and the [Translational Genomics Group](https://the-tgg.org). As we continue to expand the size of our datasets, we are constantly seeking to find ways to reduce the complexity of our workflows and to make these functions more generic. As a result, the interface for many of these functions will change over time as we generalize their implementation for more flexible use within our scripts. We are also continously adapting our code to regular changes in the Hail interface. These repos thus represent only a snapshot of the gnomAD code base and are shared without guarantees or warranties.

We therefore encourage users to browse through the [API reference](https://broadinstitute.github.io/gnomad_methods/api_reference/) to identify modules and functions that will be useful in their own pipelines, and to edit and reconfigure relevant code to suit their particular analysis and QC needs.

## Repo layout

| Directory | Purpose |
|-----------|---------|
| `gnomad/utils/` | General-purpose utilities: annotations, filtering, Ensembl VEP, constraint, sparse MTs, liftover, VCF export, etc. One module per topic. |
| `gnomad/resources/` | Resource classes (`resource_utils.py`), resource source config (`config.py`), and per-build resource definitions (`grch37/`, `grch38/`). |
| `gnomad/sample_qc/` | Sample QC: genetic ancestry, relatedness, sex inference, platform inference, filtering. |
| `gnomad/variant_qc/` | Variant QC: random forest, training, evaluation, LD. |
| `gnomad/assessment/` | Release assessment: summary stats, validity checks. |
| `tests/` | Pytest suite, mirroring the `gnomad/` layout. |
| `docs/` | Sphinx docs; the API reference is auto-generated from docstrings. |

This is a subpackage-level map only. Individual functions are not listed here
because they move and get renamed — use the
[API reference](https://broadinstitute.github.io/gnomad_methods/api_reference/)
or search the source.

## Related gnomAD repositories

- [gnomad_qc](https://github.com/broadinstitute/gnomad_qc) — the gnomAD
  production QC pipelines. The largest consumer of this library, and the best
  place to see these functions used on real data.
- [gnomad-constraint](https://github.com/broadinstitute/gnomad-constraint) —
  the gnomAD constraint pipeline.
- [gnomad-toolbox](https://github.com/broadinstitute/gnomad-toolbox) —
  higher-level helpers for working with released gnomAD Hail Tables.
- [gnomad-browser](https://github.com/broadinstitute/gnomad-browser) — the
  gnomAD browser.
