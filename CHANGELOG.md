# Changes

## Unreleased

## Version 0.2.0 - April 3rd, 2020

### Added 

* Function to subset a `MatrixTable` based on a list of samples [(#196)](https://github.com/broadinstitute/gnomad_methods/pull/196)
* Function to get file size and MD5 hash [(#186)](https://github.com/broadinstitute/gnomad_methods/pull/186)
* Developer documentation [(#185)](https://github.com/broadinstitute/gnomad_methods/pull/185)
* Include `RAW_MQ` and `AS_VQSLOD` metrics in `get_annotations_hists` [(#181)](https://github.com/broadinstitute/gnomad_methods/pull/181)
* Functions to compute coverage stats from sparse MT [(#173)](https://github.com/broadinstitute/gnomad_methods/pull/173)

### Changed

* Repo restructured - imports may need to be updated [(#207)] (https://github.com/broadinstitute/gnomad_methods/pull/207)
* Make some arguments to `get_qc_mt` optional [(#200)](https://github.com/broadinstitute/gnomad_methods/pull/200)
* Fetch VEP configuration from new Hail requestor pays buckets [(#197)](https://github.com/broadinstitute/gnomad_methods/pull/197)
* Hail must be installed separately [(#194)](https://github.com/broadinstitute/gnomad_methods/pull/194)

### Fixed

* Father/mother assignments now correct (were swapped before) `infer_families` [(#203)](https://github.com/broadinstitute/gnomad_methods/pull/203)
* Attribute assignments for `VersionedPedigreeResource` [(#198)](https://github.com/broadinstitute/gnomad_methods/pull/198)
* Field references in `get_annotations_hists` [(#181)](https://github.com/broadinstitute/gnomad_methods/pull/181)
* Use before assignment error in `default_compute_info` [(#195)](https://github.com/broadinstitute/gnomad_methods/pull/195)

## Version 0.1.0 - March 4th, 2020

Initial release
