# Changes

## Unreleased

* Added `median_impute_features` to variant QC random forest module [(224)](https://github.com/broadinstitute/gnomad_methods/pull/224)
* Created `training.py` in variant QC and added `sample_training_examples` [(224)](https://github.com/broadinstitute/gnomad_methods/pull/224)
* Added variant QC pipeline function `train_rf_model` [(224)](https://github.com/broadinstitute/gnomad_methods/pull/224)

## Version 0.3.0 - April 28th, 2020

### Changed

* Updated capitalization of ambiguous sex annotation [(#208)](https://github.com/broadinstitute/gnomad_methods/pull/208)
* Updated usage of included intervals in imputing sex ploidy, also updated interval parameter names [(#209)](https://github.com/broadinstitute/gnomad_methods/pull/209)
* Updated capitalization in relatedness constants [(#217)](https://github.com/broadinstitute/gnomad_methods/pull/217)
* Changed interface for Slack notifications [(#219)](https://github.com/broadinstitute/gnomad_methods/pull/219)

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
