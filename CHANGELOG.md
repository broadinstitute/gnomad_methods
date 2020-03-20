# Changes

## Unreleased

### Added 

* Function to subset a `MatrixTable` based on a list of samples (#196)
* Function to get file size and MD5 hash (#186)
* Developer documentation (#185)
* Include `RAW_MQ` and `AS_VQSLOD` metrics in `get_annotations_hists` (#181)
* Functions to compute coverage stats from sparse MT (#173)

### Changed

* Make some arguments to `get_qc_mt` optional (#200)
* Fetch VEP configuration from new Hail requestor pays buckets (#197)
* Hail must be installed separately (#194)

### Fixed

* Attribute assignments for `VersionedPedigreeResource` (#198)
* Field references in `get_annotations_hists` (#181)
* Use before assignment error in `default_compute_info` (#195)

## Version 0.1.0 - March 4th, 2020

Initial release
