# Changes

## Unreleased
* Added function `region_flag_expr` to flag problematic regions [(#349)](https://github.com/broadinstitute/gnomad_methods/pull/349/files)
* Added function `missing_callstats_expr` to create a Hail Struct with missing values that is inserted into frequency annotation arrays when data is missing [(#349)](https://github.com/broadinstitute/gnomad_methods/pull/349/files)
* Added function `set_female_y_metrics_to_na_expr` to set Y-variant frequency callstats for female-specific metrics to missing [(#349)](https://github.com/broadinstitute/gnomad_methods/pull/349/files)
* Added function `make_faf_index_dict` to create a look-up Dictionary for entries contained in the filter allele frequency annotation array [(#349)](https://github.com/broadinstitute/gnomad_methods/pull/349/files)
* Added function `make_freq_index_dict` to create a look-up Dictionary for entries contained in the frequency annotation array [(#349)](https://github.com/broadinstitute/gnomad_methods/pull/349/files)
* VersionedResource objects are no longer subclasses of BaseResource [(#359)](https://github.com/broadinstitute/gnomad_methods/pull/359)
* gnomAD resources can now be imported from different sources [(#373)](https://github.com/broadinstitute/gnomad_methods/pull/373)
* Replaced `ht_to_vcf_mt` with `adjust_vcf_incompatible_types` which maintains all functionality except turning the ht into a mt because it is no longer needed for use of the Hail module `export_vcf` [(#365)](https://github.com/broadinstitute/gnomad_methods/pull/365/files)
* Added function `remove_fields_from_constant` to remove fields from a list and notify which requested fields to remove were missing [(#381)](https://github.com/broadinstitute/gnomad_methods/pull/381)
* Added function `create_label_groups` to generate a list of label group dictionaries needed to populate the info dictionary for vcf export [(#381)](https://github.com/broadinstitute/gnomad_methods/pull/381)
* Added function `build_vcf_export_reference` to create a subset reference based on an existing reference genome [(#381)](https://github.com/broadinstitute/gnomad_methods/pull/381)
* Added function `rekey_new_reference` to re-key a Table or MatrixTable with a new reference genome [(#381)](https://github.com/broadinstitute/gnomad_methods/pull/381)
* Modified `SEXES` in utils/vcf to be 'XX' and 'XY' instead of 'female' and 'male' [(#381)](https://github.com/broadinstitute/gnomad_methods/pull/381)
* Fix `annotation_type_is_numeric` and `annotation_type_in_vcf_info` [(#379)](https://github.com/broadinstitute/gnomad_methods/pull/379)
* Added function `parallel_file_exists` to check whether a large number of files exist [(#394)](https://github.com/broadinstitute/gnomad_methods/pull/394)
* Changed module `sanity_checks` to `validity_checks`, modified functions `generic_field_check`, `make_filters_expr_dict` (previously `make_filters_sanity_check_expr`), and `make_group_sum_expr_dict` (previously `sample_sum_check`), and added functions `summarize_variant_filters`, `generic_field_check_loop`, `compare_subset_freqs`, `sum_group_callstats`, `summarize_variants`, `check_raw_and_adj_callstats`, `check_sex_chr_metrics`, `compute_missingness`, `vcf_field_check`, and `validate_release_t` [(#395)](https://github.com/broadinstitute/gnomad_methods/pull/389)

## Version 0.5.0 - April 22nd, 2021

### Fixed

* Fix for error in `generate_trio_stats_expr` that led to an incorrect untransmitted count. [(#238)](https://github.com/broadinstitute/gnomad_methods/pull/238)
* Fix for error in `compute_quantile_bin` that caused incorrect binning when a single score overlapped multiple bins [(#238)](https://github.com/broadinstitute/gnomad_methods/pull/238)
* Fixed `create_binned_ht` because it produced a "Cannot combine expressions from different source objects error" [(#238)](https://github.com/broadinstitute/gnomad_methods/pull/238)
* Fixed handling of missing entries (not within a ref block / alt site) when computing `coverage_stats` in `sparse_mt.py` [[#242]](https://github.com/broadinstitute/gnomad_methods/pull/242)
* Fix for error in `compute_stratified_sample_qc` where `gt_expr` caused error [(#259)](https://github.com/broadinstitute/gnomad_methods/pull/259)
* Fix for error in `default_lift_data` caused by missing `results` field in `new_locus` [(#270)](https://github.com/broadinstitute/gnomad_methods/pull/270)
* Fix to dbSNP b154 resource (resources.grch38.reference_data) import to allow for multiple rsIDs per variant [(#345)](https://github.com/broadinstitute/gnomad_methods/pull/345)
* Fix to `set_female_metrics_to_na` to correctly update chrY metrics to be missing [(#347)](https://github.com/broadinstitute/gnomad_methods/pull/347)
* Fixed available versions for gnomAD v2 `coverage` and `liftover` resources [(#352)](https://github.com/broadinstitute/gnomad_methods/pull/352)
* Removed side effect of accessing gnomAD v2 `coverage` and `liftover` exome resources that would edit available versions for other resources [(#352)](https://github.com/broadinstitute/gnomad_methods/pull/352)
* Use `overwrite` argument for importing a BlockMatrixResource [(#342)](https://github.com/broadinstitute/gnomad_methods/pull/342)

### Changed

* Removed assumption of `snv` annotation from `compute_quantile_bin`. [(#238)](https://github.com/broadinstitute/gnomad_methods/pull/238)
* Modified `compute_binned_truth_sample_concordance` to handle additional binning for subsets of variants. [(#240)](https://github.com/broadinstitute/gnomad_methods/pull/240)
* Updated liftover functions to be more generic [(#246)](https://github.com/broadinstitute/gnomad_methods/pull/246)
* Changed quality histograms to label histograms calculated on raw and not adj data [(#247)](https://github.com/broadinstitute/gnomad_methods/pull/247)
* Updated some VCF export constants [(#249)](https://github.com/broadinstitute/gnomad_methods/pull/249)
* Changed default DP threshold to 5 for hemi genotype calls in `annotate_adj` and `get_adj_expr` [(#252)](https://github.com/broadinstitute/gnomad_methods/pull/252)
* Updated coverage resources to version 3.0.1 [[#242]](https://github.com/broadinstitute/gnomad_methods/pull/242)
* Update to `compute_last_ref_block_end`, removing assumption that sparse MatrixTables are keyed only by `locus` by default [(#279)](https://github.com/broadinstitute/gnomad_methods/pull/279)
* Update `generic_field_check` to have option to show percentage of sites that fail checks. [(#284)](https://github.com/broadinstitute/gnomad_methods/pull/284)
* Modified `vep_or_lookup_vep` to support the use of different VEP versions [(#282)](https://github.com/broadinstitute/gnomad_methods/pull/282)
* Modified `create_truth_sample_ht` to add adj annotation information in the returned Table if present in the supplied MatrixTables [(#300)](https://github.com/broadinstitute/gnomad_methods/pull/300)

### Added

* Added constants and functions relevant to VCF export [(#241)](https://github.com/broadinstitute/gnomad_methods/pull/241)
* Add reference genome to call of `has_liftover` in `get_liftover_genome` [(#259)](https://github.com/broadinstitute/gnomad_methods/pull/259)
* Added fix for MQ calculation in `_get_info_agg_expr`, switched `RAW_MQ` and `MQ_DP` in calculation [(#262)](https://github.com/broadinstitute/gnomad_methods/pull/262)
* Add importable method for filtering clinvar to pathogenic sites [(#257)](https://github.com/broadinstitute/gnomad_methods/pull/257)
* Added common variant QC functions `get_rf_runs` and `get_run_data` to `random_forest.py` [(#278)](https://github.com/broadinstitute/gnomad_methods/pull/278)
* Add calculation for the strand odds ratio (SOR) to `get_site_info_expr` and `get_as_info_expr` [(#281)](https://github.com/broadinstitute/gnomad_methods/pull/281)
* Added VEPed context HT to resource files and included support for versioning [(#282)](https://github.com/broadinstitute/gnomad_methods/pull/282)
* Added code to generate summary statistics (total number of variants, number of LoF variants, LOFTEE summaries) [(#285)](https://github.com/broadinstitute/gnomad_methods/pull/285)
* Added additional counts to summary statistics (added autosome/sex chromosome counts, allele counts, counts for missense and synomymous variants) [(#289)](https://github.com/broadinstitute/gnomad_methods/pull/289)
* Added function, `default_generate_gene_lof_matrix`, to generate gene matrix [(#290)](https://github.com/broadinstitute/gnomad_methods/pull/290)
* Added function `default_generate_gene_lof_summary` to summarize gene matrix results [(#292)](https://github.com/broadinstitute/gnomad_methods/pull/292)
* Add resource for v3.1.1 release [(#364)](https://github.com/broadinstitute/gnomad_methods/pull/364)

### Removed

* Removed `rep_on_read`; this function is no longer necessary, as MatrixTables/Tables can be repartitioned on read with `_n_partitions` added by this [hail update](https://github.com/hail-is/hail/pull/9887) [(#283)](https://github.com/broadinstitute/gnomad_methods/pull/283)
* Removed `compute_quantile_bin` and added `compute_ranked_bin` as an alternative that provides more even binning. This is now used by `create_binned_ht` instead. [(#288)](https://github.com/broadinstitute/gnomad_methods/pull/288)
* Removed `prefix` parameter from  to `make_combo_header_text`, as this was only used to check if samples were from gnomAD [(#348)](https://github.com/broadinstitute/gnomad_methods/pull/348)

## Version 0.4.0 - July 9th, 2020

**Note** gnomAD resources have been moved to a [requester pays bucket](https://cloud.google.com/storage/docs/requester-pays).
Dataproc clusters must be [configured to allow reading from it](https://hail.is/docs/0.2/cloud/google_cloud.html#requester-pays).

* Added `VEP_CSQ_HEADER` to generate vep description necessary for VCF export. [(#230)](https://github.com/broadinstitute/gnomad_methods/pull/230)
* Modified variant QC pipeline functions `generate_trio_stats` and `generate_sib_stats` to add filter parameter for autosomes and bi-allelic sites [(#223)](https://github.com/broadinstitute/gnomad_methods/pull/223)
* `score_bin_agg` now requires additional annotations `ac` and `ac_qc_samples_unrelated_raw` and no longer needs `tdt` [(#223)](https://github.com/broadinstitute/gnomad_methods/pull/223)
* Changed `score_bin_agg` to use `ac_qc_samples_unrelated_raw` annotation instead of `unrelated_qc_callstats` [(#223)](https://github.com/broadinstitute/gnomad_methods/pull/223)
* Added singleton de novo counts to variant QC pipeline function `score_bin_agg` [(#223)](https://github.com/broadinstitute/gnomad_methods/pull/223)
* Modified `filter_mt_to_trios` to no longer filter to autosomes as this should be handled during the variant QC pipeline [(#223)](https://github.com/broadinstitute/gnomad_methods/pull/223)
* Updated `annotate_sex` to add globals to `sex_ht` [(#227)](https://github.com/broadinstitute/gnomad_methods/pull/227)
* Document `slack_notifications` function [(#228)](https://github.com/broadinstitute/gnomad_methods/pull/228)
* Added `median_impute_features` to variant QC random forest module [(224)](https://github.com/broadinstitute/gnomad_methods/pull/224)
* Created `training.py` in variant QC and added `sample_training_examples` [(224)](https://github.com/broadinstitute/gnomad_methods/pull/224)
* Added variant QC pipeline function `train_rf_model` [(224)](https://github.com/broadinstitute/gnomad_methods/pull/224)
* Use local copy of VEP config instead of reading from bucket [(#231)](https://github.com/broadinstitute/gnomad_methods/pull/231)
* Updated gnomAD resources paths for hail tables to requester pays buckets [(#233)](https://github.com/broadinstitute/gnomad_methods/pull/233)

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
