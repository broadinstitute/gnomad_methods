# Constraint Field Descriptions

Descriptions of columns in the gnomAD v4.1.1 gene constraint metrics tsv. Descriptions also apply to rows in Hail Tables, where a "." in the field name indicates a struct.

## General Fields

- `gene`: Gene name
- `gene_id`: Ensembl gene ID
- `transcript`: Ensembl or RefSeq transcript ID (GENCODE v39)
- `transcript_version`: Ensembl or RefSeq transcript version
- `canonical`: Boolean indicator as to whether the transcript is the canonical transcript for the gene
- `mane_select`: Boolean indicator as to whether the transcript is the MANE Select transcript for the gene

## Loss-of-Function (High and Low Confidence)

- `lof_hc_lc.obs`: Number of observed high and low confidence predicted loss-of-function (pLoF) variants in transcript
- `lof_hc_lc.exp`: Number of expected high and low confidence pLoF variants in transcript
- `lof_hc_lc.possible`: Number of possible high and low confidence pLoF variants in transcript
- `lof_hc_lc.mu`: Mutation rate summed across all possible high and low confidence pLoF variants in transcript
- `lof_hc_lc.oe_ci.lower`: Lower bound of 90% confidence interval (CI) for observed to expected (`oe`) ratio for high and low confidence pLoF variants
- `lof_hc_lc.oe_ci.upper`: Upper bound of 90% confidence interval for `oe` ratio for high and low confidence pLoF variants (lower values indicate more constrained)
- `lof_hc_lc.oe_ci.upper_rank`: Transcript's rank of upper bound `oe` CI value compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.oe_ci.upper_bin_decile`: Decile bin of upper bound of 90% CI `oe` for given transcript (lower values indicate more constrained).  This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.oe_ci.upper_bin_sextile`: Sextile bin of upper bound of 90% CI `oe` for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof_hc_lc.z_score`: Z-score for high and low confidence pLoF variants in transcript. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained)

## Loss-of-Function (High Confidence Only)

- `lof.obs`: Number of observed high confidence predicted loss-of-function (pLoF) variants in transcript
- `lof.exp`: Number of expected high confidence pLoF variants in transcript
- `lof.possible`: Number of possible high confidence pLoF variants in transcript
- `lof.oe`: Observed over expected ratio for high confidence pLoF variants in transcript (`lof.obs` divided by `lof.exp`)
- `lof.mu`: Mutation rate summed across all possible high confidence pLoF variants in transcript
- `lof.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for high confidence pLoF variants
- `lof.oe_ci.upper`: LOEUF: upper bound of 90% confidence interval for `oe` ratio for high confidence pLoF variants (lower values indicate more constrained)
- `lof.oe_ci.upper_rank`: Transcript's rank of LOEUF value compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.oe_ci.upper_bin_decile`: Decile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.oe_ci.upper_bin_sextile`: Sextile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `lof.z_score`: Z-score for pLoF variants in transcript. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained)
- `pLI`: Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% `oe` pLoF ratio;  computed from high confidence pLoF gnomAD data)
- `pRec`: Probability that transcript falls into distribution of recessive genes (~71% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)
- `pNull`: Probability that transcript falls into distribution of unconstrained genes (~100% `oe` pLoF ratio; computed from high confidence pLoF gnomAD data)

## Missense

- `mis.obs`: Number of observed missense variants in transcript
- `mis.exp`: Number of expected missense variants in transcript
- `mis.possible`: Number of possible missense variants in transcript
- `mis.oe`: Observed over expected ratio for missense variants in transcript (`mis.obs` divided by `mis.exp`)
- `mis.mu`: Mutation rate summed across all possible missense variants in transcript
- `mis.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for missense variants
- `mis.oe_ci.upper`: Upper bound of 90% confidence interval for `oe` ratio for missense variants
- `mis.oe_ci.upper_rank`: Transcript's rank of upper bound `oe` CI value for missense variants compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `mis.oe_ci.upper_bin_decile`: Decile bin of upper bound `oe` CI for missense variants for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `mis.oe_ci.upper_bin_sextile`: Sextile bin of upper bound `oe` CI for missense variants for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `mis.z_score`: Z-score for missense variants in transcript. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained)

## Synonymous

- `syn.obs`: Number of observed synonymous variants in transcript
- `syn.exp`: Number of expected synonymous variants in transcript
- `syn.possible`: Number of possible synonymous variants in transcript
- `syn.oe`: Observed over expected ratio for synonymous variants in transcript (syn.obs divided by syn.exp)
- `syn.mu`: Mutation rate summed across all synonymous variants in transcript
- `syn.oe_ci.lower`: Lower bound of 90% confidence interval for `oe` ratio for synonymous variants
- `syn.oe_ci.upper`: Upper bound of 90% confidence interval for `oe` ratio for synonymous variants
- `syn.oe_ci.upper_rank`: Transcript's rank of upper bound `oe` CI value for synonymous variants compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `syn.oe_ci.upper_bin_decile`: Decile bin of upper bound `oe` CI for synonymous variants for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `syn.oe_ci.upper_bin_sextile`: Sextile bin of upper bound `oe` CI for synonymous variants for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript will be used instead if available.
- `syn.z_score`: Z-score for synonymous variants in transcript. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained). Extreme values of syn.z_score indicate likely data quality issues

## Constraint Flags

- `constraint_flags`: Reason transcript does not have constraint metrics. One of:
  - `no_variants`: Zero observed synonymous, missense, pLoF variants
  - `no_exp_lof`: Zero expected pLoF variants
  - `outlier_lof`: Number of pLoF variants is significantly different than expectation
  - `no_exp_mis`: Zero expected missense variants
  - `outlier_mis`: Number of missense variants is significantly different than expectation
  - `no_exp_syn`: Zero expected synonymous variants
  - `outlier_syn`: Number of synonymous variants is significantly different than expectation

## Additional Transcript Metadata

- `chromosome`: Chromosome on which the transcript is located
- `transcript_level`: Transcript level from [Gencode](https://www.gencodegenes.org/pages/data_format.html)
- `transcript_type`: Transcript biotype from [Gencode](https://www.gencodegenes.org/pages/biotypes.html)
- `cds_length`: Length of the coding sequences (CDS) in the transcript
- `num_coding_exons`: Number of coding exons in the transcript
