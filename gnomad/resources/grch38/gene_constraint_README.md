# Constraint Field Descriptions

Descriptions of columns in the gnomAD v4 constraint metrics tsv. Descriptions also apply to rows in Hail Tables, where a "." in the field name indicates a struct.

## General Fields

- `gene`: Gene name
- `transcript`: Ensembl or RefSeq transcript ID (GENCODE v39)
- `canonical`: Boolean indicator as to whether the transcript is the canonical transcript for the gene
- `mane_select`: Boolean indicator as to whether the transcript is the MANE Select transcript for the gene

## Loss-of-Function (High and Low Confidence)

- `lof_hc_lc.obs`: Number of observed high and low confidence loss-of-function (pLoF) variants in transcript
- `lof_hc_lc.exp`: Number of expected high and low confidence pLoF variants in transcript
- `lof_hc_lc.possible`: Number of possible high and low confidence pLoF variants in transcript
- `lof_hc_lc.oe`: Observed over expected ratio for high and low confidence pLoF variants in transcript (lof_hc_lc.obs divided by lof_hc_lc.exp)
- `lof_hc_lc.mu`: Mutation rate summed across all possible high and low confidence pLoF variants in transcript
- `lof_hc_lc.pLI`: Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% o/e pLoF ratio; computed from high and low confidence pLoF gnomAD data)
- `lof_hc_lc.pRec`: Probability that transcript falls into distribution of recessive genes (~71% o/e pLoF ratio; computed from high and low confidence pLoF gnomAD data)
- `lof_hc_lc.pNull`: Probability that transcript falls into distribution of unconstrained genes (~100% o/e pLoF ratio; computed from high and low confidence pLoF gnomAD data)

## Loss-of-Function (High Confidence Only)

- `lof.obs`: Number of observed high confidence loss-of-function (pLoF) variants in transcript
- `lof.exp`: Number of expected high confidence pLoF variants in transcript
- `lof.possible`: Number of possible high confidence pLoF variants in transcript
- `lof.oe`: Observed over expected ratio for high confidence pLoF variants in transcript (lof.obs divided by lof.exp)
- `lof.mu`: Mutation rate summed across all possible high confidence pLoF variants in transcript
- `lof.pLI`: Probability of loss-of-function intolerance; probability that transcript falls into distribution of haploinsufficient genes (~21% o/e pLoF ratio; computed from high confidence pLoF gnomAD data)
- `lof.pRec`: Probability that transcript falls into distribution of recessive genes (~71% o/e pLoF ratio; computed from high confidence pLoF gnomAD data)
- `lof.pNull`: Probability that transcript falls into distribution of unconstrained genes (~100% o/e pLoF ratio; computed from high confidence pLoF gnomAD data)
- `lof.oe_ci.lower`: Lower bound of 90% confidence interval for o/e ratio for high confidence pLoF variants
- `lof.oe_ci.upper`: LOEUF: upper bound of 90% confidence interval for o/e ratio for high confidence pLoF variants (lower values indicate more constrained)
- `lof.oe_ci.upper_rank`: Transcript's rank of LOEUF value compared to other transcripts (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript (if available) will be used instead.
- `lof.oe_ci.upper_bin_decile`: Decile bin of LOEUF for given transcript (lower values indicate more constrained). This annotation is only applied to MANE Select transcripts unless a gene does not have a MANE Select transcript, in which case the canonical transcript (if available) will be used instead.
- `lof.z_raw`: Raw Z-score for pLoF variants in transcript
- `lof.z_score`: Z-score for pLoF variants in transcript. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained)

## Missense

- `mis.obs`: Number of observed missense variants in transcript
- `mis.exp`: Number of expected missense variants in transcript
- `mis.possible`: Number of possible missense variants in transcript
- `mis.oe`: Observed over expected ratio for missense variants in transcript (mis.obs divided by mis.exp)
- `mis.mu`: Mutation rate summed across all possible missense variants in transcript
- `mis.oe_ci.lower`: Lower bound of 90% confidence interval for o/e ratio for missense variants
- `mis.oe_ci.upper`: Upper bound of 90% confidence interval for o/e ratio for missense variants
- `mis.z_raw`: Raw Z-score for missense variants in transcript
- `mis.z_score`: Z-score for missense variants in transcript. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained)

## Synonymous

- `syn.obs`: Number of observed synonymous variants in transcript
- `syn.exp`: Number of expected synonymous variants in transcript
- `syn.possible`: Number of possible synonymous variants in transcript
- `syn.oe`: Observed over expected ratio for synonymous variants in transcript (syn.obs divided by syn.exp)
- `syn.mu`: Mutation rate summed across all synonymous variants in transcript
- `syn.oe_ci.lower`: Lower bound of 90% confidence interval for o/e ratio for synonymous variants
- `syn.oe_ci.upper`: Upper bound of 90% confidence interval for o/e ratio for synonymous variants
- `syn.z_raw`: Raw Z-score for synonymous variants in transcript
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
- `level`: Transcript level from [Gencode](https://www.gencodegenes.org/pages/data_format.html)
- `transcript_type`: Transcript biotype from [Gencode](https://www.gencodegenes.org/pages/biotypes.html)
- `cds_length`: Length of the coding sequences (CDS) in the transcript
- `num_coding_exons`: Number of coding exons in the transcript
