# gnomAD Methylation Sites Data

This directory contains methylation site annotations for the GRCh38 reference genome, used in gnomAD constraint calculations and variant annotation pipelines.

## Files Description

### BED Files
- **methylation.bed**: Methylation site annotations for autosomes (chr1-22)
- **methylation_chrX.bed**: Methylation site annotations for chromosome X

### Hail Table Files
- **methylation.ht/**: Hail Table containing methylation annotations for autosomes only
- **methylation_chrX.ht/**: Hail Table containing methylation annotations for chromosome X only
- **methylation_all.ht/**: Merged Hail Table containing methylation annotations for all chromosomes (autosomes, chrX, and chrY)

## Recommended Usage

**We recommend using `methylation_all.ht/` for all applications** as it provides comprehensive methylation annotations across all chromosomes in a single file. This eliminates the need to handle multiple files and ensures consistent annotation coverage across the entire genome.

The individual chromosome-specific files (`methylation.ht/` and `methylation_chrX.ht/`) are legacy files that were created when complete methylation data was not available for all chromosomes. These files are maintained for backward compatibility but should not be used for new analyses.

## Methylation Score Scales

The methylation scores use different scales depending on the genomic region:

- **Autosomes (chr1-22)**: 0-15 scale
- **chrX PAR regions**: 0-15 scale (same as autosomes)
- **chrX and chrY non-PAR regions**: 0-12 scale


## Data Format

### BED Files
Standard BED format with the following columns:
- Column 1: Chromosome
- Column 2: Start position (0-based)
- Column 3: End position (1-based)
- Column 4: Methylation score

### Hail Tables
The Hail Tables contain the same methylation score information in a format optimized for large-scale genomic analyses. Key annotations include:
- `locus`: Genomic position
- `methylation_level`: Methylation score (0-15 for autosomes/PAR, 0-12 for non-PAR)

## Citation

If you use these methylation annotations in your research, please cite:

Chen, S., et al. "A genomic mutational constraint map using variation in 76,156 human genomes." bioRxiv (2022). https://www.biorxiv.org/content/10.1101/2022.03.20.485034v2.full
