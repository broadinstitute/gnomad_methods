from .basics import DataException


def qc_mt_path(data_type: str, ld_pruned: bool = False) -> str:
    """
    Returns MatrixTable for sample QC purposes: can be exomes, genomes, or joint (joint dataset can also be ld_pruned=True)
    Criteria: callrate > 0.99, AF > 0.001, SNPs only, bi-allelics only
    """
    if data_type not in ('exomes', 'genomes', 'joint'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes' or 'joint'")
    if ld_pruned and data_type != 'joint':
        raise DataException("ld_pruned = True is only available for 'joint'")
    ld_pruned = '.pruned' if ld_pruned else ''
    return f'gs://gnomad/sample_qc/mt/gnomad.{data_type}.high_callrate_common_biallelic_snps{ld_pruned}.mt'


def qc_ht_path(data_type: str, part: str) -> str:
    """
    Interim sample metadata tables generated in the sample qc process.
    Generally not to be used: use tables from basics.py instead (e.g. metadata_*omes_ht_path)
    hard_filters (contains hard filter, permission information, and sex)
    platforms (contains imputed platform information for exomes only)
    pop_platform (contains final related information and population/platform-specific QC filters)
    """
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")
    if part not in ('hard_filters', 'platforms', 'pop_platform'):
        raise DataException("Select part as one of 'hard_filters', 'platforms', or 'pop_platform'")
    if data_type == 'genomes' and part == 'platforms':
        raise DataException("'platforms' only available for 'genomes'")
    return f'gs://gnomad/sample_qc/ht/gnomad.{data_type}.{part}.ht'


def rank_annotations_path(data_type: str) -> str:
    """
    Path to annotation data for ranking samples for related pruning. The 'joint' dataset is the results of the ranking.
    """
    if data_type not in ('exomes', 'genomes', 'joint'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes' or 'joint'")
    return f'gs://gnomad/sample_qc/tsv/gnomad.{data_type}.rank_list_annotations.txt.bgz'


def qc_temp_data_prefix(data_type: str) -> str:
    """
    Path to directory with intermediate files for sample QC
    """
    if data_type not in ('exomes', 'genomes', 'joint'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes' or 'joint'")
    return f'gs://gnomad/sample_qc/temp/{data_type}/gnomad.{data_type}'


def qc_meta_path(data_type: str) -> str:
    """
    Input metadata file for sample QC
    """
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.streamlined_metadata.2018-03-21.txt.bgz'
    elif data_type == 'genomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-03-21.txt.bgz'
    else:
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")


exome_callrate_scores_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.exomes.callrate_pca_scores.ht'
exome_callrate_mt_path = 'gs://gnomad/sample_qc/mt/gnomad.exomes.callrate.mt'

relatedness_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.joint.relatedness.ht'

ancestry_pca_scores_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.joint.unrelated.pca_scores.ht'
ancestry_pca_loadings_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.joint.unrelated.pca_loadings.ht'

known_population_annotations = 'gs://gnomad/sample_qc/input_meta/gnomad.pop_annots.txt'
estonian_batches = 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.estonian_samples.txt'
