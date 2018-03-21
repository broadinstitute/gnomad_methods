

# Sample QC files
def qc_mt_path(data_type: str):
    return 'gs://gnomad/sample_qc/mt/gnomad.{}.high_callrate_common_biallelic_snps.mt'.format(data_type)


def qc_ht_path(data_type: str):
    return 'gs://gnomad/sample_qc/ht/gnomad.{}.high_callrate_common_biallelic_snps.ht'.format(data_type)


def qc_temp_data_prefix(data_type: str):
    # can be joint
    return 'gs://gnomad/sample_qc/temp/{0}/gnomad.{0}'.format(data_type)


def qc_meta_path(data_type: str):
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.metadata.2018-01-31.reformatted_colnames.tsv.txt.gz'
    else:
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.metadata.2018-01-31.reformatted_colnames.tsv.txt.gz'

exome_extra_meta_path = 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.metadata_import_table.2018-01-31.one_hot_encoded.txt'
callrate_scores_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.exomes.callrate_pca_scores.ht'
exome_platform_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.exomes.platforms.ht'

qc_joined_mt_path = 'gs://gnomad/sample_qc/gnomad.joint_genomes_exomes.high_callrate_common_biallelic_snps.mt'
qc_joined_pruned_mt_path = 'gs://gnomad/sample_qc/gnomad.joint_genomes_exomes.high_callrate_common_biallelic_snps.pruned.mt'

exome_callrate_mt_path = 'gs://gnomad/sample_qc/mt/gnomad.exomes.callrate.mt'