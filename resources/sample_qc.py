

# Sample QC files
def qc_mt_path(data_type: str, pruned: bool = False):
    # can be joint, pruned only for joint
    return 'gs://gnomad/sample_qc/mt/gnomad.{}.high_callrate_common_biallelic_snps.{}mt'.format(data_type,
                                                                                                '.pruned' if pruned else '')


def qc_ht_path(data_type: str, part: str):
    # part can be hard_filters, platforms, or pop_platform
    return 'gs://gnomad/sample_qc/ht/gnomad.{}.{}.ht'.format(data_type, part)


def rank_annotations_path(data_type: str):
    # can be joint
    return 'gs://gnomad/sample_qc/tsv/gnomad.{0}.rank_list_annotations.txt.bgz'.format(data_type)


def qc_temp_data_prefix(data_type: str):
    # can be joint
    return 'gs://gnomad/sample_qc/temp/{0}/gnomad.{0}'.format(data_type)


def qc_meta_path(data_type: str):
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-03-21.txt.bgz'
    else:
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.streamlined_metadata.2018-03-21.txt.bgz'

callrate_scores_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.exomes.callrate_pca_scores.ht'

exome_callrate_mt_path = 'gs://gnomad/sample_qc/mt/gnomad.exomes.callrate.mt'

estonian_batches = 'gs://gnomad/projects/unify_sample_qc/gnomad.genomes.estonian_samples.txt'