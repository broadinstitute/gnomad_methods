import hail as hl

CURRENT_HAIL_VERSION = "0.2"
CURRENT_RELEASE = "2.0.2"
CURRENT_GENOME_META = "2017-06-02"  # YYYY-MM-DD
CURRENT_EXOME_META = "2017-06-02"

RELEASES = ["2.0.1", "2.0.2"]

GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
EXOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]


def public_exomes_vds_path(split=True, version=CURRENT_RELEASE):
    return 'gs://gnomad-public/release/{0}/vds/exomes/gnomad.exomes.r{0}.sites{1}.vds'.format(version, "" if split else ".unsplit")


def public_genomes_vds_path(split=True, version=CURRENT_RELEASE):
    return 'gs://gnomad-public/release/{0}/vds/genomes/gnomad.genomes.r{0}.sites{1}.vds'.format(version, "" if split else ".unsplit")


def get_gnomad_public_data(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get public gnomAD data as VDS.

    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Chosen VDS
    :rtype: MatrixTable
    """
    return hl.methods.read_matrix_table(get_gnomad_public_data_path(data_type, split=split, version=version))


def get_gnomad_data(data_type, adj=False, split=True, raw=False, hail_version=CURRENT_HAIL_VERSION,
                    meta_version=None, meta_root='meta', fam_root='fam', duplicate_mapping_root=None,
                    release_samples=False, release_annotations=None):
    """
    Wrapper function to get gnomAD data as VDS.

    :param str data_type: One of `exomes` or `genomes`
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool raw: Whether to return the raw (10T+) data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param str hail_version: One of the HAIL_VERSIONs
    :param str meta_version: Version of metadata (None for current)
    :param str meta_root: Where to put metadata. Set to None if no metadata is desired.
    :param str fam_root: Where to put the pedigree information. Set to None if no pedigree information is desired.
    :param str duplicate_mapping_root: Where to put the duplicate genome/exome samples ID mapping (default is None -- do not annotate)
    :param bool release_samples: When set, filters the data to release samples only
    :param str release_annotations: One of the RELEASES to add variant annotations (into va), or None for no data
    :return: gnomAD hardcalls dataset with chosen annotations
    :rtype: MatrixTable
    """
    if raw and split:
        raise DataException('No split raw data. Use of hardcalls is recommended.')

    vds = hl.methods.read_matrix_table(get_gnomad_data_path(data_type, hardcalls=not raw, split=split, hail_version=hail_version))
    if adj:
        vds = filter_to_adj(vds)

    if meta_root:
        meta_kt = get_gnomad_meta(data_type, meta_version)
        vds = vds.annotate_cols(**{meta_root: meta_kt[vds.s]})

    if duplicate_mapping_root:
        dup_kt = hl.methods.import_table(genomes_exomes_duplicate_ids_tsv_path, impute=True,
                                         key='exome_id' if data_type == "exomes" else 'genome_id')
        vds = vds.annotate_cols(**{duplicate_mapping_root: dup_kt[vds.s]})

    if fam_root:
        fam_kt = hl.methods.import_fam(exomes_fam_path if data_type == "exomes" else genomes_fam_path)
        vds = vds.annotate_cols(**{fam_root: fam_kt[vds.s]})

    pops = EXOME_POPS if data_type == 'exomes' else GENOME_POPS
    vds = vds.annotate_globals(pops=map(lambda x: x.lower(), pops))

    if release_samples:
        vds = vds.filter_cols(vds.meta.release)

    if release_annotations:
        sites_vds = get_gnomad_public_data(data_type, split, release_annotations)
        vds = vds.select_rows(release=sites_vds[vds.v, :])  # TODO: replace with ** to nuke old annotations

    return vds


def get_gnomad_meta(data_type, version=None):
    """
    Wrapper function to get gnomAD metadata as Table

    :param str data_type: One of `exomes` or `genomes`
    :param str version: Metadata version (None for current)
    :return: Metadata Table
    :rtype: Table
    """
    meta_kt = hl.methods.import_table(get_gnomad_meta_path(data_type, version), impute=True,
                                      key="sample" if data_type == "exomes" else "Sample")

    return meta_kt.annotate(
        release=meta_kt.drop_status == "keep" if data_type == 'exomes' else meta_kt.keep,  # unify_sample_qc: this is version dependent will need fixing when new metadata arrives
        population=meta_kt.population if data_type == 'exomes' else hl.functions.cond(meta_kt.final_pop == 'sas', 'oth', meta_kt.final_pop)
    )


def get_gnomad_public_data_path(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get paths to gnomAD data

    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Path to chosen VDS
    :rtype: str
    """
    if version not in RELEASES:
        return DataException("Select version as one of: {}".format(','.join(RELEASES)))

    if data_type == 'exomes':
        return public_exomes_vds_path(split, version)
    elif data_type == 'genomes':
        return public_genomes_vds_path(split, version)
    return DataException("Select data_type as one of 'genomes' or 'exomes'")


def get_gnomad_data_path(data_type, hardcalls=False, split=True, hail_version=CURRENT_HAIL_VERSION):
    """
    Wrapper function to get paths to gnomAD data

    :param str data_type: One of `exomes` or `genomes`
    :param bool hardcalls: Whether hardcalls should be returned
    :param bool split: Whether the dataset should be split (only applies to hardcalls)
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to chosen VDS
    :rtype: str
    """
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")
    if hardcalls:
        return hardcalls_vds_path(data_type, split, hail_version)
    else:
        return raw_exomes_vds_path(hail_version) if data_type == 'exomes' else raw_genomes_vds_path(hail_version)


def get_gnomad_meta_path(data_type, version=None):
    """
    Wrapper function to get paths to gnomAD metadata

    :param str data_type: One of `exomes` or `genomes`
    :param str version: String with version (date) for metadata
    :return: Path to chosen metadata file
    :rtype: str
    """
    if data_type == 'exomes':
        if version:
            return metadata_exomes_tsv_path(version)
        return metadata_exomes_tsv_path()
    elif data_type == 'genomes':
        if version:
            return metadata_genomes_tsv_path(version)
        return metadata_genomes_tsv_path()
    return DataException("Select data_type as one of 'genomes' or 'exomes'")


def raw_exomes_vds_path(hail_version=CURRENT_HAIL_VERSION):
    """
    Warning: unsplit and no special consideration on sex chromosomes
    """
    return 'gs://gnomad/raw/hail-{0}/vds/exomes/gnomad.exomes.vds'.format(hail_version)


def raw_genomes_vds_path(hail_version=CURRENT_HAIL_VERSION):
    """
    Warning: unsplit and no special consideration on sex chromosomes
    """
    return 'gs://gnomad/raw/hail-{0}/vds/genomes/gnomad.genomes.vds'.format(hail_version)


def raw_exac_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/raw/hail-{0}/vds/exac/exac.vds'.format(hail_version)


def exac_release_sites_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/raw/hail-{}/vds/exac/exac.r1.sites.vep.vds'.format(hail_version)


def hardcalls_vds_path(data_type, split=True, hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/hardcalls/hail-{0}/vds/{1}/gnomad.{1}{2}.vds'.format(hail_version, data_type,
                                                                             "" if split else ".unsplit")


def annotations_vds_path(data_type, annotation_type, hail_version=CURRENT_HAIL_VERSION):
    """
    Get sites-level annotations

    :param str data_type: One of "exomes" or "genomes"
    :param str annotation_type: One of "vep", "qc_stats", "frequencies", "rf", "concordance"
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to annotations VDS
    :rtype: str
    """
    return 'gs://gnomad/annotations/hail-{0}/vds/{1}/gnomad.{1}.{2}.vds'.format(hail_version, data_type,
                                                                                annotation_type)

gnomad_pca_vds_path = "gs://gnomad-genomes/sampleqc/gnomad.pca.vds"


def gnomad_public_pca_vds_path(version=CURRENT_RELEASE):
    """
    Returns the path for the public gnomAD VDS containing sites and loadings from the PCA

    :param str version: One of the RELEASEs
    :return: path to gnomAD public PCA VDS
    :rtype: str
    """
    return "gs://gnomad-public/release/{}/pca/gnomad_pca_loadings.vds".format(version)


def metadata_genomes_tsv_path(version=CURRENT_GENOME_META):
    return 'gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.tsv.bgz'.format(version)


def metadata_exomes_tsv_path(version=CURRENT_EXOME_META):
    return 'gs://gnomad/metadata/exomes/gnomad.exomes.metadata.{0}.tsv.bgz'.format(version)


genomes_fam_path = "gs://gnomad/metadata/genomes/gnomad.genomes.fam"
exomes_fam_path = "gs://gnomad/metadata/exomes/gnomad.exomes.fam"
genomes_exomes_duplicate_ids_tsv_path = "gs://gnomad/metadata/genomes_exomes_duplicate_ids.tsv"


def omni_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_omni2.5.b37.vds'.format(hail_version)


def mills_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/Mills_and_1000G_gold_standard.indels.b37.vds'.format(hail_version)


def hapmap_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hapmap_3.3.b37.vds'.format(hail_version)


def kgp_high_conf_snvs_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_phase1.snps.high_confidence.b37.vds'.format(hail_version)


def NA12878_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vds'.format(hail_version)


def syndip_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hybrid.m37m.vds'.format(hail_version)


def cpg_sites_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/resources/hail-{}/cpg.vds'.format(hail_version)


def methylation_sites_vds_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-resources/methylation/hail-{}/methylation.kt'.format(hail_version)


dbsnp_vcf_path = "gs://gnomad-public/truth-sets/source/All_20160601.vcf.bgz"

NA12878_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed"
NA12878_high_conf_exome_regions_bed_path = "gs://gnomad-public/truth-sets/source/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
syndip_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/hybrid.m37m.bed"
clinvar_tsv_path = "gs://gnomad-resources/annotations/clinvar_alleles.single.b37.tsv.gz"
clinvar_vds_path = "gs://gnomad-resources/annotations/clinvar_alleles.single.b37.vds"

# Useful intervals
lcr_intervals_path = "gs://gnomad-public/intervals/LCR.interval_list"
decoy_intervals_path = "gs://gnomad-public/intervals/mm-2-merged.bed.gz"
purcell5k_intervals_path = "gs://gnomad-public/intervals/purcell5k.interval_list"
segdup_intervals_path = "gs://gnomad-public/intervals/hg19_self_chain_split_both.bed.gz"

# Exome intervals
exomes_high_conf_regions_intervals_path = "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list"
exome_calling_intervals_path = 'gs://gnomad-public/intervals/exome_calling_regions.v1.interval_list'
evaluation_intervals_path = 'gs://gnomad-public/intervals/exome_evaluation_regions.v1.noheader.interval_list'
high_coverage_intervals_path = 'gs://gnomad-public/intervals/high_coverage.auto.interval_list'

# Genome intervals
genome_evaluation_intervals_path = "gs://gnomad-public/intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list"  # from Broad GP
genome_evaluation_intervals_path_hg38 = "gs://gnomad-public/intervals/hg38-v0-wgs_evaluation_regions.hg38.interval_list"
# More can be found at gs://broad-references/hg19

vep_config = "/vep/vep-gcloud.properties"

# Annotations
context_vds_path = 'gs://gnomad-resources/constraint/context_processed.vds'


class DataException(Exception):
    pass
