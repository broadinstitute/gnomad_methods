from gnomad_hail.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
    MatrixTableResource,
    VersionedMatrixTableResource,
)

from hail import Table


# Resources with no versioning needed
purcell_5k = TableResource(
    path="gs://gnomad-public/resources/grch38/purcell5k.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/purcell5k.interval_list",
    },
)

na12878_giab = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
    },
)

# Versioned resources: versions should be listed from most recent to oldest
clinvar = VersionedTableResource(
    default_version="20190923",
    versions={
        "20190923": TableResource(
            path="gs://gnomad-public/resources/grch38/clinvar_20190923.ht",
            import_sources={
                "source_path": "gs://gnomad-public/resources/grch38/clinvar_20190923.vcf.gz",
            },
        )
    },
)

dbsnp = VersionedTableResource(
    default_version="b151",
    versions={
        "b151": TableResource(
            path="gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.ht",
            import_sources={
                "source_path": "gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.vcf.bgz",
                "vcf_header_path": "gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.vcf.header",
            },
        )
    },
)

hapmap = TableResource(
    path="gs://gnomad-public/resources/grch38/hapmap_3.3.hg38.ht",
    import_sources={
        "source_path": "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz",
    },
)

kgp_omni = TableResource(
    path="gs://gnomad-public/resources/grch38/1000G_omni2.5.hg38.ht",
    import_sources={
        "source_path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
    },
)

kgp = VersionedTableResource(
    default_version="phase_1_hc",
    versions={
        "phase_1_hc": TableResource(
            path="gs://gnomad-public/resources/grch38/1000G_phase1.snps.high_confidence.hg38.ht",
            import_sources={
                "source_path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            },
        )
    },
)

mills = TableResource(
    path="gs://gnomad-public/resources/grch38/Mills_and_1000G_gold_standard.indels.hg38.ht",
    import_sources={
        "source_path": "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    },
)


def get_truth_ht() -> Table:
    """
    Returns a table with the following annotations from the latest version of the corresponding truth data:
    - hapmap
    - kgp_omni (1000 Genomes intersection Onni 2.5M array)
    - kgp_phase_1_hc (high confidence sites in 1000 genonmes)
    - mills (Mills & Devine indels)

    :return: A table with the latest version of popular truth data annotations
    """

    return hapmap.ht().select(hapmap=True).join(
        kgp_omni.ht().select(omni=True), how="outer"
    ).join(
        kgp.versions['phase_1_hc'].ht().select(kgp_phase1_hc=True), how="outer"
    ).join(
        mills.ht().select(mills=True), how="outer"
    ).repartition(200, shuffle=False).persist()
