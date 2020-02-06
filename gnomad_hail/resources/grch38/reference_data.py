from gnomad_hail.resources.resource_utils import (
    TableResource,
    VersionedTableResource,
    MatrixTableResource,
    VersionedMatrixTableResource,
    import_sites_vcf,
    NO_CHR_TO_CHR_CONTIG_RECODING
)
import hail as hl

from hail import Table


def _import_purcell_5k(path) -> hl.Table:
    p5k = hl.import_locus_intervals(path, reference_genome='GRCh37')
    p5k = p5k.annotate(
        start=hl.liftover(p5k.interval.start, 'GRCh38'),
        end=hl.liftover(p5k.interval.start, 'GRCh38')
    )
    p5k = p5k.filter(
        (p5k.start.contig == 'chr' + p5k.interval.start.contig) &
        (p5k.end.contig == 'chr' + p5k.interval.end.contig)
    )
    p5k = p5k.key_by()
    p5k = p5k.select(
        locus=p5k.start,
        locus_b37=p5k.interval.start
    )
    return p5k.key_by('locus')

# Resources with no versioning needed
purcell_5k = TableResource(
    path="gs://gnomad-public/resources/grch38/purcell5k.ht",
    import_func=_import_purcell_5k,
    import_args={
        "path": "gs://gnomad-public/intervals/purcell5k.interval_list",
    }
)

na12878_giab = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt",
    import_func=hl.import_vcf,
    import_args={
        "path": "gs://gnomad-public/truth-sets/source/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
        "force_bgz": True,
        "min_partitions": 100
    }
)

na12878_giab_hc_intervals = TableResource(
    path='gs://gnomad-public/truth-sets/source/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7_hc_regions.ht',
    import_func=hl.import_bed,
    import_args={
        "path": 'gs://gnomad-public/truth-sets/source/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed',
        "reference_genome": 'GRCh38',
        "skip_invalid_intervals": True
    }
)

syndip = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/gnomad_v3_syndip.b38.mt"
)

syndip_hc_intervals = TableResource(
    path='gs://gnomad-public/truth-sets/source/syndip_b38_hc_regions.ht',
    import_func=hl.import_bed,
    import_args={
        "path": 'gs://gnomad-public/truth-sets/source/syndip.b38.bed',
        "reference_genome": 'GRCh38',
        "skip_invalid_intervals": True,
        "min_partitions": 10
    }
)

# Versioned resources: versions should be listed from most recent to oldest
clinvar = VersionedTableResource(
    default_version="20190923",
    versions={
        "20190923": TableResource(
            path="gs://gnomad-public/resources/grch38/clinvar_20190923.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/clinvar_20190923.vcf.gz",
                "force_bgz": True,
                "contig_recoding": NO_CHR_TO_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 100
            },
        )
    }
)

dbsnp = VersionedTableResource(
    default_version="b151",
    versions={
        "b151": TableResource(
            path="gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.vcf.bgz",
                "header_file": "gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.vcf.header",
                "force_bgz": True,
                "contig_recoding": NO_CHR_TO_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 400
            },
        )
    }
)

hapmap = TableResource(
    path="gs://gnomad-public/resources/grch38/hapmap_3.3.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz",
        "force_bgz": True
    }
)

kgp_omni = TableResource(
    path="gs://gnomad-public/resources/grch38/1000G_omni2.5.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
        "force_bgz": True
    }
)

kgp = VersionedTableResource(
    default_version="phase_1_hc",
    versions={
        "phase_1_hc": TableResource(
            path="gs://gnomad-public/resources/grch38/1000G_phase1.snps.high_confidence.hg38.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "force_bgz": True
            }
        )
    }
)

mills = TableResource(
    path="gs://gnomad-public/resources/grch38/Mills_and_1000G_gold_standard.indels.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "force_bgz": True
    }
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
