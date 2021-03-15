from gnomad.resources.resource_utils import (
    DBSNP_B154_CHR_CONTIG_RECODING,
    TableResource,
    VersionedTableResource,
    MatrixTableResource,
    VersionedMatrixTableResource,
    import_sites_vcf,
    NO_CHR_TO_CHR_CONTIG_RECODING,
)
from gnomad.utils.vep import vep_or_lookup_vep
import hail as hl

from hail import Table


def _import_purcell_5k(path) -> hl.Table:
    p5k = hl.import_locus_intervals(path, reference_genome="GRCh37")
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
        )
    p5k = p5k.annotate(
        start=hl.liftover(p5k.interval.start, "GRCh38"),
        end=hl.liftover(p5k.interval.start, "GRCh38"),
    )
    p5k = p5k.filter(
        (p5k.start.contig == "chr" + p5k.interval.start.contig)
        & (p5k.end.contig == "chr" + p5k.interval.end.contig)
    )
    p5k = p5k.key_by()
    p5k = p5k.select(locus=p5k.start, locus_b37=p5k.interval.start)
    return p5k.key_by("locus")


def _import_clinvar(**kwargs) -> hl.Table:
    clinvar = import_sites_vcf(**kwargs)
    clinvar = clinvar.filter(
        hl.len(clinvar.alleles) > 1
    )  # Get around problematic single entry in alleles array in the clinvar vcf
    clinvar = vep_or_lookup_vep(clinvar, reference="GRCh38")
    return clinvar


def _import_dbsnp(**kwargs) -> hl.Table:
    dbsnp = import_sites_vcf(**kwargs)
    # Note: permit_shuffle is set because the dbsnp vcf has duplicate loci (turned into a set) so might be out of order
    dbsnp = hl.split_multi(dbsnp, permit_shuffle=True)
    dbsnp = dbsnp.group_by(dbsnp.locus, dbsnp.alleles).aggregate(
        rsid=hl.agg.collect_as_set(dbsnp.rsid)
    )

    return dbsnp


# Resources with no versioning needed
purcell_5k_intervals = TableResource(
    path="gs://gnomad-public/resources/grch38/purcell_5k_intervals/purcell5k.ht",
    import_func=_import_purcell_5k,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/purcell_5k_intervals/purcell5k.interval_list",
    },
)

na12878_giab = MatrixTableResource(
    path="gs://gnomad-public/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt",
    import_func=hl.import_vcf,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": "GRCh38",
    },
)

na12878_giab_hc_intervals = TableResource(
    path="gs://gnomad-public/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7_hc_regions.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed",
        "reference_genome": "GRCh38",
        "skip_invalid_intervals": True,
    },
)

# Versioned resources: versions should be listed from most recent to oldest
vep_context = VersionedTableResource(
    default_version="95",
    versions={
        "95": TableResource(
            path="gs://gnomad-public-requester-pays/resources/context/grch38_context_vep_annotated.ht",
        )
    },
)

syndip = VersionedMatrixTableResource(
    default_version="20180222",
    versions={
        "20180222": MatrixTableResource(
            path="gs://gnomad-public/resources/grch38/syndip/syndip.b38_20180222.mt",
            import_func=hl.import_vcf,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/syndip/full.38.20180222.vcf.gz",
                "force_bgz": True,
                "min_partitions": 100,
                "reference_genome": "GRCh38",
            },
        )
    },
)

syndip_hc_intervals = VersionedTableResource(
    default_version="20180222",
    versions={
        "20180222": TableResource(
            path="gs://gnomad-public/resources/grch38/syndip/syndip_b38_20180222_hc_regions.ht",
            import_func=hl.import_bed,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/syndip/syndip.b38_20180222.bed",
                "reference_genome": "GRCh38",
                "skip_invalid_intervals": True,
                "min_partitions": 10,
            },
        )
    },
)

clinvar = VersionedTableResource(
    default_version="20190923",
    versions={
        "20190923": TableResource(
            path="gs://gnomad-public/resources/grch38/clinvar/clinvar_20190923.ht",
            import_func=_import_clinvar,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/clinvar/clinvar_20190923.vcf.gz",
                "force_bgz": True,
                "contig_recoding": NO_CHR_TO_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 100,
                "reference_genome": "GRCh38",
            },
        )
    },
)

dbsnp = VersionedTableResource(
    default_version="b154",
    versions={
        "b154": TableResource(
            path="gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b154_grch38_all_20200514.ht",
            import_func=_import_dbsnp,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.bgz",
                "header_file": "gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.header",
                "force_bgz": True,
                "contig_recoding": DBSNP_B154_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 400,
                "reference_genome": "GRCh38",
            },
        ),
        "b151": TableResource(
            path="gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.vcf.bgz",
                "header_file": "gs://gnomad-public/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.vcf.header",
                "force_bgz": True,
                "contig_recoding": NO_CHR_TO_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 400,
                "reference_genome": "GRCh38",
            },
        ),
    },
)

hapmap = TableResource(
    path="gs://gnomad-public/resources/grch38/hapmap/hapmap_3.3.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz",
        "force_bgz": True,
        "reference_genome": "GRCh38",
    },
)

kgp_omni = TableResource(
    path="gs://gnomad-public/resources/grch38/kgp/1000G_omni2.5.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
        "force_bgz": True,
        "reference_genome": "GRCh38",
    },
)

kgp = VersionedTableResource(
    default_version="phase_1_hc",
    versions={
        "phase_1_hc": TableResource(
            path="gs://gnomad-public/resources/grch38/kgp/1000G_phase1.snps.high_confidence.hg38.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "force_bgz": True,
                "reference_genome": "GRCh38",
            },
        )
    },
)

mills = TableResource(
    path="gs://gnomad-public/resources/grch38/mills/Mills_and_1000G_gold_standard.indels.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "force_bgz": True,
        "reference_genome": "GRCh38",
    },
)

lcr_intervals = TableResource(
    path="gs://gnomad-public/resources/grch38/lcr_intervals/LCRFromHengHg38.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/lcr_intervals/LCRFromHengHg38.txt",
        "reference_genome": "GRCh38",
        "skip_invalid_intervals": True,
    },
)

seg_dup_intervals = TableResource(
    path="gs://gnomad-public/resources/grch38/seg_dup_intervals/GRCh38_segdups.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/seg_dup_intervals/GRCh38_segdups.bed",
        "reference_genome": "GRCh38",
    },
)

telomeres_and_centromeres = TableResource(
    path="gs://gnomad-public/resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.bed",
        "reference_genome": "GRCh38",
        "skip_invalid_intervals": True,
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

    return (
        hapmap.ht()
        .select(hapmap=True)
        .join(kgp_omni.ht().select(omni=True), how="outer")
        .join(kgp.versions["phase_1_hc"].ht().select(kgp_phase1_hc=True), how="outer")
        .join(mills.ht().select(mills=True), how="outer")
        .repartition(200, shuffle=False)
        .persist()
    )
