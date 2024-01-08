# noqa: D100

import hail as hl
from hail import Table

from gnomad.resources.resource_utils import (
    DBSNP_B154_CHR_CONTIG_RECODING,
    NO_CHR_TO_CHR_CONTIG_RECODING,
    GnomadPublicMatrixTableResource,
    GnomadPublicTableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
    import_sites_vcf,
)
from gnomad.utils.vep import vep_or_lookup_vep


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
    # Note: permit_shuffle is set because the dbsnp vcf has duplicate loci
    # (turned into a set) so might be out of order
    dbsnp = hl.split_multi(dbsnp, permit_shuffle=True)
    dbsnp = dbsnp.group_by(dbsnp.locus, dbsnp.alleles).aggregate(
        rsid=hl.agg.collect_as_set(dbsnp.rsid)
    )

    return dbsnp


def _import_methylation_sites(path) -> hl.Table:
    """
    Import methylation data from bed file.

    :param path: Path to bed file containing methylation scores.
    :return: Table with methylation data.
    """
    ht = hl.import_bed(path, min_partitions=100, reference_genome="GRCh38")
    ht = ht.select(
        locus=ht.interval.start,
        methylation_level=hl.int32(ht.target),
    )

    return ht.key_by("locus").drop("interval")


def _import_ensembl_interval(path) -> hl.Table:
    """
    Import and parse Ensembl intervals of protein-coding genes to a Hail Table.

    File is expected to include only the following fields: gene_stable_ID, chr, start, end, source_gene, gene_name, and type.

    :param path: Path to the interval Table file.
    """
    ensembl = hl.import_table(
        path,
        delimiter="\t",
        min_partitions=100,
        impute=True,
    )

    ensembl = ensembl.key_by(
        interval=hl.locus_interval(
            "chr" + ensembl.chr,
            ensembl.start,
            ensembl.end,
            reference_genome="GRCh38",
        )
    )
    return ensembl


def _import_gencode_cds(gtf_path: str) -> hl.Table:
    """
    Get CDS intervals from GENCODE GTF file.

    :param gtf_path: Path to GENCODE GTF file.
    :return: Table with CDS intervals.
    """
    ht = hl.experimental.import_gtf(
        gtf_path,
        "GRCh38",
        force_bgz=True,
        min_partitions=12,
    )
    ht = ht.annotate(
        gene_id=ht.gene_id.split("\\.")[0],
        transcript_id=ht.transcript_id.split("\\.")[0],
    )
    ht = (
        ht.filter((ht.feature == "CDS") & (ht.transcript_type == "protein_coding"))
        .select("gene_id", "transcript_id")
        .distinct()
    )
    return ht


def _import_gtex_rsem(gtex_path: str, meta_path: str) -> hl.MatrixTable:
    """
    Import GTEx RSEM data from expression data and sample attributes file.

    .. note::
        Files are downloaded from https://www.gtexportal.org/home/downloads/adult-gtex.
        We get the transcript TPM under Bulk tissue expression and sample attributes
        under Metadata. The transcript TPM file is expected to have transcript
        expression data, with transcript IDs as the first column and gene IDs as the
        second column.
    :param gtex_path: Path to the GTEx RSEM file.
    :param meta_path: Path to the GTEx sample attributes file.
    :return: Matrix Table with GTEx RSEM data with tissue information.
    """
    meta_ht = hl.import_table(meta_path, force_bgz=True, impute=True)
    meta_ht = meta_ht.key_by("SAMPID")

    mt = hl.import_matrix_table(
        gtex_path,
        row_key="transcript_id",
        row_fields={"transcript_id": hl.tstr, "gene_id": hl.tstr},
        entry_type=hl.tfloat64,
        force_bgz=True,
        min_partitions=1000,
    )

    mt = mt.rename({"x": "transcript_tpm"})

    # GTEx data has gene IDs and transcript IDs with version numbers, we need
    # to remove the version numbers so that it can later be joined with the
    # variant Table
    mt = mt.annotate_cols(
        tissue=meta_ht[mt.col_id]
        .SMTSD.replace(" ", "")
        .replace("-", "_")
        .replace("\\(", "_")
        .replace("\\)", "")
    )
    mt = mt.key_rows_by()
    mt = mt.annotate_rows(
        transcript_id=mt.transcript_id.split("\\.")[0],
        gene_id=mt.gene_id.split("\\.")[0],
    )
    mt = mt.key_rows_by("transcript_id")
    return mt


# Resources with no versioning needed
purcell_5k_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/purcell_5k_intervals/purcell5k.ht",
    import_func=_import_purcell_5k,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch38/purcell_5k_intervals/purcell5k.interval_list",
    },
)

na12878_giab = GnomadPublicMatrixTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt",
    import_func=hl.import_vcf,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": "GRCh38",
    },
)

na12878_giab_hc_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7_hc_regions.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed",
        "reference_genome": "GRCh38",
        "skip_invalid_intervals": True,
    },
)

# Versioned resources: versions should be listed from most recent to oldest
vep_context = VersionedTableResource(
    default_version="95",
    versions={
        "95": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/context/grch38_context_vep_annotated.ht",
        ),
        "101": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/context/grch38_context_vep_annotated.v101.ht",
        ),
        "105": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/context/grch38_context_vep_annotated.v105.ht",
        ),
    },
)

syndip = VersionedMatrixTableResource(
    default_version="20180222",
    versions={
        "20180222": GnomadPublicMatrixTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/syndip/syndip.b38_20180222.mt",
            import_func=hl.import_vcf,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/syndip/full.38.20180222.vcf.gz",
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
        "20180222": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/syndip/syndip_b38_20180222_hc_regions.ht",
            import_func=hl.import_bed,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/syndip/syndip.b38_20180222.bed",
                "reference_genome": "GRCh38",
                "skip_invalid_intervals": True,
                "min_partitions": 10,
            },
        )
    },
)

# These Ensembl Interval Tables are focused on protein-coding genes on chr1-22,X,Y.
# Downloaded from the biomart of Ensembl Archive (https://useast.ensembl.org/info/website/archives/index.html)
# Ensembl 101 & 105 are included, since 101 was used to annotate gnomAD v3 and 105 to gnomAD v4.
# Basic stats: 19924 protein-coding genes in Ensembl 101, and1 19951
# protein-coding genes in Ensembl 105.
ensembl_interval = VersionedTableResource(
    default_version="105",
    versions={
        "105": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/ensembl/ensembl_105_pc_genes_grch38.ht",
            import_func=_import_ensembl_interval,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/ensembl/ensembl_105_pc_genes_grch38.tsv",
            },
        ),
        "101": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/ensembl/ensembl_101_pc_genes_grch38.ht",
            import_func=_import_ensembl_interval,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/ensembl/ensembl_101_pc_genes_grch38.tsv",
            },
        ),
    },
)

clinvar = VersionedTableResource(
    default_version="20190923",
    versions={
        "20190923": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/clinvar/clinvar_20190923.ht",
            import_func=_import_clinvar,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/clinvar/clinvar_20190923.vcf.gz",
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
    default_version="b156",
    versions={
        "b156": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b156_grch38_all_20221116.ht",
            import_func=_import_dbsnp,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/dbsnp/dbsnp_b156_grch38_all_GCF_000001405.40_20221116.vcf.bgz",
                "header_file": "gs://gcp-public-data--gnomad/resources/grch38/dbsnp/dbsnp_b156_grch38_all_GCF_000001405.40_20221116.vcf.header",
                "force_bgz": True,
                "contig_recoding": DBSNP_B154_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 400,
                "reference_genome": "GRCh38",
            },
        ),
        "b154": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b154_grch38_all_20200514.ht",
            import_func=_import_dbsnp,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.bgz",
                "header_file": "gs://gcp-public-data--gnomad/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.header",
                "force_bgz": True,
                "contig_recoding": DBSNP_B154_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 400,
                "reference_genome": "GRCh38",
            },
        ),
        "b151": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.vcf.bgz",
                "header_file": "gs://gcp-public-data--gnomad/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.vcf.header",
                "force_bgz": True,
                "contig_recoding": NO_CHR_TO_CHR_CONTIG_RECODING,
                "skip_invalid_loci": True,
                "min_partitions": 400,
                "reference_genome": "GRCh38",
            },
        ),
    },
)

hapmap = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/hapmap/hapmap_3.3.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": (
            "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz"
        ),
        "force_bgz": True,
        "reference_genome": "GRCh38",
    },
)

kgp_omni = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/kgp/1000G_omni2.5.hg38.ht",
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
        "phase_1_hc": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/kgp/1000G_phase1.snps.high_confidence.hg38.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "force_bgz": True,
                "reference_genome": "GRCh38",
            },
        )
    },
)

mills = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/mills/Mills_and_1000G_gold_standard.indels.hg38.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        "force_bgz": True,
        "reference_genome": "GRCh38",
    },
)

# Methylation scores range from 0-15 and are described in Chen et al
# (https://www.biorxiv.org/content/10.1101/2022.03.20.485034v2.full).
methylation_sites = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/methylation_sites/methylation.ht",
    import_func=_import_methylation_sites,
    import_args={
        "path": "gs://gnomad-public-requester-pays/resources/grch38/methylation_sites/methylation.bed",
    },
)

# Methylation scores for chromosome X range from 0-12 and are described in Chen et al
# (https://www.biorxiv.org/content/10.1101/2022.03.20.485034v2.full).
methylation_sites_chrx = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/methylation_sites/methylation_chrX.ht",
    import_func=_import_methylation_sites,
    import_args={
        "path": "gs://gnomad-public-requester-pays/resources/grch38/methylation_sites/methylation_chrX.bed",
    },
)

lcr_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/lcr_intervals/LCRFromHengHg38.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch38/lcr_intervals/LCRFromHengHg38.txt",
        "reference_genome": "GRCh38",
        "skip_invalid_intervals": True,
    },
)

seg_dup_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/seg_dup_intervals/GRCh38_segdups.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch38/seg_dup_intervals/GRCh38_segdups.bed",
        "reference_genome": "GRCh38",
    },
)

telomeres_and_centromeres = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.bed",
        "reference_genome": "GRCh38",
        "skip_invalid_intervals": True,
    },
)

gtex_rsem = VersionedTableResource(
    default_version="v10",
    versions={
        "v10": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/gtex_rsem/gtex_rsem_v10.mt",
            import_func=_import_gtex_rsem,
            import_args={  # TODO: update path and name
                "gtex_path": "gs://gnomad-qin/bulk-gex_v7_rna-seq_GTEx_Analysis_20231130_v10_RSEMv1.2.22_transcript_tpm.txt.gz",
                "meta_path": "gs://gnomad-qin/annotations_v10_GTEx_v10_Annotations_SampleAttributesDS.txt.gz",
            },
        ),
    },
)

gencode_cds = VersionedTableResource(
    default_version="v39",
    versions={
        "v39": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch38/gencode_cds/gencode_v39_cds.ht",
            import_func=_import_gencode_cds,
            import_args={
                "gtf_path": "gs://gcp-public-data--gnomad/resources/grch38/gencode/gencode.v39.annotation.gtf.gz",
            },
        ),
    },
)


def get_truth_ht() -> Table:
    """
    Return a table with annotations from the latest version of the corresponding truth data.

    The following annotations are included:
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
