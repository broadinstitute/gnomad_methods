# noqa: D100

from typing import Optional, Union

import hail as hl

from gnomad.resources.resource_utils import (
    GnomadPublicMatrixTableResource,
    GnomadPublicTableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
    import_gencode,
    import_sites_vcf,
)
from gnomad.utils.constraint import transform_methylation_level


def _import_gtex_rsem(gtex_path: str, meta_path: str, **kwargs) -> hl.MatrixTable:
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
    :param kwargs: Any additional parameters to be passed to Hail's `import_matrix_table`.
    :return: Matrix Table with GTEx RSEM data with tissue information.
    """
    meta_ht = hl.import_table(meta_path, force_bgz=True, impute=True)
    meta_ht = meta_ht.key_by("SAMPID")

    mt = hl.import_matrix_table(
        gtex_path,
        row_fields={"transcript_id": hl.tstr, "gene_id": hl.tstr},
        entry_type=hl.tfloat64,
        force_bgz=True,
        **kwargs,
    )

    mt = mt.rename({"x": "transcript_tpm", "col_id": "s"})

    # GTEx data has gene IDs and transcript IDs with version numbers, we need
    # to remove the version numbers so that it can later be joined with VEP
    # transcript consequences transcript_id.
    mt = mt.annotate_cols(
        tissue=meta_ht[mt.s]
        .SMTSD.replace(" ", "")
        .replace("-", "_")
        .replace("\\(", "_")
        .replace("\\)", "")
    )
    mt = mt.annotate_rows(
        transcript_id=mt.transcript_id.split("\\.")[0],
        gene_id=mt.gene_id.split("\\.")[0],
    )
    mt = mt.key_rows_by("transcript_id").drop("row_id")

    return mt


na12878_giab = GnomadPublicMatrixTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/na12878/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.mt",
    import_func=hl.import_vcf,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/na12878/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.bgz",
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": "GRCh37",
    },
)

hapmap = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/hapmap/hapmap_3.3.b37.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/hapmap/hapmap_3.3.b37.vcf.bgz",
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": "GRCh37",
    },
)

kgp_omni = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/kgp/1000G_omni2.5.b37.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/kgp/1000G_omni2.5.b37.vcf.bgz",
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": "GRCh37",
    },
)

mills = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/mills/Mills_and_1000G_gold_standard.indels.b37.ht",
    import_func=import_sites_vcf,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/mills/Mills_and_1000G_gold_standard.indels.b37.vcf.bgz",
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": "GRCh37",
    },
)

syndip = GnomadPublicMatrixTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/syndip/hybrid.m37m.mt",
    import_func=hl.import_vcf,
    import_args={
        "path": (
            "gs://gcp-public-data--gnomad/resources/grch37/syndip/hybrid.m37m.vcf.bgz"
        ),
        "min_partitions": 100,
        "reference_genome": "GRCh37",
    },
)

# Versioned resources: versions should be listed from most recent to oldest
vep_context = VersionedTableResource(
    default_version="85",
    versions={
        "85": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/context/grch37_context_vep_annotated.ht",
        )
    },
)

dbsnp = VersionedTableResource(
    default_version="20180423",
    versions={
        "20180423": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/dbsnp/All_20180423.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch37/dbsnp/All_20180423.vcf.bgz",
                "force_bgz": True,
                "skip_invalid_loci": True,
                "min_partitions": 100,
                "reference_genome": "GRCh37",
            },
        )
    },
)

clinvar = VersionedTableResource(
    default_version="20181028",
    versions={
        "20181028": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/clinvar/clinvar_20181028.vep.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch37/clinvar/clinvar_20181028.vcf.bgz",
                "force_bgz": True,
                "skip_invalid_loci": True,
                "min_partitions": 100,
                "reference_genome": "GRCh37",
            },
        )
    },
)

kgp_phase_3 = VersionedMatrixTableResource(
    default_version="phase_3_split",
    versions={
        "phase_3_split": GnomadPublicMatrixTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/kgp/1000Genomes_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.split.mt",
            import_func=hl.import_vcf,
            import_args={
                "path": "gs://genomics-public-data/1000-genomes-phase-3/vcf-20150220/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf",
                "force_bgz": True,
                "skip_invalid_loci": True,
                "min_partitions": 300,
                "reference_genome": "GRCh37",
            },
        ),
        "phase_3": GnomadPublicMatrixTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/kgp/1000Genomes_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.mt",
            import_func=hl.import_vcf,
            import_args={
                "path": "gs://genomics-public-data/1000-genomes-phase-3/vcf-20150220/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf",
                "force_bgz": True,
                "skip_invalid_loci": True,
                "min_partitions": 300,
                "reference_genome": "GRCh37",
            },
        ),
    },
)

kgp = VersionedTableResource(
    default_version="phase_1_hc",
    versions={
        "phase_1_hc": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/kgp/1000G_phase1.snps.high_confidence.b37.ht",
            import_func=import_sites_vcf,
            import_args={
                "path": "gs://gcp-public-data--gnomad/resources/grch37/kgp/1000G_phase1.snps.high_confidence.b37.vcf.bgz",
                "force_bgz": True,
                "skip_invalid_loci": True,
                "min_partitions": 100,
                "reference_genome": "GRCh37",
            },
        ),
    },
)

cpg_sites = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/cpg_sites/cpg.ht"
)

methylation_sites = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/methylation_sites/methylation.ht"
)

lcr_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/lcr_intervals/LCR.GRCh37_compliant.interval_list.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/lcr_intervals/LCR.GRCh37_compliant.interval_list",
        "reference_genome": "GRCh37",
    },
)

decoy_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/decoy_intervals/mm-2-merged.GRCh37_compliant.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/decoy_intervals/mm-2-merged.GRCh37_compliant.bed",
        "reference_genome": "GRCh37",
    },
)

purcell_5k_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/purcell_5k_intervals/purcell5k.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/purcell_5k_intervals/purcell5k.interval_list",
        "reference_genome": "GRCh37",
    },
)

seg_dup_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/seg_dup_intervals/hg19_self_chain_split_both.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/seg_dup_intervals/hg19_self_chain_split_both.bed",
        "reference_genome": "GRCh37",
    },
)

exome_hc_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/broad_intervals/exomes_high_coverage.auto.interval_list.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/broad_intervals/exomes_high_coverage.auto.interval_list",
        "reference_genome": "GRCh37",
    },
)

high_coverage_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/broad_intervals/high_coverage.auto.interval_list.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/broad_intervals/high_coverage.auto.interval_list",
        "reference_genome": "GRCh37",
    },
)

exome_calling_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/broad_intervals/exome_calling_regions.v1.interval_list.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/broad_intervals/exome_calling_regions.v1.interval_list",
        "reference_genome": "GRCh37",
    },
)

exome_evaluation_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/broad_intervals/exome_evaluation_regions.v1.noheader.interval_list.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/broad_intervals/exome_evaluation_regions.v1.noheader.interval_list",
        "reference_genome": "GRCh37",
    },
)

genome_evaluation_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/broad_intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/broad_intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list",
        "reference_genome": "GRCh37",
    },
)

na12878_hc_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/na12878/NA12878_GIAB_highconf_intervals.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/na12878/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed",
        "reference_genome": "GRCh37",
    },
)

syndip_hc_intervals = GnomadPublicTableResource(
    path="gs://gnomad-public-requester-pays/resources/grch37/syndip/syndip_highconf_genome_intervals.ht",
    import_func=hl.import_bed,
    import_args={
        "path": "gs://gcp-public-data--gnomad/resources/grch37/syndip/hybrid.m37m.bed",
        "reference_genome": "GRCh37",
    },
)


def get_truth_ht() -> hl.Table:
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
        .join(
            kgp.versions["phase_1_hc"].mt().rows().select(kgp_phase1_hc=True),
            how="outer",
        )
        .join(mills.ht().select(mills=True), how="outer")
        .repartition(200, shuffle=False)
        .persist()
    )


gtex_rsem = VersionedMatrixTableResource(
    default_version="v7",
    versions={
        "v7": GnomadPublicMatrixTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/gtex_rsem/gtex_rsem_v7.mt",
            import_func=_import_gtex_rsem,
            import_args={
                "gtex_path": "gs://gcp-public-data--gnomad/resources/grch37/gtex/bulk-gex_v7_rna-seq_GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz",
                "meta_path": "gs://gcp-public-data--gnomad/resources/grch37/gtex/annotations_v7_GTEx_v7_Annotations_SampleAttributesDS.txt.gz",
                "min_partitions": 1000,
            },
        ),
    },
)

gencode = VersionedTableResource(
    default_version="v19",
    versions={
        "v19": GnomadPublicTableResource(
            path="gs://gnomad-public-requester-pays/resources/grch37/gencode/gencode.v19.annotation.ht",
            import_func=import_gencode,
            import_args={
                "gtf_path": "gs://gcp-public-data--gnomad/resources/grch37/gencode/gencode.v19.annotation.gtf.gz",
                "reference_genome": "GRCh37",
                "force_bgz": True,
                "min_partitions": 10,
            },
        ),
    },
)


def transform_grch37_methylation(
    ht: Optional[hl.Table] = None,
    methylation_expr: Optional[hl.expr.NumericExpression] = None,
) -> Union[hl.Table, hl.expr.NumericExpression]:
    """
    Transform methylation level from the GRCh37 methylation resource to a 0-2 scale.

    .. note::

        One of ht or methylation_expr must be provided.

    The GRCh37 methylation resource provides a MEAN score ranging from 0 to 1. We transform
    this to a 0-2 scale by converting any value greater than 0.6 to 2, any value less
    than or equal to 0.2 to 0, and any value in between to 1.

    :param ht: Optional Hail Table with methylation data. Default is None.
    :param methylation_expr: Optional methylation level expression. Default is None.
    :return: Transformed methylation level expression or annotated Hail Table.
    """
    if ht is None and methylation_expr is None:
        raise ValueError("One of 'ht' or 'methylation_expr' must be provided")

    return transform_methylation_level(
        methylation_expr="MEAN" if methylation_expr is None else methylation_expr,
        methylation_cutoffs=(0.2, 0.6),
        ht=ht,
    )
