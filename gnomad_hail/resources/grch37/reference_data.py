from gnomad_hail.resources.resource_utils import (
    MatrixTableResource,
    TableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
)

from hail import Table


na12878_giab = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.mt",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.bgz",
    },
)

hapmap = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/hapmap_3.3.b37.mt",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/hapmap_3.3.b37.vcf.bgz",
    },
)

kgp_omni = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/1000G_omni2.5.b37.mt",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/1000G_omni2.5.b37.vcf.bgz",
    },
)

mills = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/Mills_and_1000G_gold_standard.indels.b37.mt",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/Mills_and_1000G_gold_standard.indels.b37.vcf.bgz",
    },
)

syndip = MatrixTableResource(
    path="gs://gnomad-public/truth-sets/hail-0.2/hybrid.m37m.mt",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/hybrid.m37m.vcf.bgz",
    },
)

# Versioned resources: versions should be listed from most recent to oldest
dbsnp = VersionedTableResource(
    default_version="2018-04-23",
    versions={
        "2018-04-23": TableResource(
            path="gs://gnomad-public/truth-sets/source/All_20180423.ht",
            import_sources={
                "source_path": "gs://gnomad-public/truth-sets/source/All_20180423.vcf.bgz",
            },
        )
    },
)

clinvar = VersionedTableResource(
    default_version="20181028",
    versions={
        "20181028": TableResource(
            path="gs://gnomad-resources/clinvar/hail-0.2/clinvar_20181028.vep.ht",
            import_sources={
                "source_path": "gs://gnomad-resources/clinvar/source/clinvar_20181028.vcf.bgz",
            },
        )
    },
)

kgp = VersionedMatrixTableResource(
    default_version="phase_3_split",
    versions={
        "phase_3_split": MatrixTableResource(
            path="gs://gnomad-public/truth-sets/hail-0.2/1000Genomes_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.split.mt",
            import_sources={
                "source_path": "gs://genomics-public-data/1000-genomes-phase-3/vcf-20150220/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf",
            },
        ),
        "phase_3": MatrixTableResource(
            path="gs://gnomad-public/truth-sets/hail-0.2/1000Genomes_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.mt",
            import_sources={
                "source_path": "gs://genomics-public-data/1000-genomes-phase-3/vcf-20150220/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf",
            },
        ),
        "phase_1_hc": MatrixTableResource(
            path="gs://gnomad-public/truth-sets/hail-0.2/1000G_phase1.snps.high_confidence.b37.mt",
            import_sources={
                "source_path": "gs://gnomad-public/truth-sets/source/1000G_phase1.snps.high_confidence.b37.vcf.bgz",
            },
        ),
    },
)

cpg_sites = TableResource(path="gs://gnomad-public/resources/methylation/cpg.ht")

methylation_sites = TableResource(
    path="gs://gnomad-public/resources/methylation/methylation.ht"
)

lcr_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/LCR.GRCh37_compliant.interval_list.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/LCR.GRCh37_compliant.interval_list",
    },
)

decoy_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/mm-2-merged.GRCh37_compliant.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/mm-2-merged.GRCh37_compliant.bed",
    },
)

purcell_5k_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/purcell5k.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/purcell5k.interval_list",
    },
)

seg_dup_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/hg19_self_chain_split_both.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/hg19_self_chain_split_both.bed",
    },
)

exome_hc_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/exomes_high_coverage.auto.interval_list.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list",
    },
)

exome_calling_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/exome_calling_regions.v1.interval_list.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/exome_calling_regions.v1.interval_list",
    },
)

exome_evaluation_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/exome_evaluation_regions.v1.noheader.interval_list.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/exome_evaluation_regions.v1.noheader.interval_list",
    },
)

genome_evaluation_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list",
    },
)

high_coverage_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/high_coverage.auto.interval_list.ht",
    import_sources={
        "source_path": "gs://gnomad-public/intervals/high_coverage.auto.interval_list",
    },
)

na12878_hc_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/NA12878_GIAB_highconf_intervals.ht",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed"
    },
)

na12878_hc_exome_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/NA12878_GIAB_highconf_exome_intervals.ht",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed",
    },
)

syndip_hc_intervals = TableResource(
    path="gs://gnomad-public/resources/intervals/syndip_highconf_exome_intervals.ht",
    import_sources={
        "source_path": "gs://gnomad-public/truth-sets/source/hybrid.m37m.bed"
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
