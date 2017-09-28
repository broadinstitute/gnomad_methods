#Evaluation data
truth_dir = 'gs://gnomad-public/truth-sets/hail-0.2'
omni_vds_path = "%s/1000G_omni2.5.b37.vds" % truth_dir
mills_vds_path = "%s/Mills_and_1000G_gold_standard.indels.b37.vds" % truth_dir
hapmap_vds_path = "%s/hapmap_3.3.b37.vds" % truth_dir
dbsnp_vcf_path = "%s/vcf/All_20160601.vcf.bgz" % truth_dir
kgp_high_conf_snvs_vds_path = "%s/1000G_phase1.snps.high_confidence.b37.vds" % truth_dir
NA12878_vds_path = "%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vds" % truth_dir
NA12878_high_conf_regions_bed_path = "%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed" % truth_dir
NA12878_high_conf_exome_regions_bed_path = "gs://gnomad-exomes/intervals/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
syndip_vds_path = "%s/hybrid.m37m.vds" % truth_dir
syndip_high_conf_regions_bed_path = "%s/hybrid.m37m.bed" % truth_dir
clinvar_tsv_path = "gs://gnomad-resources/annotations/clinvar_alleles.single.b37.tsv.gz"
clinvar_vds_path = "gs://gnomad-resources/annotations/clinvar_alleles.single.b37.vds"

#Exome/genome duplicate samples
exomes_to_combined_IDs_tsv_path = "gs://gnomad-resources/exomes_to_combined.IDs.txt"
exomes_qc_pass_samples_list_path = "gs://gnomad-resources/exomes_qc_pass_samples.txt.gz" #Samples that didn't fail QC-metric filters (contains relateds and non-releasable samples)
genomes_to_combined_IDs_tsv_path = "gs://gnomad-resources/genomes_to_combined.IDs.txt"
genomes_qc_pass_samples_list_path = "gs://gnomad-resources/genomes_qc_pass_samples.txt.gz" #Samples that didn't fail QC-metric filters (contains relateds and non-releasable samples)

#Usefult intervals
lcr_intervals_path = "gs://gnomad-public/intervals/LCR.interval_list"
decoy_intervals_path = "gs://gnomad-public/intervals/mm-2-merged.bed.gz"
purcell5k_intervals_path = "gs://gnomad-public/intervals/purcell5k.interval_list"

#Exome intervals
exomes_high_conf_regions_intervals_path = "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list"
exome_calling_intervals_path = 'gs://gnomad-exomes/intervals/exome_calling_regions.v1.interval_list'
exome_calling_noheader_intervals_path = 'gs://gnomad-exomes/intervals/exome_calling_regions.v1.nohead.interval_list'
evaluation_intervals_path = 'gs://gnomad-exomes/intervals/exome_evaluation_regions.v1.intervals'
high_coverage_intervals_path = 'gs://gnomad-exomes/intervals/high_coverage.auto.interval_list'

additional_vcf_header_path = "gs://gnomad-resources/gnomad.extra_header_fields.vcf"

vep_config = "/vep/vep-gcloud.properties"

# Full VDSs
full_exac_v1_vds_path = 'gs://gnomad-exomes-raw/exacv1/exac.all.vds'
full_exome_vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
full_exome_hardcalls_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw_hardcalls.vds'
full_exome_hardcalls_split_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw_hardcalls.split.vds'
full_genome_vds_path = 'gs://gnomad-genomes-raw/full/gnomad.genomes.all.vds'
full_genome_hardcalls_vds_path = "gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.vds"
full_genome_hardcalls_split_vds_path = "gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.split.vds"
full_genomes_vep_split_vds_path = "gs://gnomad-genomes-raw/gnomad.genomes.vep.split.vds"
full_exomes_vep_split_vds_path = "gs://gnomad-exomes-raw/gnomad.exomes.vep.split.vds"

# Adj hardcalls
full_exome_hardcalls_adj_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.adj_hardcalls.vds'
full_genome_hardcalls_adj_vds_path = 'gs://gnomad-genomes-raw/hardcalls/gnomad.adj_hardcalls.vds'

#Filtering VDSs
vqsr_vds_path = 'gs://gnomad-exomes/variantqc/gnomad.exomes.vqsr.unsplit.vds'
genomes_rf_vds_path = "gs://gnomad-genomes/variantqc/RF/gnomad.sites.RF.newStats24.vds"
exomes_rf_vds_path = "gs://gnomad-exomes/variantqc/gnomad.exomes.rf.vds"

# Release Sites VDSs
final_exac_sites_vds_path = 'gs://gnomad-exomes-raw/exacv1/exac.sites.vds'
final_exome_autosomes_vds_path = 'gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.autosomes.vds'
final_genome_autosomes_vds_path = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.autosomes.vds'
final_exome_vds_path = 'gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.vds'
final_exome_split_vds_path = 'gs://gnomad-public/release/2.0.1/vds/exomes/gnomad.exomes.r2.0.1.sites.split.vds'
final_genome_vds_path = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.vds'
final_genome_split_vds_path = 'gs://gnomad-public/release/2.0.1/vds/genomes/gnomad.genomes.r2.0.1.sites.split.vds'

# Metadata
genomes_meta_tsv_path = "gs://gnomad-genomes-raw/gnomad.final.all_meta.txt.bgz"
genomes_fam_path = "gs://gnomad-genomes-raw/gnomad.final.goodTrios.fam"
exomes_meta_tsv_path = 'gs://gnomad-exomes-raw/super_meta_june_02_2017.txt.gz'
exomes_fam_path = "gs://gnomad-exomes/variantqc/gnomad_exomes.qctrios.fam"

# PCA
gnomad_pca_vds_path = "gs://gnomad-genomes/sampleqc/gnomad.pca.vds"

# Annotations
methylation_kt_path = "gs://gnomad-resources/methylation.kt"


def check_resources():
    import subprocess
    for f, g in globals().items():
        if isinstance(g, str) and g.startswith('gs://'):
            try:
                subprocess.check_output(['gsutil', 'ls', g])
            except subprocess.CalledProcessError:
                print('WARNING: Missing {} ({})'.format(f, g))