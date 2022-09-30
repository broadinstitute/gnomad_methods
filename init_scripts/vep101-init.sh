#!/bin/bash

export PROJECT="$(gcloud config get-value project)"
export VEP_CONFIG_PATH="$(/usr/share/google/get_metadata_value attributes/VEP_CONFIG_PATH)"
export VEP_REPLICATE="$(/usr/share/google/get_metadata_value attributes/VEP_REPLICATE)"
export VEP_BUCKET=hail-${VEP_REPLICATE}-vep
export ASSEMBLY=GRCh38
export VEP_DOCKER_IMAGE=gcr.io/broad-mpg-gnomad/vep_101

mkdir -p /vep_data/loftee_data
mkdir -p /vep_data/homo_sapiens

# Install docker
apt-get update
apt-get -y install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg2 \
    software-properties-common \
    tabix
curl -fsSL https://download.docker.com/linux/debian/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"
apt-get update
apt-get install -y --allow-unauthenticated docker-ce


gsutil -u $PROJECT cat gs://${VEP_BUCKET}/loftee-beta/${ASSEMBLY}.tar | tar -xf - -C /vep_data/ &
gsutil -u $PROJECT cat gs://${VEP_BUCKET}/Plugins.tar /vep_data/Plugins.tar | tar -xf - -C /vep_data
docker pull ${VEP_DOCKER_IMAGE} &
wait

################################################################
# Added stuff
################################################################

# GCS copy of ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_merged_vep_101_GRCh38.tar.gz
gsutil -u $PROJECT cat gs://gcp-public-data--gnomad/resources/vep/v101/homo_sapiens_merged_vep_101_GRCh38.tar.gz | tar -xzf - -C /vep_data

wait

# FASTA file from Hail's VEP 101 data
gsutil -u $PROJECT cp gs://gcp-public-data--gnomad/resources/vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /vep_data/
gsutil -u $PROJECT cp gs://gcp-public-data--gnomad/resources/vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai /vep_data/
gsutil -u $PROJECT cp gs://gcp-public-data--gnomad/resources/vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi /vep_data/

# Create config file.

cat > /vep_data/vep101-GRCh38.json <<EOF
{"command": [
    "/vep",
    "--format", "vcf",
    "__OUTPUT_FORMAT_FLAG__",
    "--everything",
    "--allele_number",
    "--no_stats",
    "--cache", "--offline",
    "--minimal",
    "--assembly", "GRCh38",
    "--merged",
    "--fasta", "/opt/vep/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    "--plugin", "LoF,loftee_path:/opt/vep/Plugins/,gerp_bigwig:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee.sql",
    "--dir_plugins", "/opt/vep/Plugins/",
    "-o", "STDOUT"
],
 "env": {},
 "vep_json_schema": "Struct{assembly_name:String,allele_string:String,ancestral:String,colocated_variants:Array[Struct{aa_allele:String,aa_maf:Float64,afr_allele:String,afr_maf:Float64,allele_string:String,amr_allele:String,amr_maf:Float64,clin_sig:Array[String],end:Int32,eas_allele:String,eas_maf:Float64,ea_allele:String,ea_maf:Float64,eur_allele:String,eur_maf:Float64,exac_adj_allele:String,exac_adj_maf:Float64,exac_allele:String,exac_afr_allele:String,exac_afr_maf:Float64,exac_amr_allele:String,exac_amr_maf:Float64,exac_eas_allele:String,exac_eas_maf:Float64,exac_fin_allele:String,exac_fin_maf:Float64,exac_maf:Float64,exac_nfe_allele:String,exac_nfe_maf:Float64,exac_oth_allele:String,exac_oth_maf:Float64,exac_sas_allele:String,exac_sas_maf:Float64,id:String,minor_allele:String,minor_allele_freq:Float64,phenotype_or_disease:Int32,pubmed:Array[Int32],sas_allele:String,sas_maf:Float64,somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,appris:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,minimised:Int32,polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,strand:Int32,swissprot:String,transcript_id:String,trembl:String,tsl:Int32,uniparc:String,variant_allele:String}],variant_class:String}"
}
EOF

ln -s /vep_data/vep101-GRCh38.json $VEP_CONFIG_PATH


################################################################
# End added stuff
################################################################

cat >/vep.c <<EOF
#include <unistd.h>
#include <stdio.h>

int
main(int argc, char *const argv[]) {
  if (setuid(geteuid()))
    perror( "setuid" );

  execv("/vep.sh", argv);
  return 0;
}
EOF
gcc -Wall -Werror -O2 /vep.c -o /vep
chmod u+s /vep

cat >/vep.sh <<EOF
#!/bin/bash

docker run -i -v /vep_data/:/opt/vep/.vep/:ro ${VEP_DOCKER_IMAGE} \
  /opt/vep/src/ensembl-vep/vep "\$@"
EOF
chmod +x /vep.sh
