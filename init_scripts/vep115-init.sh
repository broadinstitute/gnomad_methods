#!/bin/bash

export PROJECT="$(gcloud config get-value project)"
export VEP_CONFIG_PATH="$(/usr/share/google/get_metadata_value attributes/VEP_CONFIG_PATH)"
export VEP_REPLICATE="$(/usr/share/google/get_metadata_value attributes/VEP_REPLICATE)"
export VEP_BUCKET=hail-${VEP_REPLICATE}-vep
export ASSEMBLY=GRCh38
export VEP_DOCKER_IMAGE=us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vep115

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


gsutil -u ${PROJECT} cat gs://${VEP_BUCKET}/loftee-beta/${ASSEMBLY}.tar | tar -xf - -C /vep_data/

docker pull ${VEP_DOCKER_IMAGE}


################################################################
# Added stuff
################################################################

# GCS copy of https://ftp.ebi.ac.uk/ensemblorg/pub/release-115/variation/indexed_vep_cache/homo_sapiens_merged_vep_115_GRCh38.tar.gz
gsutil -u $PROJECT cat gs://gcp-public-data--gnomad/resources/vep/v115/homo_sapiens_merged_vep_115_GRCh38.tar.gz | tar -xzf - -C /vep_data


# FASTA file from Hail's VEP 95 data
gsutil -u "$PROJECT" cp gs://gcp-public-data--gnomad/resources/vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /vep_data/
gsutil -u $PROJECT cp gs://gcp-public-data--gnomad/resources/vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai /vep_data/
gsutil -u $PROJECT cp gs://gcp-public-data--gnomad/resources/vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi /vep_data/

# Create config file.
cat > /vep_data/vep115-GRCh38.json <<EOF
{"command": [
    "/vep",
    "--format", "vcf",
    "__OUTPUT_FORMAT_FLAG__",
    "--everything",
    "--allele_number",
    "--no_stats",
    "--regulatory",
    "--cache", "--offline",
    "--minimal",
    "--assembly", "GRCh38",
    "--merged",
    "--fasta", "/opt/vep/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    "--plugin", "LoF,loftee_path:/opt/vep/Plugins,gerp_bigwig:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee.sql",
    "--plugin", "context",
    "--plugin", "ancestral,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz",
    "--dir_plugins", "/opt/vep/Plugins/",
    "-o", "STDOUT"
],
 "env": {},
"vep_json_schema": "Struct{allele_string:String,assembly_name:String,colocated_variants:Array[Struct{allele_string:String,clin_sig:Array[String],clin_sig_allele:String,end:Int32,evidence:Array[String],id:String,minor_allele:String,minor_allele_freq:Float64,phenotype_or_disease:Int32,phenotypes:Array[String],pubmed:Array[String],somatic:Int32,seq_region_name:String,start:Int32,strand:Int32}],end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,ancestral:String,consequence_terms:Array[String],context:String,impact:String,variant_allele:String}],minimised:Int32,most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],high_inf_pos:String,impact:String,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,transcription_factors:Array[String],strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,ancestral:String,biotype:String,consequence_terms:Array[String],context:String,impact:String,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,ancestral:String,appris:String,bam_edit:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],context:String,distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,given_ref:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_filter:String,lof_flags:String,lof_info:String,mane:Array[String],mane_select:String,mane_plus_clinical:String,polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_id:String,protein_start:Int32,sift_prediction:String,sift_score:Float64,source:String,strand:Int32,swissprot:Array[String],transcript_id:String,trembl:Array[String],tsl:String,uniparc:Array[String],uniprot_isoform:Array[String],used_ref:String,variant_allele:String}],variant_class:String}"
}
EOF

ln -s /vep_data/vep115-GRCh38.json $VEP_CONFIG_PATH


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
