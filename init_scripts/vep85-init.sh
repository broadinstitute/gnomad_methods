#!/bin/bash


# Copy VEP
mkdir -p /vep/homo_sapiens
gsutil -m cp -r gs://hail-common/vep/vep/loftee /vep
gsutil -m cp -r gs://hail-common/vep/vep/ensembl-tools-release-85 /vep
gsutil -m cp -r gs://hail-common/vep/vep/loftee_data /vep
gsutil -m cp -r gs://hail-common/vep/vep/Plugins /vep
gsutil -m cp -r gs://hail-common/vep/vep/homo_sapiens/85_GRCh37 /vep/homo_sapiens/
gsutil cp gs://hail-common/vep/vep/vep-gcloud.properties /vep

#Create symlink to vep
ln -s /vep/ensembl-tools-release-85/scripts/variant_effect_predictor /vep

# Hack since Hail VEP version is hardcoded
ln -s /vep/homo_sapiens/85_GRCh37/ /vep/homo_sapiens/81_GRCh37

#Give perms
chmod -R 777 /vep

# Copy perl JSON module
gsutil -m cp -r gs://hail-common/vep/perl-JSON/* /usr/share/perl/5.20/

#Copy perl DBD::SQLite module
gsutil -m cp -r gs://hail-common/vep/perl-SQLITE/* /usr/share/perl/5.20/


# Copy htslib and samtools
gsutil cp gs://hail-common/vep/htslib/* /usr/bin/
gsutil cp gs://hail-common/vep/samtools /usr/bin/
chmod a+rx  /usr/bin/tabix
chmod a+rx  /usr/bin/bgzip
chmod a+rx  /usr/bin/htsfile
chmod a+rx  /usr/bin/samtools

#Run VEP on the 1-variant VCF to create fasta.index file -- caution do not make fasta.index file writeable afterwards!
gsutil cp gs://hail-common/vep/vep/1var.vcf /vep
gsutil cp gs://hail-common/vep/vep/run_hail_vep_vcf.sh /vep
chmod a+rx /vep/run_hail_vep_vcf.sh

/vep/run_hail_vep_vcf.sh /vep/1var.vcf