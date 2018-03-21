# Streamline genomes metadata file

args=commandArgs(TRUE)

# Load data
genomes_meta_input=args[1] # gs://gnomad/sample_qc/input_meta/gnomad.genomes.metadata.2017-06-02.tsv.bgz
new_releasables_input=args[2] # gs://gnomad/sample_qc/input_meta/gnomad.genomes.2.1_additional_releasables.txt
colmap_path=args[3] # gs://gnomad/sample_qc/input_meta/gnomad.genomes.column_map.2018-03-21.txt
streamlined_meta_output=args[4] # NOTE: output is not gzipped


# Update releasables for 2.1
genomes=read.delim(genomes_meta_input, header=T)
new_releasables=read.delim(new_releasables_input, header=T)

genomes$releasable_updated=genomes$releasable
genomes[genomes$Sample %in% new_releasables$sample,]$releasable_updated="TRUE"
genomes$releasable_updated=as.factor(genomes$releasable_updated)
table(genomes$releasable, genomes$releasable_updated)
#         FALSE  TRUE
#   FALSE  3154   817
#   TRUE      0 16343

#------------------
# Reformat column names
colmap=read.delim(colmap_path, header=T)
colmap=subset(colmap, selected_for_condensed_meta=="TRUE")

genomes_reformat=genomes[,names(genomes) %in% colmap$existing_meta_annotation]
names(genomes_reformat)=colmap$condensed_meta_colname[match(names(genomes_reformat), colmap$existing_meta_annotation)]

#------------------
# Transform column values

# Flip Boolean: non_psych -> neuro
summary(genomes_reformat$neuro)
#    Mode   FALSE    TRUE    NA's 
# logical    5619   14695       0 
genomes_reformat$neuro[genomes_reformat$neuro=="TRUE"]="False"
genomes_reformat$neuro[genomes_reformat$neuro=="FALSE"]="True"
genomes_reformat$neuro=as.factor(genomes_reformat$neuro)

# Transform to Boolean: control
summary(genomes_reformat$control)
#    case control    NA's 
#   10901    5203    4210 
genomes_reformat$control=as.character(genomes_reformat$control)
genomes_reformat$control[genomes_reformat$control=="control"]="True"
genomes_reformat$control[genomes_reformat$control=="case"]="False"
genomes_reformat$control=as.factor(genomes_reformat$control)


# Add columns for parity with exomes
genomes_reformat$internal="True"
genomes_reformat$internal=as.factor(genomes_reformat$internal)

genomes_reformat$cloudable="True"
genomes_reformat$cloudable=as.factor(genomes_reformat$cloudable)

#------------------
# Create detailed platforms annotation ("qc_platform")

genomes_reformat$qc_platform=as.character(genomes_reformat$product_simplified)
genomes_reformat$qc_platform[genomes_reformat$product_simplified=="Standard High Coverage Whole Genome Sequencing (30x)" & genomes_reformat$mean_read_length<125]="Standard High Coverage Whole Genome Sequencing (30x) | Low Read Length"
genomes_reformat$qc_platform[genomes_reformat$product_simplified=="Standard High Coverage Whole Genome Sequencing (30x)" & genomes_reformat$mean_read_length>=125]="Standard High Coverage Whole Genome Sequencing (30x) | High Read Length"
decoy_projects=c("G89396", "G89400", "G89401", "G89402", "G89403", "G90053", "G90699", "G90700", "G90701")
genomes_reformat$qc_platform[genomes_reformat$project_id %in% decoy_projects]="PCR-Free Human WGS (Standard Coverage) + Decoy"
genomes_reformat$qc_platform=as.factor(genomes_reformat$qc_platform)

#------------------
# Create syndip annotation

genomes_reformat$syndip="False"
genomes_reformat$syndip[genomes_reformat$s=="CHMI_CHMI3_WGS1"]="True"
genomes_reformat$syndip=as.factor(genomes_reformat$syndip)

write.table(genomes_reformat, streamlined_meta_output, col.names=T, row.names=F, quote=F, sep="\t")

