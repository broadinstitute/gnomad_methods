# Streamline exomes metadata file

args=commandArgs(TRUE)

# Load data
exomes_meta_input=args[1] # gs://gnomad/sample_qc/input_meta/gnomad.exomes.metadata_import_table.2.1.txt
colmap_path=args[2] # gs://gnomad/sample_qc/input_meta/gnomad.exomes.column_map.2018-03-21.txt
old_releasables_input=args[3] # gs://gnomad/sample_qc/input_meta/gnomad.exomes.metadata.2017-06-02.reformatted_colnames.tsv.txt.gz
streamlined_meta_output=args[4] # NOTE: output is not gzipped


exomes=read.delim(exomes_meta_input, header=T)
# NOTE: Both permissions and metrics have been updated to 2.1 in this file by Kristen

#---------------
# Reformat column names

colmap=read.delim(colmap_path, header=T)
colmap=subset(colmap, selected_for_condensed_meta=="TRUE")

exomes_reformat=exomes[,names(exomes) %in% colmap$existing_meta_annotation]
names(exomes_reformat)=colmap$condensed_meta_colname[match(names(exomes_reformat), colmap$existing_meta_annotation)]

#----------------
# Add one-hot-encoding for "exclusion_reasons"

exomes_reformat$non_releasable="False"
exomes_reformat$non_releasable[grep("NR", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$non_releasable[grep("DNE", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$non_releasable[grep("NO_PERM", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$non_releasable=as.factor(exomes_reformat$non_releasable)

exomes_reformat$esp="False"
exomes_reformat$esp[grep("ESP", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$esp=as.factor(exomes_reformat$esp)

exomes_reformat$tcga_below_30="False"
exomes_reformat$tcga_below_30[grep("TCGA_below_30", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$tcga_below_30=as.factor(exomes_reformat$tcga_below_30)

exomes_reformat$tcga_tumor="False"
exomes_reformat$tcga_tumor[grep("TCGA_tumor", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$tcga_tumor=as.factor(exomes_reformat$tcga_tumor)

exomes_reformat$tcga_weird_barcode="False"
exomes_reformat$tcga_weird_barcode[grep("TCGA_weird_barcode", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$tcga_weird_barcode=as.factor(exomes_reformat$tcga_weird_barcode)

exomes_reformat$specific_exclusion="False"
exomes_reformat$specific_exclusion[grep("specific_exclusion", exomes_reformat$exclusion_reasons)]="True"
exomes_reformat$specific_exclusion=as.factor(exomes_reformat$specific_exclusion)

exomes_reformat$syndip="False"
exomes_reformat$syndip[exomes_reformat$s == "CHMI_CHMI3_Nex1"]="True"
exomes_reformat$syndip=as.factor(exomes_reformat$syndip)

# Checks
# table(exomes_reformat$exclusion_reasons, exomes_reformat$non_releasable)
# table(exomes_reformat$exclusion_reasons, exomes_reformat$esp)
# table(exomes_reformat$exclusion_reasons, exomes_reformat$tcga_below_30)
# table(exomes_reformat$exclusion_reasons, exomes_reformat$tcga_tumor)
# table(exomes_reformat$exclusion_reasons, exomes_reformat$tcga_weird_barcode)
# table(exomes_reformat$exclusion_reasons, exomes_reformat$specific_exclusion)
# table(exomes_reformat$exclusion_reasons, exomes_reformat$syndip)
# table(exomes_reformat$exclusion_reasons, exomes_reformat[,"releasable_2_1"])

#--------------------
# Add old release status

old_meta=read.delim(old_releasables_input, header=T)
old_meta=old_meta[,c("sample", "releasable")]
names(old_meta)=c("s", "releasable_2_0_2")
exomes_reformat=merge(exomes_reformat, old_meta, by="s", all.x=T)

table(exomes_reformat[,"releasable_2_0_2"], exomes_reformat[,"releasable_2_1"])
#             NO    YES
#   False  52842   2040
#   True       1 144675

#------------------
# Transform column values

exomes_reformat$internal=as.character(exomes_reformat$internal)
exomes_reformat$internal[exomes_reformat$internal=="internal"]="True"
exomes_reformat$internal[exomes_reformat$internal=="external"]="False"
exomes_reformat$internal=as.factor(exomes_reformat$internal)

exomes_reformat$exac_joint=as.character(exomes_reformat$exac_joint)
exomes_reformat$exac_joint[exomes_reformat$exac_joint=="ExACv1"]="True"
exomes_reformat$exac_joint[exomes_reformat$exac_joint=="ExACv2"]="False"
exomes_reformat$exac_joint=as.factor(exomes_reformat$exac_joint)

exomes_reformat$control=as.character(exomes_reformat$control)
exomes_reformat$control[exomes_reformat$control=="case"]="False"
exomes_reformat$control[exomes_reformat$control=="control"]="True"
exomes_reformat$control=as.factor(exomes_reformat$control)

exomes_reformat$neuro=as.character(exomes_reformat$neuro)
exomes_reformat$neuro[exomes_reformat$neuro=="yes"]="True"
exomes_reformat$neuro[exomes_reformat$neuro=="no"]="False"
exomes_reformat$neuro[exomes_reformat$neuro=="partial" & exomes_reformat$control=="False"]="True"
exomes_reformat$neuro[exomes_reformat$neuro=="partial" & exomes_reformat$control=="True"]="False"
exomes_reformat$neuro[exomes_reformat$neuro=="partial" & is.na(exomes_reformat$control)]="NA"
exomes_reformat$neuro=as.factor(exomes_reformat$neuro)

exomes_reformat[,"releasable_2_1"]=as.character(exomes_reformat[,"releasable_2_1"])
exomes_reformat[,"releasable_2_1"][exomes_reformat[,"releasable_2_1"]=="NO"]="False"
exomes_reformat[,"releasable_2_1"][exomes_reformat[,"releasable_2_1"]=="YES"]="True"
exomes_reformat[,"releasable_2_1"]=as.factor(exomes_reformat[,"releasable_2_1"])


# Drop exclusion_reasons on write:
write.table(exomes_reformat[,names(exomes_reformat)!="exclusion_reasons"], streamlined_meta_output, col.names=T, row.names=F, quote=F, sep="\t")
