import argparse
import logging
from typing import List, Optional, Dict
import hailtop.batch as hb
from hailtop.batch.job import Job
import json


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def _add_split_intervals_job(
        b: hb.Batch,
        utils: Dict,
) -> Job:
    """
    Split genome into intervals to parallelise VQSR for large sample sizes
    :param b: Batch object to add jobs to
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :return: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f"""Make {utils['NUMBER_OF_GENOMICS_DB_INTERVALS']} intervals""")
    j.image(utils['GATK_IMAGE'])
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(utils['NUMBER_OF_GENOMICS_DB_INTERVALS'])
        }
    )

    j.command(
        f"""set -e
    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the NUMBER_OF_GENOMICS_DB_INTERVALS number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
      -L {utils['UNPADDED_INTERVALS']} \\
      -O {j.intervals} \\
      -scatter {utils['NUMBER_OF_GENOMICS_DB_INTERVALS']} \\
      -R {utils['ref_fasta']} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


# SNPs
def SNPsVariantRecalibratorCreateModel(
        b: hb.Batch,
        sites_only_vcf: hb.ResourceGroup,
        hapmap_resource_vcf: hb.ResourceGroup,
        omni_resource_vcf: hb.ResourceGroup,
        one_thousand_genomes_resource_vcf: hb.ResourceGroup,
        dbsnp_resource_vcf: hb.ResourceGroup,
        utils: Dict,
        disk_size: int,
        use_as_annotations: bool,
        transmitted_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        sibling_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        out_bucket: str = None,
        is_small_callset: bool = False,
        is_huge_callset: bool = False,
        max_gaussians: int = 6,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches
    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param hapmap_resource_vcf: HapMap ResourceGroup VCF file to be used in building the model
    :param omni_resource_vcf: OMNI ResourceGroup VCF file to be used in building the model
    :param one_thousand_genomes_resource_vcf: 1KG ResourceGroup VCF file to be used in building the model
    :param dbsnp_resource_vcf: DBSNP ResourceGroup VCF file to be used in building the model
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: If supplied, Transmitted Singletons VCF will be used in building the model
    :param sibling_singletons_resource_vcf: Sibling Singletons VCF will be used in building the model
    :param out_bucket: full path to output bucket to write model and plots to
    :param is_small_callset: whether or not the dataset is small. Used to set number of CPUs for the job
    :param is_huge_callset: whether or not the dataset is huge. Used to set number of CPUs for the job
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.model_file and j.snp_rscript_file.
    """
    j = b.new_job('VQSR: SNPsVariantRecalibratorCreateModel')
    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G
    else:
        ncpu = 16  # ~ 8G/core ~ 128G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage(f'{disk_size}G')

    downsample_factor = 75 if is_huge_callset else 10

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils['SNP_RECALIBRATION_TRANCHE_VALUES']])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
            utils['SNP_RECALIBRATION_ANNOTATION_VALUES_AS']
            if use_as_annotations
            else utils['SNP_RECALIBRATION_ANNOTATION_VALUES']
        )
        ]
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options -Xms{java_mem}g \\
          VariantRecalibrator \\
          -V {sites_only_vcf['vcf.gz']} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          {"--use-allele-specific-annotations" if use_as_annotations else ""} \\
          --sample-every-Nth-variant {downsample_factor} \\
          --output-model {j.model_file} \\
          --max-gaussians {max_gaussians} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf.base}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf.base}' if sibling_singletons_resource_vcf else ''} \\
          --rscript-file {j.snp_rscript}
          ls $(dirname {j.snp_rscript})
          ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
          ln {j.tranches}.pdf {j.tranches_pdf}
          """
    )

    if out_bucket:
        b.write_output(
            j.snp_rscript, f'{out_bucket}model/SNPS/recalibration-snps-features-build.RScript'
        )
        b.write_output(
            j.snp_rscript_pdf, f'{out_bucket}model/SNPS/recalibration-snps-features-build.pdf'
        )
        b.write_output(
            j.tranches_pdf, f'{out_bucket}model/SNPS/recalibration-snps-tranches-build.pdf'
        )
        b.write_output(
            j.model_file, f'{out_bucket}model/SNPS/recalibration-snps-model-file.recal'
        )
    return j


def SNPsVariantRecalibrator(
        b: hb.Batch,
        sites_only_vcf: hb.ResourceGroup,
        hapmap_resource_vcf: hb.ResourceGroup,
        omni_resource_vcf: hb.ResourceGroup,
        one_thousand_genomes_resource_vcf: hb.ResourceGroup,
        dbsnp_resource_vcf: hb.ResourceGroup,
        utils: Dict,
        out_bucket: str,
        disk_size: int,
        use_as_annotations: bool,
        transmitted_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        sibling_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        max_gaussians: int = 6,
) -> Job:
    """
    Recalibrate SNPs in one run (alternative to scatter-gather approach)
    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param hapmap_resource_vcf: HapMap ResourceGroup VCF file to be used in building the model
    :param omni_resource_vcf: OMNI ResourceGroup VCF file to be used in building the model
    :param one_thousand_genomes_resource_vcf: 1KG ResourceGroup VCF file to be used in building the model
    :param dbsnp_resource_vcf: DBSNP ResourceGroup VCF file to be used in building the model
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param out_bucket: full path to output bucket to write model and plots to
    :param disk_size: disk size to be used for the job
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: If supplied, Transmitted Singletons VCF will be used in building the model
    :param sibling_singletons_resource_vcf: Sibling Singletons VCF will be used in building the model
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job('VQSR: SNPsVariantRecalibrator')

    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    ncpu = 8  # ~ 8G/core ~ 64G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 8
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils['SNP_RECALIBRATION_TRANCHE_VALUES']])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
            utils['SNP_RECALIBRATION_ANNOTATION_VALUES_AS']
            if use_as_annotations
            else utils['SNP_RECALIBRATION_ANNOTATION_VALUES']
        )
        ]
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options -Xms{java_mem}g \\
          VariantRecalibrator \\
          -V {sites_only_vcf['vcf.gz']} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
          --max-gaussians {max_gaussians} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf.base}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf.base}' if sibling_singletons_resource_vcf else ''} \\
          --rscript-file {j.snp_rscript}
          ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
          ln {j.tranches}.pdf {j.tranches_pdf}
          """
    )

    if out_bucket:
        b.write_output(
            j.snp_rscript, f'{out_bucket}model/SNPS/recalibration-snps-features-apply.RScript'
        )
        b.write_output(
            j.snp_rscript_pdf, f'{out_bucket}model/SNPS/recalibration-snps-features-apply.pdf'
        )
        b.write_output(
            j.tranches_pdf, f'{out_bucket}model/recalibration-snps-tranches-apply.pdf'
        )
    return j


def SNPsVariantRecalibratorScattered(
        b: hb.Batch,
        sites_only_vcf: hb.ResourceGroup,
        model_file: hb.ResourceFile,
        hapmap_resource_vcf: hb.ResourceGroup,
        omni_resource_vcf: hb.ResourceGroup,
        one_thousand_genomes_resource_vcf: hb.ResourceGroup,
        dbsnp_resource_vcf: hb.ResourceGroup,
        utils: Dict,
        disk_size: int,
        out_bucket: str,
        tranche_idx: int,
        use_as_annotations: bool,
        transmitted_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        sibling_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        interval: Optional[hb.ResourceGroup] = None,
        max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is broken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches
    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param model_file: model file to be applied
    :param hapmap_resource_vcf: HapMap ResourceGroup VCF file to be used in building the model
    :param omni_resource_vcf: OMNI ResourceGroup VCF file to be used in building the model
    :param one_thousand_genomes_resource_vcf: 1KG ResourceGroup VCF file to be used in building the model
    :param dbsnp_resource_vcf: DBSNP ResourceGroup VCF file to be used in building the model
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param out_bucket: full path to output bucket to write model and plots to
    :param tranche_idx: index for the tranches file
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: If supplied, Transmitted Singletons VCF will be used in building the model
    :param sibling_singletons_resource_vcf: Sibling Singletons VCF will be used in building the model
    :param interval: genomic interval to apply the model to
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job('VQSR: SNPsVariantRecalibratorScattered')

    j.image(utils['GATK_IMAGE'])
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils['SNP_RECALIBRATION_TRANCHE_VALUES']])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
            utils['SNP_RECALIBRATION_ANNOTATION_VALUES_AS']
            if use_as_annotations
            else utils['SNP_RECALIBRATION_ANNOTATION_VALUES']
        )
        ]
    )
    j.command(
        f"""set -euo pipefail
        MODEL_REPORT={model_file}
        gatk --java-options -Xms{mem_gb - 1}g \\
          VariantRecalibrator \\
          -V {sites_only_vcf['vcf.gz']} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          {f'-L {interval} ' if interval else ''} \\
          {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
          --input-model {model_file} --output-tranches-for-scatter \\
          --max-gaussians {max_gaussians} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf.base}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf.base}' if sibling_singletons_resource_vcf else ''}"""
    )
    if out_bucket:
        b.write_output(j.tranches, f'{out_bucket}model/SNPS/tranches/recalibration-snps-tranches-{tranche_idx}')
    return j


# INDELs
def IndelsVariantRecalibratorCreateModel(
        b: hb.Batch,
        sites_only_vcf: hb.ResourceGroup,
        mills_resource_vcf: hb.ResourceGroup,
        axiomPoly_resource_vcf: hb.ResourceGroup,
        dbsnp_resource_vcf: hb.ResourceGroup,
        utils: Dict,
        disk_size: int,
        use_as_annotations: bool,
        transmitted_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        sibling_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        out_bucket: str = None,
        is_small_callset: bool = False,
        max_gaussians: int = 4,
) -> Job:
    """
    First step of VQSR for INDELs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.
    To support cohorts with more than 10,000 WGS samples, the INDEL recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches
    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value

    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param mills_resource_vcf: MILLS ResourceGroup VCF file to be used in building the model
    :param axiomPoly_resource_vcf: Axiom Poly ResourceGroup VCF file to be used in building the model
    :param dbsnp_resource_vcf: DBSNP ResourceGroup VCF file to be used in building the model
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: If supplied, Transmitted Singletons VCF will be used in building the model
    :param sibling_singletons_resource_vcf: Sibling Singletons VCF will be used in building the model
    :param out_bucket: full path to output bucket to write model and plots to
    :param is_small_callset: whether or not the dataset is small. Used to set number of CPUs for the job
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.model_file and j.indel_rscript_file.
    The latter is useful to produce the optional tranche plot.
    """
    j = b.new_job('VQSR: INDELsVariantRecalibratorCreateModel')
    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G
    else:
        ncpu = 16  # ~ 8G/core ~ 128G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage(f'{disk_size}G')

    # downsample_factor = 75 if is_huge_callset else 10
    downsample_factor = 10

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils['INDEL_RECALIBRATION_TRANCHE_VALUES']])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
            utils['INDEL_RECALIBRATION_ANNOTATION_VALUES_AS']
            if use_as_annotations
            else utils['INDEL_RECALIBRATION_ANNOTATION_VALUES']
        )
        ]
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options -Xms{java_mem}g \\
          VariantRecalibrator \\
          -V {sites_only_vcf['vcf.gz']} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
          --sample-every-Nth-variant {downsample_factor} \\
          --output-model {j.model_file} \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiomPoly_resource_vcf.base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf.base}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf.base}' if sibling_singletons_resource_vcf else ''} 
          """
    )

    if out_bucket:
        b.write_output(
            j.model_file, f'{out_bucket}model/INDELS/recalibration-indels-model-file.recal'
        )
    return j


def IndelsVariantRecalibrator(
        b: hb.Batch,
        sites_only_vcf: hb.ResourceGroup,
        mills_resource_vcf: hb.ResourceGroup,
        axiomPoly_resource_vcf: hb.ResourceGroup,
        dbsnp_resource_vcf: hb.ResourceGroup,
        utils: Dict,
        disk_size: int,
        use_as_annotations: bool,
        transmitted_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        sibling_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        out_bucket: str = None,
        max_gaussians: int = 4,
) -> Job:
    """
    Run VariantRecalibrator to calculate VQSLOD tranches for indels
    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value. 4 is a
    reasonable default for indels, as their number is smaller than SNPs

    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param mills_resource_vcf: MILLS ResourceGroup VCF file to be used in building the model
    :param axiomPoly_resource_vcf: Axiom Poly ResourceGroup VCF file to be used in building the model
    :param dbsnp_resource_vcf: DBSNP ResourceGroup VCF file to be used in building the model
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: If supplied, Transmitted Singletons VCF will be used in building the model
    :param sibling_singletons_resource_vcf: Sibling Singletons VCF will be used in building the model
    :param out_bucket: full path to output bucket to write model and plots to
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 3 outputs: j.recalibration (ResourceGroup), j.tranches,
    and j.indel_rscript_file. The latter is usedful to produce the optional tranche plot.
    """
    j = b.new_job('VQSR: INDELsVariantRecalibrator')
    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    ncpu = 4  # ~ 8G/core ~ 32G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 4
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join(
        [f'-tranche {v}' for v in utils['INDEL_RECALIBRATION_TRANCHE_VALUES']]
    )
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
            utils['INDEL_RECALIBRATION_ANNOTATION_VALUES_AS']
            if use_as_annotations
            else utils['INDEL_RECALIBRATION_ANNOTATION_VALUES']
        )
        ]
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options -Xms{java_mem}g \\
          VariantRecalibrator \\
          -V {sites_only_vcf['vcf.gz']} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiomPoly_resource_vcf.base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf.base}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf.base}' if sibling_singletons_resource_vcf else ''} \\
          --rscript-file {j.indel_rscript_file}

          ls $(dirname {j.indel_rscript_file})
          ln {j.indel_rscript}.pdf {j.indel_rscript_pdf}
          """
    )
    if out_bucket:
        b.write_output(
            j.indel_rscript_file, f'{out_bucket}model/INDELS/recalibration-indels-features-apply.Rscript'
        )
        b.write_output(
            j.indel_rscript_pdf, f'{out_bucket}model/INDELS/recalibration-indels-features-apply.pdf'
        )

    return j


def IndelsVariantRecalibratorScattered(
        b: hb.Batch,
        sites_only_vcf: hb.ResourceGroup,
        model_file: hb.ResourceGroup,
        mills_resource_vcf: hb.ResourceGroup,
        axiomPoly_resource_vcf: hb.ResourceGroup,
        dbsnp_resource_vcf: hb.ResourceGroup,
        utils: Dict,
        disk_size: int,
        out_bucket: str,
        tranche_idx: int,
        use_as_annotations: bool,
        transmitted_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        sibling_singletons_resource_vcf: Optional[hb.ResourceGroup] = None,
        interval: Optional[hb.ResourceGroup] = None,
        max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for INDELs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches
    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param model_file: model file to be applied
    :param mills_resource_vcf: MILLS ResourceGroup VCF file to be used in building the model
    :param axiomPoly_resource_vcf: Axiom Poly ResourceGroup VCF file to be used in building the model
    :param dbsnp_resource_vcf: DBSNP ResourceGroup VCF file to be used in building the model
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param out_bucket: full path to output bucket to write model and plots to
    :param tranche_idx: index for the tranches file
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: If supplied, Transmitted Singletons VCF will be used in building the model
    :param sibling_singletons_resource_vcf: Sibling Singletons VCF will be used in building the model
    :param interval: genomic interval to apply the model to
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job('VQSR: INDELsVariantRecalibratorScattered')

    j.image(utils.GATK_IMAGE)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in utils['INDEL_RECALIBRATION_TRANCHE_VALUES']])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
            utils['INDEL_RECALIBRATION_ANNOTATION_VALUES_AS']
            if use_as_annotations
            else utils['INDEL_RECALIBRATION_TRANCHE_VALUES']
        )
        ]
    )
    j.command(
        f"""set -euo pipefail
        MODEL_REPORT={model_file}
        gatk --java-options -Xms{mem_gb - 1}g \\
          VariantRecalibrator \\
          -V {sites_only_vcf['vcf.gz']} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          {f'-L {interval} ' if interval else ''} \\
          {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
          --input-model {model_file} --output-tranches-for-scatter \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiomPoly_resource_vcf.base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf.base}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf.base}' if sibling_singletons_resource_vcf else ''}"""
    )
    if out_bucket:
        b.write_output(j.tranches, f'{out_bucket}model/INDELS/tranches/recalibration-indels-tranches-{tranche_idx}')
    return j


# other
def GatherTranches(
        b: hb.Batch,
        tranches: List[hb.ResourceFile],
        mode: str,
        disk_size: int,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval
    tranches outputs.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is broken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches
    Returns: a Job object with one output j.out_tranches.

    :param b: Batch object to add jobs to
    :param tranches: index for the tranches file
    :param mode: Recalibration mode to employ, either SNP or INDEL
    :param disk_size: disk size to be used for the job
    :return: a Job object with one output j.out_tranches
    """
    j = b.new_job(f'VQSR: {mode}GatherTranches')
    j.image('us.gcr.io/broad-dsde-methods/gatk-for-ccdg@sha256:9e9f105ecf3534fbda91a4f2c2816ec3edf775882917813337a8d6e18092c959')
    j.memory('8G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    inputs_cmdl = ' '.join([f'--input {t}' for t in tranches])
    j.command(
        f"""set -euo pipefail
        gatk --java-options -Xms6g \\
          GatherTranches \\
          --mode {mode} \\
          {inputs_cmdl} \\
          --output {j.out_tranches}"""
    )
    return j


def ApplyRecalibration(
        b: hb.Batch,
        input_vcf: hb.ResourceGroup,
        out_vcf_name: str,
        indels_recalibration: hb.ResourceGroup,
        indels_tranches: hb.ResourceFile,
        snps_recalibration: hb.ResourceGroup,
        snps_tranches: hb.ResourceFile,
        utils: Dict,
        disk_size: int,
        use_as_annotations: bool,
        indel_filter_level: float,
        snp_filter_level: float,
        scatter: Optional[int] = None,
        interval: Optional[hb.ResourceGroup] = None,
        out_bucket: Optional[str] = None,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.
    Runs ApplyVQSR twice to apply first indel, then SNP recalibrations.
    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.
    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truthset. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    :param b: Batch object to add jobs to
    :param input_vcf: sites only VCF file to be used
    :param out_vcf_name: output vcf filename
    :param indels_recalibration: input recal file (ResourceGroup) for INDELs
    :param indels_tranches: input tranches file (ResourceFile) for INDELs
    :param snps_recalibration: nput recal file (ResourceGroup) for SNPs
    :param snps_tranches: input tranches file (ResourceFile) for SNPs
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param indel_filter_level: the truth sensitivity level at which to start filtering for INDELs
    :param snp_filter_level: the truth sensitivity level at which to start filtering for SNPs
    :param scatter: scatter index to be used in output VCF filename if running in scattered mode
    :param interval: genomic interval to apply the model to
    :param out_bucket: full path to output bucket to write output(s) to
    :return: a Job object with one ResourceGroup output j.output_vcf, corresponding
    to a VCF with tranche annotated in the FILTER field
    """
    if scatter:
        filename = f'{out_vcf_name}_vqsr_recalibrated_{scatter}'
        outpath = f'{out_bucket}scatter/'
    else:
        filename = f'{out_vcf_name}_vqsr_recalibrated'
        outpath = out_bucket

    j = b.new_job('VQSR: ApplyRecalibration')
    j.image(utils['GATK_IMAGE'])
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': f'{filename}.vcf.gz', 'vcf.gz.tbi': f'{filename}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail
        df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
        gatk --java-options -Xms5g \\
          ApplyVQSR \\
          -O tmp.indel.recalibrated.vcf \\
          -V {input_vcf['vcf.gz']} \\
          --recal-file {indels_recalibration} \\
          --tranches-file {indels_tranches} \\
          --truth-sensitivity-filter-level {indel_filter_level} \\
          --create-output-variant-index true \\
          {f'-L {interval} ' if interval else ''} \\
          {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
          -mode INDEL

        df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
        rm {input_vcf['vcf.gz']} {indels_recalibration} {indels_tranches}
        df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
        gatk --java-options -Xms5g \\
          ApplyVQSR \\
          -O {j.output_vcf['vcf.gz']} \\
          -V tmp.indel.recalibrated.vcf \\
          --recal-file {snps_recalibration} \\
          --tranches-file {snps_tranches} \\
          --truth-sensitivity-filter-level {snp_filter_level} \\
          --create-output-variant-index true \\
          {f'-L {interval} ' if interval else ''} \\
          {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
          -mode SNP
        df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
          """
    )

    if out_bucket:
        b.write_output(j.output_vcf, f'{outpath}{filename}')
    return j


def GatherVcfs(
        b: hb.Batch,
        input_vcfs: List[hb.ResourceGroup],
        out_vcf_name: str,
        utils: Dict,
        disk_size: int,
        out_bucket: str = None,
) -> Job:
    """
    Combines recalibrated VCFs into a single VCF.
    Saves the output VCF to a bucket `out_bucket`

    :param b: Batch object to add jobs to
    :param input_vcfs: list of VCFs to be gathered
    :param out_vcf_name: output vcf filename
    :param utils: a dictionary containing paths to resource files to be used to split genome
    :param disk_size: disk size to be used for the job
    :param out_bucket: full path to output bucket to write the gathered VCF to
    :return: a Job object with one ResourceGroup output j.output_vcf
    """
    filename = f'{out_vcf_name}_vqsr_recalibrated'
    j = b.new_job('VQSR: FinalGatherVcf')
    j.image(utils['GATK_IMAGE'])
    j.memory(f'16G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': f'{filename}.vcf.gz', 'vcf.gz.tbi': f'{filename}vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail
        # --ignore-safety-checks makes a big performance difference so we include it in 
        # our invocation. This argument disables expensive checks that the file headers 
        # contain the same set of genotyped samples and that files are in order 
        # by position of first record.
        gatk --java-options -Xms6g \\
          GatherVcfsCloud \\
          --gather-type BLOCK \\
          {input_cmdl} \\
          --output {j.output_vcf['vcf.gz']}
        tabix {j.output_vcf['vcf.gz']}"""
    )
    if out_bucket:
        b.write_output(j.output_vcf, f'{out_bucket}{filename}')
    return j


def make_vqsr_jobs(
        b: hb.Batch,
        sites_only_vcf: str,
        is_small_callset: bool,
        is_huge_callset: bool,
        output_vcf_name: str,
        utils: Dict,
        out_bucket: str,
        intervals: Dict,
        use_as_annotations: bool,
        transmitted_singletons: Optional[str] = None,
        sibling_singletons: Optional[str] = None,
):
    """
    Add jobs that perform the allele-specific VQSR variant QC
    :param b: Batch object to add jobs to
    :param sites_only_vcf: path to a sites only VCF created using gnomAD default_compute_info()
    :param is_small_callset: for small callsets, we gather the VCF shards and collect
        QC metrics directly. For anything larger, we need to keep the VCF sharded and
        gather metrics collected from them
    :param is_huge_callset: For huge callsets, we allocate more memory for the SNPs
        Create Model step
    :param output_vcf_name: name, without extension, to use for the output VCF file(s)
    :param utils: a dictionary containing resource files and parameters to be used in VQSR
    :param out_bucket: path to write, plots, evaluation results, and recalibrated VCF to
    :param intervals: ResourceGroup object with intervals to scatter
    :param use_as_annotations: use allele-specific annotation for VQSR
    :param transmitted_singletons: full path to transmitted singletons VCF file and its index
    :param sibling_singletons: full path to sibling singletons VCF file and its index
    :return: a final Job, and a path to the VCF with VQSR annotations
    """
    # Reference files. All options have defaults.
    dbsnp_vcf = b.read_input_group(base=utils['dbsnp_vcf'], index=utils['dbsnp_vcf_index'])
    hapmap_resource_vcf = b.read_input_group(
        base=utils['hapmap_resource_vcf'], index=utils['hapmap_resource_vcf_index']
    )
    omni_resource_vcf = b.read_input_group(
        base=utils['omni_resource_vcf'], index=utils['omni_resource_vcf_index']
    )
    one_thousand_genomes_resource_vcf = b.read_input_group(
        base=utils['one_thousand_genomes_resource_vcf'],
        index=utils['one_thousand_genomes_resource_vcf_index'],
    )
    mills_resource_vcf = b.read_input_group(
        base=utils['mills_resource_vcf'], index=utils['mills_resource_vcf_index']
    )
    axiom_poly_resource_vcf = b.read_input_group(
        base=utils['axiom_poly_resource_vcf'], index=utils['axiom_poly_resource_vcf_index']
    )
    dbsnp_resource_vcf = dbsnp_vcf

    # To fit only a sites-only VCF
    if is_small_callset:
        small_disk = 50
    elif not is_huge_callset:
        small_disk = 100
    else:
        small_disk = 200

    if is_small_callset:
        huge_disk = 200
    elif not is_huge_callset:
        huge_disk = 500
    else:
        huge_disk = 2000

    gathered_vcf = b.read_input_group(
        **{
            'vcf.gz': sites_only_vcf,
            'vcf.gz.tbi': sites_only_vcf + '.tbi',
        }
    )

    if transmitted_singletons:
        transmitted_singletons_resource_vcf = b.read_input_group(
            base=transmitted_singletons, index=f'{transmitted_singletons}.tbi'
        )
    else:
        transmitted_singletons_resource_vcf = None

    if sibling_singletons:
        sibling_singletons_resource_vcf = b.read_input_group(
            base=sibling_singletons, index=f'{sibling_singletons}.tbi'
        )
    else:
        sibling_singletons_resource_vcf = None

    snp_max_gaussians = 6
    indel_max_gaussians = 4

    if is_huge_callset:
        # 1. Run SNP recalibrator in a scattered mode
        # file exists:
        if hl.hadoop_exists(f'{out_bucket}model/SNPS/recalibration-snps-model-file.recal'):
            snps_model_file = b.read_input(f'{out_bucket}model/SNPS/recalibration-snps-model-file.recal')
        else:
            snps_model_file = SNPsVariantRecalibratorCreateModel(
                b=b,
                sites_only_vcf=gathered_vcf,
                transmitted_singletons_resource_vcf=transmitted_singletons_resource_vcf,
                sibling_singletons_resource_vcf=sibling_singletons_resource_vcf,
                hapmap_resource_vcf=hapmap_resource_vcf,
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                utils=utils,
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                out_bucket=out_bucket,
                is_small_callset=is_small_callset,
                is_huge_callset=is_huge_callset,
                max_gaussians=snp_max_gaussians,
            ).model_file

        snps_recalibrator_jobs = [
            SNPsVariantRecalibratorScattered(
                b=b,
                sites_only_vcf=gathered_vcf,
                interval=intervals[f'interval_{idx}'],
                model_file=snps_model_file,
                hapmap_resource_vcf=hapmap_resource_vcf,
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                utils=utils,
                disk_size=small_disk,
                out_bucket=out_bucket,
                tranche_idx=idx,
                use_as_annotations=use_as_annotations,
                transmitted_singletons_resource_vcf=transmitted_singletons_resource_vcf,
                sibling_singletons_resource_vcf=sibling_singletons_resource_vcf,
                max_gaussians=snp_max_gaussians,
            )
            for idx in range(utils['NUMBER_OF_GENOMICS_DB_INTERVALS'])
        ]

        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        snps_gathered_tranches = GatherTranches(
            b=b,
            tranches=snps_tranches,
            mode='SNP',
            disk_size=small_disk,
        ).out_tranches

        # 2. Run INDEL recalibrator in a scattered mode
        if hl.hadoop_exists(f'{out_bucket}model/INDELS/recalibration-indels-model-file.recal'):
            indels_model_file = b.read_input(f'{out_bucket}model/INDELS/recalibration-indels-model-file.recal')
        else:
            indels_model_file = IndelsVariantRecalibratorCreateModel(
                b=b,
                sites_only_vcf=gathered_vcf,
                transmitted_singletons_resource_vcf=transmitted_singletons_resource_vcf,
                sibling_singletons_resource_vcf=sibling_singletons_resource_vcf,
                mills_resource_vcf=hapmap_resource_vcf,
                axiomPoly_resource_vcf=axiom_poly_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                utils=utils,
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                out_bucket=out_bucket,
                is_small_callset=is_small_callset,
                max_gaussians=indel_max_gaussians,
            ).model_file

        indels_recalibrator_jobs = [
            IndelsVariantRecalibratorScattered(
                b=b,
                sites_only_vcf=gathered_vcf,
                interval=intervals[f'interval_{idx}'],
                model_file=indels_model_file,
                mills_resource_vcf=hapmap_resource_vcf,
                axiomPoly_resource_vcf=axiom_poly_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                utils=utils,
                disk_size=small_disk,
                out_bucket=out_bucket,
                tranche_idx=idx,
                use_as_annotations=use_as_annotations,
                transmitted_singletons_resource_vcf=transmitted_singletons_resource_vcf,
                sibling_singletons_resource_vcf=sibling_singletons_resource_vcf,
                max_gaussians=indel_max_gaussians,
            )
            for idx in range(utils['NUMBER_OF_GENOMICS_DB_INTERVALS'])
        ]

        indels_recalibrations = [j.recalibration for j in indels_recalibrator_jobs]
        indels_tranches = [j.tranches for j in indels_recalibrator_jobs]
        indels_gathered_tranches = GatherTranches(
            b=b,
            tranches=indels_tranches,
            mode='INDEL',
            disk_size=small_disk,
        ).out_tranches

        # 3. Apply recalibration
        scattered_vcfs = [
            ApplyRecalibration(
                b=b,
                input_vcf=gathered_vcf,
                out_vcf_name=output_vcf_name,
                indels_recalibration=indels_recalibrations[idx],
                indels_tranches=indels_gathered_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                utils=utils,
                disk_size=huge_disk,
                use_as_annotations=use_as_annotations,
                indel_filter_level=utils.INDEL_HARD_FILTER_LEVEL,
                snp_filter_level=utils.SNP_HARD_FILTER_LEVEL,
                scatter=idx,
                interval=intervals[f'interval_{idx}'],
                out_bucket = out_bucket,
            ).output_vcf
            for idx in range(utils['NUMBER_OF_GENOMICS_DB_INTERVALS'])
        ]

        # 4. Gather VCFs
        recalibrated_gathered_vcf_job = GatherVcfs(
            b=b,
            input_vcfs=scattered_vcfs,
            out_vcf_name=output_vcf_name,
            utils=utils,
            disk_size=huge_disk,
            out_bucket=out_bucket,
        )

    else:
        snps_recalibrator_job = SNPsVariantRecalibrator(
            b=b,
            sites_only_vcf=gathered_vcf,
            transmitted_singletons_resource_vcf=transmitted_singletons_resource_vcf,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            utils=utils,
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            max_gaussians=snp_max_gaussians,
            out_bucket=out_bucket
        )
        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        indels_variant_recalibrator_job = IndelsVariantRecalibrator(
            b=b,
            sites_only_vcf=gathered_vcf,
            transmitted_singletons_resource_vcf=transmitted_singletons_resource_vcf,
            mills_resource_vcf=mills_resource_vcf,
            axiomPoly_resource_vcf=axiom_poly_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            utils=utils,
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            max_gaussians=indel_max_gaussians,
            out_bucket=out_bucket,
        )
        indels_recalibration = indels_variant_recalibrator_job.recalibration
        indels_tranches = indels_variant_recalibrator_job.tranches

        recalibrated_gathered_vcf_job = ApplyRecalibration(
            b=b,
            input_vcf=gathered_vcf,
            out_vcf_name=output_vcf_name,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            utils=utils,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            indel_filter_level=utils.INDEL_HARD_FILTER_LEVEL,
            snp_filter_level=utils.SNP_HARD_FILTER_LEVEL,
            out_bucket=out_bucket,
        )

    # return recalibrated_gathered_vcf_job


# VQSR workflow including scatter step
def vqsr_workflow(
        sites_only_vcf: str,
        output_vcf_filename: str,
        transmitted_singletons: str,
        sibling_singletons: str,
        resources: str,
        out_bucket: str,
        billing_project: str,
        n_samples: int,
        use_as_annotations: bool = True
):
    """
    Wraps all the functions into a workflow
    :param sites_only_vcf: path to a sites only VCF created using gnomAD default_compute_info()
    :param output_vcf_filename: name, without extension, to use for the output VCF file(s)
    :param transmitted_singletons: full path to transmitted singletons VCF file and its index
    :param sibling_singletons: full path to sibling singletons VCF file and its index
    :param resources: json file (vqsr_resources.json) with paths to the resource files and parameters to be used in VQSR
    :param out_bucket: path to write, plots, evaluation results, and recalibrated VCF to
    :param billing_project: billing project to be used for the workflow
    :param n_samples: number of samples in the VCF to determine if we should run in scatter or full VCF mode
    :param use_as_annotations: use allele-specific annotation for VQSR

    :return:
    """

    tmp_vqsr_bucket = f'{out_bucket}/vqsr/'

    with open(resources, 'r') as f:
        utils = json.load(f)

    logger.info(
        f'Starting hail Batch with the project {billing_project}, '
        f'bucket {tmp_vqsr_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=tmp_vqsr_bucket,
    )

    b = hb.Batch(
        'VQSR pipeline',
        backend=backend,
    )

    intervals_j = _add_split_intervals_job(
        b=b,
        utils=utils
    )

    is_small_callset = n_samples < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = n_samples >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    make_vqsr_jobs(
        b=b,
        sites_only_vcf=sites_only_vcf,
        is_small_callset=is_small_callset,
        is_huge_callset=is_huge_callset,
        output_vcf_name=output_vcf_filename,
        utils=utils,
        out_bucket=tmp_vqsr_bucket,
        intervals=intervals_j.intervals,
        use_as_annotations=use_as_annotations,
        transmitted_singletons=transmitted_singletons,
        sibling_singletons=sibling_singletons)

    b.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcf', type=str, required=True)
    parser.add_argument('--out-bucket', type=str, required=True)
    parser.add_argument('--out-vcf-name', type=str, required=True)
    parser.add_argument('--resources', type=str, required=True)
    parser.add_argument('--n-samples', type=int, required=True)
    parser.add_argument('--transmitted-singletons', type=str, required=False)
    parser.add_argument('--sibling-singletons', type=str, required=False)

    args = parser.parse_args()

    vqsr_workflow(sites_only_vcf=args.input_vcf,
                  output_vcf_filename=args.out_vcf_name,
                  transmitted_singletons=args.transmitted_singletons,
                  sibling_singletons=args.sibling_singletons,
                  resources=args.resources,
                  out_bucket=args.out_bucket,
                  billing_project=args.billing_project,
                  n_samples=args.n_samples)


if __name__ == '__main__':
    main()
