"""
Script to run VQSR on an AS-annotated Sites VCF
Example input:
python3 ~/Documents/GitHub/gnomad_methods/gnomad/variant_qc/vqsr.py \\
    --input-vcf gs://gnomad-tmp/gnomad_v4.0_testing/annotations/exomes/gnomad.exomes.v4.0.info.AS.chr22.vcf.bgz \\
    --out-bucket gs://gnomad-marten/vqsr_20230921/results --out-vcf-name chr22tester \\
    --resources ~/Documents/GitHub/gnomad_methods/gnomad/variant_qc/vqsr_resources.json --billing-project marten-trial \\
    --transmitted-singletons gs://gnomad-tmp/gnomad_v4.0_testing/annotations/exomes/gnomad.exomes.v4.0.transmitted_singleton.raw.chr22.vcf.bgz \\
    --batch-suffix chr22testlarge --n-samples 700000
"""
import argparse
import json
import logging
from typing import Dict, List, Optional

import hail as hl
import hailtop.batch as hb
from hailtop.batch.job import Job

logger = logging.getLogger(__file__)
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger.setLevel(logging.INFO)


def split_intervals(
    b: hb.Batch,
    utils: Dict,
    gcp_billing_project: str,
) -> Job:
    """
    Split genome into intervals to parallelize VQSR for large sample sizes
    :param b: Batch object to add jobs to
    :param utils: a dictionary containing resources (file paths and arguments) to be used to split genome
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :return: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f"""Make {utils['NUMBER_OF_GENOMICS_DB_INTERVALS']} intervals""")
    j.image(utils["GATK_IMAGE"])
    java_mem = 3
    j.memory("standard")  # ~ 4G/core ~ 4G
    j.storage("16G")
    j.declare_resource_group(
        intervals={
            f"interval_{idx}": f"{{root}}/{str(idx).zfill(4)}-scattered.interval_list"
            for idx in range(utils["NUMBER_OF_GENOMICS_DB_INTERVALS"])
        }
    )

    j.command(f"""set -e
    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number
    # of intervals. But we have to expect exactly the NUMBER_OF_GENOMICS_DB_INTERVALS_CHR22 number of
    # output files because our workflow is not dynamic.
    gatk --java-options "-Xms{java_mem}g" SplitIntervals \\
      -L {utils['CALLING_INTERVALS']} \\
      --interval-padding 150 \\
      -O {j.intervals} \\
      -scatter {utils['NUMBER_OF_GENOMICS_DB_INTERVALS']} \\
      -R {utils['ref_fasta']} \\
      --gcs-project-for-requester-pays {gcp_billing_project} \\
      -mode INTERVAL_SUBDIVISION
      """)
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


# SNPs
def snps_variant_recalibrator_create_model(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    use_as_annotations: bool,
    gcp_billing_project: str,
    transmitted_singletons_resource_vcf: Optional[str] = None,
    sibling_singletons_resource_vcf: Optional[str] = None,
    out_bucket: str = None,
    is_small_callset: bool = False,
    is_large_callset: bool = False,
    max_gaussians: int = 6,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
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
    :param utils: a dictionary containing resources (file paths and arguments) to be used to create the model
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :param transmitted_singletons_resource_vcf: Optional transmitted singletons VCF to be used in building the model
    :param sibling_singletons_resource_vcf: Optional sibling singletons VCF to be used in building the model
    :param out_bucket: full path to output bucket to write model and plots to
    :param is_small_callset: whether the dataset is small. Used to set number of CPUs for the job
    :param is_large_callset: whether the dataset is huge. Used to set number of CPUs for the job
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.model_file and j.snp_rscript_file.
    """
    j = b.new_job("VQSR: SNPsVariantRecalibratorCreateModel")
    j.image(utils["GATK_IMAGE"])
    j.memory("highmem")
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G
    else:
        ncpu = 16  # ~ 8G/core ~ 128G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage("50G")

    downsample_factor = 75  # if is_large_callset else 10

    tranche_cmdl = " ".join(
        [f"-tranche {v}" for v in utils["SNP_RECALIBRATION_TRANCHE_VALUES"]]
    )
    an_cmdl = " ".join(
        [
            f"-an {v}"
            for v in (
                utils["SNP_RECALIBRATION_ANNOTATION_VALUES_AS"]
                if use_as_annotations
                else utils["SNP_RECALIBRATION_ANNOTATION_VALUES"]
            )
        ]
    )

    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms{java_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads={ncpu-2}" \\
          VariantRecalibrator \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          -L {utils['EVALUATION_INTERVALS']} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          {"--use-allele-specific-annotations" if use_as_annotations else ""} \\
          --sample-every-Nth-variant {downsample_factor} \\
          --output-model {j.model_file} \\
          --max-gaussians {max_gaussians} \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {utils['hapmap_resource_vcf']} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {utils['omni_resource_vcf']} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {utils['one_thousand_genomes_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {utils['dbsnp_resource_vcf']} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf}' if sibling_singletons_resource_vcf else ''} \\
          --rscript-file {j.snp_rscript}
          ls $(dirname {j.snp_rscript})
          ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
          ln {j.tranches}.pdf {j.tranches_pdf}
          """)

    if out_bucket:
        b.write_output(
            j.snp_rscript, f"{out_bucket}model/SNPS/snps.features.build.RScript"
        )
        b.write_output(
            j.snp_rscript_pdf, f"{out_bucket}model/SNPS/snps.features.build.pdf"
        )
        b.write_output(
            j.tranches_pdf, f"{out_bucket}model/SNPS/snps.tranches.build.pdf"
        )
        b.write_output(j.model_file, f"{out_bucket}model/SNPS/snps.model.report")

    return j


def snps_variant_recalibrator(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    out_bucket: str,
    use_as_annotations: bool,
    gcp_billing_project,
    interval: Optional[hb.ResourceGroup] = None,
    tranche_idx: Optional[int] = None,
    model_file: Optional[hb.ResourceFile] = None,
    transmitted_singletons_resource_vcf: Optional[str] = None,
    sibling_singletons_resource_vcf: Optional[str] = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
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
    :param utils: a dictionary containing resources (file paths and arguments)
    :param out_bucket: full path to output bucket to write model and plots to
    :param tranche_idx: index for the tranches file
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param transmitted_singletons_resource_vcf: Optional transmitted singletons VCF to include in VariantRecalibrator
    :param sibling_singletons_resource_vcf: Optional sibling singletons VCF to include in VariantRecalibrator
    :param interval: genomic interval to apply the model to
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job("VQSR: SNPsVariantRecalibratorScattered")

    j.image(utils["GATK_IMAGE"])
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f"{mem_gb}G")
    j.cpu(4)
    j.storage("20G")

    j.declare_resource_group(recalibration={"index": "{root}.idx", "base": "{root}"})

    tranche_cmdl = " ".join(
        [f"-tranche {v}" for v in utils["SNP_RECALIBRATION_TRANCHE_VALUES"]]
    )
    an_cmdl = " ".join(
        [
            f"-an {v}"
            for v in (
                utils["SNP_RECALIBRATION_ANNOTATION_VALUES_AS"]
                if use_as_annotations
                else utils["SNP_RECALIBRATION_ANNOTATION_VALUES"]
            )
        ]
    )

    cmd = f"""set -euo pipefail
        gatk --java-options "-Xms{mem_gb-1}g -XX:+UseParallelGC -XX:ParallelGCThreads=3" \\
          VariantRecalibrator \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          --max-gaussians {max_gaussians} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {utils['hapmap_resource_vcf']} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {utils['omni_resource_vcf']} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {utils['one_thousand_genomes_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {utils['dbsnp_resource_vcf']} \\
        """

    if interval:
        cmd += f" -L {interval}"
    if use_as_annotations:
        cmd += " --use-allele-specific-annotations"
    if model_file:
        cmd += f" --input-model {model_file} --output-tranches-for-scatter"
    if transmitted_singletons_resource_vcf:
        cmd += (
            " -resource:singletons,known=true,training=true,truth=true,prior=10"
            f" {transmitted_singletons_resource_vcf}"
        )
    if sibling_singletons_resource_vcf:
        cmd += (
            " -resource:singletons,known=true,training=true,truth=true,prior=10"
            f" {sibling_singletons_resource_vcf}"
        )

    j.command(cmd)

    if out_bucket:
        if tranche_idx is not None:
            b.write_output(
                j.tranches,
                f"{out_bucket}recalibration/SNPS/snps.{tranche_idx}.tranches",
            )
        else:
            b.write_output(j.tranches, f"{out_bucket}recalibration/SNPS/snps.tranches")

    return j


# INDELs
def indels_variant_recalibrator_create_model(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    use_as_annotations: bool,
    gcp_billing_project: str,
    transmitted_singletons_resource_vcf: str = None,
    sibling_singletons_resource_vcf: str = None,
    out_bucket: str = None,
    is_small_callset: bool = False,
    max_gaussians: int = 4,
) -> Job:
    """
    First step of VQSR for INDELs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.
    To support cohorts with more than 10,000 WGS samples, the INDEL recalibration process
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
    message. In this case, try decrementing the --max-gaussians value

    :param b: Batch object to add jobs to
    :param sites_only_vcf: sites only VCF file to be used to build the model
    :param utils: a dictionary containing resources (file paths and arguments)
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :param transmitted_singletons_resource_vcf: Optional transmitted singletons VCF to be used in building the model
    :param sibling_singletons_resource_vcf: Optional sibling singletons VCF to be used in building the model
    :param out_bucket: full path to output bucket to write model and plots to
    :param is_small_callset: whether the dataset is small. Used to set number of CPUs for the job
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.model_file and j.indel_rscript_file.
    The latter is useful to produce the optional tranche plot.
    """
    j = b.new_job("VQSR: INDELsVariantRecalibratorCreateModel")
    j.image(utils["GATK_IMAGE"])
    j.memory("highmem")
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G
    else:
        ncpu = 16  # ~ 8G/core ~ 128G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage("50G")

    downsample_factor = 75 if not is_small_callset else 10

    tranche_cmdl = " ".join(
        [f"-tranche {v}" for v in utils["INDEL_RECALIBRATION_TRANCHE_VALUES"]]
    )
    an_cmdl = " ".join(
        [
            f"-an {v}"
            for v in (
                utils["INDEL_RECALIBRATION_ANNOTATION_VALUES_AS"]
                if use_as_annotations
                else utils["INDEL_RECALIBRATION_ANNOTATION_VALUES"]
            )
        ]
    )

    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms{java_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads={ncpu-2}" \\
          VariantRecalibrator \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          -L {utils['EVALUATION_INTERVALS']} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          {"--use-allele-specific-annotations" if use_as_annotations else ""} \\
          --sample-every-Nth-variant {downsample_factor} \\
          --output-model {j.model_file} \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {utils['mills_resource_vcf']} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {utils['axiom_poly_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {utils['dbsnp_resource_vcf']} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {transmitted_singletons_resource_vcf}' if transmitted_singletons_resource_vcf else ''} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {sibling_singletons_resource_vcf}' if sibling_singletons_resource_vcf else ''} \\
          --rscript-file {j.indel_rscript}
          ls $(dirname {j.indel_rscript})
          ln {j.indel_rscript}.pdf {j.indel_rscript_pdf}
        """)

    if out_bucket:
        b.write_output(
            j.indel_rscript, f"{out_bucket}model/INDELS/indels.features.build.RScript"
        )
        b.write_output(
            j.indel_rscript_pdf, f"{out_bucket}model/INDELS/indels.features.build.pdf"
        )
        b.write_output(j.model_file, f"{out_bucket}model/INDELS/indels.model.report")

    return j


def indels_variant_recalibrator(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    out_bucket: str,
    use_as_annotations: bool,
    gcp_billing_project: str,
    interval: Optional[hb.ResourceGroup] = None,
    tranche_idx: Optional[int] = None,
    model_file: Optional[hb.ResourceFile] = None,
    transmitted_singletons_resource_vcf: Optional[str] = None,
    sibling_singletons_resource_vcf: Optional[str] = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for INDELs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
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
    :param utils: a dictionary containing resources (file paths and arguments)
    :param out_bucket: full path to output bucket to write model and plots to
    :param tranche_idx: index for the tranches file
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :param transmitted_singletons_resource_vcf: Optional transmitted singletons VCF to include in VariantRecalibrator
    :param sibling_singletons_resource_vcf: Optional sibling singletons VCF to include in VariantRecalibrator
    :param interval: genomic interval to apply the model to
    :param max_gaussians: maximum number of Gaussians for the positive model
    :return: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job("VQSR: INDELsVariantRecalibratorScattered")

    j.image(utils["GATK_IMAGE"])
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f"{mem_gb}G")
    j.cpu(4)
    j.storage("20G")

    j.declare_resource_group(recalibration={"index": "{root}.idx", "base": "{root}"})

    tranche_cmdl = " ".join(
        [f"-tranche {v}" for v in utils["INDEL_RECALIBRATION_TRANCHE_VALUES"]]
    )
    an_cmdl = " ".join(
        [
            f"-an {v}"
            for v in (
                utils["INDEL_RECALIBRATION_ANNOTATION_VALUES_AS"]
                if use_as_annotations
                else utils["INDEL_RECALIBRATION_TRANCHE_VALUES"]
            )
        ]
    )

    cmd = f"""set -euo pipefail
        gatk --java-options "-Xms{mem_gb-1}g -XX:+UseParallelGC -XX:ParallelGCThreads=3" \\
          VariantRecalibrator \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {utils['mills_resource_vcf']} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {utils['axiom_poly_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {utils['dbsnp_resource_vcf']} \\
        """

    if interval:
        cmd += f" -L {interval}"
    if use_as_annotations:
        cmd += " --use-allele-specific-annotations"
    if model_file:
        cmd += f" --input-model {model_file} --output-tranches-for-scatter"
    if transmitted_singletons_resource_vcf:
        cmd += (
            " -resource:singletons,known=true,training=true,truth=true,prior=10"
            f" {transmitted_singletons_resource_vcf}"
        )
    if sibling_singletons_resource_vcf:
        cmd += (
            " -resource:singletons,known=true,training=true,truth=true,prior=10"
            f" {sibling_singletons_resource_vcf}"
        )

    j.command(cmd)

    if out_bucket:
        if tranche_idx is not None:
            b.write_output(
                j.tranches,
                f"{out_bucket}recalibration/INDELS/indels.{tranche_idx}.tranches",
            )
        else:
            b.write_output(
                j.tranches, f"{out_bucket}recalibration/INDELS/indels.tranches"
            )
    return j


# other
def gather_tranches(
    b: hb.Batch,
    tranches: List[hb.ResourceFile],
    mode: str,
    disk_size: int,
    gcp_billing_project: str,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval
    tranches outputs.
    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
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
    j = b.new_job(f"VQSR: {mode}GatherTranches")
    j.image(
        "us.gcr.io/broad-dsde-methods/gatk-for-ccdg@sha256:9e9f105ecf3534fbda91a4f2c2816ec3edf775882917813337a8d6e18092c959"
    )
    j.memory("8G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    inputs_cmdl = " ".join([f"--input {t}" for t in tranches])
    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms6g" \\
          GatherTranches \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --mode {mode} \\
          {inputs_cmdl} \\
          --output {j.out_tranches}""")

    return j


def apply_recalibration(
    b: hb.Batch,
    input_vcf: str,
    out_vcf_name: str,
    indels_recalibration: hb.ResourceGroup,
    indels_tranches: hb.ResourceFile,
    snps_recalibration: hb.ResourceGroup,
    snps_tranches: hb.ResourceFile,
    utils: Dict,
    disk_size: int,
    use_as_annotations: bool,
    gcp_billing_project: str,
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
    :param snps_recalibration: input recal file (ResourceGroup) for SNPs
    :param snps_tranches: input tranches file (ResourceFile) for SNPs
    :param utils: a dictionary containing resources (file paths and arguments)
    :param disk_size: disk size to be used for the job
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used
    :param scatter: scatter index to be used in output VCF filename if running in scattered mode
    :param interval: genomic interval to apply the model to
    :param out_bucket: full path to output bucket to write output(s) to
    :return: a Job object with one ResourceGroup output j.output_vcf, corresponding
    to a VCF with tranche annotated in the FILTER field
    """
    if scatter is not None:
        filename = f"{out_vcf_name}_vqsr_recalibrated_{scatter}"
        outpath = f"{out_bucket}apply_recalibration/scatter/"
    else:
        filename = f"{out_vcf_name}_vqsr_recalibrated"
        outpath = out_bucket

    j = b.new_job("VQSR: ApplyRecalibration")
    # couldn't find a public image with both gatk and bcftools installed
    j.image("docker.io/lindonkambule/vqsr_gatk_bcftools_img:latest")
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(
        output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"},
        intermediate_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"},
    )

    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms5g" \\
          ApplyVQSR \\
          -O tmp.indel.recalibrated.vcf \\
          -V {input_vcf} \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --recal-file {indels_recalibration} \\
          --tranches-file {indels_tranches} \\
          --truth-sensitivity-filter-level {utils['INDEL_HARD_FILTER_LEVEL']} \\
          --create-output-variant-index true \\
          {f'-L {interval} ' if interval else ''} \\
          {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
          -mode INDEL

        gatk --java-options "-Xms5g" \\
          ApplyVQSR \\
          -O {j.output_vcf['vcf.gz']} \\
          -V tmp.indel.recalibrated.vcf \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --recal-file {snps_recalibration} \\
          --tranches-file {snps_tranches} \\
          --truth-sensitivity-filter-level {utils['SNP_HARD_FILTER_LEVEL']} \\
          --create-output-variant-index true \\
          {f'-L {interval} ' if interval else ''} \\
          {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
          -mode SNP""")

    # An INDEL at the beginning of a chunk will overlap with the previous chunk and will cause issues when trying to
    # merge. This makes sure the INDEL is ONLY in ONE of two consecutive
    # chunks (not both)
    if interval:
        j.command(
            f"""bcftools query -f '%CHROM\n' {j.output_vcf['vcf.gz']} | cut -f1 | uniq -c """
        )
        j.command(f"""
            interval=$(cat {interval} | tail -n1 | awk '{{print $1":"$2"-"$3}}')
            echo $interval
        """)
        #  overwrite VCF with overlap issue addressed
        j.command(f"""
                interval=$(cat {interval} | tail -n1 | awk '{{print $1":"$2"-"$3}}')
                bcftools view -t $interval {j.output_vcf['vcf.gz']} --output-file {j.output_vcf['vcf.gz']} --output-type z
                tabix -f {j.output_vcf['vcf.gz']}
                df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
            """)

    if out_bucket:
        b.write_output(j.output_vcf, f"{outpath}{filename}")

    return j


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    out_vcf_name: str,
    utils: Dict,
    disk_size: int,
    gcp_billing_project: str,
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
    filename = f"{out_vcf_name}_vqsr_recalibrated"
    j = b.new_job("VQSR: FinalGatherVcf")
    j.image(utils["GATK_IMAGE"])
    j.memory(f"16G")
    j.storage(f"{disk_size}G")

    j.declare_resource_group(
        output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"}
    )

    input_cmdl = " ".join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(f"""set -euo pipefail
        # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
        # This argument disables expensive checks that the file headers contain the same set of
        # genotyped samples and that files are in order by position of first record.
        cd /io
        mkdir tmp/
        gatk --java-options "-Xms6g -Djava.io.tmpdir=`pwd`/tmp" \\
          GatherVcfsCloud \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --ignore-safety-checks \\
          --gather-type BLOCK \\
          {input_cmdl} \\
          --output {j.output_vcf['vcf.gz']} \\
          --tmp-dir `pwd`/tmp

        tabix {j.output_vcf['vcf.gz']}""")
    if out_bucket:
        b.write_output(j.output_vcf, f"{out_bucket}{filename}")
    return j


def make_vqsr_jobs(
    b: hb.Batch,
    sites_only_vcf: str,
    is_small_callset: bool,
    is_large_callset: bool,
    output_vcf_name: str,
    utils: Dict,
    out_bucket: str,
    intervals: Dict,
    use_as_annotations: bool,
    gcp_billing_project: str,
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
    :param is_large_callset: For huge callsets, we allocate more memory for the SNPs
        Create Model step
    :param output_vcf_name: name, without extension, to use for the output VCF file(s)
    :param utils: a dictionary containing resource files and parameters to be used in VQSR
    :param out_bucket: path to write, plots, evaluation results, and recalibrated VCF to
    :param intervals: ResourceGroup object with intervals to scatter
    :param use_as_annotations: use allele-specific annotation for VQSR
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :param transmitted_singletons: full path to transmitted singletons VCF file and its index
    :param sibling_singletons: full path to sibling singletons VCF file and its index
    :return: a final Job, and a path to the VCF with VQSR annotations
    """
    # To fit only a sites-only VCF
    if is_small_callset:
        small_disk = 50
    elif not is_large_callset:
        small_disk = 100
    else:
        small_disk = 200

    if is_small_callset:
        huge_disk = 200
    elif not is_large_callset:
        huge_disk = 500
    else:
        huge_disk = 2000

    snp_max_gaussians = 6
    indel_max_gaussians = 4

    # Iif it is a large callset, run in scatter mode
    if is_large_callset:
        # 1. Run SNP recalibrator in a scattered mode
        # file exists:
        if hl.hadoop_exists(f"{out_bucket}model/SNPS/snps.model.report"):
            print(
                "Found existing model for SNPs:"
                f" {out_bucket}model/SNPS/snps.model.report"
            )
            snps_model_file = b.read_input(f"{out_bucket}model/SNPS/snps.model.report")
        else:
            snps_model_file = snps_variant_recalibrator_create_model(
                b=b,
                sites_only_vcf=sites_only_vcf,
                transmitted_singletons_resource_vcf=transmitted_singletons,
                sibling_singletons_resource_vcf=sibling_singletons,
                utils=utils,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                out_bucket=out_bucket,
                is_small_callset=is_small_callset,
                is_large_callset=is_large_callset,
                max_gaussians=snp_max_gaussians,
            ).model_file

        snps_recalibrator_jobs = [
            snps_variant_recalibrator(
                b=b,
                sites_only_vcf=sites_only_vcf,
                utils=utils,
                out_bucket=out_bucket,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                interval=intervals[f"interval_{idx}"],
                tranche_idx=idx,
                model_file=snps_model_file,
                transmitted_singletons_resource_vcf=transmitted_singletons,
                sibling_singletons_resource_vcf=sibling_singletons,
                max_gaussians=snp_max_gaussians,
            )
            for idx in range(utils["NUMBER_OF_GENOMICS_DB_INTERVALS"])
        ]

        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        snps_gathered_tranches = gather_tranches(
            b=b,
            tranches=snps_tranches,
            mode="SNP",
            disk_size=small_disk,
            gcp_billing_project=gcp_billing_project,
        ).out_tranches

        # 2. Run INDEL recalibrator in a scattered mode
        if hl.hadoop_exists(f"{out_bucket}model/INDELS/indels.model.report"):
            print(
                "Found existing model for INDELs:"
                f" {out_bucket}model/INDELS/indels.model.report"
            )
            indels_model_file = b.read_input(
                f"{out_bucket}model/INDELS/indels.model.report"
            )
        else:
            indels_model_file = indels_variant_recalibrator_create_model(
                b=b,
                sites_only_vcf=sites_only_vcf,
                transmitted_singletons_resource_vcf=transmitted_singletons,
                sibling_singletons_resource_vcf=sibling_singletons,
                utils=utils,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                out_bucket=out_bucket,
                is_small_callset=is_small_callset,
                max_gaussians=indel_max_gaussians,
            ).model_file

        indels_recalibrator_jobs = [
            indels_variant_recalibrator(
                b=b,
                sites_only_vcf=sites_only_vcf,
                utils=utils,
                out_bucket=out_bucket,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                interval=intervals[f"interval_{idx}"],
                tranche_idx=idx,
                model_file=indels_model_file,
                transmitted_singletons_resource_vcf=transmitted_singletons,
                sibling_singletons_resource_vcf=sibling_singletons,
                max_gaussians=indel_max_gaussians,
            )
            for idx in range(utils["NUMBER_OF_GENOMICS_DB_INTERVALS"])
        ]

        indels_recalibrations = [j.recalibration for j in indels_recalibrator_jobs]
        indels_tranches = [j.tranches for j in indels_recalibrator_jobs]
        indels_gathered_tranches = gather_tranches(
            b=b,
            tranches=indels_tranches,
            mode="INDEL",
            disk_size=small_disk,
            gcp_billing_project=gcp_billing_project,
        ).out_tranches

        # 3. Apply recalibration
        # Doesn't require too much storage in scatter mode (<500MB on gnomad VCF
        # for each scatter), use small_disk
        scattered_vcfs = [
            apply_recalibration(
                b=b,
                input_vcf=sites_only_vcf,
                out_vcf_name=output_vcf_name,
                indels_recalibration=indels_recalibrations[idx],
                indels_tranches=indels_gathered_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                utils=utils,
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                scatter=idx,
                interval=intervals[f"interval_{idx}"],
                out_bucket=out_bucket,
            ).output_vcf
            for idx in range(utils["NUMBER_OF_GENOMICS_DB_INTERVALS"])
        ]

        # 4. Gather VCFs
        gathered_vcf_job = gather_vcfs(
            b=b,
            input_vcfs=scattered_vcfs,
            out_vcf_name=output_vcf_name,
            utils=utils,
            disk_size=huge_disk,
            out_bucket=out_bucket,
            gcp_billing_project=gcp_billing_project,
        )

    else:
        snps_recalibrator_job = snps_variant_recalibrator(
            b=b,
            sites_only_vcf=sites_only_vcf,
            utils=utils,
            out_bucket=out_bucket,
            use_as_annotations=use_as_annotations,
            gcp_billing_project=gcp_billing_project,
            transmitted_singletons_resource_vcf=transmitted_singletons,
            sibling_singletons_resource_vcf=sibling_singletons,
            max_gaussians=snp_max_gaussians,
        )
        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        indels_variant_recalibrator_job = indels_variant_recalibrator(
            b=b,
            sites_only_vcf=sites_only_vcf,
            transmitted_singletons_resource_vcf=transmitted_singletons,
            sibling_singletons_resource_vcf=sibling_singletons,
            utils=utils,
            out_bucket=out_bucket,
            use_as_annotations=use_as_annotations,
            gcp_billing_project=gcp_billing_project,
            max_gaussians=indel_max_gaussians,
        )
        indels_recalibration = indels_variant_recalibrator_job.recalibration
        indels_tranches = indels_variant_recalibrator_job.tranches

        recalibrated_gathered_vcf_job = apply_recalibration(
            b=b,
            input_vcf=sites_only_vcf,
            out_vcf_name=output_vcf_name,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            utils=utils,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            gcp_billing_project=gcp_billing_project,
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
    batch_billing_project: str,
    gcp_billing_project: str,
    run_mode: str,
    batch_suffix: str,
    use_as_annotations: bool = True,
):
    """
    Wraps all the functions into a workflow
    :param sites_only_vcf: path to a sites only VCF created using gnomAD default_compute_info()
    :param output_vcf_filename: name, without extension, to use for the output VCF file(s)
    :param transmitted_singletons: full path to transmitted singletons VCF file and its index
    :param sibling_singletons: full path to sibling singletons VCF file and its index
    :param resources: json file (vqsr_resources.json) with paths to the resource files and parameters to be used in VQSR
    :param out_bucket: path to write, plots, evaluation results, and recalibrated VCF to
    :param batch_billing_project: Batch billing project to be used for the workflow
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :param use_as_annotations: use allele-specific annotation for VQSR

    :return:
    """
    hl.init(
        backend="batch",
        tmp_dir="gs://gnomad-tmp-4day/",
        gcs_requester_pays_configuration=gcp_billing_project,
        regions=["us-central1"],
    )

    tmp_vqsr_bucket = f"{out_bucket}/"

    backend = hb.ServiceBackend(
        billing_project=batch_billing_project,
        remote_tmpdir="gs://gnomad-tmp-4day/",
    )

    with open(resources, "r") as f:
        utils = json.load(f)

    logger.info(
        f"Starting hail Batch with the project {batch_billing_project}, "
        f"bucket {tmp_vqsr_bucket}"
    )

    b = hb.Batch(
        f"VQSR pipeline{batch_suffix}",
        backend=backend,
    )

    intervals_j = split_intervals(
        b=b, utils=utils, gcp_billing_project=gcp_billing_project
    )

    is_small_callset = False
    is_large_callset = False
    if run_mode == "small":
        is_small_callset = True
    elif run_mode == "large":
        is_large_callset = True

    make_vqsr_jobs(
        b=b,
        sites_only_vcf=sites_only_vcf,
        is_small_callset=is_small_callset,
        is_large_callset=is_large_callset,
        output_vcf_name=output_vcf_filename,
        utils=utils,
        out_bucket=tmp_vqsr_bucket,
        intervals=intervals_j.intervals,
        use_as_annotations=use_as_annotations,
        gcp_billing_project=gcp_billing_project,
        transmitted_singletons=transmitted_singletons,
        sibling_singletons=sibling_singletons,
    )

    b.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-vcf",
        type=str,
        required=True,
        help="Input sites VCF containing AS annotations.",
    )
    parser.add_argument(
        "--out-bucket", type=str, required=True, help="Bucket to store VQSR outputs in."
    )
    parser.add_argument(
        "--out-vcf-name",
        type=str,
        required=True,
        help="Required prefix for VQSR outputs.",
    )
    parser.add_argument(
        "--resources",
        type=str,
        required=True,
        help=(
            "Path to .json file containing paths and information for the VQSR pipeline."
        ),
    )
    parser.add_argument(
        "--run-mode",
        type=str,
        default="standard",
        choices=["small", "standard", "large"],
        help=(
            "Option to pass so that the mode/resources fit the size of the database."
            " This affects the size of the clusters and, if --large is set, will run in"
            " a scattered mode (one job for each partition)."
        ),
    )
    parser.add_argument(
        "--batch-billing-project",
        type=str,
        required=True,
        help="Hail Batch billing project.",
    )
    parser.add_argument(
        "--gcp-billing-project",
        type=str,
        required=True,
        help="Google Cloud billing project for reading requester pays buckets.",
    )
    parser.add_argument(
        "--transmitted-singletons",
        type=str,
        required=False,
        help="Path to transmitted singletons or first singleton truth set VCF.",
    )
    parser.add_argument(
        "--sibling-singletons",
        type=str,
        required=False,
        help="Path to sibling singletons or second singleton truth set VCF.",
    )

    parser.add_argument(
        "--no-as-annotations",
        action="store_true",
        help="Option to pass if you do not want to use AS annotation.",
    )
    parser.add_argument(
        "--batch-suffix",
        type=str,
        default="",
        help="String to add to end of batch name.",
    )
    parser.add_argument(
        "--test-on-chr22",
        action='store_true',
        help='If passed, will search resource file for _CHR22 versions of some files'

    )

    args = parser.parse_args()

    use_as_annotations = False if args.no_as_annotations else True

    # a smarter man with more time would have an idea for implementing this!
    # arg_suffix = ""

    # if args.test_on_chr22:
    #     arg_uffix="_CHR22"


    print("billing project as: ", args.batch_billing_project)

    vqsr_workflow(
        sites_only_vcf=args.input_vcf,
        output_vcf_filename=args.out_vcf_name,
        transmitted_singletons=args.transmitted_singletons,
        sibling_singletons=args.sibling_singletons,
        resources=args.resources,
        out_bucket=args.out_bucket,
        batch_billing_project=args.batch_billing_project,
        gcp_billing_project=args.gcp_billing_project,
        run_mode=args.run_mode,
        use_as_annotations=use_as_annotations,
        batch_suffix=args.batch_suffix,
    )


if __name__ == "__main__":
    main()
