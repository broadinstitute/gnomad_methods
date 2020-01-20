import argparse
import hail as hl
from .reference_data import (
    purcell_5k_intervals,
    clinvar,
    na12878_giab,
    dbsnp,
    hapmap,
    kgp_omni,
    kgp,
    mills,
    syndip,
    lcr_intervals,
    decoy_intervals,
    seg_dup_intervals,
    exome_hc_intervals,
    exome_calling_intervals,
    exome_evaluation_intervals,
    genome_evaluation_intervals,
    high_coverage_intervals,
    na12878_hc_exome_intervals,
    na12878_hc_intervals,
    syndip_hc_intervals,
)


def main(args):

    hl.init(log="/load_resources.hail.log", default_reference="GRCh37")
    if args.purcell_5k:
        hl.import_locus_intervals(
            purcell_5k_intervals.import_sources["source_path"],
            reference_genome="GRCh37",
        ).write(purcell_5k_intervals.path, overwrite=args.overwrite)

    if args.clinvar:
        hl.import_vcf(
            clinvar.import_sources["source_path"],
            force_bgz=True,
            skip_invalid_loci=True,
            min_partitions=100,
        ).rows().write(clinvar.path, overwrite=args.overwrite)

    if args.na12878:
        hl.import_vcf(
            na12878_giab.import_sources["source_path"],
            force_bgz=True,
            min_partitions=100,
        ).write(na12878_giab.path, overwrite=args.overwrite)

    if args.dbsnp:
        hl.import_vcf(
            dbsnp.import_sources["source_path"],
            force_bgz=True,
            skip_invalid_loci=True,
            min_partitions=400,
        ).rows().write(dbsnp.path, overwrite=args.overwrite)

    if args.hapmap:
        hl.import_vcf(hapmap.import_sources["source_path"], force_bgz=True).write(
            hapmap.path, overwrite=args.overwrite
        )

    if args.kgp_omni:
        hl.import_vcf(kgp_omni.import_sources["source_path"], force_bgz=True).write(
            kgp_omni.path, overwrite=args.overwrite
        )

    if args.kgp:
        hl.import_vcf(kgp.import_sources["source_path"], force_bgz=True).rows().write(
            kgp.path, overwrite=args.overwrite
        )

    if args.mills:
        hl.import_vcf(mills.import_sources["source_path"], force_bgz=True).write(
            mills.path, overwrite=args.overwrite
        )

    if args.syndip:
        hl.import_vcf(syndip.import_sources["source_path"], force_bgz=True).write(
            syndip.path, overwrite=args.overwrite
        )

    if args.lcr_intervals:
        hl.import_locus_intervals(syndip.import_sources["source_path"]).write(
            lcr_intervals.path, overwrite=args.overwrite
        )

    if args.decoy_intervals:
        hl.import_bed(decoy_intervals.import_sources["source_path"]).write(
            decoy_intervals.path, overwrite=args.overwrite
        )

    if args.seg_dup_intervals:
        hl.import_bed(seg_dup_intervals.import_sources["source_path"]).write(
            seg_dup_intervals.path, overwrite=args.overwrite
        )

    if args.exome_hc_intervals:
        hl.import_locus_intervals(
            exome_hc_intervals.import_sources["source_path"]
        ).write(exome_hc_intervals.path, overwrite=args.overwrite)

    if args.exome_calling_intervals:
        hl.import_locus_intervals(
            exome_calling_intervals.import_sources["source_path"]
        ).write(exome_calling_intervals.path, overwrite=args.overwrite)

    if args.exome_evaluation_intervals:
        hl.import_locus_intervals(
            exome_evaluation_intervals.import_sources["source_path"]
        ).write(exome_evaluation_intervals.path, overwrite=args.overwrite)

    if args.genome_evaluation_intervals:
        hl.import_locus_intervals(
            genome_evaluation_intervals.import_sources["source_path"]
        ).write(genome_evaluation_intervals.path, overwrite=args.overwrite)

    if args.high_coverage_intervals:
        hl.import_locus_intervals(
            high_coverage_intervals.import_sources["source_path"]
        ).write(high_coverage_intervals.path, overwrite=args.overwrite)

    if args.na12878_hc_intervals:
        hl.import_bed(na12878_hc_intervals.import_sources["source_path"]).write(
            na12878_hc_intervals.path, overwrite=args.overwrite
        )

    if args.na12878_exome_hc_intervals:
        hl.import_bed(na12878_hc_exome_intervals.import_sources["source_path"]).write(
            na12878_hc_exome_intervals.path, overwrite=args.overwrite
        )

    if args.syndyp_hc_intervals:
        hl.import_bed(syndip_hc_intervals.import_sources["source_path"]).write(
            syndip_hc_intervals.path, overwrite=args.overwrite
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--purcell_5k", help="Purcell5k sites (from intervals)", action="store_true"
    )
    parser.add_argument("--clinvar", help="Import clinvar VCF", action="store_true")
    parser.add_argument("--na12878", help="Imports GiaB NA12878", action="store_true")
    parser.add_argument("--dbsnp", help="Imports DBSNP", action="store_true")
    parser.add_argument("--hapmap", help="Imports HapMap", action="store_true")
    parser.add_argument( "--kgp_omni", help="Imports Omni / KGP sites VCF", action="store_true")
    parser.add_argument("--kgp", help="Imports 1000 Genomes sites VCF", action="store_true")
    parser.add_argument("--mills", help="Imports Mills Devine indel sites VCF", action="store_true")
    parser.add_argument( "--syndip", help="Imports synthetic duplicates VCF", action="store_true")
    parser.add_argument("--lcr_intervals", help="Imports low confidence regions' intervals", action="store_true")
    parser.add_argument("--decoy_intervals", help="Imports decoy intervals", action="store_true")
    parser.add_argument("--seg_dup_intervals", help="Imports segmental duplication intervals", action="store_true")
    parser.add_argument("--exomes_hc_intervals", help="Imports high confidence exome intervals", action="store_true")
    parser.add_argument("--exomes_calling_intervals", help="Imports exome calling intervals", action="store_true")
    parser.add_argument("--exomes_evaluation_intervals", help="Imports exome evaluation intervals", action="store_true")
    parser.add_argument("--genomes_evaluation_intervals", help="Imports genome evaluation intervals", action="store_true")
    parser.add_argument("--high_coverage_intervals", help="Imports genome evaluation intervals", action="store_true")
    parser.add_argument("--na12878_hc_intervals", help="Imports NA12878 high confidence regions", actions="store_true")
    parser.add_argument("--na12878_hc_exome_intervals", help="Imports NA12878 ecxome high confidence regions", actions="store_true")
    parser.add_argument("--syndyp_hc_intervals", help="Imports high confidence synthetic duplicate intervals", actions="store_true")
    parser.add_argument("--overwrite", help="Overwrites existing files", action="store_true")
    main(parser.parse_args())
