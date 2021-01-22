import logging
from typing import Dict

import hail as hl

from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.utils.vep import (
    filter_vep_to_canonical_transcripts,
    get_worst_consequence_for_summary,
)


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def freq_criteria(
    freq_expr: hl.expr.ArrayExpression, index: int = 0
) -> hl.expr.builders.CaseBuilder:
    """
	Returns case statement adding frequency string annotations based on input AC or AF.

	:param freq_expr: Array of structs containing frequency information.
	:param index: Which index of freq_expr to use for annotation. Default is 0 (all pops calculated on adj genotypes only).
	:return: Case statement adding string expressions based on input AC or AF.
	"""
    return (
        hl.case(missing_false=True)
        .when(freq_expr[index].AC == 0, "Not found")
        .when(freq_expr[index].AC == 1, "Singleton")
        .when(freq_expr[index].AC == 2, "Doubleton")
        .when(freq_expr[index].AC <= 5, "AC 3 - 5")
        .when(freq_expr[index].AF < 1e-4, "AC 6 - 0.01%")
        .when(freq_expr[index].AF < 1e-3, "0.01% - 0.1%")
        .when(freq_expr[index].AF < 1e-2, "0.1% - 1%")
        .when(freq_expr[index].AF < 1e-1, "1% - 10%")
        .when(freq_expr[index].AF > 0.95, ">95%")
        .default("10% - 95%")
    )


def get_summary_counts_dict(
    allele_expr: hl.expr.ArrayExpression,
    lof_expr: hl.expr.StringExpression,
    no_lof_flags_expr: hl.expr.BooleanExpression,
    prefix_str: str = "",
) -> Dict[str, hl.expr.Int64Expression]:
    """
	Returns dictionary containing desired categories and their counts.

	Assumes allele_expr contains only two variants (multi-allelics have been split).

	:param allele_expr: ArrayExpression containing alleles.
	:param lof_expr: StringExpression containing LOFTEE annotation.
	:param no_lof_flags_expr: BooleanExpression indicating whether LoF variant has any flags.
	:param prefix_str: Desired prefix string for category names. Default is empty str.
	:return: Dict of categories and counts per category.
	"""
    return {
        f"{prefix_str}num_variants": hl.agg.count(),
        f"{prefix_str}indels": hl.agg.count_where(
            hl.is_indel(allele_expr[0], allele_expr[1])
        ),
        f"{prefix_str}snps": hl.agg.count_where(
            hl.is_snp(allele_expr[0], allele_expr[1])
        ),
        f"{prefix_str}LOF": hl.agg.count_where(hl.is_defined(lof_expr)),
        f"{prefix_str}pass_loftee": hl.agg.count_where(lof_expr == "HC"),
        f"{prefix_str}pass_loftee_no_flag": hl.agg.count_where(
            (lof_expr == "HC") & (no_lof_flags_expr)
        ),
        f"{prefix_str}loftee_os": hl.agg.count_where(lof_expr == "OS"),
        f"{prefix_str}fail_loftee": hl.agg.count_where(lof_expr == "LC"),
    }


def get_summary_per_variant(
    ht: hl.Table, freq_str: str = "freq", filter_field: str = "filters"
) -> hl.Table:
    """
	Generates a struct with summary counts per variant.

	Annotates Table's globals with total variant counts.

	Summary counts:
		- Number of variants
		- Number of indels
		- Number of SNPs
		- Number of LoF variants
		- Number of LoF variants that pass LOFTEE (including with LoF flags)
		- Number of LoF variants that pass LOFTEE without LoF flags
		- Number of OS (other splice) variants annotated by LOFTEE
		- Number of loF variants that fail LOFTEE filters

	Assumes that:
		- Input HT is annotated with VEP.
		- Multiallelic variants have been split and/or input HT contains bi-allelic variants only.

	:param ht: Input Table.
	:param freq: Name of field in HT containing frequency annotation (array of structs). Default is "freq".
	:param filter_field: Name of field in HT containing variant filter information. Default is "filters".
	:return: Table grouped by frequency bin and aggregated across summary count categories. 
	"""
    logger.info("Filtering to PASS variants in high confidence regions...")
    ht = ht.filter((hl.len(ht[filter_field]) == 0))
    ht = filter_low_conf_regions(ht)

    logger.info(
        "Filtering to canonical transcripts and getting VEP summary annotations..."
    )
    ht = filter_vep_to_canonical_transcripts(ht)
    ht = get_worst_consequence_for_summary(ht)

    logger.info("Annotating with frequency bin information...")
    ht = ht.annotate(freq_bin=freq_criteria(ht[freq_str]))

    logger.info("Annotating HT globals with total counts...")
    summary_per_variant = ht.aggregate(
        hl.struct(
            **get_summary_counts_dict(
                ht.alleles, ht.lof, ht.no_lof_flags, prefix_str="total_"
            )
        )
    )
    ht = ht.annotate_globals(summary_per_variant=summary_per_variant)
    return ht.group_by("freq_bin").aggregate(
        **ht.aggregate(
            hl.struct(**get_summary_counts_dict(ht.alleles, ht.lof, ht.no_lof_flags))
        )
    )
