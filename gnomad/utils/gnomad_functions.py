
import logging
import gzip
import os
import hail as hl
from pprint import pformat
import pandas as pd
from typing import *
import json

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)


def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter genotypes to adj criteria
    """
    if 'adj' not in list(mt.entry):
        mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)
    return mt.drop(mt.adj)


def get_adj_expr(
        gt_expr: hl.expr.CallExpression,
        gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
        dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
        ad_expr: hl.expr.ArrayNumericExpression,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.expr.BooleanExpression:
    """
    Gets adj genotype annotation.
    Defaults correspond to gnomAD values.
    """
    return (
            (gq_expr >= adj_gq) &
            hl.cond(
                gt_expr.is_haploid(),
                dp_expr >= haploid_adj_dp,
                dp_expr >= adj_dp
            ) &
            (
                hl.case()
                .when(~gt_expr.is_het(), True)
                .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
                .default((ad_expr[gt_expr[0]] / dp_expr >= adj_ab ) & (ad_expr[gt_expr[1]] / dp_expr >= adj_ab ))
            )
    )


def annotate_adj(
        mt: hl.MatrixTable,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid)
    Defaults correspond to gnomAD values.
    """
    return mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD, adj_gq, adj_dp, adj_ab, haploid_adj_dp))


def add_variant_type(alt_alleles: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """
    Get Struct of variant_type and n_alt_alleles from ArrayExpression of Strings (all alleles)
    """
    ref = alt_alleles[0]
    alts = alt_alleles[1:]
    non_star_alleles = hl.filter(lambda a: a != '*', alts)
    return hl.struct(variant_type=hl.cond(
        hl.all(lambda a: hl.is_snp(ref, a), non_star_alleles),
        hl.cond(hl.len(non_star_alleles) > 1, "multi-snv", "snv"),
        hl.cond(
            hl.all(lambda a: hl.is_indel(ref, a), non_star_alleles),
            hl.cond(hl.len(non_star_alleles) > 1, "multi-indel", "indel"),
            "mixed")
    ), n_alt_alleles=hl.len(non_star_alleles))


def adjusted_sex_ploidy_expr(
        locus_expr: hl.expr.LocusExpression,
        gt_expr: hl.expr.CallExpression,
        karyotype_expr: hl.expr.StringExpression,
        xy_karyotype_str: str = 'XY',
        xx_karyotype_str: str = 'XX'
) -> hl.expr.CallExpression:
    """
    Creates an entry expression to convert males to haploid on non-PAR X/Y and females to missing on Y

    :param locus_expr: Locus
    :param gt_expr: Genotype
    :param karyotype_expr: Karyotype
    :param xy_karyotype_str: Male sex karyotype representation
    :param xx_karyotype_str: Female sex karyotype representation
    :return: Genotype adjusted for sex ploidy
    """
    male = karyotype_expr == xy_karyotype_str
    female = karyotype_expr == xx_karyotype_str
    x_nonpar = locus_expr.in_x_nonpar()
    y_par = locus_expr.in_y_par()
    y_nonpar = locus_expr.in_y_nonpar()
    return (
        hl.case(missing_false=True)
            .when(female & (y_par | y_nonpar), hl.null(hl.tcall))
            .when(male & (x_nonpar | y_nonpar) & gt_expr.is_het(), hl.null(hl.tcall))
            .when(male & (x_nonpar | y_nonpar), hl.call(gt_expr[0], phased=False))
            .default(gt_expr)
    )


def adjust_sex_ploidy(mt: hl.MatrixTable, sex_expr: hl.expr.StringExpression,
                      male_str: str = 'male', female_str: str = 'female') -> hl.MatrixTable:
    """
    Converts males to haploid on non-PAR X/Y, sets females to missing on Y

    :param mt: Input MatrixTable
    :param sex_expr: Expression pointing to sex in MT (if not male_str or female_str, no change)
    :param male_str: String for males (default 'male')
    :param female_str: String for females (default 'female')
    :return: MatrixTable with fixed ploidy for sex chromosomes
    """
    return mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus,
            mt.GT,
            sex_expr,
            male_str,
            female_str
        )
    )


def read_list_data(input_file_path: str) -> List[str]:
    """
    Reads a file input into a python list (each line will be an element).
    Supports Google storage paths and .gz compression.

    :param input_file_path: File path
    :return: List of lines
    """
    if input_file_path.startswith('gs://'):
        hl.hadoop_copy(input_file_path, 'file:///' + input_file_path.split("/")[-1])
        f = gzip.open("/" + os.path.basename(input_file_path)) if input_file_path.endswith('gz') else open("/" + os.path.basename(input_file_path))
    else:
        f = gzip.open(input_file_path) if input_file_path.endswith('gz') else open(input_file_path)
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def liftover_using_gnomad_map(ht, data_type):
    """
    Liftover a gnomAD table using already-established liftover file. Warning: shuffles!

    :param ht: Input Hail table
    :param data_type: one of "exomes" or "genomes" which to map across
    :return: Lifted over table
    """
    from gnomad.resources.grch37.gnomad import liftover
    lift_ht = liftover(data_type).ht()
    ht = ht.key_by(original_locus=ht.locus, original_alleles=ht.alleles).drop('locus', 'alleles')
    return lift_ht.annotate(**ht[(lift_ht.original_locus, lift_ht.original_alleles)]).key_by('locus', 'alleles')


def filter_by_frequency(t: Union[hl.MatrixTable, hl.Table], direction: str,
                        frequency: float = None, allele_count: int = None,
                        population: str = None, subpop: str = None, downsampling: int = None,
                        keep: bool = True, adj: bool = True) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter MatrixTable or Table with gnomAD-format frequency data (assumed bi-allelic/split)
    (i.e. Array[Struct(Array[AC], Array[AF], AN, homozygote_count, meta)])

    At least one of frequency or allele_count is required.

    Subpop can be specified without a population if desired.

    :param t: Input MatrixTable or Table
    :param direction: One of "above", "below", and "equal" (how to apply the filter)
    :param frequency: Frequency to filter by (one of frequency or allele_count is required)
    :param allele_count: Allele count to filter by (one of frequency or allele_count is required)
    :param population: Population in which to filter frequency
    :param subpop: Sub-population in which to filter frequency
    :param downsampling: Downsampling in which to filter frequency
    :param keep: Whether to keep rows passing this frequency (passed to filter_rows)
    :param adj: Whether to use adj frequency
    :return: Filtered MatrixTable or Table
    """
    if frequency is None and allele_count is None:
        raise ValueError('At least one of frequency or allele_count must be specified')
    if direction not in ('above', 'below', 'equal'):
        raise ValueError('direction needs to be one of "above", "below", or "equal"')
    group = 'adj' if adj else 'raw'
    criteria = [lambda f: f.meta.get('group', '') == group]
    if frequency is not None:
        if direction == 'above':
            criteria.append(lambda f: f.AF[1] > frequency)
        elif direction == 'below':
            criteria.append(lambda f: f.AF[1] < frequency)
        else:
            criteria.append(lambda f: f.AF[1] == frequency)
    if allele_count is not None:
        if direction == 'above':
            criteria.append(lambda f: f.AC[1] > allele_count)
        elif direction == 'below':
            criteria.append(lambda f: f.AC[1] < allele_count)
        else:
            criteria.append(lambda f: f.AC[1] == allele_count)
    size = 1
    if population:
        criteria.append(lambda f: f.meta.get('pop', '') == population)
        size += 1
    if subpop:
        criteria.append(lambda f: f.meta.get('subpop', '') == subpop)
        size += 1
        # If one supplies a subpop but not a population, this will ensure this gets it right
        if not population: size += 1
    if downsampling:
        criteria.append(lambda f: f.meta.get('downsampling', '') == str(downsampling))
        size += 1
        if not population:
            size += 1
            criteria.append(lambda f: f.meta.get('pop', '') == 'global')
        if subpop:
            raise Exception('No downsampling data for subpopulations implemented')
    criteria.append(lambda f: f.meta.size() == size)

    def combine_functions(func_list, x):
        cond = func_list[0](x)
        for c in func_list[1:]:
            cond &= c(x)
        return cond

    filt = lambda x: combine_functions(criteria, x)
    criteria = hl.any(filt, t.freq)
    return t.filter_rows(criteria, keep=keep) if isinstance(t, hl.MatrixTable) else t.filter(criteria, keep=keep)


def pretty_print_runs(runs: Dict, label_col: str = 'rf_label', prediction_col_name: str = 'rf_prediction') -> None:
    """
    Prints the information for the RF runs loaded from the json file storing the RF run hashes -> info

    :param runs: Dictionary containing JSON input loaded from RF run file
    :param label_col: Name of the RF label column
    :param prediction_col_name: Name of the RF prediction column
    :return: Nothing -- only prints information
    """

    for run_hash, run_data in runs.items():
        print(f"\n=== {run_hash} ===")
        testing_results = run_data.pop('test_results') if 'test_results' in run_data else None
        # print(testing_results)
        print(json.dumps(run_data, sort_keys=True, indent=4, separators=(',', ': ')))
        if testing_results is not None:
            # Print results
            res_pd = pd.DataFrame(testing_results)
            res_pd = res_pd.pivot(index=label_col, columns=prediction_col_name, values='n')
            logger.info("Testing results:\n{}".format(pformat(res_pd)))


def add_rank(ht: hl.Table,
             score_expr: hl.expr.NumericExpression,
             subrank_expr: Optional[Dict[str, hl.expr.BooleanExpression]] = None) -> hl.Table:
    """
    Adds rank based on the `score_expr`. Rank is added for snvs and indels separately.
    If one or more `subrank_expr` are provided, then subrank is added based on all sites for which the boolean expression is true.

    In addition, variant counts (snv, indel separately) is added as a global (`rank_variant_counts`).

    :param ht: input Hail Table containing variants (with QC annotations) to be ranked
    :param score_expr: the Table annotation by which ranking should be scored
    :param subrank_expr: Any subranking to be added in the form name_of_subrank: subrank_filtering_expr
    :return: Table with rankings added
    """

    key = ht.key
    if subrank_expr is None:
        subrank_expr = {}

    temp_expr = {'_score': score_expr}
    temp_expr.update({f'_{name}': expr for name, expr in subrank_expr.items()})
    rank_ht = ht.select(**temp_expr, is_snv=hl.is_snp(ht.alleles[0], ht.alleles[1]))

    rank_ht = rank_ht.key_by('_score').persist()
    scan_expr = {'rank': hl.cond(rank_ht.is_snv, hl.scan.count_where(rank_ht.is_snv), hl.scan.count_where(~rank_ht.is_snv))}
    scan_expr.update({name: hl.or_missing(rank_ht[f'_{name}'],
                                          hl.cond(rank_ht.is_snv,
                                                  hl.scan.count_where(rank_ht.is_snv & rank_ht[f'_{name}']),
                                                  hl.scan.count_where(~rank_ht.is_snv & rank_ht[f'_{name}'])))
                      for name in subrank_expr})
    rank_ht = rank_ht.annotate(**scan_expr)

    rank_ht = rank_ht.key_by(*key).persist()
    rank_ht = rank_ht.select(*scan_expr.keys())

    ht = ht.annotate(**rank_ht[key])
    return ht
