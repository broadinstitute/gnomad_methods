
import re
import sys
import logging
import gzip
import os

import hail as hl
from hail.expr.expressions import *
from collections import defaultdict, namedtuple, OrderedDict
from pprint import pprint, pformat
import argparse
from typing import *

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


def annotate_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid)
    """
    adj_gq = 20
    adj_dp = 10
    adj_ab = 0.2

    return mt.annotate_entries(adj=(mt.GQ >= adj_gq) & (mt.DP >= adj_dp) & (
                                   ~mt.GT.is_het() |
                                   ((mt.GT[0] == 0) & (mt.AD[mt.GT[1]] / mt.DP >= adj_ab)) |
                                   ((mt.GT[0] > 0) & (mt.AD[mt.GT[0]] / mt.DP >= adj_ab) &
                                    (mt.AD[mt.GT[1]] / mt.DP >= adj_ab))
                               ))


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


def adjust_sex_ploidy(mt: hl.MatrixTable, sex_expr: hl.expr.StringExpression,
                      male_str: str = 'male', female_str: str = 'female') -> hl.MatrixTable:
    """
    Converts males to haploid on non-PAR X/Y, sets females to missing on Y

    :param MatrixTable mt: Input MatrixTable
    :param StringExpression sex_expr: Expression pointing to sex in MT (if not male_str or female_str, no change)
    :param str male_str: String for males (default 'male')
    :param str female_str: String for females (default 'female')
    :return: MatrixTable with fixed ploidy for sex chromosomes
    :rtype: MatrixTable
    """
    male = sex_expr == male_str
    female = sex_expr == female_str
    x_nonpar = mt.locus.in_x_nonpar()
    y_par = mt.locus.in_y_par()
    y_nonpar = mt.locus.in_y_nonpar()
    return mt.annotate_entries(
        GT=hl.case(missing_false=True)
        .when(female & (y_par | y_nonpar), hl.null(hl.tcall))
        .when(male & (x_nonpar | y_nonpar) & mt.GT.is_het(), hl.null(hl.tcall))
        .when(male & (x_nonpar | y_nonpar), hl.call(mt.GT[0], phased=False))
        .default(mt.GT)
    )


def add_popmax_expr(freq: hl.expr.ArrayExpression) -> hl.expr.ArrayExpression:
    """
    Calculates popmax (add an additional entry into freq with popmax: pop)

    :param ArrayExpression freq: ArrayExpression of Structs with ['ac', 'an', 'hom', 'meta']
    :return: Frequency data with annotated popmax
    :rtype: ArrayExpression
    """
    freq_filtered = hl.filter(lambda x: (x.meta.keys() == ['population']) & (x.meta['population'] != 'oth'), freq)
    sorted_freqs = hl.sorted(freq_filtered, key=lambda x: x.AC[1] / x.AN, reverse=True)
    return hl.cond(hl.len(sorted_freqs) > 0, freq.append(
        hl.struct(AC=sorted_freqs[0].AC, AF=sorted_freqs[0].AF, AN=sorted_freqs[0].AN,
                  homozygote_count=sorted_freqs[0].homozygote_count,
                  meta={'popmax': sorted_freqs[0].meta['population']})), freq)


def get_projectmax(mt: hl.MatrixTable, loc: hl.expr.StringExpression) -> hl.MatrixTable:
    """
    First pass of projectmax (returns aggregated MT with project_max field)

    :param MatrixTable mt: Input MT
    :param StringExpression loc: Column expression location of project ID (e.g. mt.meta.pid)
    :return: Frequency data with annotated project_max
    :rtype: MatrixTable
    """
    mt = mt.annotate_cols(project=loc)
    agg_mt = mt.group_cols_by(mt.project).aggregate(callstats=hl.agg.call_stats(mt.GT, mt.alleles))
    return agg_mt.annotate_rows(project_max=hl.agg.take(hl.struct(**agg_mt.callstats, project=agg_mt.project),
                                                        5, -agg_mt.callstats.AF[1]))


def read_list_data(input_file: str) -> List[str]:
    if input_file.startswith('gs://'):
        hl.hadoop_copy(input_file, 'file:///' + input_file.split("/")[-1])
        f = gzip.open("/" + os.path.basename(input_file)) if input_file.endswith('gz') else open("/" + os.path.basename(input_file))
    else:
        f = gzip.open(input_file) if input_file.endswith('gz') else open(input_file)
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def melt_kt(kt, columns_to_melt, key_column_name='variable', value_column_name='value'):
    """
    Go from wide to long, or from:

    +---------+---------+---------+
    | Variant | AC_NFE  | AC_AFR  |
    +=========+=========+=========+
    | 1:1:A:G |      1  |      8  |
    +---------+---------+---------+
    | 1:2:A:G |     10  |    100  |
    +---------+---------+---------+

    to:

    +---------+----------+--------+
    | Variant | variable | value  |
    +=========+==========+========+
    | 1:1:A:G |   AC_NFE |     1  |
    +---------+----------+--------+
    | 1:1:A:G |   AC_AFR |     8  |
    +---------+----------+--------+
    | 1:2:A:G |   AC_NFE |    10  |
    +---------+----------+--------+
    | 1:2:A:G |   AC_AFR |   100  |
    +---------+----------+--------+

    :param KeyTable kt: Input KeyTable
    :param list of str columns_to_melt: Which columns to spread out
    :param str key_column_name: What to call the key column
    :param str value_column_name: What to call the value column
    :return: melted Key Table
    :rtype: KeyTable
    return (kt
            .annotate('comb = [{}]'.format(', '.join(['{{k: "{0}", value: {0}}}'.format(x) for x in columns_to_melt])))
            .drop(columns_to_melt)
            .explode('comb')
            .annotate('{} = comb.k, {} = comb.value'.format(key_column_name, value_column_name))
            .drop('comb'))
    """
    raise NotImplementedError


def melt_kt_grouped(kt, columns_to_melt, value_column_names, key_column_name='variable'):
    """
    Go from wide to long for a group of variables, or from:

    +---------+---------+---------+---------+---------+
    | Variant | AC_NFE  | AC_AFR  | Hom_NFE | Hom_AFR |
    +=========+=========+=========+=========+=========+
    | 1:1:A:G |      1  |      8  |       0 |       0 |
    +---------+---------+---------+---------+---------+
    | 1:2:A:G |     10  |    100  |       1 |      10 |
    +---------+---------+---------+---------+---------+

    to:

    +---------+----------+--------+--------+
    | Variant |      pop |    AC  |   Hom  |
    +=========+==========+========+========+
    | 1:1:A:G |      NFE |     1  |     0  |
    +---------+----------+--------+--------+
    | 1:1:A:G |      AFR |     8  |     0  |
    +---------+----------+--------+--------+
    | 1:2:A:G |      NFE |    10  |     1  |
    +---------+----------+--------+--------+
    | 1:2:A:G |      AFR |   100  |    10  |
    +---------+----------+--------+--------+

    This is done with:

    columns_to_melt = {
        'NFE': ['AC_NFE', 'Hom_NFE'],
        'AFR': ['AC_AFR', 'Hom_AFR']
    }
    value_column_names = ['AC', 'Hom']
    key_column_name = 'pop'

    Note that len(value_column_names) == len(columns_to_melt[i]) for all in columns_to_melt

    :param KeyTable kt: Input KeyTable
    :param dict of list of str columns_to_melt: Which columns to spread out
    :param list of str value_column_names: What to call the value columns
    :param str key_column_name: What to call the key column
    :return: melted Key Table
    :rtype: KeyTable

    if any([len(value_column_names) != len(v) for v in columns_to_melt.values()]):
        logger.warning('Length of columns_to_melt sublist is not equal to length of value_column_names')
        logger.warning('value_column_names = %s', value_column_names)
        logger.warning('columns_to_melt = %s', columns_to_melt)

    # I think this goes something like this:
    fields = []
    for k, v in columns_to_melt.items():
        subfields = [': '.join(x) for x in zip(value_column_names, v)]
        field = '{{k: "{0}", {1}}}'.format(k, ', '.join(subfields))
        fields.append(field)

    split_text = ', '.join(['{0} = comb.{0}'.format(x) for x in value_column_names])

    return (kt
            .annotate('comb = [{}]'.format(', '.join(fields)))
            .drop([y for x in columns_to_melt.values() for y in x])
            .explode('comb')
            .annotate('{} = comb.k, {}'.format(key_column_name, split_text))
            .drop('comb'))
    """
    raise NotImplementedError


def get_duplicated_samples(
        kin_ht: hl.Table,
        i_col: str = 'i',
        j_col: str = 'j',
        kin_col: str = 'kin',
        duplicate_threshold: float = 0.4
) -> List[Set[str]]:
    """
    Given a pc_relate output Table, extract the list of duplicate samples. Returns a list of set of samples that are duplicates.


    :param Table kin_ht: pc_relate output table
    :param str i_col: Column containing the 1st sample
    :param str j_col: Column containing the 2nd sample
    :param str kin_col: Column containing the kinship value
    :param float duplicate_threshold: Kinship threshold to consider two samples duplicated
    :return: List of samples that are duplicates
    :rtype: list of set of str
    """

    def get_all_dups(s, dups, samples_duplicates):
        if s in samples_duplicates:
            dups.add(s)
            s_dups = samples_duplicates.pop(s)
            for s_dup in s_dups:
                if s_dup not in dups:
                    dups = get_all_dups(s_dup, dups, samples_duplicates)
        return dups

    dup_rows = kin_ht.filter(kin_ht[kin_col] > duplicate_threshold).collect()

    samples_duplicates = defaultdict(set)
    for row in dup_rows:
        samples_duplicates[row[i_col]].add(row[j_col])
        samples_duplicates[row[j_col]].add(row[i_col])

    duplicated_samples = []
    while len(samples_duplicates) > 0:
        duplicated_samples.append(get_all_dups(list(samples_duplicates)[0], set(), samples_duplicates))

    return duplicated_samples


def infer_families(kin_ht: hl.Table,
                   sex: Dict[str, bool],
                   duplicated_samples: Set[str],
                   i_col: str = 'i',
                   j_col: str = 'j',
                   kin_col: str = 'kin',
                   ibd2_col: str = 'ibd2',
                   first_degree_threshold: Tuple[float, float] = (0.2, 0.4),
                   second_degree_threshold: Tuple[float, float] = (0.05, 0.16),
                   ibd2_parent_offspring_threshold: float = 0.2
                   ) -> hl.Pedigree:
    """

    Infers familial relationships from the results of pc_relate and sex information.
    Note that both kinship and ibd2 are needed in the pc_relate output.

    This function returns a pedigree containing trios inferred from the data. Family ID can be the same for multiple
    trios if one or more members of the trios are related (e.g. sibs, multi-generational family). Trios are ordered by family ID.

    Note that this function only returns complete trios defined as:
    one child, one father and one mother (sex is required for both parents)

    :param Table kin_ht: pc_relate output table
    :param dict of str -> bool sex: A dict containing the sex for each sample. True = female, False = male, None = unknown
    :param set of str duplicated_samples: Duplicated samples to remove (If not provided, this function won't work as it assumes that each child has exactly two parents)
    :param str i_col: Column containing the 1st sample id in the pc_relate table
    :param str j_col: Column containing the 2nd sample id in the pc_relate table
    :param str kin_col: Column containing the kinship in the pc_relate table
    :param str ibd2_col: Column containing ibd2 in the pc_relate table
    :param (float, float) first_degree_threshold: Lower/upper bounds for kin for 1st degree relatives
    :param (float, float) second_degree_threshold: Lower/upper bounds for kin for 2nd degree relatives
    :param float ibd2_parent_offspring_threshold: Upper bound on ibd2 for a parent/offspring
    :return: Pedigree containing all trios in the data
    :rtype: Pedigree
    """

    def get_fam_samples(sample: str,
                        fam: Set[str],
                        samples_rel: Dict[str, Set[str]],
                        ) -> Set[str]:
        """
        Given a sample, its known family and a dict that links samples with their relatives, outputs the set of
        samples that constitute this sample family.

        :param str sample: sample
        :param dict of str -> set of str samples_rel: dict(sample -> set(sample_relatives))
        :param set of str fam: sample known family
        :return: Family including the sample
        :rtype: set of str
        """
        fam.add(sample)
        for s2 in samples_rel[sample]:
            if s2 not in fam:
                fam = get_fam_samples(s2, fam, samples_rel)
        return fam

    def get_indexed_ibd2(
            pc_relate_rows: List[hl.Struct]
    ) -> Dict[Tuple[str, str], float]:
        """
        Given rows from a pc_relate table, creates a dict with:
        keys: Pairs of individuals, lexically ordered
        values: ibd2

        :param list of hl.Struct pc_relate_rows: Rows from a pc_relate table
        :return: Dict of lexically ordered pairs of individuals -> kinship
        :rtype: dict of (str, str) -> float
        """
        ibd2 = dict()
        for row in pc_relate_rows:
            ibd2[tuple(sorted((row[i_col], row[j_col])))] = row[ibd2_col]
        return ibd2

    def get_parents(
            possible_parents: List[str],
            indexed_kinship: Dict[Tuple[str, str], Tuple[float, float]],
            sex: Dict[str, bool]
    ) -> Tuple[str, str]:
        """
        Given a list of possible parents for a sample (first degree relatives with low ibd2),
        looks for a single pair of samples that are unrelated with different sexes.
        If a single pair is found, return the pair (father, mother)

        :param list of str possible_parents: Possible parents
        :param dict of (str, str) -> (float, float)) indexed_kinship: Dict mapping pairs of individuals to their kinship and ibd2 coefficients
        :param dict of str -> bool sex: Dict mapping samples to their sex (True = female, False = male, None or missing = unknown)
        :return: (father, mother)
        :rtype: (str, str)
        """

        parents = []
        while len(possible_parents) > 1:
            p1 = possible_parents.pop()
            for p2 in possible_parents:
                if tuple(sorted((p1,p2))) not in indexed_kinship:
                    if sex.get(p1) is False and sex.get(p2):
                        parents.append((p1,p2))
                    elif sex.get(p1) and sex.get(p2) is False:
                        parents.append((p2,p1))

        if len(parents) == 1:
            return parents[0]

        return None

    # Get first degree relatives - exclude duplicate samples
    dups = hl.literal(duplicated_samples)
    first_degree_pairs = kin_ht.filter(
        (kin_ht[kin_col] > first_degree_threshold[0]) &
        (kin_ht[kin_col] < first_degree_threshold[1]) &
        ~dups.contains(kin_ht[i_col]) &
        ~dups.contains(kin_ht[j_col])
    ).collect()
    first_degree_relatives = defaultdict(set)
    for row in first_degree_pairs:
        first_degree_relatives[row[i_col]].add(row[j_col])
        first_degree_relatives[row[j_col]].add(row[i_col])

    # Add second degree relatives for those samples
    # This is needed to distinguish grandparent - child - parent from child - mother, father down the line
    first_degree_samples = hl.literal(set(first_degree_relatives.keys()))
    second_degree_samples = kin_ht.filter(
        (first_degree_samples.contains(kin_ht[i_col]) | first_degree_samples.contains(kin_ht[j_col])) &
        (kin_ht[kin_col] > second_degree_threshold[0]) &
        (kin_ht[kin_col] < first_degree_threshold[1])
    ).collect()

    ibd2 = get_indexed_ibd2(second_degree_samples)

    fam_id = 1
    trios = []
    while len(first_degree_relatives) > 0:
        s_fam = get_fam_samples(list(first_degree_relatives)[0], set(),
                                first_degree_relatives)
        for s in s_fam:
            s_rel = first_degree_relatives.pop(s)
            possible_parents = []
            for rel in s_rel:
                if ibd2[tuple(sorted((s, rel)))] < ibd2_parent_offspring_threshold:
                    possible_parents.append(rel)

            parents = get_parents(possible_parents, ibd2, sex)

            if parents is not None:
                trios.append(hl.Trio(s=s,
                                     fam_id=str(fam_id),
                                     pat_id=parents[0],
                                     mat_id=parents[1],
                                     is_female=sex.get(s)))

        fam_id += 1

    return hl.Pedigree(trios)
