
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


POP_NAMES = {
    'AFR': "African/African American",
    'AMR': "Admixed American",
    'ASJ': "Ashkenazi Jewish",
    'EAS': "East Asian",
    'FIN': "Finnish",
    'NFE': "Non-Finnish European",
    'OTH': "Other (population not assigned)",
    'SAS': "South Asian"
}

SEXES = {
    'Male': 'Male',
    'Female': 'Female'
}

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost"]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant"
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant"]

CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant"
]

CSQ_ORDER = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING


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


def split_multi_dynamic(mt: hl.MatrixTable, keep_star: bool = False, left_aligned: bool = True) -> hl.MatrixTable:
    """
    Splits MatrixTable based on entry fields found. Downcodes whatever it can. Supported so far:
    GT, DP, AD, PL, GQ
    PGT, PID
    ADALL

    :param MatrixTable mt: Input MatrixTable
    :param bool keep_star: whether to keep star alleles (passed to SplitMulti)
    :param bool left_aligned: whether matrix table is already left_aligned (passed to SplitMulti)
    :return: Split MatrixTable
    :rtype: MatrixTable
    """
    fields = list(mt.entry)
    sm = hl.SplitMulti(mt, keep_star=keep_star, left_aligned=left_aligned)
    sm.update_rows(a_index=sm.a_index(), was_split=sm.was_split())
    expression = {}

    # HTS/standard
    if 'GT' in fields:
        expression['GT'] = hl.downcode(mt.GT, sm.a_index())
    if 'DP' in fields:
        expression['DP'] = mt.DP
    if 'AD' in fields:
        expression['AD'] = hl.or_missing(hl.is_defined(mt.AD),
                                         [hl.sum(mt.AD) - mt.AD[sm.a_index()], mt.AD[sm.a_index()]])
    if 'PL' in fields:
        pl = hl.or_missing(
            hl.is_defined(mt.PL),
            (hl.range(0, 3).map(lambda i:
                                hl.min((hl.range(0, hl.triangle(mt.alleles.length()))
                                        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j),
                                                                      sm.a_index()) == hl.unphased_diploid_gt_index_call(i)
                                                ).map(lambda j: mt.PL[j]))))))
        expression['PL'] = pl
        if 'GQ' in fields:
            expression['GQ'] = hl.gq_from_pl(pl)
    else:
        if 'GQ' in fields:
            expression['GQ'] = mt.GQ

    # Phased data
    if 'PGT' in fields:
        expression['PGT'] = hl.downcode(mt.PGT, sm.a_index())
    if 'PID' in fields:
        expression['PGT'] = mt.PID

    # Custom data
    if 'ADALL' in fields:  # found in NA12878
        expression['ADALL'] = hl.or_missing(hl.is_defined(mt.ADALL),
                                            [hl.sum(mt.ADALL) - mt.ADALL[sm.a_index()], mt.ADALL[sm.a_index()]])

    sm.update_entries(**expression)
    return sm.result()


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


def get_sample_data(mt: hl.MatrixTable, fields: List[hl.expr.StringExpression], sep: str = '\t', delim: str = '|'):
    """
    Hail devs hate this one simple py4j trick to speed up sample queries

    :param MatrixTable or Table mt: MT
    :param list of StringExpression fields: fields
    :param sep: Separator to use (tab usually fine)
    :param delim: Delimiter to use (pipe usually fine)
    :return: Sample data
    :rtype: list of list of str
    """
    field_expr = fields[0]
    for field in fields[1:]:
        field_expr = field_expr + '|' + field
    if isinstance(mt, hl.MatrixTable):
        mt_agg = mt.aggregate_cols
    else:
        mt_agg = mt.aggregate
    return [x.split(delim) for x in mt_agg(hl.delimit(hl.agg.collect(field_expr), sep)).split(sep) if x != 'null']


def add_popmax_expr(freq: hl.expr.ArrayExpression) -> hl.expr.ArrayExpression:
    """
    Calculates popmax (add an additional entry into freq with popmax: pop)

    :param ArrayExpression freq: ArrayExpression of Structs with ['ac', 'an', 'hom', 'meta']
    :return: Frequency data with annotated popmax
    :rtype: ArrayExpression
    """
    freq_filtered = hl.filter(lambda x: (x.meta.keys() == ['population']) & (x.meta['population'] != 'oth'), freq)
    sorted_freqs = hl.sorted(freq_filtered, key=lambda x: x.ac / x.an, reverse=True)
    return hl.cond(hl.len(sorted_freqs) > 0, freq.append(
        hl.struct(ac=sorted_freqs[0].ac, an=sorted_freqs[0].an, hom=sorted_freqs[0].hom,
                  meta={'popmax': sorted_freqs[0].meta['population']})), freq)


def get_projectmax(mt: hl.MatrixTable, loc: hl.expr.StringExpression) -> hl.MatrixTable:
    """
    First pass of projectmax (returns aggregated MT with project_max field)

    :param MatrixTable mt: Input MT
    :param StringExpression loc: Column expression location of project ID (e.g. mt.meta.pid)
    :return: Frequency data with annotated project_max
    :rtype: MatrixTable
    """
    # TODO: add hom count
    mt = mt.annotate_cols(project=loc)
    agg_mt = mt.group_cols_by(mt.project).aggregate(ac=hl.agg.sum(mt.GT.n_alt_alleles()),
                                                    an=2 * hl.agg.count_where(hl.is_defined(mt.GT)),
                                                    hom=hl.agg.count_where(mt.GT.n_alt_alleles() == 2))
    agg_mt = agg_mt.annotate_entries(af=agg_mt.ac / agg_mt.an)
    return agg_mt.annotate_rows(project_max=hl.agg.take(hl.struct(project=agg_mt.project, ac=agg_mt.ac,
                                                                  af=agg_mt.af, an=agg_mt.an, hom=agg_mt.hom),
                                                        5, -agg_mt.af))


def annotation_type_is_numeric(t: Any) -> bool:
    """
    Given an annotation type, returns whether it is a numerical type or not.

    :param Type t: Type to test
    :return: If the input type is numeric
    :rtype: bool
    """
    return (isinstance(t, hl.tint32) or
            isinstance(t, hl.tint64) or
            isinstance(t, hl.tfloat32) or
            isinstance(t, hl.tfloat64)
            )


def annotation_type_in_vcf_info(t: Any) -> bool:
    """
    Given an annotation type, returns whether that type can be natively exported to a VCF INFO field.
    Note types that aren't natively exportable to VCF will be converted to String on export.

    :param Type t: Type to test
    :return: If the input type can be exported to VCF
    :rtype: bool
    """
    return (annotation_type_is_numeric(t) or
            isinstance(t, hl.tstr) or
            isinstance(t, hl.tarray) or
            isinstance(t, hl.tset) or
            isinstance(t, hl.tbool)
            )


def pc_project(mt: hl.MatrixTable, pc_loadings: hl.Table,
               loading_location: str = "loadings", af_location: str = "pca_af") -> hl.MatrixTable:
    """
    Projects samples in `mt` on PCs computed in `pc_mt`
    :param MatrixTable mt: MT containing the samples to project
    :param Table pc_loadings: MT containing the PC loadings for the variants
    :param str loading_location: Location of expression for loadings in `pc_loadings`
    :param str af_location: Location of expression for allele frequency in `pc_loadings`
    :return: MT with scores calculated from loadings
    """
    n_variants = mt.count_rows()

    mt = mt.annotate_rows(**pc_loadings[mt.locus, mt.alleles])
    mt = mt.filter_rows(hl.is_defined(mt[loading_location]) & hl.is_defined(mt[af_location]) &
                        (mt[af_location] > 0) & (mt[af_location] < 1))

    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt[af_location]) / hl.sqrt(n_variants * 2 * mt[af_location] * (1 - mt[af_location]))
    return mt.annotate_cols(pca_scores=hl.agg.array_sum(mt[loading_location] * gt_norm))


def filter_to_autosomes(mt: hl.MatrixTable) -> hl.MatrixTable:
    return hl.filter_intervals(mt, [hl.parse_locus_interval('1-22')])


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


def filter_low_conf_regions(mt: hl.MatrixTable, filter_lcr: bool = True, filter_decoy: bool = True,
                            filter_segdup: bool = True, high_conf_regions: Optional[List[str]] = None) -> hl.MatrixTable:
    """
    Filters low-confidence regions

    :param MatrixTable mt: MT to filter
    :param bool filter_lcr: Whether to filter LCR regions
    :param bool filter_decoy: Whether to filter decoy regions
    :param bool filter_segdup: Whether to filter Segdup regions
    :param list of str high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MT with low confidence regions removed
    :rtype: MatrixTable
    """
    from gnomad_hail.resources import lcr_intervals_path, decoy_intervals_path, segdup_intervals_path

    if filter_lcr:
        lcr = hl.import_locus_intervals(lcr_intervals_path)
        mt = mt.filter_rows(lcr[mt.locus], keep=False)

    if filter_decoy:
        decoy = hl.import_locus_intervals(decoy_intervals_path)
        mt = mt.filter_rows(decoy[mt.locus], keep=False)

    if filter_segdup:
        segdup = hl.import_locus_intervals(segdup_intervals_path)
        mt = mt.filter_rows(segdup[mt.locus], keep=False)

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            mt = mt.filter_rows(region[mt.locus], keep=True)

    return mt


def process_consequences(mt: hl.MatrixTable, vep_root: str = 'vep', penalize_flags: bool = True) -> hl.MatrixTable:
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq_by_gene, any_lof into [vep_root]

    :param MatrixTable mt: Input MT
    :param str vep_root: Root for vep annotation (probably vep)
    :param bool penalize_flags: Whether to penalize LOFTEE flagged variants, or treat them as equal to HC
    :return: MT with better formatted consequences
    :rtype: MatrixTable
    """
    csqs = hl.literal(CSQ_ORDER)
    csq_dict = hl.literal(dict(zip(CSQ_ORDER, range(len(CSQ_ORDER)))))

    def add_most_severe_consequence(tc: hl.expr.StructExpression) -> hl.expr.StructExpression:
        """
        Add most_severe_consequence annotation to transcript consequences
        This is for a given transcript, as there are often multiple annotations for a single transcript:
        e.g. splice_region_variant&intron_variant -> splice_region_variant
        """
        return tc.annotate(
            most_severe_consequence=csqs.find(lambda c: tc.consequence_terms.contains(c))
        )

    def find_worst_transcript_consequence(tcl: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
        """
        Gets worst transcript_consequence from an array of em
        """
        flag_score = 500
        no_flag_score = flag_score * (1 + penalize_flags)

        def csq_score(tc):
            return csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]
        tcl = tcl.map(lambda tc: tc.annotate(
            csq_score=hl.case(missing_false=True)
            .when((tc.lof == 'HC') & (tc.lof_flags == ''), csq_score(tc) - no_flag_score)
            .when((tc.lof == 'HC') & (tc.lof_flags != ''), csq_score(tc) - flag_score)
            .when(tc.lof == 'LC', csq_score(tc) - 10)
            .when(tc.polyphen_prediction == 'probably_damaging', csq_score(tc) - 0.5)
            .when(tc.polyphen_prediction == 'possibly_damaging', csq_score(tc) - 0.25)
            .when(tc.polyphen_prediction == 'benign', csq_score(tc) - 0.1)
            .default(csq_score(tc))
        ))
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])

    transcript_csqs = mt[vep_root].transcript_consequences.map(add_most_severe_consequence)

    gene_dict = transcript_csqs.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene = gene_dict.map_values(find_worst_transcript_consequence)
    sorted_scores = hl.sorted(worst_csq_gene.values(), key=lambda tc: tc.csq_score)
    lowest_score = hl.or_missing(hl.len(sorted_scores) > 0, sorted_scores[0].csq_score)
    gene_with_worst_csq = sorted_scores.filter(lambda tc: tc.csq_score == lowest_score).map(lambda tc: tc.gene_symbol)
    ensg_with_worst_csq = sorted_scores.filter(lambda tc: tc.csq_score == lowest_score).map(lambda tc: tc.gene_id)

    vep_data = mt[vep_root].annotate(transcript_consequences=transcript_csqs,
                                     worst_csq_by_gene=worst_csq_gene,
                                     any_lof=hl.any(lambda x: x.lof == 'HC', worst_csq_gene.values()),
                                     gene_with_most_severe_csq=gene_with_worst_csq,
                                     ensg_with_most_severe_csq=ensg_with_worst_csq)

    return mt.annotate_rows(**{vep_root: vep_data})


def filter_vep_to_canonical_transcripts(mt: hl.MatrixTable, vep_root: str = 'vep') -> hl.MatrixTable:
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.canonical == 1)
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return mt.annotate_rows(**{vep_root: vep_data})


def filter_vep_to_synonymous_variants(mt: hl.MatrixTable, vep_root: str = 'vep') -> hl.MatrixTable:
    synonymous = mt[vep_root].transcript_consequences.filter(lambda csq: csq.most_severe_consequence == "synonymous_variant")
    vep_data = mt[vep_root].annotate(transcript_consequences=synonymous)
    return mt.annotate_rows(**{vep_root: vep_data})


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
