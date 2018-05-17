
import hail as hl
from hail.expr.expressions import *
from collections import defaultdict, namedtuple, OrderedDict
from typing import *


def unphase_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Generate unphased version of MatrixTable (assumes call is in mt.GT and is diploid or haploid only)
    """
    return mt.annotate_entries(GT=hl.case()
                               .when(mt.GT.is_diploid(), hl.call(mt.GT[0], mt.GT[1], phased=False))
                               .when(mt.GT.is_haploid(), hl.call(mt.GT[0], phased=False))
                               .default(hl.null(hl.tcall))
    )


def filter_to_autosomes(mt: hl.MatrixTable) -> hl.MatrixTable:
    return hl.filter_intervals(mt, [hl.parse_locus_interval('1-22')])


def write_temp_gcs(t: Union[hl.MatrixTable, hl.Table], gcs_path: str,
                   overwrite: bool = False, temp_path: str = '/tmp.h') -> None:
    t.write(temp_path, overwrite=True)
    t = hl.read_matrix_table(temp_path) if isinstance(t, hl.MatrixTable) else hl.read_table(temp_path)
    t.write(gcs_path, overwrite=overwrite)


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


def split_multi_dynamic(t: Union[hl.MatrixTable, hl.Table], keep_star: bool = False,
                        left_aligned: bool = True, vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    """
    Splits MatrixTable based on entry fields found. Downcodes whatever it can. Supported so far:
    GT, DP, AD, PL, GQ
    PGT, PID
    ADALL

    :param MatrixTable t: Input MatrixTable
    :param bool keep_star: whether to keep star alleles (passed to SplitMulti)
    :param bool left_aligned: whether matrix table is already left_aligned (passed to SplitMulti)
    :param str vep_root: If provided and exists in t, splits multi-allelics in VEP field properly (default "vep")
    :return: Split MatrixTable
    :rtype: MatrixTable
    """
    rows = list(t.row)

    if isinstance(t, hl.Table):
        t = t.annotate(a_index=hl.range(0, hl.len(t.alleles) - 1)).explode('a_index')
        update_rows_expr = {'alleles': [t.alleles[0], t.alleles[t.a_index]]}
        if vep_root in rows:
            update_rows_expr[vep_root] = t[vep_root].annotate(
                intergenic_consequences=t[vep_root].intergenic_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index),
                motif_feature_consequences=t[vep_root].motif_feature_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index),
                regulatory_feature_consequences=t[vep_root].motif_feature_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index),
                transcript_consequences=t[vep_root].transcript_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index))
        return t.annotate(**update_rows_expr)  # Note: does not minrep at the moment

    fields = list(t.entry)
    sm = hl.SplitMulti(t, keep_star=keep_star, left_aligned=left_aligned)
    update_rows_expr = {'a_index': sm.a_index(), 'was_split': sm.was_split()}
    if vep_root in rows:
        update_rows_expr[vep_root] = t[vep_root].annotate(
            intergenic_consequences=t[vep_root].intergenic_consequences.filter(
                lambda csq: csq.allele_num == sm.a_index()),
            motif_feature_consequences=t[vep_root].motif_feature_consequences.filter(
                lambda csq: csq.allele_num == sm.a_index()),
            regulatory_feature_consequences=t[vep_root].motif_feature_consequences.filter(
                lambda csq: csq.allele_num == sm.a_index()),
            transcript_consequences=t[vep_root].transcript_consequences.filter(
                lambda csq: csq.allele_num == sm.a_index()))
    sm.update_rows(**update_rows_expr)
    expression = {}

    # HTS/standard
    if 'GT' in fields:
        expression['GT'] = hl.downcode(t.GT, sm.a_index())
    if 'DP' in fields:
        expression['DP'] = t.DP
    if 'AD' in fields:
        expression['AD'] = hl.or_missing(hl.is_defined(t.AD),
                                         [hl.sum(t.AD) - t.AD[sm.a_index()], t.AD[sm.a_index()]])
    if 'PL' in fields:
        pl = hl.or_missing(
            hl.is_defined(t.PL),
            (hl.range(0, 3).map(lambda i:
                                hl.min((hl.range(0, hl.triangle(t.alleles.length()))
                                        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j),
                                                                      sm.a_index()) == hl.unphased_diploid_gt_index_call(i)
                                                ).map(lambda j: t.PL[j]))))))
        expression['PL'] = pl
        if 'GQ' in fields:
            expression['GQ'] = hl.gq_from_pl(pl)
    else:
        if 'GQ' in fields:
            expression['GQ'] = t.GQ

    # Phased data
    if 'PGT' in fields:
        expression['PGT'] = hl.downcode(t.PGT, sm.a_index())
    if 'PID' in fields:
        expression['PID'] = t.PID

    # Custom data
    if 'ADALL' in fields:  # found in NA12878
        expression['ADALL'] = hl.or_missing(hl.is_defined(t.ADALL),
                                            [hl.sum(t.ADALL) - t.ADALL[sm.a_index()], t.ADALL[sm.a_index()]])

    sm.update_entries(**expression)
    return sm.result()


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
        mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)

    if filter_decoy:
        decoy = hl.import_bed(decoy_intervals_path)
        mt = mt.filter_rows(hl.is_defined(decoy[mt.locus]), keep=False)

    if filter_segdup:
        segdup = hl.import_bed(segdup_intervals_path)
        mt = mt.filter_rows(hl.is_defined(segdup[mt.locus]), keep=False)

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            mt = mt.filter_rows(hl.is_defined(region[mt.locus]), keep=True)

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
    from .constants import CSQ_ORDER

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
                                     worst_consequence_term=csqs.find(lambda c: transcript_csqs.map(lambda csq: csq.most_severe_consequence).contains(c)),
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
