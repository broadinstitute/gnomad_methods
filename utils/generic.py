
import hail as hl
from hail.expr.expressions import *
from collections import defaultdict, namedtuple, OrderedDict
from typing import *
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import random
import warnings
import uuid


def unphase_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Generate unphased version of MatrixTable (assumes call is in mt.GT and is diploid or haploid only)
    """
    return mt.annotate_entries(GT=hl.case()
                               .when(mt.GT.is_diploid(), hl.call(mt.GT[0], mt.GT[1], phased=False))
                               .when(mt.GT.is_haploid(), hl.call(mt.GT[0], phased=False))
                               .default(hl.null(hl.tcall))
    )


def filter_to_autosomes(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    autosomes = hl.parse_locus_interval('1-22')
    if isinstance(t, hl.MatrixTable):
        return hl.filter_intervals(t, [autosomes])
    else:
        return t.filter(autosomes.contains(t.locus))


def write_temp_gcs(t: Union[hl.MatrixTable, hl.Table], gcs_path: str,
                   overwrite: bool = False, temp_path: Optional[str] = None) -> None:
    if not temp_path:
        temp_path = f'/tmp_{uuid.uuid4()}.h'
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
    Splits MatrixTable based on entry fields found. All `Call` fields (e.g. `GT`, `PGT`) are downcoded whatever it can.
    In addition to all `Call` fields, the following entry fields are split correctly:
    DP, AD, PL, GQ
    PID
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
        if {'locus', 'alleles'} != set(t.key):
            raise Exception('Table not keyed by locus and alleles - cannot split_multi')

        a_index_expr = hl.range(1, hl.len(t.alleles))
        if not keep_star:
            a_index_expr = a_index_expr.filter(lambda i: t.alleles[i] != '*')

        t = t.annotate(a_index=a_index_expr, was_split=hl.len(t.alleles) > 2)
        t = t.explode('a_index')

        new_locus_alleles = hl.min_rep(t.locus, hl.array([t.alleles[0], t.alleles[t.a_index]]))
        t = t.key_by(**{'locus': new_locus_alleles[0], 'alleles': new_locus_alleles[1]})

        if vep_root in rows:
            t = t.annotate(**{vep_root: t[vep_root].annotate(
                intergenic_consequences=t[vep_root].intergenic_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index),
                motif_feature_consequences=t[vep_root].motif_feature_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index),
                regulatory_feature_consequences=t[vep_root].motif_feature_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index),
                transcript_consequences=t[vep_root].transcript_consequences.filter(
                    lambda csq: csq.allele_num == t.a_index)
            )})

        return t

    call_fields = [f[0] for f in t.entry.items() if f[1].dtype == hl.tcall]
    non_call_fields = [f for f in list(t.entry) if f not in call_fields]

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

    # Split all Call fields
    for f in call_fields:
        expression[f] = hl.downcode(t[f], sm.a_index())

    # HTS/standard
    if 'DP' in non_call_fields:
        expression['DP'] = t.DP
    if 'AD' in non_call_fields:
        expression['AD'] = hl.or_missing(hl.is_defined(t.AD),
                                         [hl.sum(t.AD) - t.AD[sm.a_index()], t.AD[sm.a_index()]])
    if 'PL' in non_call_fields:
        pl = hl.or_missing(
            hl.is_defined(t.PL),
            (hl.range(0, 3).map(lambda i:
                                hl.min((hl.range(0, hl.triangle(t.alleles.length()))
                                        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j),
                                                                      sm.a_index()) == hl.unphased_diploid_gt_index_call(i)
                                                ).map(lambda j: t.PL[j]))))))
        expression['PL'] = pl
        if 'GQ' in non_call_fields:
            expression['GQ'] = hl.gq_from_pl(pl)
    else:
        if 'GQ' in non_call_fields:
            expression['GQ'] = t.GQ

    # Phased data
    if 'PID' in non_call_fields:
        expression['PID'] = t.PID

    # Custom data
    if 'ADALL' in non_call_fields:  # found in NA12878
        expression['ADALL'] = hl.or_missing(hl.is_defined(t.ADALL),
                                            [hl.sum(t.ADALL) - t.ADALL[sm.a_index()], t.ADALL[sm.a_index()]])

    sm.update_entries(**expression)
    return sm.result()


def pc_project(
        mt: hl.MatrixTable,
        loadings_ht: hl.Table,
        loading_location: str = "loadings",
        af_location: str = "pca_af"
) -> hl.Table:
    """
    Projects samples in `mt` on pre-computed PCs.

    :param MatrixTable mt: MT containing the samples to project
    :param Table loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param str loading_location: Location of expression for loadings in `loadings_ht`
    :param str af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    :rtype: Table
    """

    n_variants = loadings_ht.count()

    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location]
    )

    mt = mt.filter_rows(hl.is_defined(mt.pca_loadings) & hl.is_defined(mt.pca_af) &
                        (mt.pca_af > 0) & (mt.pca_af < 1))

    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(n_variants * 2 * mt.pca_af * (1 - mt.pca_af))

    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))

    return mt.cols().select('scores')


def sample_pcs_uniformly(scores_table: hl.Table, num_pcs: int = 5, num_bins: int = 10, num_per_bin: int = 20) -> hl.Table:
    """
    Sample somewhat uniformly in num_pcs-dimensional PC space, by:
    1. Binning each PC axis into num_bins bins, creating an array of num_pcs with num_bins possible values (total of num_bins ^ num_pcs sectors)
    2. For each k-dimensional sector, take up to num_per_bin samples
    Max number of samples return is num_per_bin * num_bins ^ num_pcs, but in practice, typically much fewer (corners of PC space are sparse)

    Assumes your scores are in scores_table.scores (and sample stored in `s`)
    """
    ranges = scores_table.aggregate([hl.agg.stats(scores_table.scores[i]) for i in range(num_pcs)])
    ranges = [x.annotate(r=x.max - x.min) for x in ranges]
    ranges = hl.literal([x.annotate(step=x.r / num_bins) for x in ranges])

    scores_table = scores_table.annotate(
        scores_bin=hl.range(0, num_pcs).map(lambda i:
                                            hl.int((scores_table.scores[i] - ranges[i].min) / ranges[i].step))
    )
    per_bin = scores_table.group_by(scores_table.scores_bin).aggregate(s=hl.agg.take(scores_table.s, num_per_bin))
    return per_bin.explode(per_bin.s)


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


def add_most_severe_consequence_to_consequence(tc: hl.expr.StructExpression) -> hl.expr.StructExpression:

    """
    Add most_severe_consequence annotation to transcript consequences
    This is for a given transcript, as there are often multiple annotations for a single transcript:
    e.g. splice_region_variant&intron_variant -> splice_region_variant
    """
    from .constants import CSQ_ORDER

    csqs = hl.literal(CSQ_ORDER)

    return tc.annotate(
        most_severe_consequence=csqs.find(lambda c: tc.consequence_terms.contains(c))
    )


def process_consequences(mt: Union[hl.MatrixTable, hl.Table], vep_root: str = 'vep',
                         penalize_flags: bool = True) -> Union[hl.MatrixTable, hl.Table]:
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

    transcript_csqs = mt[vep_root].transcript_consequences.map(add_most_severe_consequence_to_consequence)

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

    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})


def filter_vep_to_canonical_transcripts(mt: Union[hl.MatrixTable, hl.Table],
                                        vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.canonical == 1)
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})


def filter_vep_to_synonymous_variants(mt: Union[hl.MatrixTable, hl.Table],
                                      vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    synonymous = mt[vep_root].transcript_consequences.filter(lambda csq: csq.most_severe_consequence == "synonymous_variant")
    vep_data = mt[vep_root].annotate(transcript_consequences=synonymous)
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})


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


def phase_by_transmission(
        locus: hl.expr.LocusExpression,
        alleles: hl.expr.ArrayExpression,
        proband_call: hl.expr.CallExpression,
        father_call: hl.expr.CallExpression,
        mother_call: hl.expr.CallExpression
) -> hl.expr.ArrayExpression:
    """
    Phases genotype calls in a trio based allele transmission.

    In the phased calls returned, the order is as follows:
    * Proband: father_allele | mother_allele
    * Parents: transmitted_allele | untransmitted_allele

    Phasing of sex chromosomes:
    Sex chromosomes of male individuals should be haploid to be phased correctly.
    If `proband_call` is diploid on non-par regions of the sex chromosomes, it is assumed to be female.

    Returns `NA` when genotype calls cannot be phased.
    The following genotype calls combinations cannot be phased by transmission:
    1. One of the calls in the trio is missing
    2. The proband genotype cannot be obtained from the parents alleles (Mendelian violation)
    3. All individuals of the trio are heterozygous for the same two alleles
    4. Father is diploid on non-PAR region of X or Y
    5. Proband is diploid on non-PAR region of Y

    In addition, individual phased genotype calls are returned as missing in the following situations:
    1. All mother genotype calls non-PAR region of Y
    2. Diploid father genotype calls on non-PAR region of X for a male proband (proband and mother are still phased as father doesn't participate in allele transmission)

    :param LocusExpression locus: Locus in the trio MatrixTable
    :param ArrayExpression alleles: Alleles in the trio MatrixTable
    :param CallExpression proband_call: Input proband genotype call
    :param CallExpression father_call: Input father genotype call
    :param CallExpression mother_call: Input mother genotype call
    :return: Array containing: phased proband call, phased father call, phased mother call
    :rtype: ArrayExpression
    """

    def call_to_one_hot_alleles_array(call: hl.expr.CallExpression, alleles: hl.expr.ArrayExpression) -> hl.expr.ArrayExpression:
        """
        Get the set of all different one-hot-encoded allele-vectors in a genotype call.
        It is returned as an ordered array where the first vector corresponds to the first allele,
        and the second vector (only present if het) the second allele.

        :param CallExpression call: genotype
        :param ArrayExpression alleles: Alleles at the site
        :return: Array of one-hot-encoded alleles
        :rtype: ArrayExpression
        """
        return hl.cond(
            call.is_het(),
            hl.array([
                hl.call(call[0]).one_hot_alleles(alleles),
                hl.call(call[1]).one_hot_alleles(alleles),
            ]),
            hl.array([hl.call(call[0]).one_hot_alleles(alleles)])
        )

    def phase_parent_call(call: hl.expr.CallExpression, transmitted_allele_index: int):
        """
        Given a genotype and which allele was transmitted to the offspring, returns the parent phased genotype.

        :param CallExpression call: Parent genotype
        :param int transmitted_allele_index: index of transmitted allele (0 or 1)
        :return: Phased parent genotype
        :rtype: CallExpression
        """
        return hl.call(
            call[transmitted_allele_index],
            call[hl.int(transmitted_allele_index == 0)],
            phased=True
        )

    def phase_diploid_proband(
            locus: hl.expr.LocusExpression,
            alleles: hl.expr.ArrayExpression,
            proband_call: hl.expr.CallExpression,
            father_call: hl.expr.CallExpression,
            mother_call: hl.expr.CallExpression
    ) -> hl.expr.ArrayExpression:
        """
        Returns phased genotype calls in the case of a diploid proband
        (autosomes, PAR regions of sex chromosomes or non-PAR regions of a female proband)

        :param LocusExpression locus: Locus in the trio MatrixTable
        :param ArrayExpression alleles: Alleles in the trio MatrixTable
        :param CallExpression proband_call: Input proband genotype call
        :param CallExpression father_call: Input father genotype call
        :param CallExpression mother_call: Input mother genotype call
        :return: Array containing: phased proband call, phased father call, phased mother call
        :rtype: ArrayExpression
        """

        proband_v = proband_call.one_hot_alleles(alleles)
        father_v = hl.cond(
            locus.in_x_nonpar() | locus.in_y_nonpar(),
            hl.or_missing(father_call.is_haploid(), hl.array([father_call.one_hot_alleles(alleles)])),
            call_to_one_hot_alleles_array(father_call, alleles)
        )
        mother_v = call_to_one_hot_alleles_array(mother_call, alleles)

        combinations = hl.flatmap(
            lambda f:
            hl.zip_with_index(mother_v)
                .filter(lambda m: m[1] + f[1] == proband_v)
                .map(lambda m: hl.struct(m=m[0], f=f[0])),
            hl.zip_with_index(father_v)
        )

        return (
            hl.cond(
                hl.is_defined(combinations) & (hl.len(combinations) == 1),
                hl.array([
                    hl.call(father_call[combinations[0].f], mother_call[combinations[0].m], phased=True),
                    hl.cond(father_call.is_haploid(), hl.call(father_call[0], phased=True), phase_parent_call(father_call, combinations[0].f)),
                    phase_parent_call(mother_call, combinations[0].m)
                ]),
                hl.null(hl.tarray(hl.tcall))
            )
        )

    def phase_haploid_proband_x_nonpar(
            proband_call: hl.expr.CallExpression,
            father_call: hl.expr.CallExpression,
            mother_call: hl.expr.CallExpression
    ) -> hl.expr.ArrayExpression:
        """
        Returns phased genotype calls in the case of a haploid proband in the non-PAR region of X

        :param CallExpression proband_call: Input proband genotype call
        :param CallExpression father_call: Input father genotype call
        :param CallExpression mother_call: Input mother genotype call
        :return: Array containing: phased proband call, phased father call, phased mother call
        :rtype: ArrayExpression
        """

        transmitted_allele = hl.zip_with_index(hl.array([mother_call[0], mother_call[1]])).find(lambda m: m[1] == proband_call[0])
        return hl.cond(
            hl.is_defined(transmitted_allele),
            hl.array([
                hl.call(proband_call[0], phased=True),
                hl.or_missing(father_call.is_haploid(), hl.call(father_call[0], phased=True)),
                phase_parent_call(mother_call, transmitted_allele[0])
            ]),
            hl.null(hl.tarray(hl.tcall))
        )

    def phase_y_nonpar(
            proband_call: hl.expr.CallExpression,
            father_call: hl.expr.CallExpression,
    ) -> hl.expr.ArrayExpression:
        """
        Returns phased genotype calls in the non-PAR region of Y (requires both father and proband to be haploid to return phase)

        :param CallExpression proband_call: Input proband genotype call
        :param CallExpression father_call: Input father genotype call
        :return: Array containing: phased proband call, phased father call, phased mother call
        :rtype: ArrayExpression
        """
        return hl.cond(
            proband_call.is_haploid() & father_call.is_haploid() & (father_call[0] == proband_call[0]),
            hl.array([
                hl.call(proband_call[0], phased=True),
                hl.call(father_call[0], phased=True),
                hl.null(hl.tcall)
            ]),
            hl.null(hl.tarray(hl.tcall))
        )

    return (
        hl.case()
            .when(locus.in_x_nonpar() & proband_call.is_haploid(), phase_haploid_proband_x_nonpar(proband_call, father_call, mother_call))
            .when(locus.in_y_nonpar(), phase_y_nonpar(proband_call, father_call))
            .when(proband_call.is_diploid(), phase_diploid_proband(locus, alleles, proband_call, father_call, mother_call))
            .default(hl.null(hl.tarray(hl.tcall)))
    )


def phase_trio_matrix_by_transmission(tm: hl.MatrixTable, call_field: str = 'GT' ,phased_call_field: str = 'PBT_GT') -> hl.MatrixTable:
    """
        Adds a phased genoype entry to a trio MatrixTable based allele transmission in the trio.
        Uses only a `Call` field to phase and only phases when all 3 members of the trio are present and have a call.

        In the phased genotypes, the order is as follows:
        * Proband: father_allele | mother_allele
        * Parents: transmitted_allele | untransmitted_allele

        Phasing of sex chromosomes:
        Sex chromosomes of male individuals should be haploid to be phased correctly.
        If a proband is diploid on non-par regions of the sex chromosomes, it is assumed to be female.

        Genotypes that cannot be phased are set to `NA`.
        The following genotype calls combinations cannot be phased by transmission (all trio members phased calls set to missing):
        1. One of the calls in the trio is missing
        2. The proband genotype cannot be obtained from the parents alleles (Mendelian violation)
        3. All individuals of the trio are heterozygous for the same two alleles
        4. Father is diploid on non-PAR region of X or Y
        5. Proband is diploid on non-PAR region of Y

        In addition, individual phased genotype calls are returned as missing in the following situations:
        1. All mother genotype calls non-PAR region of Y
        2. Diploid father genotype calls on non-PAR region of X for a male proband (proband and mother are still phased as father doesn't participate in allele transmission)


        Typical usage:
        ```
            trio_matrix = hl.trio_matrix(mt, ped)
            phased_trio_matrix = phase_trio_matrix_by_transmission(trio_matrix)
        ```

        :param MatrixTable tm: Trio MatrixTable (entries should be a Struct with `proband_entry`, `mother_entry` and `father_entry` present)
        :param str call_field: genotype field name to phase
        :param str phased_call_field: name for the phased genotype field
        :return: trio MatrixTable entry with additional phased genotype field for each individual
        :rtype: MatrixTable
        """

    tm = tm.annotate_entries(
        __phased_GT=phase_by_transmission(
            tm.locus,
            tm.alleles,
            tm.proband_entry[call_field],
            tm.father_entry[call_field],
            tm.mother_entry[call_field]
        )
    )

    return tm.select_entries(
        proband_entry=hl.struct(
            **tm.proband_entry,
            **{phased_call_field: tm.__phased_GT[0]}
        ),
        father_entry=hl.struct(
            **tm.father_entry,
            **{phased_call_field: tm.__phased_GT[1]}
        ),
        mother_entry=hl.struct(
            **tm.mother_entry,
            **{phased_call_field: tm.__phased_GT[2]}
        )
    )


def explode_trio_matrix(tm: hl.MatrixTable, col_keys: List[str] = ['s']) -> hl.MatrixTable:
    """

    Splits a trio MatrixTable back into a sample MatrixTable.
    It assumes that the input MatrixTable schema

    :param MatrixTable tm: Input trio MatrixTable
    :param list of str col_keys: Column keys for the sample MatrixTable
    :return: Sample MatrixTable
    :rtype: MatrixTable
    """
    tm = tm.select_entries(
        __trio_entries=hl.array([tm.proband_entry, tm.father_entry, tm.mother_entry])
    )

    tm = tm.select_cols(
        __trio_members=hl.zip_with_index(hl.array([tm.proband, tm.father, tm.mother]))
    )
    mt = tm.explode_cols(tm.__trio_members)

    mt = mt.select_entries(
        **mt.__trio_entries[mt.__trio_members[0]]
    )

    mt = mt.key_cols_by()
    mt = mt.select_cols(**mt.__trio_members[1])

    if col_keys:
        mt = mt.key_cols_by(*col_keys)

    return mt


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
            relative_pairs: List[Tuple[str, str]],
            sex: Dict[str, bool]
    ) -> Union[Tuple[str, str], None]:
        """
        Given a list of possible parents for a sample (first degree relatives with low ibd2),
        looks for a single pair of samples that are unrelated with different sexes.
        If a single pair is found, return the pair (father, mother)

        :param list of str possible_parents: Possible parents
        :param list of (str, str) relative_pairs: Pairs of relatives, used to check that parents aren't related with each other
        :param dict of str -> bool sex: Dict mapping samples to their sex (True = female, False = male, None or missing = unknown)
        :return: (father, mother) if found, `None` otherwise
        :rtype: (str, str) or None
        """

        parents = []
        while len(possible_parents) > 1:
            p1 = possible_parents.pop()
            for p2 in possible_parents:
                if tuple(sorted((p1,p2))) not in relative_pairs:
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

            parents = get_parents(possible_parents, list(ibd2.keys()), sex)

            if parents is not None:
                trios.append(hl.Trio(s=s,
                                     fam_id=str(fam_id),
                                     pat_id=parents[0],
                                     mat_id=parents[1],
                                     is_female=sex.get(s)))

        fam_id += 1

    return hl.Pedigree(trios)


def expand_pd_array_col(
        df: pd.DataFrame,
        array_col: str,
        num_out_cols: int = 0,
        out_cols_prefix=None,
        out_1based_indexing: bool = True
) -> pd.DataFrame:
    """
    Expands a Dataframe column containing an array into multiple columns.

    :param DataFrame df: input dataframe
    :param str array_col: Column containing the array
    :param int num_out_cols: Number of output columns. If set, only the `n_out_cols` first elements of the array column are output.
                             If <1, the number of output columns is equal to the length of the shortest array in `array_col`
    :param out_cols_prefix: Prefix for the output columns (uses `array_col` as the prefix unless set)
    :param bool out_1based_indexing: If set, the output column names indexes start at 1. Otherwise they start at 0.
    :return: dataframe with expanded columns
    :rtype: DataFrame
    """

    if out_cols_prefix is None:
        out_cols_prefix = array_col

    if num_out_cols < 1:
        num_out_cols = min([len(x) for x in df[array_col].values.tolist()])

    cols = ['{}{}'.format(out_cols_prefix, i + out_1based_indexing) for i in range(num_out_cols)]
    df[cols] = pd.DataFrame(df[array_col].values.tolist())[list(range(num_out_cols))]

    return df


def assign_population_pcs(
        pop_pc_pd: pd.DataFrame,
        pc_cols: List[str],
        known_col: str = 'known_pop',
        fit: RandomForestClassifier = None,
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth'
) -> Tuple[pd.DataFrame, RandomForestClassifier]:
    """

    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.

    Note that if PCs come from Hail, they will be stored in a single column named `scores` by default.
    The `expand_pd_array_col` can be used to expand this `scores` column into multiple `PC` columns:

    ```
        data = pca_ht.to_pandas() #Load PCA results
        data = expand_pd_array_col(data, 'scores', 10 , 'PC') # Expand `scores` column for 10 PCs into columns `PC1` ... `PC10`
        results, rf_model = assign_population_pcs(
            data,
            ['PC{}'.format(i + 1) for i in range(10)]),
            ...
        )

    ```


    :param Table pop_pc_pd: Pandas dataframe containing population PCs as well as a column with population labels
    :param list of str pc_cols: Columns storing the PCs to use
    :param str known_col: Column storing the known population labels
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Dataframe containing sample IDs and imputed population labels, trained random forest model
    :rtype: DataFrame, RandomForestClassifier
    """

    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

    N = len(train_data)

    # Split training data into subsamples for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(list(range(0, N)), int(N * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit['s']]
        evaluate_fit = train_data.loc[~train_data['s'].isin(fit_samples)]

        # Train RF
        training_set_known_labels = train_fit[known_col].as_matrix()
        training_set_pcs = train_fit[pc_cols].as_matrix()
        evaluation_set_pcs = evaluate_fit[pc_cols].as_matrix()

        pop_clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print('Random forest feature importances are as follows: {}'.format(pop_clf.feature_importances_))

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(len(predictions))
        print('Estimated error rate for RF model is {}'.format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].as_matrix())
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].as_matrix())
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])
    pop_pc_pd = pd.concat([pop_pc_pd, probs], axis=1)
    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label
    pop_pc_pd = pop_pc_pd.drop(pc_cols, axis='columns')

    return pop_pc_pd, pop_clf
