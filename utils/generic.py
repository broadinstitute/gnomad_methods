import hail as hl
from hail.expr.expressions import *
from collections import defaultdict, namedtuple, OrderedDict, Counter
from typing import *
import pandas as pd
import random
import warnings
import uuid
import operator
import functools
from hail.utils.misc import divide_null
from .gnomad_functions import logger
import os

INFO_VCF_AS_PIPE_DELIMITED_FIELDS = ['AS_QUALapprox', 'AS_VarDP', 'AS_MQ_DP', 'AS_RAW_MQ', 'AS_SB_TABLE']

def file_exists(fname: str) -> bool:
    """
    Check whether a file exists.
    Supports either local or Google cloud (gs://) paths.
    If the file is a Hail file (.ht, .mt extensions), it checks that _SUCCESS is present.

    :param str fname: File name
    :return: Whether the file exists
    :rtype: bool
    """
    fext = os.path.splitext(fname)[1]
    if fext in ['.ht', '.mt']:
        fname += '/_SUCCESS'
    if fname.startswith('gs://'):
        return hl.hadoop_exists(fname)
    else:
        return os.path.isfile(fname)


def unphase_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Generate unphased version of MatrixTable (assumes call is in mt.GT and is diploid or haploid only)
    """
    return mt.annotate_entries(GT=hl.case()
                               .when(mt.GT.is_diploid(), hl.call(mt.GT[0], mt.GT[1], phased=False))
                               .when(mt.GT.is_haploid(), hl.call(mt.GT[0], phased=False))
                               .default(hl.null(hl.tcall))
                               )


def get_reference_genome(locus: Union[hl.expr.LocusExpression, hl.expr.IntervalExpression]) -> hl.ReferenceGenome:
    """
    Returns the reference genome associated with the input Locus expression

    :param LocusExpression or IntervalExpression locus: Input locus
    :return: Reference genome
    :rtype: ReferenceGenome
    """
    if isinstance(locus, hl.expr.LocusExpression):
        return locus.dtype.reference_genome
    assert (isinstance(locus, hl.expr.IntervalExpression))
    return locus.dtype.point_type.reference_genome


def flip_base(base: str) -> str:
    """
    Returns the complement of a base

    :param str base: Base to be flipped
    :return: Complement of input base
    :rtype: str
    """
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(base))


def reverse_complement_bases(bases: hl.expr.StringExpression) -> hl.expr.StringExpression:
    """
    Returns the reverse complement of a sequence

    :param StringExpression bases: Sequence to be flipped
    :return: Reverse complement of input sequence
    :rtype: StringExpression
    """
    return hl.delimit(hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), '')


def filter_to_autosomes(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters the Table or MatrixTable to autosomes only.
    This assumes that the input contains a field named `locus` of type Locus

    :param MatrixTable or Table t: Input MT/HT
    :return:  MT/HT autosomes
    :rtype: MatrixTable or Table
    """
    reference = get_reference_genome(t.locus)
    autosomes = hl.parse_locus_interval(f'{reference.contigs[0]}-{reference.contigs[21]}', reference_genome=reference)
    return hl.filter_intervals(t, [autosomes])


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


def filter_low_conf_regions(mt: Union[hl.MatrixTable, hl.Table], filter_lcr: bool = True, filter_decoy: bool = True,
                            filter_segdup: bool = True, filter_exome_low_coverage_regions: bool = False,
                            high_conf_regions: Optional[List[str]] = None) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters low-confidence regions

    :param MatrixTable or Table mt: MatrixTable or Table to filter
    :param bool filter_lcr: Whether to filter LCR regions
    :param bool filter_decoy: Whether to filter decoy regions
    :param bool filter_segdup: Whether to filter Segdup regions
    :param bool filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
    :param list of str high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MatrixTable or Table with low confidence regions removed
    :rtype: MatrixTable or Table
    """
    from gnomad_hail.resources import lcr_intervals_path, decoy_intervals_path, segdup_intervals_path, high_coverage_intervals_path

    criteria = []
    if filter_lcr:
        lcr = hl.import_locus_intervals(lcr_intervals_path)
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_decoy:
        decoy = hl.import_bed(decoy_intervals_path)
        criteria.append(hl.is_missing(decoy[mt.locus]))

    if filter_segdup:
        segdup = hl.import_bed(segdup_intervals_path)
        criteria.append(hl.is_missing(segdup[mt.locus]))

    if filter_exome_low_coverage_regions:
        high_cov = hl.import_locus_intervals(high_coverage_intervals_path)
        criteria.append(hl.is_missing(high_cov[mt.locus]))

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            criteria.append(hl.is_defined(region[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


def vep_or_lookup_vep(ht, reference_vep_ht=None, reference=None, vep_config=None):
    """
    VEP a table, or lookup variants in a reference database

    :param Table ht: Input Table
    :param Table reference_vep_ht: A reference database with VEP annotations (must be in top-level `vep`)
    :param ReferenceGenome reference: If reference_vep_ht is not specified, find a suitable one in reference (if None, grabs from hl.default_reference)
    :param str vep_config: vep_config to pass to hl.vep (if None, a suitable one for `reference` is chosen)
    :return: VEPped Table
    :rtype: Table
    """
    from gnomad_hail.resources.basics import vep_config_path, context_ht_path

    VEP_REFERENCES = {
        'GRCh37': context_ht_path(),
        'GRCh38': '',
    }

    if reference_vep_ht is None:
        if reference is None:
            reference = hl.default_reference().name

        if reference not in VEP_REFERENCES:
            raise ValueError(f'vep_or_lookup_vep got {reference}. Expected one of {", ".join(VEP_REFERENCES.keys())}')

        reference_vep_ht = hl.read_table(VEP_REFERENCES[reference])

    ht = ht.annotate(vep=reference_vep_ht[ht.key].vep)

    vep_ht = ht.filter(hl.is_defined(ht.vep))
    revep_ht = ht.filter(hl.is_missing(ht.vep))

    if vep_config is None:
        vep_config = vep_config_path(reference)

    revep_ht = hl.vep(revep_ht, vep_config)

    return vep_ht.union(revep_ht)


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
            .when(tc.lof == 'OS', csq_score(tc) - 20)
            .when(tc.lof == 'LC', csq_score(tc) - 10)
            .when(tc.polyphen_prediction == 'probably_damaging', csq_score(tc) - 0.5)
            .when(tc.polyphen_prediction == 'possibly_damaging', csq_score(tc) - 0.25)
            .when(tc.polyphen_prediction == 'benign', csq_score(tc) - 0.1)
            .default(csq_score(tc))
        ))
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])

    transcript_csqs = mt[vep_root].transcript_consequences.map(add_most_severe_consequence_to_consequence)

    gene_dict = transcript_csqs.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene = gene_dict.map_values(find_worst_transcript_consequence).values()
    sorted_scores = hl.sorted(worst_csq_gene, key=lambda tc: tc.csq_score)

    canonical = transcript_csqs.filter(lambda csq: csq.canonical == 1)
    gene_canonical_dict = canonical.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene_canonical = gene_canonical_dict.map_values(find_worst_transcript_consequence).values()
    sorted_canonical_scores = hl.sorted(worst_csq_gene_canonical, key=lambda tc: tc.csq_score)

    vep_data = mt[vep_root].annotate(transcript_consequences=transcript_csqs,
                                     worst_consequence_term=csqs.find(lambda c: transcript_csqs.map(lambda csq: csq.most_severe_consequence).contains(c)),
                                     worst_csq_by_gene=sorted_scores,
                                     worst_csq_for_variant=hl.or_missing(hl.len(sorted_scores) > 0, sorted_scores[0]),
                                     worst_csq_by_gene_canonical=sorted_canonical_scores,
                                     worst_csq_for_variant_canonical=hl.or_missing(hl.len(sorted_canonical_scores) > 0, sorted_canonical_scores[0])
                                     )

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


def select_primitives_from_ht(ht: hl.Table) -> hl.Table:
    """
    Select only primitive types (string, int, float, bool) from a Table.
    Particularly useful for exporting a Table.

    :param Table ht: Input Table
    :return: Table with only primitive types selected
    :rtype: Table
    """
    return ht.select(**{x: v for x, v in ht.row_value.items() if
                        v.dtype in {hl.tstr, hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64, hl.tbool}})


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


def phase_trio_matrix_by_transmission(tm: hl.MatrixTable, call_field: str = 'GT', phased_call_field: str = 'PBT_GT') -> hl.MatrixTable:
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


def infer_families(relatedness_ht: hl.Table,
                   sex: Union[hl.Table, Dict[str, bool]],
                   duplicated_samples: Set[str],
                   i_col: str = 'i',
                   j_col: str = 'j',
                   kin_col: str = 'kin',
                   ibd2_col: str = 'ibd2',
                   first_degree_threshold: Tuple[float, float] = (0.2, 0.4),
                   second_degree_threshold: float = 0.05,
                   ibd2_parent_offspring_threshold: float = 0.2
                   ) -> hl.Pedigree:
    """

    Infers familial relationships from the results of pc_relate and sex information.
    Note that both kinship and ibd2 are needed in the pc_relate output.

    This function returns a pedigree containing trios inferred from the data. Family ID can be the same for multiple
    trios if one or more members of the trios are related (e.g. sibs, multi-generational family). Trios are ordered by family ID.

    Note that this function only returns complete trios defined as:
    one child, one father and one mother (sex is required for both parents)

    :param Table relatedness_ht: pc_relate output table
    :param Table or dict of str -> bool sex: Either a Table with the standard hail `is_female` annotation, or a dict containing the sex for each sample. True = female, False = male, None = unknown
    :param set of str duplicated_samples: Duplicated samples to remove (If not provided, this function won't work as it assumes that each child has exactly two parents)
    :param str i_col: Column containing the 1st sample id in the pc_relate table
    :param str j_col: Column containing the 2nd sample id in the pc_relate table
    :param str kin_col: Column containing the kinship in the pc_relate table
    :param str ibd2_col: Column containing ibd2 in the pc_relate table
    :param (float, float) first_degree_threshold: Lower/upper bounds for kin for 1st degree relatives
    :param float second_degree_threshold: Lower bound for kin for 2nd degree relatives
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
                if tuple(sorted((p1, p2))) not in relative_pairs:
                    if sex.get(p1) is False and sex.get(p2):
                        parents.append((p1, p2))
                    elif sex.get(p1) and sex.get(p2) is False:
                        parents.append((p2, p1))

        if len(parents) == 1:
            return parents[0]

        return None

    # Check relatedness table format
    if not relatedness_ht[i_col].dtype == relatedness_ht[j_col].dtype:
        logger.error("i_col and j_col of the relatedness table need to be of the same type.")

    # If `relatedness_ht` is keyed by a Struct, attempt to extract `s` instead
    if isinstance(relatedness_ht[i_col], hl.expr.StructExpression):
        if 's' in relatedness_ht[i_col]:
            logger.warn(f"Pedigrees can only be constructed from string IDs, but your relatedness_ht ID column is of type: {relatedness_ht[i_col].dtype}. Attempting to use relatedness_ht[{i_col}].s")
            relatedness_ht = relatedness_ht.key_by(  # Needed for familial inference at this point -- should be generalized
                **{
                    i_col: relatedness_ht[i_col].s,
                    j_col: relatedness_ht[j_col].s
                }
            )
        else:
            logger.error(
                f"Pedigrees can only be constructed from string IDs, but your relatedness_ht ID column is of type: {relatedness_ht[i_col].dtype}.")

    # If sex is a Table, extract sex information as a Dict
    if isinstance(sex, hl.Table):
        sex = dict(hl.tuple([sex.s, sex.is_female]).collect())

    # Get first degree relatives - exclude duplicate samples
    dups = hl.literal(duplicated_samples)
    first_degree_pairs = relatedness_ht.filter(
        (relatedness_ht[kin_col] >= first_degree_threshold[0]) &
        (relatedness_ht[kin_col] <= first_degree_threshold[1]) &  # We do not want to include twins/duplicate samples
        ~dups.contains(relatedness_ht[i_col]) &
        ~dups.contains(relatedness_ht[j_col])
    ).collect()
    first_degree_relatives = defaultdict(set)
    for row in first_degree_pairs:
        first_degree_relatives[row[i_col]].add(row[j_col])
        first_degree_relatives[row[j_col]].add(row[i_col])

    # Add second degree relatives for those samples
    # This is needed to distinguish grandparent - child - parent from child - mother, father down the line
    first_degree_samples = hl.literal(set(first_degree_relatives.keys()))
    second_degree_samples = relatedness_ht.filter(
        (first_degree_samples.contains(relatedness_ht[i_col]) | first_degree_samples.contains(relatedness_ht[j_col])) &
        (relatedness_ht[kin_col] >= second_degree_threshold) &
        (relatedness_ht[kin_col] <= first_degree_threshold[1])
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
        pop_pca_scores: Union[hl.Table, pd.DataFrame],
        pc_cols: Union[hl.expr.ArrayExpression, List[str]],
        known_col: str = 'known_pop',
        fit: Any = None, # Type should be RandomForestClassifier but we do not want to import sklearn.RandomForestClassifier outside
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth'
) -> Tuple[Union[hl.Table, pd.DataFrame], Any]: # 2nd element of the tuple should be RandomForestClassifier but we do not want to import sklearn.RandomForestClassifier outside
    """

    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.

    As input, this function can either take:
    - A Hail Table (typically the output of `hwe_normalized_pca`). In this case,
        - `pc_cols` should be an ArrayExpression of Floats where each element is one of the PCs to use.
        - A Hail Table will be returned as output
    - A Pandas DataFrame. In this case:
        - Each PC should be in a separate column and `pc_cols` is the list of all the columns containing the PCs to use.
        - A pandas DataFrame is returned as output

    Note
    ----
    If you have a Pandas Dataframe and have all PCs as an array in a single column, the
    `expand_pd_array_col` can be used to expand this column into multiple `PC` columns.

    :param Table or DataFrame pop_pc_pd: Input Hail Table or Pandas Dataframe
    :param ArrayExpression or list of str pc_cols: Columns storing the PCs to use
    :param str known_col: Column storing the known population labels
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Hail Table or Pandas Dataframe (depending on input) containing sample IDs and imputed population labels, trained random forest model
    :rtype: (Table or DataFrame, RandomForestClassifier)
    """
    from sklearn.ensemble import RandomForestClassifier
    hail_input = isinstance(pop_pca_scores, hl.Table)
    if hail_input:
        pop_pc_pd = pop_pca_scores.select(
            known_col,
            pca_scores=pc_cols
        ).to_pandas()
        pop_pc_pd = expand_pd_array_col(pop_pc_pd, 'pca_scores', out_cols_prefix='PC')
        pc_cols = [col for col in pop_pc_pd if col.startswith('PC')]
    else:
        pop_pc_pd = pop_pca_scores

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
        training_set_known_labels = train_fit[known_col].values
        training_set_pcs = train_fit[pc_cols].values
        evaluation_set_pcs = evaluate_fit[pc_cols].values

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
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].values)
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].values)
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])
    pop_pc_pd = pd.concat([pop_pc_pd, probs], axis=1)
    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label
    pop_pc_pd = pop_pc_pd.drop(pc_cols, axis='columns')

    logger.info("Found the following sample count after population assignment: {}".format(
        ", ".join(f'{pop}: {count}' for pop, count in Counter(pop_pc_pd[output_col]).items())
    ))

    if hail_input:
        pops_ht = hl.Table.from_pandas(pop_pc_pd, key=list(pop_pca_scores.key))
        pops_ht.annotate_globals(
            assign_pops_from_pc_params=hl.struct(
                min_assignment_prob=min_prob
            )
        )
        return pops_ht, pop_clf
    else:
        return pop_pc_pd, pop_clf


def merge_stats_counters_expr(stats: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """
    Merges multiple stats counters, assuming that they were computed on non-overlapping data.
    Examples:
    - Merge stats computed on indel and snv separately
    - Merge stats computed on bi-allelic and multi-allelic variants separately
    - Merge stats computed on autosomes and sex chromosomes separately

    :param hl.expr.ArrayExpression stats: An array of stats counters to merge
    :return: Merged stats Struct
    :rtype: hl.expr.StructExpression
    """

    def add_stats(i: hl.expr.StructExpression, j: hl.expr.StructExpression) -> hl.expr.StructExpression:
        """
        This merges two stast counters together. It assumes that all stats counter fields are present in the struct.

        :param hl.expr.Tuple i: accumulator: struct with mean, n and variance
        :param hl.expr.Tuple j: new element: stats_struct -- needs to contain mean, n and variance
        :return: Accumulation over all elements: struct with mean, n and variance
        :rtype: hl.expr.StructExpression
        """
        delta = j.mean - i.mean
        n_tot = i.n + j.n
        return hl.struct(
            min=hl.min(i.min, j.min),
            max=hl.max(i.max, j.max),
            mean=(i.mean * i.n + j.mean * j.n) / n_tot,
            variance=i.variance + j.variance + (delta * delta * i.n * j.n) / n_tot,
            n=n_tot,
            sum=i.sum + j.sum
        )

    # Gather all metrics present in all stats counters
    metrics = set(stats[0])
    dropped_metrics = set()
    for stat_expr in stats[1:]:
        stat_expr_metrics = set(stat_expr)
        dropped_metrics = dropped_metrics.union(stat_expr_metrics.difference(metrics))
        metrics = metrics.intersection(stat_expr_metrics)
    if dropped_metrics:
        logger.warn("The following metrics will be dropped during stats counter merging as they do not appear in all counters:".format(",".join(dropped_metrics)))

    # Because merging standard deviation requires having the mean and n,
    # check that they are also present if `stdev` is. Otherwise remove stdev
    if 'stdev' in metrics:
        missing_fields = [x for x in ['n', 'mean'] if x not in metrics]
        if missing_fields:
            logger.warn(f'Cannot merge `stdev` from given stats counters since they are missing the following fields: {",".join(missing_fields)}')
            metrics.remove('stdev')

    # Create a struct with all possible stats for merging.
    # This step helps when folding because we can rely on the struct schema
    # Note that for intermediate merging, we compute the variance rather than the stdev
    all_stats = hl.array(stats).map(lambda x: hl.struct(
        min=x.min if 'min' in metrics else hl.null(hl.tfloat64),
        max=x.max if 'max' in metrics else hl.null(hl.tfloat64),
        mean=x.mean if 'mean' in metrics else hl.null(hl.tfloat64),
        variance=x.stdev * x.stdev if 'stdev' in metrics else hl.null(hl.tfloat64),
        n=x.n if 'n' in metrics else hl.null(hl.tfloat64),
        sum=x.sum if 'sum' in metrics else hl.null(hl.tfloat64)
    ))

    # Merge the stats
    agg_stats = all_stats[1:].fold(add_stats, all_stats[0])

    # Return only the metrics that were present in all independent stats counters
    # If `stdev` is present, then compute it from the variance
    return agg_stats.select(
        **{metric: agg_stats[metric] if metric != 'stdev' else hl.sqrt(agg_stats.variance) for metric in metrics}
    )


def merge_sample_qc_expr(sample_qc_exprs: List[hl.expr.StructExpression]) -> hl.expr.StructExpression:
    """
    Creates an expression that merges results from non-overlapping strata of hail.sample_qc
    E.g.:
    - Compute autosomes and sex chromosomes metrics separately, then merge results
    - Compute bi-allelic and multi-allelic metrics separately, then merge results

    Note regarding the merging of ``dp_stats`` and ``gq_stats``:
    Because ``n`` is needed to aggregate ``stdev``, ``n_called`` is used for this purpose.
    This should work very well on a standard GATK VCF and it essentially assumes that:
    - samples that are called have `DP` and `GQ` fields
    - samples that are not called do not have `DP` and `GQ` fields
    Even if these assumptions are broken for some genotypes, it shouldn't matter too much.

    :param list of StructExpression sample_qc_exprs: List of sample QC struct expressions for each stratification
    :return: Combined sample QC results
    :rtype: StructExpression
    """

    # List of metrics that can be aggregated by summing
    additive_metrics = [
        'n_called', 'n_not_called', 'n_hom_ref', 'n_het', 'n_hom_var', 'n_snp',
        'n_insertion', 'n_deletion', 'n_singleton', 'n_transition', 'n_transversion', 'n_star'
    ]

    # List of metrics that are ratio of summed metrics (name, nominator, denominator)
    ratio_metrics = [
        ('call_rate', 'n_called', 'n_not_called'),
        ('r_ti_tv', 'n_transition', 'n_transversion'),
        ('r_het_hom_var', 'n_het', 'n_hom_var'),
        ('r_insertion_deletion', 'n_insertion', 'n_deletion')
    ]

    # List of metrics that are struct generated by a stats counter
    stats_metrics = ['gq_stats', 'dp_stats']

    # Gather metrics present in sample qc fields
    sample_qc_fields = set(sample_qc_exprs[0])
    for sample_qc_expr in sample_qc_exprs[1:]:
        sample_qc_fields = sample_qc_fields.union(set(sample_qc_expr))

    # Merge additive metrics in sample qc fields
    merged_exprs = {
        metric: hl.sum([sample_qc_expr[metric] for sample_qc_expr in sample_qc_exprs])
        for metric in additive_metrics if metric in sample_qc_fields
    }

    # Merge ratio metrics in sample qc fields
    merged_exprs.update({
        metric: hl.float64(divide_null(merged_exprs[nom], merged_exprs[denom]))
        for metric, nom, denom in ratio_metrics if nom in sample_qc_fields and denom in sample_qc_fields
    })

    # Merge stats counter metrics in sample qc fields
    # Use n_called as n for DP and GQ stats
    if 'n_called' in sample_qc_fields:
        merged_exprs.update({
            metric: merge_stats_counters_expr([
                sample_qc_expr[metric].annotate(n=sample_qc_expr.n_called)
                for sample_qc_expr in sample_qc_exprs
            ]).drop('n')
            for metric in stats_metrics
        })

    return hl.struct(**merged_exprs)


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression selecting bi-allelic sites only,
    accounting for whether the input MT/HT was split.

    :param Table or MatrixTable t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    :rtype: BooleanExpression
    """
    return (~t.was_split if 'was_split' in t.row else (hl.len(t.alleles) == 2))


def bi_allelic_site_inbreeding_expr(call: hl.expr.CallExpression) -> hl.expr.Float32Expression:
    """
    Return the site inbreeding coefficient as an expression to be computed on a MatrixTable.
    This is implemented based on the GATK InbreedingCoeff metric:
    https://software.broadinstitute.org/gatk/documentation/article.php?id=8032

    Note
    ----
    The computation is run based on the counts of alternate alleles and thus should only be run on bi-allelic sites.

    :param CallExpression call: Expression giving the calls in the MT
    :return: Site inbreeding coefficient expression
    :rtype: Float32Expression
    """

    def inbreeding_coeff(gt_counts: hl.expr.DictExpression) -> hl.expr.Float32Expression:
        n = gt_counts.get(0, 0) + gt_counts.get(1, 0) + gt_counts.get(2, 0)
        p = (2 * gt_counts.get(0, 0) + gt_counts.get(1, 0)) / (2 * n)
        q = (2 * gt_counts.get(2, 0) + gt_counts.get(1, 0)) / (2 * n)
        return 1 - (gt_counts.get(1, 0) / (2 * p * q * n))

    return hl.bind(
        inbreeding_coeff,
        hl.agg.counter(call.n_alt_alleles())
    )


def to_phred(linear_expr: hl.expr.NumericExpression) -> hl.expr.Float64Expression:
    """
    Computes the phred-scaled value of the linear-scale input

    :param NumericExpression linear_expr: input
    :return: Phred-scaled value
    :rtype: Float64Expression
    """
    return -10 * hl.log10(linear_expr)


def from_phred(phred_score_expr: hl.expr.NumericExpression) -> hl.expr.Float64Expression:
    """
    Computes the linear-scale value of the phred-scaled input.

    :param NumericExpression phred_score_expr: phred-scaled value
    :return: linear-scale value of the phred-scaled input.
    :rtype: Float64Expression
    """
    return 10 ** -(phred_score_expr / 10)


def fs_from_sb(
        sb: Union[hl.expr.ArrayNumericExpression, hl.expr.ArrayExpression],
        normalize: bool = True,
        min_cell_count: int = 200,
        min_count: int = 4,
        min_p_value: float = 1e-320
) -> hl.expr.Int64Expression:
    """
    Computes `FS` (Fisher strand balance) annotation from  the `SB` (strand balance table) field.
    `FS` is the phred-scaled value of the double-sided Fisher exact test on strand balance.

    Using default values will have the same behavior as the GATK implementation, that is:
    - If sum(counts) > 2*`min_cell_count` (default to GATK value of 200), they are normalized
    - If sum(counts) < `min_count` (default to GATK value of 4), returns missing
    - Any p-value < `min_p_value` (default to GATK value of 1e-320) is truncated to that value

    In addition to the default GATK behavior, setting `normalize` to `False` will perform a chi-squared test
    for large counts (> `min_cell_count`) instead of normalizing the cell values.

    Note
    ----
    This function can either take
    - an array of length containing the table counts: [ref fwd, ref rev, alt fwd, alt rev]
    - an array containig 2 arrays of length 2, containing the counts: [[ref fwd, ref rev], [alt fwd, alt rev]]

    GATK code here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/FisherStrand.java

    :param ArrayNumericExpression or ArrayExpression sb: Count of ref/alt reads on each strand
    :param bool normalize: Whether to normalize counts is sum(counts) > min_cell_count (`normalize`=True), or use a chi sq instead of FET (`normalize`=False)
    :param int min_cell_count: Maximum count for performing a FET
    :param int min_count: Minimum total count to output `FS` (otherwise `null` it output)
    :return: FS value
    :rtype: Int64Expression
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(
            lambda x: hl.flatten(x),
            sb
        )

    sb_sum = hl.bind(
        lambda x: hl.sum(x),
        sb
    )

    # Normalize table if counts get too large
    if normalize:
        fs_expr = hl.bind(
            lambda sb, sb_sum: hl.cond(
                sb_sum <= 2 * min_cell_count,
                sb,
                sb.map(lambda x: hl.int(x / (sb_sum / min_cell_count)))
            ),
            sb, sb_sum
        )

        # FET
        fs_expr = to_phred(
            hl.max(
                hl.fisher_exact_test(
                    fs_expr[0], fs_expr[1], fs_expr[2], fs_expr[3]
                ).p_value,
                min_p_value
            )
        )
    else:
        fs_expr = to_phred(
            hl.max(
                hl.contingency_table_test(
                    sb[0], sb[1], sb[2], sb[3], min_cell_count=min_cell_count
                ).p_value,
                min_p_value
            )
        )

    # Return null if counts <= `min_count`
    return hl.or_missing(
        sb_sum > min_count,
        hl.max(0, fs_expr) # Needed to avoid -0.0 values
    )


def vep_struct_to_csq(
        vep_expr: hl.expr.StructExpression,
        csq_fields: str = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
                          "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
                          "ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|"
                          "TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|"
                          "HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|"
                          "LoF_flags|LoF_info"
) -> hl.expr.ArrayExpression:
    """
    Given a VEP Struct, returns and array of VEP VCF CSQ strings (one per consequence in the struct).
    The fields and their order will correspond to those passed in `csq_fields`, which corresponds to the
    VCF header that is required to interpret the VCF CSQ INFO field.

    Note that the order is flexible and that all fields that are in the default value are supported.
    These fields will be formatted in the same way that their VEP CSQ counterparts are.

    While other fields can be added if their name are the same as those in the struct. Their value will be the result of calling
    hl.str(), so it may differ from their usual VEP CSQ representation.

    :param vep_expr: The input VEP Struct
    :param csq_fields: The | delimited list of fields to include in the CSQ (in that order)
    :return: The corresponding CSQ strings
    :rtype: ArrayExpression
    """

    _csq_fields = [f.lower() for f in csq_fields.split("|")]

    def get_csq_from_struct(element: hl.expr.StructExpression, feature_type: str) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update({
            'allele': element.variant_allele,
            'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
            'feature_type': feature_type,
            'feature': (
                element.transcript_id if 'transcript_id' in element else
                element.regulatory_feature_id if 'regulatory_feature_id' in element else
                element.motif_feature_id if 'motif_feature_id' in element else ''
            ),
            'variant_class': vep_expr.variant_class
        })

        # Add exception for transcripts
        if feature_type == 'Transcript':
            fields.update({
                'canonical': hl.cond(element.canonical == 1, 'YES', ''),
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
                'symbol_source': element.gene_symbol_source,
                'cdna_position': hl.str(element.cdna_start) + hl.cond(element.cdna_start == element.cdna_end, '', "-" + hl.str(element.cdna_end)),
                'cds_position': hl.str(element.cds_start) + hl.cond(element.cds_start == element.cds_end, '', "-" + hl.str(element.cds_end)),
                'protein_position': hl.str(element.protein_start) + hl.cond(element.protein_start == element.protein_end, '', "-" + hl.str(element.protein_end)),
                'sift': element.sift_prediction + "(" + hl.format('%.3f', element.sift_score) + ")",
                'polyphen': element.polyphen_prediction + "(" + hl.format('%.3f', element.polyphen_score) + ")",
                'domains': hl.delimit(element.domains.map(lambda d: d.db + ":" + d.name), "&")
            })
        elif feature_type == 'MotifFeature':
            fields['motif_score_change'] = hl.format('%.3f', element.motif_score_change)

        return hl.delimit([hl.or_else(hl.str(fields.get(f, '')), '') for f in _csq_fields], "|")

    csq = hl.empty_array(hl.tstr)
    for feature_field, feature_type in [
        ('transcript_consequences', 'Transcript'),
        ('regulatory_feature_consequences', 'RegulatoryFeature'),
        ('motif_feature_consequences', 'MotifFeature'),
        ('intergenic_consequences', 'Intergenic')]:
        csq = csq.extend(
            hl.or_else(
                vep_expr[feature_field].map(
                    lambda x: get_csq_from_struct(x, feature_type=feature_type)
                ),
                hl.empty_array(hl.tstr)
            )
        )

    return hl.or_missing(hl.len(csq) > 0, csq)


def get_array_element_type(array_expr: hl.expr.ArrayExpression) -> hl.HailType:
    """
    Returns the type of an array element.

    :param ArrayExpression array_expr: The array expression to get the element type
    :return: Hail type
    :rtype: HailType
    """
    return array_expr.dtype.element_type


def ht_to_vcf_mt(
        info_ht: hl.Table,
        pipe_delimited_annotations  : List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS
 ) -> hl.MatrixTable:
    """
    Creates a MT ready for vcf export from a HT. In particular, the following conversions are done:
    - All int64 are coerced to int32
    - Fields specified by `pipe_delimited_annotations` will be converted from arrays to pipe-delimited strings

    Note: The MT returned has no cols.

    :param Table info_ht: Input HT
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct)
    :return: MatrixTable ready for VCF export
    :rtype: MatrixTable
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    info_ht = info_ht.key_by('locus', 'alleles')

    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in info_ht.info.dtype.items():
        if ft == hl.dtype('int64'):
            logger.warn(f"Coercing field info.{f} from int64 to int32 for VCF output. Value will be capped at int32 max value.")
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(**{f: hl.int32(hl.min(2**31 - 1, info_ht.info[f]))})
            )
        elif ft == hl.dtype('array<int64>'):
            logger.warn(f"Coercing field info.{f} from array<int64> to array<int32> for VCF output. Array values will be capped at int32 max value.")
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(**{f: info_ht.info[f].map(lambda x: hl.int32(hl.min(2**31 - 1, x)))})
            )

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, but it is overwritten below
    for f in pipe_delimited_annotations:
        if f in info_ht.info:
            info_expr[f] = "|" + get_pipe_expr(info_ht.info[f])

    # Take care of strand balance field
    if 'SB' in info_ht.info:
        info_expr['SB'] = info_ht.info.SB[0].extend(info_ht.info.SB[1])

    if 'AS_SB_TABLE' in info_ht.info:
        info_expr['AS_SB_TABLE'] = get_pipe_expr(info_ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ",")))

    # Annotate with new expression and add 's' empty string field required to cast HT to MT
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(**info_expr),
        s=hl.null(hl.tstr)
    )

    # Create an MT with no cols so that we acn export to VCF
    info_mt = info_ht.to_matrix_table_row_major(columns=['s'], entry_field_name='s')
    return info_mt.filter_cols(False)
