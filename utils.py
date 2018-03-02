
import re
import sys
import logging
import gzip
import os

import hail as hl
from hail.expr import Field
from hail.expr.expression import *
from collections import defaultdict, namedtuple, OrderedDict
from pprint import pprint, pformat
import argparse
from typing import Dict, Tuple, List

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)


POP_NAMES = {'AFR': "African/African American",
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
CSQ_CODING_HIGH_IMPACT = ["transcript_ablation",
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


def filter_to_adj(vds):
    """
    Filter genotypes to adj criteria

    :param MatrixTable vds: VDS
    :return: MT
    :rtype: MatrixTable
    """
    try:
        vds = vds.filter_entries(vds.adj)
    except AttributeError:
        vds = annotate_adj(vds)
        vds = vds.filter_entries(vds.adj)
    return vds.drop(vds.adj)


def annotate_adj(vds):
    """
    Annotate genotypes with adj criteria (assumes diploid)

    :param MatrixTable vds: MT
    :return: MT
    :rtype: MatrixTable
    """
    adj_gq = 20
    adj_dp = 10
    adj_ab = 0.2

    return vds.annotate_entries(adj=
                                (vds.GQ >= adj_gq) & (vds.DP >= adj_dp) & (
                                    ~vds.GT.is_het() |
                                    ((vds.GT[0] == 0) & (vds.AD[vds.GT[1]] / vds.DP >= adj_ab)) |
                                    ((vds.GT[0] > 0) & (vds.AD[vds.GT[0]] / vds.DP >= adj_ab) &
                                     (vds.AD[vds.GT[1]] / vds.DP >= adj_ab))
                                )
    )


def add_variant_type(alt_alleles):
    """
    Get Struct of variant_type and n_alt_alleles from ArrayExpression of AltAlleles

    :param ArrayExpression alt_alleles: Input ArrayExpression of Strings
    :return: Struct with variant_type and n_alt_alleles
    :rtype: Struct
    """
    ref = alt_alleles[0]
    alts = alt_alleles[1:]
    non_star_alleles = hl.filter(lambda a: a != '*', alts)
    return Struct(variant_type=
                  hl.cond(
                      hl.all(lambda a: hl.is_snp(ref, a), non_star_alleles),
                      hl.cond(hl.len(non_star_alleles) > 1, "multi-snv", "snv"),
                      hl.cond(
                          hl.all(lambda a: hl.is_indel(ref, a), non_star_alleles),
                          hl.cond(hl.len(non_star_alleles) > 1, "multi-indel", "indel"),
                          "mixed")
                  ),
                  n_alt_alleles=hl.len(non_star_alleles))


def split_multi_dynamic(vds, keep_star=False):
    """
    Splits MatrixTable based on entry fields found. Downcodes whatever it can. Supported so far:
    GT, DP, AD, PL, GQ
    PGT, PID
    ADALL

    :param MatrixTable vds: Input MatrixTable
    :return: Split MatrixTable
    :rtype: MatrixTable
    """
    fields = set(map(lambda x: x.name, vds.entry_schema.fields))
    sm = hl.SplitMulti(vds, keep_star=keep_star)
    sm.update_rows(a_index=sm.a_index(), was_split=sm.was_split())
    expression = {}

    # HTS/standard
    if 'GT' in fields:
        expression['GT'] = hl.downcode(vds.GT, sm.a_index())
    if 'DP' in fields:
        expression['DP'] = vds.DP
    if 'AD' in fields:
        expression['AD'] = hl.or_missing(hl.is_defined(vds.AD),
                                         [hl.sum(vds.AD) - vds.AD[sm.a_index()], vds.AD[sm.a_index()]])
    if 'PL' in fields:
        pl = hl.or_missing(
            hl.is_defined(vds.PL),
            (hl.range(0, 3).map(lambda i:
                                hl.min((hl.range(0, hl.triangle(vds.alleles.length()))
                                        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j),
                                                                      sm.a_index()) == hl.unphased_diploid_gt_index_call(i)
                                ).map(lambda j: vds.PL[j]))))))
        expression['PL'] = pl
        if 'GQ' in fields:
            expression['GQ'] = hl.gq_from_pl(pl)
    else:
        if 'GQ' in fields:
            expression['GQ'] = vds.GQ

    # Phased data
    if 'PGT' in fields:
        expression['PGT'] = hl.downcode(vds.PGT, sm.a_index())
    if 'PID' in fields:
        expression['PGT'] = vds.PID

    # Custom data
    if 'ADALL' in fields:  # found in NA12878
        expression['ADALL'] = hl.or_missing(hl.is_defined(vds.ADALL),
                                            [hl.sum(vds.ADALL) - vds.ADALL[sm.a_index()], vds.ADALL[sm.a_index()]])

    sm.update_entries(**expression)
    return sm.result()


def adjust_sex_ploidy(vds, sex_expr):
    """
    Converts males to haploid on non-PAR X/Y, sets females to missing on Y

    :param MatrixTable vds: VDS
    :param StringExpression sex_expr: Expression pointing to sex in VDS (must be "male" and "female", otherwise no change)
    :return: MatrixTable with fixed ploidy for sex chromosomes
    :rtype: MatrixTable
    """
    male = sex_expr == 'male'
    female = sex_expr == 'female'
    x_nonpar = vds.locus.in_x_nonpar()
    y_par = vds.locus.in_y_par()
    y_nonpar = vds.locus.in_y_nonpar()
    return vds.annotate_entries(
        GT=hl.case()
        .when(female & (y_par | y_nonpar), hl.null(hl.TCall()))
        .when(male & (x_nonpar | y_nonpar) & vds.GT.is_het(), hl.null(hl.TCall()))
        .when(male & (x_nonpar | y_nonpar), hl.call(vds.GT[0], phased=False))
        .default(vds.GT)
    )


def get_sample_data(vds, fields, sep='\t', delim='|'):
    """
    Hail devs hate this one simple py4j trick to speed up sample queries

    :param MatrixTable or Table vds: MT
    :param list of StringExpression fields: fields
    :param sep: Separator to use (tab usually fine)
    :param delim: Delimiter to use (pipe usually fine)
    :return: Sample data
    :rtype: list of list of str
    """
    field_expr = fields[0]
    for field in fields[1:]:
        field_expr = field_expr + '|' + field
    if isinstance(vds, hl.MatrixTable):
        vds_agg = vds.aggregate_cols
    else:
        vds_agg = vds.aggregate
    return [x.split(delim) for x in vds_agg(hl.delimit(hl.agg.collect(field_expr), sep)).split(sep) if x != 'null']


def add_popmax_expr(freq):
    """
    First pass of popmax (add an additional entry into freq with popmax: pop)
    TODO: update dict instead?

    :param ArrayStructExpression freq: Array of StructExpression with ['ac', 'an', 'hom', 'meta']
    :return: Frequency data with annotated popmax
    :rtype: ArrayStructExpression
    """
    freq_filtered = hl.filter(lambda x: (x.meta.keys() == ['population']) & (x.meta['population'] != 'oth'), freq)
    popmax_entry = hl.sorted(freq_filtered, key=lambda x: x.ac / x.an, reverse=True)[0]  # TODO: check for missing
    return freq.append(Struct(ac=popmax_entry.ac, an=popmax_entry.an, hom=popmax_entry.hom,
                              meta={'popmax': popmax_entry.meta['population']}))


def get_projectmax(vds, loc):
    """
    First pass of projectmax (returns aggregated VDS with project_max field)

    :param MatrixTable vds: Input VDS
    :param StringExpression loc: Column expression location of project ID (e.g. vds.meta.pid)
    :return: Frequency data with annotated project_max
    :rtype: MatrixTable
    """
    # TODO: add hom count
    vds = vds.annotate_cols(project=loc)
    agg_vds = vds.group_cols_by(vds.project).aggregate(ac=hl.agg.sum(vds.GT.num_alt_alleles()),
                                                       an=2 * hl.agg.count_where(hl.is_defined(vds.GT)))
    agg_vds = agg_vds.annotate_entries(af=agg_vds.ac / agg_vds.an)
    return agg_vds.annotate_rows(project_max=hl.agg.take(Struct(project=agg_vds.project, ac=agg_vds.ac,
                                                                af=agg_vds.af, an=agg_vds.an), 5, -agg_vds.af))


def flatten_struct(struct, root='va', leaf_only=True, recursive=True):
    """
    Given a `TStruct` and its `root` path, creates an `OrderedDict` of each path -> Field by flattening the `TStruct` tree.
    The order of the fields is the same as the input `Struct` fields, using a depth-first approach.
    When `leaf_only=False`, `Struct`s roots are printed as they are traversed (i.e. before their leaves).
    The following TStruct at root 'va', for example
    Struct{
     rsid: String,
     qual: Double,
     filters: Set[String],
     info: Struct{
         AC: Array[Int],
         AF: Array[Double],
         AN: Int
         }
    }

    Would give the following dict:
    {
        'va.rsid': Field(rsid),
        'va.qual': Field(qual),
        'va.filters': Field(filters),
        'va.info.AC': Field(AC),
        'va.info.AF': Field(AF),
        'va.info.AN': Field(AN)
    }

    :param TStruct struct: The struct to flatten
    :param str root: The root path of the struct to flatten (added at the beginning of all dict keys)
    :param bool leaf_only: When set to `True`, only leaf nodes in the tree are output in the output
    :param bool recursive: When set to `True`, internal `Struct`s are flatten
    :return: Dictionary of path : Field
    :rtype: OrderedDict of str:Field
    """
    result = OrderedDict()
    for f in struct.fields:
        path = '{}.{}'.format(root, f.name) if root else f.name
        if isinstance(f.typ, TStruct) and recursive:
            if not leaf_only:
                result[path] = f
            result.update(flatten_struct(f.typ, path, leaf_only))
        else:
            result[path] = f
    return result


def ann_exists(annotation, schema, root='va'):
    """
    Tests whether an annotation (given by its full path) exists in a given schema and its root.

    :param str annotation: The annotation to find (given by its full path in the schema tree)
    :param TStruct schema: The schema to find the annotation in
    :param str root: The root of the schema (or struct)
    :return: Whether the annotation was found
    :rtype: bool
    """
    anns = flatten_struct(schema, root, leaf_only=False)
    return annotation in anns


def get_ann_field(annotation, schema, root='va'):
    """
    Given an annotation path and a schema, return that annotation field.

    :param str annotation: annotation path to fetch
    :param TStruct schema: schema (or struct) in which to search
    :param str root: root of the schema (or struct)
    :return: The Field corresponding to the input annotation
    :rtype: Field
    """
    anns = flatten_struct(schema, root, leaf_only=False)
    if not annotation in anns:
        logger.error("%s missing from schema.", annotation)
        sys.exit(1)
    return anns[annotation]


def get_ann_type(annotation, schema, root='va'):
    """
     Given an annotation path and a schema, return the type of the annotation.

    :param str annotation: annotation path to fetch
    :param TStruct schema: schema (or struct) in which to search
    :param str root: root of the schema (or struct)
    :return: The type of the input annotation
    :rtype: Type
    """
    return get_ann_field(annotation, schema, root).typ


def annotation_type_is_numeric(t):
    """
    Given an annotation type, returns whether it is a numerical type or not.

    :param Type t: Type to test
    :return: If the input type is numeric
    :rtype: bool
    """
    return (isinstance(t, TInt32) or
            isinstance(t, TInt64) or
            isinstance(t, TFloat32) or
            isinstance(t, TFloat64)
            )


def annotation_type_in_vcf_info(t):
    """
    Given an annotation type, returns whether that type can be natively exported to a VCF INFO field.
    Note types that aren't natively exportable to VCF will be converted to String on export.

    :param Type t: Type to test
    :return: If the input type can be exported to VCF
    :rtype: bool
    """
    return (annotation_type_is_numeric(t) or
            isinstance(t, TString) or
            isinstance(t, TArray) or
            isinstance(t, TSet) or
            isinstance(t, TBoolean)
            )


def run_samples_sanity_checks(vds, reference_vds, n_samples=10, verbose=True):
    """
    logger.info("Running samples sanity checks on %d samples" % n_samples)

    comparison_metrics = ['nHomVar',
                          'nSNP',
                          'nTransition',
                          'nTransversion',
                          'nInsertion',
                          'nDeletion',
                          'nNonRef',
                          'nHet'
                          ]

    samples = vds.sample_ids[:n_samples]

    def get_samples_metrics(vds, samples):
        metrics = (vds.filter_samples_expr('["%s"].toSet.contains(s)' % '","'.join(samples))
                   .sample_qc()
                   .query_samples('samples.map(s => {sample: s, metrics: sa.qc }).collect()')
                   )
        return {x.sample: x.metrics for x in metrics}

    test_metrics = get_samples_metrics(vds, samples)
    ref_metrics = get_samples_metrics(reference_vds, samples)

    output = ''

    for s, m in test_metrics.iteritems():
        if s not in ref_metrics:
            output += "WARN: Sample %s not found in reference data.\n" % s
        else:
            rm = ref_metrics[s]
            for metric in comparison_metrics:
                if m[metric] == rm[metric]:
                    if verbose:
                        output += "SUCCESS: Sample %s %s matches (N = %d).\n" % (s, metric, m[metric])
                else:
                    output += "FAILURE: Sample %s, %s differs: Data: %s, Reference: %s.\n" % (
                        s, metric, m[metric], rm[metric])

    logger.info(output)
    return output
    """
    raise NotImplementedError


def filter_annotations_regex(annotation_fields, ignore_list):
    def ann_in(name, lst):
        # `list` is a list of regexes to ignore
        return any([x for x in lst if re.search('^%s$' % x, name)])

    return [x for x in annotation_fields if not ann_in(x.name, ignore_list)]


def pc_project(vds, pc_vds, pca_loadings_root='va.pca_loadings'):
    """
    Projects samples in `vds` on PCs computed in `pc_vds`
    :param vds: VDS containing the samples to project
    :param pc_vds: VDS containing the PC loadings for the variants
    :param pca_loadings_root: Annotation root for the loadings. Can be either an Array[Double] or a Struct{ PC1: Double, PC2: Double, ...}
    :return: VDS with

    pc_vds = pc_vds.annotate_variants_expr('va.pca.calldata = gs.callStats(g => v)')

    pcs_struct_to_array = ",".join(['vds.pca_loadings.PC%d' % x for x in range(1, 21)])
    arr_to_struct_expr = ",".join(['PC%d: sa.pca[%d - 1]' % (x, x) for x in range(1, 21)])

    vds = (vds.filter_multi()
           .annotate_variants_vds(pc_vds, expr = 'va.pca_loadings = [%s], va.pca_af = vds.pca.calldata.AF[1]' % pcs_struct_to_array)
           .filter_variants_expr('!isMissing(va.pca_loadings) && !isMissing(va.pca_af)')
     )

    n_variants = vds.query_variants(['variants.count()'])[0]

    return(vds
           .annotate_samples_expr('sa.pca = gs.filter(g => g.isCalled && va.pca_af > 0.0 && va.pca_af < 1.0).map(g => let p = va.pca_af in (g.gt - 2 * p) / sqrt(%d * 2 * p * (1 - p)) * va.pca_loadings).sum()' % n_variants)
           .annotate_samples_expr('sa.pca = {%s}' % arr_to_struct_expr)
    )
    """
    raise NotImplementedError


def read_list_data(input_file):
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


def rename_samples(vds, input_file, filter_to_samples_in_file=False):
    """
    names = {old: new for old, new in [x.split("\t") for x in read_list_data(input_file)]}
    logger.info("Found %d samples for renaming in input file %s." % (len(names.keys()), input_file))
    logger.info("Renaming %d samples found in VDS" % len(set(names.keys()).intersection(set(vds.sample_ids)) ))

    if filter_to_samples_in_file:
        vds = vds.filter_samples_list(names.keys())
    return vds.rename_samples(names)
    """
    raise NotImplementedError


def filter_low_conf_regions(vds, filter_lcr=True, filter_decoy=True, filter_segdup=True, high_conf_regions=None):
    """
    Filters low-confidence regions

    :param MatrixTable vds: VDS to filter
    :param bool filter_lcr: Whether to filter LCR regions
    :param bool filter_decoy: Whether to filter decoy regions
    :param bool filter_segdup: Whether to filter Segdup regions
    :param list of str high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: VDS with low confidence regions removed
    :rtype: MatrixTable
    """
    from gnomad_hail.resources import lcr_intervals_path, decoy_intervals_path, segdup_intervals_path

    if filter_lcr:
        lcr = hl.import_interval_list(lcr_intervals_path)
        vds = vds.filter_rows(lcr[vds.locus], keep=False)

    if filter_decoy:
        decoy = hl.import_interval_list(decoy_intervals_path)
        vds = vds.filter_rows(decoy[vds.locus], keep=False)

    if filter_segdup:
        segdup = hl.import_interval_list(segdup_intervals_path)
        vds = vds.filter_rows(segdup[vds.locus], keep=False)

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_interval_list(region)
            vds = vds.filter_rows(region[vds.locus], keep=True)

    return vds


def process_consequences(vds, vep_root='vep'):
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq_by_gene, canonical_csq_by_gene, any_lof into [vep_root]

    :param MatrixTable vds: Input VDS
    :param str vep_root: Root for vep annotation (probably vep)
    :return: VDS with better formatted consequences
    :rtype: MatrixTable
    """
    csqs = hl.literal(CSQ_ORDER)
    csq_dict = hl.literal(dict(zip(CSQ_ORDER, range(len(CSQ_ORDER)))))

    def add_most_severe_consequence(tc):
        """
        Add most_severe_consequence annotation to transcript consequences
        This is for a given transcript, as there are often multiple annotations for a single transcript:
        e.g. splice_region_variant&intron_variant -> splice_region_variant

        :param StructExpression tc: Transcript consequences expression
        :return: Transcript consequences expression with most_severe_consequence
        :rtype StructExpression
        """
        return tc.annotate(
            most_severe_consequence=csqs.find(lambda c: tc.consequence_terms.contains(c))
        )

    def find_worst_transcript_consequence(tcl):
        """
        Gets worst transcript_consequence from an array of em

        :param ArrayStructExpression tcl: Array of Structs, one for each consequence
        :return: Worst transcript consequence among an array
        :rtype: StructExpression
        """
        if tcl.length() == 0: return tcl
        csq_score = lambda tc: csq_dict[csqs.find(tc)]
        tcl = tcl.map(lambda tc: tc.annotate(csq_score=hl.case()
            .when((tc.lof == 'HC') & (tc.lof_flags == ''), csq_score(tc) - 1000)
            .when((tc.lof == 'HC') & (tc.lof_flags != ''), csq_score(tc) - 500)
            .when(tc.lof == 'LC', csq_score(tc) - 10)
            .when(tc.polyphen_prediction == 'probably_damaging', csq_score(tc) - 0.5)
            .when(tc.polyphen_prediction == 'possibly_damaging', csq_score(tc) - 0.25)
            .when(tc.polyphen_prediction == 'benign', csq_score(tc) - 0.1)
            .default(csq_score(tc))
        ))
        return tcl.sort_by(lambda x: x.csq_score)[0]

    transcript_csqs = vds[vep_root].transcript_consequences.map(add_most_severe_consequence)

    gene_dict = transcript_csqs.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene = gene_dict.map_values(find_worst_transcript_consequence)
    canonical_csq_gene = gene_dict.map_values(lambda tcl: tcl.filter(lambda tc: tc.canonical == 1)[0])
    worst_csq = csqs.find(lambda c: transcript_csqs.map(lambda tc: tc.most_severe_consequence).contains(c))

    vep_data = vds[vep_root].annotate(transcript_consequences=transcript_csqs,
                                      worst_csq_by_gene=worst_csq_gene,
                                      canonical_csq_by_gene=canonical_csq_gene,
                                      any_lof=worst_csq_gene.values().exists(lambda x: x.lof == 'HC'),
                                      worst_csq_overall=worst_csq)

    return vds.annotate_rows(**{vep_root: vep_data})


def filter_vep_to_canonical_transcripts(vds, vep_root='vep'):
    """

    :param MatrixTable vds: VDS
    :param str vep_root: Location of VEP data
    :return: MT
    :rtype: MatrixTable
    """
    canonical = vds[vep_root].transcript_consequences.filter(lambda csq: csq.canonical == 1)
    vep_data = vds[vep_root].annotate(transcript_consequences=canonical)
    return vds.annotate_rows(**{vep_root: vep_data})


def filter_vep_to_synonymous_variants(vds, vep_root='vep'):
    """

    :param MatrixTable vds: Input VDS
    :param str vep_root: Location of VEP data
    :return: MT
    :rtype: MatrixTable
    """
    synonymous = vds[vep_root].transcript_consequences.filter(lambda csq: csq.most_severe_consequence == "synonymous_variant")
    vep_data = vds[vep_root].annotate(transcript_consequences=synonymous)
    return vds.annotate_rows(**{vep_root: vep_data})


def toSSQL(s):
    """
        Replaces `.` with `___`, since Spark ML doesn't support column names with `.`

    :param str s: The string in which the replacement should be done
    :return: string with `___`
    :rtype: str
    """
    return s.replace('.', '___')


def fromSSQL(s):
    """
        Replaces `___` with `.`, to go back from SSQL to hail annotations

    :param str s: The string in which the replacement should be done
    :return: string with `.`
    :rtype: str
    """
    return s.replace('___', '.')


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
