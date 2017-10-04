
import re
import sys
import json
import copy
import logging
import gzip
import os
import random

from resources import *
from hail import *
from slack_utils import *
from pyspark.sql.functions import bround
from pprint import pprint, pformat

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

ADJ_GQ = 20
ADJ_DP = 10
ADJ_AB = 0.2

ADJ_CRITERIA = 'g.gq >= %(gq)s && g.dp >= %(dp)s && (' \
               '!g.isHet || ' \
               '(g.gtj == 0 && g.ad[g.gtk]/g.dp >= %(ab)s) || ' \
               '(g.gtj > 0 && g.ad[g.gtj]/g.dp >= %(ab)s && g.ad[g.gtk]/g.dp >= %(ab)s)' \
               ')' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}

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


def cut_allele_from_g_array(target, destination=None):
    if destination is None: destination = target
    return ('%s = let removed_alleles = range(1, v.nAltAlleles + 1).filter(i => !aIndices.toSet.contains(i)).toSet in\n'
            'range(%s.size).filter(i => !removed_alleles.contains(gtj(i)) && !removed_alleles.contains(gtk(i)))\n'
            '.map(i => %s[i])' % (destination, target, target))


def index_into_arrays(a_based_annotations=None, r_based_annotations=None, vep_root=None, drop_ref_ann = False):
    """

    Creates annotation expressions to get the correct values when splitting multi-allelics

    :param list of str a_based_annotations: A-based annotations
    :param list of str r_based_annotations: R-based annotations
    :param str vep_root: Root of the vep annotation
    :param bool drop_ref_ann: If set to True, then the reference value of R-based annotations is removed (effectively converting them in A-based annotations)
    :return: Annotation expressions
    :rtype: list of str
    """
    annotations = []
    if a_based_annotations:
        for ann in a_based_annotations:
            annotations.append('{0} = {0}[va.aIndex - 1]'.format(ann))
    if r_based_annotations:
        expr = '{0} = {0}[va.aIndex]' if drop_ref_ann else '{0} = [{0}[0], {0}[va.aIndex]]'
        for ann in r_based_annotations:
            annotations.append(expr.format(ann))
    if vep_root:
        sub_fields = ['transcript_consequences', 'intergenic_consequences', 'motif_feature_consequences', 'regulatory_feature_consequences']
        annotations.extend(['{0}.{1} = {0}.{1}.filter(x => x.allele_num == va.aIndex)'.format(vep_root, sub_field) for sub_field in sub_fields])

    return annotations


def unfurl_filter_alleles_annotation(a_based=None, r_based=None, g_based=None, additional_annotations=None):

    annotations = []
    if r_based:
        for ann in r_based:
            annotations.append('%s = aIndices.map(i => %s[i])' % (ann, ann))

    if a_based:
        for ann in a_based:
            annotations.append('%s = aIndices[1:].map(i => %s[i - 1])' % (ann, ann))

    if g_based:
        for ann in g_based:
            annotations.append(cut_allele_from_g_array(ann))

    if additional_annotations:
        if isinstance(additional_annotations, str):
            annotations.append(additional_annotations)
        else:
            annotations.extend(additional_annotations)

    return ',\n'.join(annotations)


def filter_to_adj(vds):
    return vds.filter_genotypes(ADJ_CRITERIA)


def filter_star(vds, a_based=None, r_based=None, g_based=None, additional_annotations=None):
    annotation = unfurl_filter_alleles_annotation(a_based=a_based, r_based=r_based, g_based=g_based,
                                                  additional_annotations=additional_annotations)
    return vds.filter_alleles('v.altAlleles[aIndex - 1].alt == "*"', annotation=annotation, keep=False)


def flatten_struct(struct, root='', leaf_only=True):
    result = {}
    for f in struct.fields:
        path = '%s.%s' % (root, f.name)
        if isinstance(f.typ, TStruct):
            result.update(flatten_struct(f.typ, path))
            if not leaf_only:
                result[path] = f
        else:
            result[path] = f
    return result


def ann_exists(annotation, schema, root='va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    return annotation in anns


def get_ann_field(annotation, schema, root='va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    if not annotation in anns:
        logger.error("%s missing from schema.", annotation)
        sys.exit(1)
    return anns[annotation]


def get_ann_type(annotation, schema, root='va'):
    return get_ann_field(annotation, schema, root).typ


def annotation_type_is_numeric(t):
    return (isinstance(t, TInt) or
            isinstance(t, TLong) or
            isinstance(t, TFloat) or
            isinstance(t, TDouble)
            )


def get_variant_type_expr(root="va.variantType"):
    return '''%s =
    let non_star = v.altAlleles.filter(a => a.alt != "*") in
        if (non_star.forall(a => a.isSNP))
            if (non_star.length > 1)
                "multi-snv"
            else
                "snv"
        else if (non_star.forall(a => a.isIndel))
            if (non_star.length > 1)
                "multi-indel"
            else
                "indel"
        else
            "mixed"''' % root


def get_allele_stats_expr(root="va.stats", medians=False, samples_filter_expr=''):
    """

    Gets allele-specific stats expression: GQ, DP, NRQ, AB, Best AB, p(AB), NRDP, QUAL, combined p(AB)

    :param str root: annotations root
    :param bool medians: Calculate medians for GQ, DP, NRQ, AB and p(AB)
    :param str samples_filter_expr: Expression for filtering samples (e.g. "sa.keep")
    :return: List of expressions for `annotate_alleles_expr`
    :rtype: list of str
    """

    if samples_filter_expr:
        samples_filter_expr = "&& " + samples_filter_expr

    stats = ['%s.gq = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).stats()',
             '%s.dp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).stats()',
             '%s.nrq = gs.filter(g => g.isCalledNonRef %s).map(g => -log10(g.gp[0])).stats()',
             '%s.ab = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).stats()',
             '%s.best_ab = gs.filter(g => g.isHet %s).map(g => abs((g.ad[1]/g.dp) - 0.5)).min()',
             '%s.pab = gs.filter(g => g.isHet %s).map(g => g.pAB()).stats()',
             '%s.nrdp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).sum()',
             '%s.qual = -10*gs.filter(g => g.isCalledNonRef %s).map(g => if(g.pl[0] > 3000) -300 else log10(g.gp[0])).sum()',
             '%s.combined_pAB = let hetSamples = gs.filter(g => g.isHet %s).map(g => log(g.pAB())).collect() in orMissing(!hetSamples.isEmpty, -10*log10(pchisqtail(-2*hetSamples.sum(),2*hetSamples.length)))']

    if medians:
        stats.extend(['%s.gq_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).collect().median',
                    '%s.dp_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).collect().median',
                    '%s.nrq_median = gs.filter(g => g.isCalledNonRef %s).map(g => -log10(g.gp[0])).collect().median',
                    '%s.ab_median = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).collect().median',
                    '%s.pab_median = gs.filter(g => g.isHet %s).map(g => g.pAB()).collect().median'])

    stats_expr = [x % (root, samples_filter_expr) for x in stats]

    return stats_expr


def run_samples_sanity_checks(vds, reference_vds, n_samples=10, verbose=True):
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


def merge_schemas(vdses):

    vds_schemas = [vds.variant_schema for vds in vdses]

    for s in vds_schemas[1:]:
        if not isinstance(vds_schemas[0], type(s)):
            logger.fatal("Cannot merge schemas as the root (va) is of different type: %s and %s", vds_schemas[0], s)
            sys.exit(1)

    if not isinstance(vds_schemas[0], TStruct):
        return vdses

    anns = [flatten_struct(s, root='va') for s in vds_schemas]

    all_anns = {}
    for i in reversed(range(len(vds_schemas))):
        common_keys = set(all_anns.keys()).intersection(anns[i].keys())
        for k in common_keys:
            if not isinstance(all_anns[k].typ, type(anns[i][k].typ)):
                logger.fatal(
                    "Cannot merge schemas as annotation %s type %s found in VDS %d is not the same as previously existing type %s"
                    % (k, anns[i][k].typ, i, all_anns[k].typ))
                sys.exit(1)
        all_anns.update(anns[i])

    for i, vds in enumerate(vdses):
        vds = vds.annotate_variants_expr(["%s = NA: %s" % (k, str(v.typ)) for k, v in
                                          all_anns.iteritems() if k not in anns[i]])
        for ann, f in all_anns.iteritems():
            vds = vds.set_va_attributes(ann, f.attributes)

    return vdses


def copy_schema_attributes(vds1, vds2):
    anns1 = flatten_struct(vds1.variant_schema, root='va')
    anns2 = flatten_struct(vds2.variant_schema, root='va')
    for ann in anns1.keys():
        if ann in anns2:
            vds1 = vds1.set_va_attributes(ann, anns2[ann].attributes)

    return vds1


def print_attributes(vds, path=None):
    anns = flatten_struct(vds.variant_schema, root='va')
    if path is not None:
        print "%s attributes: %s" % (path, anns[path].attributes)
    else:
        for ann, f in anns.iteritems():
            print "%s attributes: %s" % (ann, f.attributes)


def get_numbered_annotations(vds, root='va.info'):
    """
    Get all 1-, A-, G- numbered annotations from a VDS based on the Number va attributes.
    In addition returns arrays with no Number or Number=. va attribute separately
    :param vds: Input VDS
    :param root: Place to find annotations (defaults to va.info)
    :return: annotations, a_annotations, g_annotations, dot_annotations as list[Field]
    """
    a_annotations = []
    g_annotations = []
    dot_annotations = []
    annotations = []

    release_info = get_ann_type(root, vds.variant_schema)
    for field in release_info.fields:
        if isinstance(field.typ, TArray):
            if 'Number' in field.attributes:
                number = field.attributes['Number']
                if number == "A":
                    a_annotations.append(field)
                elif number == "G":
                    g_annotations.append(field)
                else:
                    dot_annotations.append(field)
        else:
            annotations.append(field)

    logger.info("Found the following fields:")
    logger.info("1-based annotations: " + ",".join([x.name for x in annotations]))
    logger.info("A-based annotations: " + ",".join([x.name for x in a_annotations]))
    logger.info("G-based annotations: " + ",".join([x.name for x in g_annotations]))
    logger.info("dot annotations: " + ",".join([x.name for x in dot_annotations]))

    return annotations, a_annotations, g_annotations, dot_annotations


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
    """

    pca_loadings_type = get_ann_type(pca_loadings_root, pc_vds.variant_schema)  # TODO: this isn't used?

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


def read_list_data(input_file):
    if input_file.startswith('gs://'):
        hadoop_copy(input_file, 'file:///' + input_file.split("/")[-1])
        f = gzip.open("/" + os.path.basename(input_file)) if input_file.endswith('gz') else open( "/" + os.path.basename(input_file))
    else:
        f = gzip.open(input_file) if input_file.endswith('gz') else open(input_file)
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def rename_samples(vds, input_file, filter_to_samples_in_file=False):
    names = {old: new for old, new in [x.split("\t") for x in read_list_data(input_file)]}
    logger.info("Found %d samples for renaming in input file %s." % (len(names.keys()), input_file))
    logger.info("Renaming %d samples found in VDS" % len(set(names.keys()).intersection(set(vds.sample_ids)) ))

    if filter_to_samples_in_file:
        vds = vds.filter_samples_list(names.keys())
    return vds.rename_samples(names)


def filter_low_conf_regions(vds, filter_lcr=True, filter_decoy=True, high_conf_regions=None):
    """
    Filters low-confidence regions

    :param VariantDataset vds: VDS to filter
    :param bool filter_lcr: Whether to filter LCR regions
    :param bool filter_decoy: Wheter to filter Segdup regions
    :param list of str high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return:
    """

    if filter_lcr:
        vds = vds.filter_variants_table(KeyTable.import_interval_list(lcr_intervals_path), keep=False)

    if filter_decoy:
        vds = vds.filter_variants_table(KeyTable.import_interval_list(decoy_intervals_path), keep=False)

    if high_conf_regions is not None:
        for region in high_conf_regions:
            vds = vds.filter_variants_table(KeyTable.import_interval_list(region), keep=True)

    return vds


def process_consequences(vds, vep_root='va.vep', genes_to_string=True):
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq and worst_csq_suffix (worst consequence across transcripts) into [vep_root]

    :param VariantDataset vds: Input VDS
    :param str vep_root: Root for vep annotation (probably va.vep)
    :return: VDS with better formatted consequences
    :rtype: VariantDataset
    """
    if vep_root + '.worst_csq' in flatten_struct(vds.variant_schema, root='va'):
        vds = (vds.annotate_variants_expr('%(vep)s.transcript_consequences = '
                                          ' %(vep)s.transcript_consequences.map('
                                          '     csq => drop(csq, most_severe_consequence)'
                                          ')' % {'vep': vep_root}))
    vds = (vds.annotate_global('global.csqs', CSQ_ORDER, TArray(TString()))
           .annotate_variants_expr(
        '%(vep)s.transcript_consequences = '
        '   %(vep)s.transcript_consequences.map(csq => '
        '   let worst_csq = global.csqs.find(c => csq.consequence_terms.toSet().contains(c)) in'
        # '   let worst_csq_suffix = if (csq.filter(x => x.lof == "HC").length > 0)'
        # '       worst_csq + "-HC" '
        # '   else '
        # '       if (csq.filter(x => x.lof == "LC").length > 0)'
        # '           worst_csq + "-LC" '
        # '       else '
        # '           if (csq.filter(x => x.polyphen_prediction == "probably_damaging").length > 0)'
        # '               worst_csq + "-probably_damaging"'
        # '           else'
        # '               if (csq.filter(x => x.polyphen_prediction == "possibly_damaging").length > 0)'
        # '                   worst_csq + "-possibly_damaging"'
        # '               else'
        # '                   worst_csq in'
        '   merge(csq, {most_severe_consequence: worst_csq'
        # ', most_severe_consequence_suffix: worst_csq_suffix'
        '})'
        ')' % {'vep': vep_root}
    ).annotate_variants_expr(
        '%(vep)s.worst_csq = global.csqs.find(c => %(vep)s.transcript_consequences.map(x => x.most_severe_consequence).toSet().contains(c)),'
        '%(vep)s.worst_csq_suffix = '
        'let csq = global.csqs.find(c => %(vep)s.transcript_consequences.map(x => x.most_severe_consequence).toSet().contains(c)) in '
        'if (%(vep)s.transcript_consequences.filter(x => x.lof == "HC" && x.lof_flags == "").length > 0)'
        '   csq + "-HC" '
        'else '
        '   if (%(vep)s.transcript_consequences.filter(x => x.lof == "HC").length > 0)'
        '       csq + "-HC-flag" '
        '   else '
        '       if (%(vep)s.transcript_consequences.filter(x => x.lof == "LC").length > 0)'
        '           csq + "-LC" '
        '       else '
        '           if (%(vep)s.transcript_consequences.filter(x => x.polyphen_prediction == "probably_damaging").length > 0)'
        '               csq + "-probably_damaging"'
        '           else'
        '               if (%(vep)s.transcript_consequences.filter(x => x.polyphen_prediction == "possibly_damaging").length > 0)'
        '                   csq + "-possibly_damaging"'
        '               else'
        '                   csq' % {'vep': vep_root}
    ).annotate_variants_expr(
        '{vep}.lof = "-HC" ~ {vep}.worst_csq_suffix, '
        '{vep}.worst_csq_genes = {vep}.transcript_consequences'
        '.filter(x => x.most_severe_consequence == {vep}.worst_csq).map(x => x.gene_symbol).toSet(){genes_to_string}'.format(
            vep=vep_root, genes_to_string='.mkString("|")' if genes_to_string else '')
    ))
    return vds


def filter_vep_to_canonical_transcripts(vds, vep_root='va.vep'):
    return vds.annotate_variants_expr(
        '{vep}.transcript_consequences = '
        '   {vep}.transcript_consequences.filter(csq => csq.canonical == 1)'.format(vep=vep_root))




def filter_vep(vds, vep_root='va.vep', canonical=False, synonymous=False):
    """
    Fairly specific function, but used by multiple scripts


    """
    if canonical: vds = filter_vep_to_canonical_transcripts(vds, vep_root=vep_root)
    vds = process_consequences(vds)
    if synonymous: vds = filter_vep_to_synonymous_variants(vds, vep_root=vep_root)

    return (vds.filter_variants_expr('!{}.transcript_consequences.isEmpty'.format(vep_root))
            .annotate_variants_expr('{0} = select({0}, transcript_consequences)'.format(vep_root)))


def filter_vep_to_synonymous_variants(vds, vep_root='va.vep'):
    return vds.annotate_variants_expr(
        '{vep}.transcript_consequences = '
        '   {vep}.transcript_consequences.filter(csq => csq.most_severe_consequence == "synonymous_variant")'.format(vep=vep_root))


def filter_rf_variants(vds):
    """
    Does what it says

    :param VariantDataset vds: Input VDS (assumed split, but AS_FilterStatus unsplit)
    :return: vds with only RF variants removed
    :rtype: VariantDataset
    """
    return (vds
            .annotate_variants_expr(index_into_arrays(['va.info.AS_FilterStatus']))
            .filter_variants_expr('va.info.AS_FilterStatus.toArray() != ["RF"]'))


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
    """
    return (kt
            .annotate('comb = [{}]'.format(', '.join(['{{k: "{0}", value: {0}}}'.format(x) for x in columns_to_melt])))
            .drop(columns_to_melt)
            .explode('comb')
            .annotate('{} = comb.k, {} = comb.value'.format(key_column_name, value_column_name))
            .drop('comb'))


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
    """

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


def filter_samples_then_variants(vds, sample_criteria, callstats_temp_location='va.callstats_temp', min_allele_count=0):
    """
    Filter out samples, then generate callstats to filter variants, then filter out monomorphic variants
    Assumes split VDS
    TODO: add split logic

    :param VariantDataset vds: Input VDS
    :param str sample_criteria: String to be passed to `filter_samples_expr` to filter samples
    :param str callstats_temp_location: Temporary location for callstats to use to determine variants to drop
    :param int min_allele_count: minimum allele count to filter (default 0 for monomorphic variants)

    :return: Filtered VDS
    :rtype: VariantDataset
    """
    vds = vds.filter_samples_expr(sample_criteria)
    vds = vds.annotate_variants_expr('{} = gs.callStats(g => v)'.format(callstats_temp_location))
    vds = vds.filter_variants_expr('{}.AC[1] > {}'.format(callstats_temp_location, min_allele_count))
    return vds.annotate_variants_expr('va = drop(va, {})'.format(callstats_temp_location.split('.', 1)[-1]))


def set_site_filters(vds, site_filters_dict, as_filters_root):
    """

    Sets site filters using (a) site-level filters based on expr, (b) allele-specific filters from as_filters_root.
    Notes:
    1) `NA` propagates, that is if any filter expr evaluates to `NA` the site filter will be `NA`
    2) Setting `site_filters_dict` to `{}` or `None` results in no site-level filter. Setting as_filters_root to `""` or `None` results in no AS filter.
    3) as_filters_root should be either of type Set (split VDS) or Array[Set] (non-split VDS)
    4) In order to keep existing sites filter, add them to site_filters_dict. (e.g. {'InbreedingCoeff' : 'va.filters.contains("InbreedingCoeff")'})

    :param VariantDataset vds: Input VDS
    :param dict of str:str site_filters_dict: Dictionary with name of filter as key, filter expr as value (Use None if no site filter)
    :param str as_filters_root: Allele-specific filter root (use None if no AS filter)
    :return: VDS with updated filters
    :rtype: VariantDataset
    """

    if site_filters_dict:
        site_filters_expr = ",".join(['if({0}) ["{1}"] else [""][:0]'.format(filter_expr, name) for name, filter_expr in site_filters_dict.iteritems()])
        site_filters_expr = '[{}].flatten.toSet'.format(site_filters_expr)
    else:
        site_filters_expr = '[""][:0].toSet'

    if as_filters_root:
        #Find out from schema if as_filters is a Set or an Array of Set
        if isinstance(get_ann_type(as_filters_root, vds.variant_schema), TSet):
            as_filters_expr = as_filters_root
        else:
            as_filters_expr = 'orMissing({0}.forall(x => isDefined(x) && x.forall(y => isDefined(y))),' \
                              '{0}.find(x => isDefined(x) && x.isEmpty())' \
                              '.orElse({0}.toSet.flatten))'.format(as_filters_root)
    else:
        site_filters_expr = '[""][:0].toSet'


    annotate_expr = ('va.filters = let filters = {}.union({}) in '
                     'orMissing(filters.forall(x => isDefined(x)), filters)'.format(site_filters_expr, as_filters_expr))

    logger.debug(annotate_expr)

    return vds.annotate_variants_expr(annotate_expr)


def split_vds_and_annotations(vds, hard_filters, as_filters_root, extra_ann_expr=None, vep_root='va.vep'):
    """
     Split VDS and its associated va annotations properly

     :param VariantDataset vds: Input VDS
     :param list of str hard_filters: All hard filter values to keep
     :param str as_filters_root: Root of AS filters (set to `None` if no AS filters)
     :param list of str extra_ann_expr: List of additional annotation expressions to pass to annotate_variants_expr
     :return: Split VDS with properly split va fields
     :rtype: VariantDataset
     """
    annotations, a_annotations, g_annotations, dot_annotations = get_numbered_annotations(vds, "va.info")

    vds = vds.split_multi()
    vds = vds.annotate_variants_expr(
        index_into_arrays(a_based_annotations=["va.info." + a.name for a in a_annotations], vep_root=vep_root))
    vds = set_site_filters(vds,
                           { f: 'va.filters.contains("{}")'.format(f) for f in hard_filters },
                           as_filters_root)
    ann_expr = ['va.info = drop(va.info, {0})'.format(",".join([a.name for a in g_annotations]))]
    if extra_ann_expr:
        ann_expr.extend(extra_ann_expr)
    vds = vds.annotate_variants_expr(ann_expr)
    return vds


def quote_field_name(f):
    """
    Given a field name, returns the name quote if necessary for Hail columns access.
    E.g.
    - The name contains a `.`
    - The name starts with a numeric

    :param str f: The field name
    :return: Quoted (or not) field name
    :rtype: str
    """

    return '`{}`'.format(f) if re.search('^\d|\.', f) else f



