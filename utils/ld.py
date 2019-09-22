from gnomad_hail.resources import *
from hail.linalg import BlockMatrix


def ld_matrix_path(data_type: str, pop: str, common_only: bool = True, adj: bool = True, version: str = CURRENT_RELEASE):
    if data_type == 'genomes_snv_sv':
        return f'gs://gnomad-resources/snv_sv_ld/gnomad.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.bm'
    else:
        return f'gs://gnomad-public/release/{version}/ld/gnomad.{data_type}.r{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.bm'


def ld_index_path(data_type: str, pop: str, common_only: bool = True, adj: bool = True, version: str = CURRENT_RELEASE):
    if data_type == 'genomes_snv_sv':
        return f'gs://gnomad-resources/snv_sv_ld/gnomad.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.variant_indices.ht'
    else:
        return f'gs://gnomad-public/release/{version}/ld/gnomad.{data_type}.r{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.variant_indices.ht'


def ld_scores_path(data_type: str, pop: str, adj: bool = True, version: str = CURRENT_RELEASE):
    return f'gs://gnomad-public/release/{version}/ld/scores/gnomad.{data_type}.r{version}.{pop}.{"adj." if adj else ""}ld_scores.ht'


def get_ld_matrix(pop: str):
    return BlockMatrix.read(ld_matrix_path('genomes', pop))


def get_ld_index(pop: str):
    return hl.read_table(ld_index_path('genomes', pop))


def get_ld_scores(pop: str):
    return hl.read_table(ld_scores_path('genomes', pop))


def get_r_human_readable(pop: str, var1: str, var2: str, ref_genome: str = 'GRCh37'):
    bm = get_ld_matrix(pop)
    ht = get_ld_index(pop)
    chrom, pos, ref, alt = var1.split('-')
    var1 = (hl.parse_locus(f'{chrom}:{pos}', ref_genome), [ref, alt])
    chrom, pos, ref, alt = var2.split('-')
    var2 = (hl.parse_locus(f'{chrom}:{pos}', ref_genome), [ref, alt])
    return get_r_for_pair_of_variants(bm, ht, var1, var2)


# TODO: find LD proxies


def get_r_for_pair_of_variants(bm: BlockMatrix, ld_index: hl.Table,
                               var1: (hl.tlocus, hl.tarray(hl.tstr)),
                               var2: (hl.tlocus, hl.tarray(hl.tstr))):
    """
    Get `r` value (LD) for pair of variants `var1` and `var2`.

    bm = get_ld_matrix('nfe')
    ld_index = get_ld_index('nfe')
    var1 = (hl.parse_locus('1:10146', 'GRCh37'), ['AC', 'A'])
    var2 = (hl.parse_locus('1:10151', 'GRCh37'), ['TA', 'T'])
    get_r_for_pair_of_variants(bm, ld_index, var1, var2)
    # 0.01789767935482124

    :param BlockMatrix bm: Input BlockMatrix
    :param Table ld_index: Corresponding index table
    :param tuple var1: Tuple of locus and alleles
    :param tuple var2: Tuple of locus and alleles
    :return: Correlation (r) between two variants
    :rtype: float
    """
    idx1 = ld_index.filter((ld_index.locus == var1[0]) & (ld_index.alleles == var1[1])).idx.collect()[0]
    idx2 = ld_index.filter((ld_index.locus == var2[0]) & (ld_index.alleles == var2[1])).idx.collect()[0]

    if idx1 > idx2:
        temp = idx1
        idx1 = idx2
        idx2 = temp

    return bm[idx1, idx2]


def get_r_within_gene_in_pop(pop: str, gene: str):
    """
    Gets LD information (`r`) for all pairs of variants within `gene` for a given `pop`.

    Warning: this returns a table quadratic in number of variants. Exercise caution with large genes.

    :param str pop: Population for which to get LD information
    :param str gene: Gene symbol as string
    :return: Table with pairs of variants
    :rtype: Table
    """
    return get_r_within_gene(get_ld_matrix(pop), get_ld_index(pop), gene, None, 'GRCh37')


def get_r_within_gene(bm: BlockMatrix, ld_index: hl.Table, gene: str, vep_ht: hl.Table = None, reference_genome: str = None):
    """
    Gets LD information (`r`) for all pairs of variants within `gene`.

    Warning: this returns a table quadratic in number of variants. Exercise caution with large genes.

    :param BlockMatrix bm: Input Block Matrix
    :param Table ld_index: Corresponding index table
    :param str gene: Gene symbol as string
    :param Table vep_ht: Table with VEP annotations (if None, gets from get_gnomad_public_data())
    :param str reference_genome: Reference genome to pass to get_gene_intervals for fast filtering to gene
    :return: Table with pairs of variants
    :rtype: Table
    """
    if vep_ht is None:
        vep_ht = get_gnomad_public_data('exomes')
    if reference_genome is None:
        reference_genome = hl.default_reference().name
    intervals = hl.experimental.get_gene_intervals(gene_symbols=[gene], reference_genome=reference_genome)
    ld_index = hl.filter_intervals(ld_index, intervals)
    ld_index = ld_index.annotate(vep=vep_ht[ld_index.key].vep)
    ld_index = ld_index.filter(hl.any(lambda tc: tc.gene_symbol == gene, ld_index.vep.transcript_consequences))

    indices_to_keep = ld_index.idx.collect()
    filt_bm = bm.filter(indices_to_keep, indices_to_keep)
    ht = filt_bm.entries()
    ld_index = ld_index.add_index('new_idx').key_by('new_idx')
    return ht.transmute(r=ht.entry, i_variant=ld_index[ht.i], j_variant=ld_index[ht.j])

