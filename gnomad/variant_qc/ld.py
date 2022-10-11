# noqa: D100

import hail as hl
from hail.linalg import BlockMatrix

from gnomad.resources.grch37.gnomad import public_release
from gnomad.resources.grch37.gnomad_ld import ld_index, ld_matrix


def get_r_human_readable(
    pop: str, var1: str, var2: str, ref_genome: str = "GRCh37"
):  # noqa: D103
    bm = ld_matrix(pop).bm()
    ht = ld_index(pop).ht()
    chrom, pos, ref, alt = var1.split("-")
    var1 = (hl.parse_locus(f"{chrom}:{pos}", ref_genome), [ref, alt])
    chrom, pos, ref, alt = var2.split("-")
    var2 = (hl.parse_locus(f"{chrom}:{pos}", ref_genome), [ref, alt])
    return get_r_for_pair_of_variants(bm, ht, var1, var2)


# TODO: find LD proxies


def get_r_for_pair_of_variants(
    bm: BlockMatrix,
    ld_index: hl.Table,
    var1: (hl.tlocus, hl.tarray(hl.tstr)),
    var2: (hl.tlocus, hl.tarray(hl.tstr)),
):
    """
    Get `r` value (LD) for pair of variants `var1` and `var2`.

    .. code-block:: python

        bm = get_ld_matrix('nfe')
        ld_index = get_ld_index('nfe')
        var1 = (hl.parse_locus('1:10146', 'GRCh37'), ['AC', 'A'])
        var2 = (hl.parse_locus('1:10151', 'GRCh37'), ['TA', 'T'])
        get_r_for_pair_of_variants(bm, ld_index, var1, var2)
        # 0.01789767935482124

    :param bm: Input BlockMatrix
    :param ld_index: Corresponding index table
    :param var1: Tuple of locus and alleles
    :param var2: Tuple of locus and alleles
    :return: Correlation (r) between two variants
    """
    idx1 = ld_index.filter(
        (ld_index.locus == var1[0]) & (ld_index.alleles == var1[1])
    ).idx.collect()[0]
    idx2 = ld_index.filter(
        (ld_index.locus == var2[0]) & (ld_index.alleles == var2[1])
    ).idx.collect()[0]

    if idx1 > idx2:
        temp = idx1
        idx1 = idx2
        idx2 = temp

    return bm[idx1, idx2]


def get_r_within_gene_in_pop(pop: str, gene: str):
    """
    Get LD information (`r`) for all pairs of variants within `gene` for a given `pop`.

    Warning: this returns a table quadratic in number of variants. Exercise caution with large genes.

    :param pop: Population for which to get LD information
    :param gene: Gene symbol as string
    :return: Table with pairs of variants
    """
    return get_r_within_gene(
        ld_matrix(pop).bm(), ld_index(pop).ht(), gene, None, "GRCh37"
    )


def get_r_within_gene(
    bm: BlockMatrix,
    ld_index: hl.Table,
    gene: str,
    vep_ht: hl.Table = None,
    reference_genome: str = None,
):
    """
    Get LD information (`r`) for all pairs of variants within `gene`.

    Warning: this returns a table quadratic in number of variants. Exercise caution with large genes.

    :param bm: Input Block Matrix
    :param ld_index: Corresponding index table
    :param gene: Gene symbol as string
    :param vep_ht: Table with VEP annotations (if None, gets from get_gnomad_public_data())
    :param reference_genome: Reference genome to pass to get_gene_intervals for fast filtering to gene
    :return: Table with pairs of variants
    """
    if vep_ht is None:
        vep_ht = public_release("exomes").ht()
    if reference_genome is None:
        reference_genome = hl.default_reference().name
    intervals = hl.experimental.get_gene_intervals(
        gene_symbols=[gene], reference_genome=reference_genome
    )
    ld_index = hl.filter_intervals(ld_index, intervals)
    ld_index = ld_index.annotate(vep=vep_ht[ld_index.key].vep)
    ld_index = ld_index.filter(
        hl.any(lambda tc: tc.gene_symbol == gene, ld_index.vep.transcript_consequences)
    )

    indices_to_keep = ld_index.idx.collect()
    filt_bm = bm.filter(indices_to_keep, indices_to_keep)
    ht = filt_bm.entries()
    ld_index = ld_index.add_index("new_idx").key_by("new_idx")
    return ht.transmute(r=ht.entry, i_variant=ld_index[ht.i], j_variant=ld_index[ht.j])
