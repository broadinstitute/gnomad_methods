import logging
from typing import List, Optional, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_reference_ht(
    ref: hl.ReferenceGenome,
    contigs: Optional[List[str]] = None,
    excluded_intervals: Optional[List[hl.Interval]] = None,
    add_all_substitutions: bool = False,
    filter_n: bool = True,
) -> hl.Table:
    """
    Creates a reference Table with locus and alleles (containing only the reference allele by default) from the given reference genome.

    .. note::

        If the `contigs` argument is not provided, all contigs (including obscure ones) will be added to the table.
        This can be slow as contigs are added one by one.

    :param ref: Input reference genome
    :param contigs: An optional list of contigs that the Table should include
    :param excluded_intervals: An optional list of intervals to exclude
    :param add_all_substitutions: If set, then all possible substitutions are added in the alleles array
    :param filter_n: If set, bases where the reference is unknown (n) are filtered.
    :return:
    """
    if not ref.has_sequence():
        add_reference_sequence(ref)

    if not contigs:
        contigs = ref.contigs

    if add_all_substitutions:
        SUBSTITUTIONS_TABLE = hl.literal(
            {
                "a": ["c", "g", "t"],
                "c": ["a", "g", "t"],
                "g": ["a", "c", "t"],
                "t": ["a", "c", "g"],
            }
        )

    context = []
    for contig in contigs:
        n_partitions = max(1, int(ref.contig_length(contig) / 5000000))
        logger.info(
            f"Creating reference contig {contig} with {n_partitions} partitions."
        )
        _context = hl.utils.range_table(
            ref.contig_length(contig), n_partitions=n_partitions
        )

        locus_expr = hl.locus(contig=contig, pos=_context.idx + 1, reference_genome=ref)
        ref_allele_expr = locus_expr.sequence_context().lower()
        if add_all_substitutions:
            alleles_expr = hl.array([ref_allele_expr]).extend(
                SUBSTITUTIONS_TABLE.get(ref_allele_expr, hl.empty_array(hl.tstr))
            )
        else:
            alleles_expr = [ref_allele_expr]

        _context = (
            _context.select(locus=locus_expr, alleles=alleles_expr)
            .key_by("locus", "alleles")
            .drop("idx")
        )

        if excluded_intervals is not None:
            _context = hl.filter_intervals(_context, excluded_intervals, keep=False)

        if filter_n:
            _context = _context.filter(_context.alleles[0] == "n", keep=False)

        context.append(_context)

    return context.pop().union(*context)


def add_reference_sequence(ref: hl.ReferenceGenome) -> hl.ReferenceGenome:
    """
    Adds the fasta sequence to a Hail reference genome.
    Only GRCh37 and GRCh38 references are supported.

    :param ref: Input reference genome.
    :return:
    """
    if not ref.has_sequence():
        if ref.name == "GRCh38":
            ref.add_sequence(
                "gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz",
                "gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai",
            )
        elif ref.name == "GRCh37":
            ref.add_sequence(
                "gs://hail-common/references/human_g1k_v37.fasta.gz",
                "gs://hail-common/references/human_g1k_v37.fasta.fai",
            )
        else:
            raise NotImplementedError(
                f"No known location for the fasta/fai files for genome {ref.name}. Only GRCh37 and GRCh38 are supported at this time."
            )
    else:
        logger.info(
            "Reference genome sequence already present. Ignoring add_reference_sequence."
        )

    return ref


def get_reference_genome(
    locus: Union[hl.expr.LocusExpression, hl.expr.IntervalExpression],
    add_sequence: bool = False,
) -> hl.ReferenceGenome:
    """
    Returns the reference genome associated with the input Locus expression

    :param locus: Input locus
    :param add_sequence: If set, the fasta sequence is added to the reference genome
    :return: Reference genome
    """
    if isinstance(locus, hl.expr.LocusExpression):
        ref = locus.dtype.reference_genome
    else:
        assert isinstance(locus, hl.expr.IntervalExpression)
        ref = locus.dtype.point_type.reference_genome
    if add_sequence:
        ref = add_reference_sequence(ref)
    return ref
