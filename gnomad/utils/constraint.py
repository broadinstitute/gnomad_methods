# noqa: D100

from typing import Union

import hail as hl

def annotate_mutation_type(
    t: Union[hl.MatrixTable, hl.Table], trimer: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotate variant types.

    The following annotations are added to the output Table:
        - cpg
        - transition
        - variant_type
        - variant_type_model

    :param t: Input Table or MatrixTable.
    :param trimer: Whether to use trimers for context. Defaults to False (uses heptamers as context).
    :return: Table with annotations.
    """
    mid_index = 1 if trimer else 3
    transition_expr = hl.is_transition(t.ref, t.alt)
    cpg_expr = (
        (t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1 : mid_index] == "C")
    ) | (
        (t.ref == "C")
        & (t.alt == "T")
        & (t.context[mid_index + 1 : mid_index + 2] == "G")
    )
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (
        hl.case()
        .when(t.cpg, "CpG")
        .when(t.transition, "non-CpG transition")
        .default("transversion")
    )
    variant_type_model_expr = hl.if_else(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(
            mutation_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )
    else:
        return t.annotate(
            mutation_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )


def trimer_from_heptamer(
    t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Trim heptamer context to create trimer context.

    :param t: Input MatrixTable or Table with context annotation.
    :return: MatrixTable or Table with trimer context annotated.
    """
    trimer_expr = hl.if_else(hl.len(t.context) == 7, t.context[2:5], t.context)
    return (
        t.annotate_rows(context=trimer_expr)
        if isinstance(t, hl.MatrixTable)
        else t.annotate(context=trimer_expr)
    )


def collapse_strand(
    t: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Return the deduplicated context by collapsing DNA strands.

    Function returns the reverse complement for 'ref, 'alt', and 'context' if the reference allele is either 'G' or 'T'.
    
    The following annotations are added to the output Table:
        - was_flipped - whether the 'ref, 'alt', and 'context' were flipped

    :param ht: Input Table.
    :return: Table with deduplicated context annotation (ref, alt, context, was_flipped).
    """
    collapse_expr = {
        "ref": hl.if_else(
            ((t.ref == "G") | (t.ref == "T")), hl.reverse_complement(t.ref), t.ref
        ),
        "alt": hl.if_else(
            ((t.ref == "G") | (t.ref == "T")), hl.reverse_complement(t.alt), t.alt
        ),
        "context": hl.if_else(
            ((t.ref == "G") | (t.ref == "T")),
            hl.reverse_complement(t.context),
            t.context,
        ),
        "was_flipped": (t.ref == "G") | (t.ref == "T"),
    }
    return (
        t.annotate(**collapse_expr)
        if isinstance(t, hl.Table)
        else t.annotate_rows(**collapse_expr)
    )
