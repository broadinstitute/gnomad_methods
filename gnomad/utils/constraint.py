from typing import Union

import hail as hl


def annotate_variant_types(
    t: Union[hl.MatrixTable, hl.Table], trimer: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotate variant types.
    
    The following annotations are added to the output Table:
        - cpg
        - transition
        - variant_type
        - variant_type_model columns

    :param t: Input Table or MatrixTable.
    :param trimer: Whether to use trimers or heptamers, defaults to False.
    :return: Table with annotations.
    """
    mid_index = 1 if trimer else 3
    transition_expr = (
        ((t.ref == "A") & (t.alt == "G"))
        | ((t.ref == "G") & (t.alt == "A"))
        | ((t.ref == "T") & (t.alt == "C"))
        | ((t.ref == "C") & (t.alt == "T"))
    )
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
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )
    else:
        return t.annotate(
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
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
    ht: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Return the deduplicated context by collapsing DNA strands.
    
    Function returns the reverse complement if the reference allele is either 'G' or 'T'.
    
    The following annotations are added to the output Table:
        - was_flipped
    
    :param ht: Input Table.
    :return: Table with deduplicated context annotation (ref, alt, context, was_flipped).
    """
    collapse_expr = {
        "ref": hl.if_else(
            ((ht.ref == "G") | (ht.ref == "T")), hl.reverse_complement(ht.ref), ht.ref
        ),
        "alt": hl.if_else(
            ((ht.ref == "G") | (ht.ref == "T")), hl.reverse_complement(ht.alt), ht.alt
        ),
        "context": hl.if_else(
            ((ht.ref == "G") | (ht.ref == "T")),
            hl.reverse_complement(ht.context),
            ht.context,
        ),
        "was_flipped": (ht.ref == "G") | (ht.ref == "T"),
    }
    return (
        ht.annotate(**collapse_expr)
        if isinstance(ht, hl.Table)
        else ht.annotate_rows(**collapse_expr)
    )
