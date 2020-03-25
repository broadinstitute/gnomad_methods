from typing import Union

import hail as hl
import pandas as pd


def expand_pd_array_col(
        df: pd.DataFrame,
        array_col: str,
        num_out_cols: int = 0,
        out_cols_prefix=None,
        out_1based_indexing: bool = True
) -> pd.DataFrame:
    """
    Expands a Dataframe column containing an array into multiple columns.

    :param df: input dataframe
    :param array_col: Column containing the array
    :param num_out_cols: Number of output columns. If set, only the `n_out_cols` first elements of the array column are output.
                             If <1, the number of output columns is equal to the length of the shortest array in `array_col`
    :param out_cols_prefix: Prefix for the output columns (uses `array_col` as the prefix unless set)
    :param out_1based_indexing: If set, the output column names indexes start at 1. Otherwise they start at 0.
    :return: dataframe with expanded columns
    """

    if out_cols_prefix is None:
        out_cols_prefix = array_col

    if num_out_cols < 1:
        num_out_cols = min([len(x) for x in df[array_col].values.tolist()])

    cols = ['{}{}'.format(out_cols_prefix, i + out_1based_indexing) for i in range(num_out_cols)]
    df[cols] = pd.DataFrame(df[array_col].values.tolist())[list(range(num_out_cols))]

    return df


def get_array_element_type(array_expr: hl.expr.ArrayExpression) -> hl.HailType:
    """
    Returns the type of an array element.

    :param array_expr: The array expression to get the element type
    :return: Hail type
    """
    return array_expr.dtype.element_type


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression selecting bi-allelic sites only,
    accounting for whether the input MT/HT was split.

    :param t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    """
    return (~t.was_split if 'was_split' in t.row else (hl.len(t.alleles) == 2))
