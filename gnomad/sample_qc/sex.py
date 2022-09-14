# noqa: D100

import logging
from typing import Tuple, Union

import hail as hl
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SEXES = {"Male": "Male", "Female": "Female"}


def adjusted_sex_ploidy_expr(
    locus_expr: hl.expr.LocusExpression,
    gt_expr: hl.expr.CallExpression,
    karyotype_expr: hl.expr.StringExpression,
    xy_karyotype_str: str = "XY",
    xx_karyotype_str: str = "XX",
) -> hl.expr.CallExpression:
    """
    Create an entry expression to convert males to haploid on non-PAR X/Y and females to missing on Y.

    :param locus_expr: Locus
    :param gt_expr: Genotype
    :param karyotype_expr: Karyotype
    :param xy_karyotype_str: Male sex karyotype representation
    :param xx_karyotype_str: Female sex karyotype representation
    :return: Genotype adjusted for sex ploidy
    """
    male = karyotype_expr == xy_karyotype_str
    female = karyotype_expr == xx_karyotype_str
    x_nonpar = locus_expr.in_x_nonpar()
    y_par = locus_expr.in_y_par()
    y_nonpar = locus_expr.in_y_nonpar()
    return (
        hl.case(missing_false=True)
        .when(female & (y_par | y_nonpar), hl.null(hl.tcall))
        .when(male & (x_nonpar | y_nonpar) & gt_expr.is_het(), hl.null(hl.tcall))
        .when(male & (x_nonpar | y_nonpar), hl.call(gt_expr[0], phased=False))
        .default(gt_expr)
    )


def adjust_sex_ploidy(
    mt: hl.MatrixTable,
    sex_expr: hl.expr.StringExpression,
    male_str: str = "male",
    female_str: str = "female",
) -> hl.MatrixTable:
    """
    Convert males to haploid on non-PAR X/Y, sets females to missing on Y.

    :param mt: Input MatrixTable
    :param sex_expr: Expression pointing to sex in MT (if not male_str or female_str, no change)
    :param male_str: String for males (default 'male')
    :param female_str: String for females (default 'female')
    :return: MatrixTable with fixed ploidy for sex chromosomes
    """
    return mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, sex_expr, male_str, female_str)
    )


def gaussian_mixture_model_karyotype_assignment(
    sex_ht: hl.Table,
    chrx_ploidy_expr: Union[hl.expr.NumericExpression, str] = "chrX_ploidy",
    chry_ploidy_expr: Union[hl.expr.NumericExpression, str] = "chrY_ploidy",
    karyotype_output_prefix: str = "gmm",
) -> hl.Table:
    """
    Annotate the input Table with an X karyotype, Y karyotype, and sex karyotype based on a gaussian mixture model.

    This function uses two component gaussian mixture models on `chrx_ploidy_expr` and `chry_ploidy_expr` to assign
    an X karyotype and a Y karyotype which are then combined into the sex karyotype.

    The following annotations are added:
        - `karyotype_output_prefix`_x_karyotype
        - `karyotype_output_prefix`_y_karyotype
        - `karyotype_output_prefix`_karyotype = `karyotype_output_prefix`_x_karyotype + `karyotype_output_prefix`_y_karyotype

    .. note::

        This uses a two component gaussian mixture model so all samples are given one of the following sex karyotypes:
        X, XX, XY, YY. It's recommended that this annotation is only used to split samples into XX and
        XY groups that can then be used in `get_ploidy_cutoffs` to determine XX and XY ploidy means and stdevs.

    :param sex_ht: Input Table with chromosome X and chromosome Y ploidy values.
    :param chrx_ploidy_expr: Expression pointing to chromosome X ploidy in `sex_ht`. Default is 'chrX_ploidy'.
    :param chry_ploidy_expr: Expression pointing to chromosome Y ploidy in `sex_ht`. Default is 'chrY_ploidy'.
    :param karyotype_output_prefix: String to use as the prefix for the gaussian mixture model karyotype output.
    :return: Input Table with gaussian mixture model karyotype annotations added.
    """
    if isinstance(chrx_ploidy_expr, str):
        chrx_ploidy_expr = sex_ht[chrx_ploidy_expr]
    if isinstance(chry_ploidy_expr, str):
        chry_ploidy_expr = sex_ht[chry_ploidy_expr]

    sex_pd = sex_ht.select(
        chrX_ploidy=chrx_ploidy_expr,
        chrY_ploidy=chry_ploidy_expr,
    ).to_pandas()

    def _run_gaussian_mixture_model(feature, karyotypes, karyotype_name):
        df = sex_pd[["s", feature]].set_index("s")
        gmm = GaussianMixture(n_components=2)
        gmm.fit(df)
        probs = gmm.predict_proba(df)
        cluster_to_karyotype = dict(
            zip(np.argsort([m[0] for m in gmm.means_]), karyotypes)
        )

        df[f"{feature}_cluster"] = gmm.predict(df)
        df[
            f"{feature if karyotype_name is None else karyotype_name}_karyotype"
        ] = df.apply(
            lambda row: cluster_to_karyotype[row[f"{feature}_cluster"]], axis=1
        )
        for i in cluster_to_karyotype:
            df[f"{feature}_prob_{cluster_to_karyotype[i]}"] = probs[:, i]

        return df

    x_df = _run_gaussian_mixture_model(
        "chrX_ploidy", ["X", "XX"], f"{karyotype_output_prefix}_x"
    )
    y_df = _run_gaussian_mixture_model(
        "chrY_ploidy", ["", "Y"], f"{karyotype_output_prefix}_y"
    )
    xy_df = pd.concat(
        [
            x_df[f"{karyotype_output_prefix}_x_karyotype"],
            y_df[f"{karyotype_output_prefix}_y_karyotype"],
        ],
        axis=1,
    )
    xy_df[f"{karyotype_output_prefix}_karyotype"] = (
        xy_df[f"{karyotype_output_prefix}_x_karyotype"]
        + xy_df[f"{karyotype_output_prefix}_y_karyotype"]
    )
    xy_ht = hl.Table.from_pandas(xy_df.reset_index(), key=["s"])

    return sex_ht.annotate(**xy_ht[sex_ht.key])


def get_ploidy_cutoffs(
    ht: hl.Table,
    f_stat_cutoff: float = None,
    normal_ploidy_cutoff: int = 5,
    aneuploidy_cutoff: int = 6,
    group_by_expr: hl.expr.StringExpression = None,
) -> Tuple[Tuple[float, Tuple[float, float], float], Tuple[Tuple[float, float], float]]:
    """
    Get chromosome X and Y ploidy cutoffs for XY and XX samples.

    .. note::

        This assumes the input hail Table has the fields chrX_ploidy, and chrY_ploidy, and f_stat if `f_stat_cutoff` is
        set.

    Return a tuple of sex chromosome ploidy cutoffs: ((x_ploidy_cutoffs), (y_ploidy_cutoffs)).
    x_ploidy_cutoffs: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for triple X)
    y_ploidy_cutoffs: ((lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)

    Uses the normal_ploidy_cutoff parameter to determine the ploidy cutoffs for XX and XY karyotypes.
    Uses the aneuploidy_cutoff parameter to determine the cutoffs for sex aneuploidies.

    .. note::

        `f_stat_cutoff` or `group_by_expr` must be supplied. If `f_stat_cutoff` is supplied then f-stat is used to
        split the samples into roughly 'XX' and 'XY'. If `group_by_expr` is supplied instead, then it must include an
        annotation grouping samples by 'XX' and 'XY'. These are both only used to divide samples into XX and XY to
        determine means and standard deviations for these categories and are not used in the final karyotype annotation.

    :param ht: Table with f_stat and sex chromosome ploidies
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :param normal_ploidy_cutoff: Number of standard deviations to use when determining sex chromosome ploidy cutoffs for XX, XY karyotypes.
    :param aneuploidy_cutoff: Number of standard deviations to use when sex chromosome ploidy cutoffs for aneuploidies.
    :param group_by_expr: Expression grouping samples into 'XX' and 'XY'. Can be used instead of and `f_stat_cutoff`.
    :return: Tuple of ploidy cutoff tuples: ((x_ploidy_cutoffs), (y_ploidy_cutoffs))
    """
    if f_stat_cutoff is None and group_by_expr is None:
        raise ValueError("One of 'f_stat_cutoff' or 'group_by_expr' must be supplied!")

    # If 'f_stat_cutoff' is supplied, group the sex chromosome ploidy table by f_stat cutoff
    if f_stat_cutoff is not None:
        group_by_expr = hl.cond(ht.f_stat < f_stat_cutoff, "xx", "xy")

    # Get mean/stdev for chrX/Y ploidies based on 'group_by_expr'
    sex_stats = ht.aggregate(
        hl.agg.group_by(
            group_by_expr,
            hl.struct(x=hl.agg.stats(ht.chrX_ploidy), y=hl.agg.stats(ht.chrY_ploidy)),
        )
    )
    if "xx" not in sex_stats:
        raise ValueError("No samples are grouped as XX!")
    if "xy" not in sex_stats:
        raise ValueError("No samples are grouped as XY!")
    logger.info("XX stats: %s", sex_stats["xx"])
    logger.info("XY stats: %s", sex_stats["xy"])

    cutoffs = (
        (
            sex_stats["xy"].x.mean + (normal_ploidy_cutoff * sex_stats["xy"].x.stdev),
            (
                sex_stats["xx"].x.mean
                - (normal_ploidy_cutoff * sex_stats["xx"].x.stdev),
                sex_stats["xx"].x.mean
                + (normal_ploidy_cutoff * sex_stats["xx"].x.stdev),
            ),
            sex_stats["xx"].x.mean + (aneuploidy_cutoff * sex_stats["xx"].x.stdev),
        ),
        (
            (
                sex_stats["xx"].y.mean
                + (normal_ploidy_cutoff * sex_stats["xx"].y.stdev),
                sex_stats["xy"].y.mean
                + (normal_ploidy_cutoff * sex_stats["xy"].y.stdev),
            ),
            sex_stats["xy"].y.mean + (aneuploidy_cutoff * sex_stats["xy"].y.stdev),
        ),
    )

    logger.info("X ploidy cutoffs: %s", cutoffs[0])
    logger.info("Y ploidy cutoffs: %s", cutoffs[1])
    return cutoffs


def get_sex_expr(
    chr_x_ploidy: hl.expr.NumericExpression,
    chr_y_ploidy: hl.expr.NumericExpression,
    x_ploidy_cutoffs: Tuple[float, Tuple[float, float], float],
    y_ploidy_cutoffs: Tuple[Tuple[float, float], float],
) -> hl.expr.StructExpression:
    """
    Create a struct with X_karyotype, Y_karyotype, and sex_karyotype.

    Note that X0 is currently returned as 'X'.

    :param chr_x_ploidy: Chromosome X ploidy (or relative ploidy)
    :param chr_y_ploidy: Chromosome Y ploidy (or relative ploidy)
    :param x_ploidy_cutoffs: Tuple of X chromosome ploidy cutoffs: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for triple X)
    :param y_ploidy_cutoffs: Tuple of Y chromosome ploidy cutoffs: ((lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)
    :return: Struct containing X_karyotype, Y_karyotype, and sex_karyotype
    """
    sex_expr = hl.struct(
        X_karyotype=(
            hl.case()
            .when(chr_x_ploidy < x_ploidy_cutoffs[0], "X")
            .when(
                (
                    (chr_x_ploidy > x_ploidy_cutoffs[1][0])
                    & (chr_x_ploidy < x_ploidy_cutoffs[1][1])
                ),
                "XX",
            )
            .when((chr_x_ploidy >= x_ploidy_cutoffs[2]), "XXX")
            .default("ambiguous")
        ),
        Y_karyotype=(
            hl.case()
            .when(chr_y_ploidy < y_ploidy_cutoffs[0][0], "")
            .when(
                (
                    (chr_y_ploidy > y_ploidy_cutoffs[0][0])
                    & (chr_y_ploidy < y_ploidy_cutoffs[0][1])
                ),
                "Y",
            )
            .when(chr_y_ploidy >= y_ploidy_cutoffs[1], "YY")
            .default("ambiguous")
        ),
    )

    return sex_expr.annotate(
        sex_karyotype=hl.if_else(
            (sex_expr.X_karyotype == "ambiguous")
            | (sex_expr.Y_karyotype == "ambiguous"),
            "ambiguous",
            sex_expr.X_karyotype + sex_expr.Y_karyotype,
        )
    )
