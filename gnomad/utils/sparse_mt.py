# noqa: D100

import logging
from typing import Callable, Dict, List, Optional, Set, Tuple, Union

import hail as hl

from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    _read_reduction_globals,
    agg_by_strata,
    annotate_adj,
    expand_strata_array_from_leaves,
    fs_from_sb,
    generate_freq_group_membership_array,
    get_adj_expr,
    get_lowqual_expr,
    merge_array_expressions,
    merge_histograms,
    pab_max_expr,
    sor_from_sb,
)
from gnomad.utils.intervals import interval_length, union_intervals
from gnomad.utils.reference_genome import get_reference_genome

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

INFO_AGG_FIELDS = {
    "sum_agg_fields": ["QUALapprox"],
    "int32_sum_agg_fields": ["VarDP"],
    "median_agg_fields": ["ReadPosRankSum", "MQRankSum"],
    "array_sum_agg_fields": ["SB", "RAW_MQandDP"],
}

AS_INFO_AGG_FIELDS = {
    "sum_agg_fields": ["AS_QUALapprox", "AS_RAW_MQ"],
    "int32_sum_agg_fields": ["AS_VarDP"],
    "median_agg_fields": ["AS_RAW_ReadPosRankSum", "AS_RAW_MQRankSum"],
    "array_sum_agg_fields": ["AS_SB_TABLE"],
}


def compute_last_ref_block_end(mt: hl.MatrixTable) -> hl.Table:
    """
    Compute the genomic position of the most upstream reference block overlapping each row on a sparse MT.

    Note that since reference blocks do not extend beyond contig boundaries, only the position is kept.

    This function returns a Table with that annotation.  (`last_END_position`).

    :param mt: Input MatrixTable
    :return: Output Table with `last_END_position` annotation
    """
    mt = mt.select_entries("END")

    # Localize entries, so that they can be viewed as an array and scanned
    # over using hl.scan.array_agg
    ht = mt._localize_entries("__entries", "__cols")

    # Compute the position by using hl.scan._prev_nonnull.
    # This was inspired by hl.experimental.densify
    # _prev_non_null is an aggregator that keeps the previous record in memory
    # and updates it with the given value at the row if it's not null (missing)
    # The following code computes the following annotation for each row:
    # 1. Keep a scan of the entries using _prev_nonnull, keeping the start (ht.locus) and end (entry.END) of each ref block  (1.1)
    # 2. For the current row locus, record the start of the block that starts the furthest away,
    # that is the minimum position in the current scan for any block that
    # overlaps the current locus (2.1)
    ht = ht.select(
        last_END_position=hl.or_else(
            hl.min(  # 2. For the current row locus, record the start of the block that starts the furthest away
                hl.scan.array_agg(
                    lambda entry: hl.scan._prev_nonnull(  # 1. Keep a scan of the entries using _prev_nonnull
                        hl.or_missing(
                            hl.is_defined(
                                entry.END
                            ),  # Update the scan whenever a new ref block is encountered
                            hl.tuple(
                                [  # 1.1 keep the start (ht.locus) and end (entry.END) of each ref block
                                    ht.locus,
                                    entry.END,
                                ]
                            ),
                        )
                    ),
                    ht.__entries,
                ).map(
                    lambda x: hl.or_missing(  # 2.1 get the start position of blocks that overlap the current locus
                        (x[1] >= ht.locus.position) & (x[0].contig == ht.locus.contig),
                        x[0].position,
                    )
                )
            ),
            ht.locus.position,
        )
    )
    return ht.select_globals().key_by("locus")


def densify_sites(
    mt: hl.MatrixTable,
    sites_ht: hl.Table,
    last_END_positions_ht: hl.Table,
    semi_join_rows: bool = True,
) -> hl.MatrixTable:
    """
    Create a dense version of the input sparse MT at the sites in `sites_ht` reading the minimal amount of data required.

    Note that only rows that appear both in `mt` and `sites_ht` are returned.

    :param mt: Input sparse MT
    :param sites_ht: Desired sites to densify
    :param last_END_positions_ht: Table storing positions of the furthest ref block (END tag)
    :param semi_join_rows: Whether to filter the MT rows based on semi-join (default, better if sites_ht is large) or based on filter_intervals (better if sites_ht only contains a few sites)
    :return: Dense MT filtered to the sites in `sites_ht`
    """
    logger.info("Computing intervals to densify from sites Table.")
    sites_ht = sites_ht.key_by("locus")
    sites_ht = sites_ht.annotate(
        interval=hl.locus_interval(
            sites_ht.locus.contig,
            last_END_positions_ht[sites_ht.key].last_END_position,
            end=sites_ht.locus.position,
            includes_end=True,
            reference_genome=sites_ht.locus.dtype.reference_genome,
        )
    )
    sites_ht = sites_ht.filter(hl.is_defined(sites_ht.interval))

    if semi_join_rows:
        mt = mt.filter_rows(hl.is_defined(sites_ht.key_by("interval")[mt.locus]))
    else:
        logger.info("Collecting intervals to densify.")
        intervals = sites_ht.interval.collect()

        print(
            "Found {0} intervals, totalling {1} bp in the dense Matrix.".format(
                len(intervals),
                sum(
                    [
                        interval_length(interval)
                        for interval in union_intervals(intervals)
                    ]
                ),
            )
        )

        mt = hl.filter_intervals(mt, intervals)

    mt = hl.experimental.densify(mt)

    return mt.filter_rows(hl.is_defined(sites_ht[mt.locus]))


def _get_info_agg_expr(
    mt: hl.MatrixTable,
    sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["sum_agg_fields"],
    int32_sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["int32_sum_agg_fields"],
    median_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["median_agg_fields"],
    array_sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.ArrayNumericExpression]
    ] = INFO_AGG_FIELDS["array_sum_agg_fields"],
    prefix: str = "",
    treat_fields_as_allele_specific: bool = False,
    retain_cdfs: bool = False,
    cdf_k: int = 200,
) -> Dict[str, hl.expr.Aggregation]:
    """
    Create Aggregators for both site or AS info expression aggregations.

    .. note::

        - If `SB` is specified in array_sum_agg_fields, it will be aggregated as
          `AS_SB_TABLE`, according to GATK standard nomenclature.
        - If `RAW_MQandDP` is specified in array_sum_agg_fields, it will be used for
          the `MQ` calculation and then dropped according to GATK recommendation.
        - If `RAW_MQ` and `MQ_DP` are given, they will be used for the `MQ` calculation
          and then dropped according to GATK recommendation.
        - If the fields to be aggregated (`sum_agg_fields`, `int32_sum_agg_fields`,
          `median_agg_fields`) are passed as list of str, then they should correspond
          to entry fields in `mt` or in mt.gvcf_info`.
        - Priority is given to entry fields in `mt` over those in `mt.gvcf_info` in
          case of a name clash.

    :param mt: Input MT
    :param sum_agg_fields: Fields to aggregate using sum.
    :param int32_sum_agg_fields: Fields to aggregate using sum using int32.
    :param median_agg_fields: Fields to aggregate using (approximate) median.
    :param array_sum_agg_fields: Fields to aggregate using element-wise summing over an
        array.
    :param prefix: Optional prefix for the fields. Used for adding 'AS_' in the AS case.
    :param treat_fields_as_allele_specific: Treat info fields as allele-specific.
        Defaults to False.
    :param retain_cdfs: If True, retains the cumulative distribution functions (CDFs)
        as an annotation for `median_agg_fields`. Keeping the CDFs is useful for
        annotations that require calculating the median across combined datasets at a
        later stage. Default is False.
    :param cdf_k: Parameter controlling the accuracy vs. memory usage tradeoff when
        retaining CDFs. A higher value of `cdf_k` results in a more accurate CDF
        approximation but increases memory usage and computation time. Default is 200.
    :return: Dictionary of expression names and their corresponding aggregation
        Expression.
    """

    def _agg_list_to_dict(
        mt: hl.MatrixTable, fields: List[str]
    ) -> Dict[str, hl.expr.NumericExpression]:
        out_fields = {}
        if "gvcf_info" in mt.entry:
            out_fields = {f: mt.gvcf_info[f] for f in fields if f in mt.gvcf_info}

        out_fields.update({f: mt[f] for f in fields if f in mt.entry})

        # Check that all fields were found.
        missing_fields = [f for f in fields if f not in out_fields]
        if missing_fields:
            raise ValueError(
                "Could not find the following field(s)in the MT entry schema (or nested"
                " under mt.gvcf_info: {}".format(",".join(missing_fields))
            )

        if treat_fields_as_allele_specific:
            # TODO: Change to use hl.vds.local_to_global when fill_value can accept
            #  missing (error in v0.2.119).
            out_fields = {
                f: hl.bind(
                    lambda x: hl.if_else(f == "AS_SB_TABLE", x, x[1:]),
                    hl.range(hl.len(mt.alleles)).map(
                        lambda i: hl.or_missing(
                            mt.LA.contains(i), out_fields[f][mt.LA.index(i)]
                        )
                    ),
                )
                for f in fields
            }

        return out_fields

    # Map str to expressions where needed.
    if isinstance(sum_agg_fields, list):
        sum_agg_fields = _agg_list_to_dict(mt, sum_agg_fields)

    if isinstance(int32_sum_agg_fields, list):
        int32_sum_agg_fields = _agg_list_to_dict(mt, int32_sum_agg_fields)

    if isinstance(median_agg_fields, list):
        median_agg_fields = _agg_list_to_dict(mt, median_agg_fields)

    if isinstance(array_sum_agg_fields, list):
        array_sum_agg_fields = _agg_list_to_dict(mt, array_sum_agg_fields)

    aggs = [
        (median_agg_fields, lambda x: hl.agg.approx_quantiles(x, 0.5)),
        (sum_agg_fields, hl.agg.sum),
        (int32_sum_agg_fields, lambda x: hl.int32(hl.agg.sum(x))),
        (array_sum_agg_fields, hl.agg.array_sum),
    ]

    if retain_cdfs:
        # Note: hl.agg.approx_cdf is a non-deterministic method and cannot be seeded.
        # Results may vary with each rerun.
        cdf_median_agg_fields = {}
        # Store values for each median agg fields in a new dictionary with "_cdf"
        # appended to the annotation name.
        for k, v in median_agg_fields.items():
            cdf_median_agg_fields[f"{k}_cdf"] = v
        # Append the cdf annotations to the aggs list. Set '_raw' to True to return
        # a representation of the internal state of the CDF, which allows for mergining
        # with other CDFs downstream.
        aggs.append(
            (cdf_median_agg_fields, lambda x: hl.agg.approx_cdf(x, k=cdf_k, _raw=True))
        )

    # Create aggregators.
    agg_expr = {}
    for agg_fields, agg_func in aggs:
        for k, expr in agg_fields.items():
            if treat_fields_as_allele_specific:
                # If annotation is of the form 'AS_RAW_*_RankSum' it has a histogram
                # representation where keys give the per-variant rank sum value to one
                # decimal place followed by a comma and the corresponding count for
                # that value, so we want to sum the rank sum value (first element).
                # Rename annotation in the form 'AS_RAW_*_RankSum' to 'AS_*_RankSum'.
                if k.startswith("AS_RAW_") and (
                    k.endswith("RankSum") or k.endswith("RankSum_cdf")
                ):
                    agg_expr[f"{prefix}{k.replace('_RAW', '')}"] = hl.agg.array_agg(
                        lambda x: agg_func(hl.or_missing(hl.is_defined(x), x[0])), expr
                    )
                else:
                    agg_expr[f"{prefix}{k}"] = hl.agg.array_agg(
                        lambda x: agg_func(x), expr
                    )
            else:
                agg_expr[f"{prefix}{k}"] = agg_func(expr)

    if treat_fields_as_allele_specific:
        prefix = "AS_"

    # Handle annotations combinations and casting for specific annotations
    # If RAW_MQandDP is in agg_expr or if both MQ_DP and RAW_MQ are, compute MQ instead
    mq_tuple = None
    if f"{prefix}RAW_MQandDP" in agg_expr:
        logger.info(
            "Computing %sMQ as sqrt(%sRAW_MQandDP[0]/%sRAW_MQandDP[1]). "
            "Note that %sMQ will be set to 0 if %sRAW_MQandDP[1] == 0.",
            *[prefix] * 5,
        )
        mq_tuple = agg_expr.pop(f"{prefix}RAW_MQandDP")
    elif "AS_RAW_MQ" in agg_expr and treat_fields_as_allele_specific:
        logger.info(
            "Computing AS_MQ as sqrt(AS_RAW_MQ[i]/AD[i+1]). "
            "Note that AS_MQ will be set to 0 if AS_RAW_MQ == 0."
        )
        ad_expr = hl.vds.local_to_global(
            mt.LAD, mt.LA, hl.len(mt.alleles), fill_value=0, number="R"
        )
        mq_tuple = hl.zip(agg_expr.pop("AS_RAW_MQ"), hl.agg.array_sum(ad_expr[1:]))
    elif f"{prefix}RAW_MQ" in agg_expr and f"{prefix}MQ_DP" in agg_expr:
        logger.info(
            "Computing %sMQ as sqrt(%sRAW_MQ/%sMQ_DP). "
            "Note that MQ will be set to 0 if %sRAW_MQ == 0.",
            *[prefix] * 4,
        )
        mq_tuple = (agg_expr.pop(f"{prefix}RAW_MQ"), agg_expr.pop(f"{prefix}MQ_DP"))

    if mq_tuple is not None:
        if treat_fields_as_allele_specific:
            agg_expr[f"{prefix}MQ"] = mq_tuple.map(
                lambda x: hl.if_else(x[1] > 0, hl.sqrt(x[0] / x[1]), 0)
            )
        else:
            agg_expr[f"{prefix}MQ"] = hl.if_else(
                mq_tuple[1] > 0, hl.sqrt(mq_tuple[0] / mq_tuple[1]), 0
            )

    # If both VarDP and QUALapprox are present, also compute QD.
    if f"{prefix}VarDP" in agg_expr and f"{prefix}QUALapprox" in agg_expr:
        logger.info(
            "Computing %sQD as %sQUALapprox/%sVarDP. "
            "Note that %sQD will be set to 0 if %sVarDP == 0.",
            *[prefix] * 5,
        )
        var_dp = agg_expr[f"{prefix}VarDP"]
        qual_approx = agg_expr[f"{prefix}QUALapprox"]
        if treat_fields_as_allele_specific:
            agg_expr[f"{prefix}QD"] = hl.map(
                lambda x: hl.if_else(x[1] > 0, x[0] / x[1], 0),
                hl.zip(qual_approx, var_dp),
            )
        else:
            agg_expr[f"{prefix}QD"] = hl.if_else(var_dp > 0, qual_approx / var_dp, 0)

    # SB needs to be cast to int32 for FS down the line.
    if f"{prefix}SB" in agg_expr:
        agg_expr[f"{prefix}SB"] = agg_expr[f"{prefix}SB"].map(lambda x: hl.int32(x))

    # SB needs to be cast to int32 for FS down the line.
    if "AS_SB_TABLE" in agg_expr:
        agg_expr["AS_SB_TABLE"] = agg_expr["AS_SB_TABLE"].map(
            lambda x: x.map(lambda y: hl.int32(y))
        )

    return agg_expr


def get_as_info_expr(
    mt: hl.MatrixTable,
    sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["sum_agg_fields"],
    int32_sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["int32_sum_agg_fields"],
    median_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["median_agg_fields"],
    array_sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.ArrayNumericExpression]
    ] = INFO_AGG_FIELDS["array_sum_agg_fields"],
    alt_alleles_range_array_field: str = "alt_alleles_range_array",
    treat_fields_as_allele_specific: bool = False,
    retain_cdfs: bool = False,
    cdf_k: int = 200,
) -> hl.expr.StructExpression:
    """
    Return an allele-specific annotation Struct containing typical VCF INFO fields from GVCF INFO fields stored in the MT entries.

    .. note::

        - If `SB` is specified in array_sum_agg_fields, it will be aggregated as
          `AS_SB_TABLE`, according to GATK standard nomenclature.
        - If `RAW_MQandDP` is specified in array_sum_agg_fields, it will be used for
          the `MQ` calculation and then dropped according to GATK recommendation.
        - If `RAW_MQ` and `MQ_DP` are given, they will be used for the `MQ` calculation
          and then dropped according to GATK recommendation.
        - If the fields to be aggregate (`sum_agg_fields`, `int32_sum_agg_fields`,
          `median_agg_fields`) are passed as list of str, then they should correspond
          to entry fields in `mt` or in `mt.gvcf_info`.
        - Priority is given to entry fields in `mt` over those in `mt.gvcf_info` in
          case of a name clash.
        - If `treat_fields_as_allele_specific` is False, it's expected that there is a
          single value for each entry field to be aggregated. Then when performing the
          aggregation per global alternate allele, that value is included in the
          aggregation if the global allele is present in the entry's list of local
          alleles. If `treat_fields_as_allele_specific` is True, it's expected that
          each entry field to be aggregated has one value per local allele, and each
          of those is mapped to a global allele for aggregation.

    :param mt: Input Matrix Table
    :param sum_agg_fields: Fields to aggregate using sum.
    :param int32_sum_agg_fields: Fields to aggregate using sum using int32.
    :param median_agg_fields: Fields to aggregate using (approximate) median.
    :param array_sum_agg_fields: Fields to aggregate using array sum.
    :param alt_alleles_range_array_field: Annotation containing an array of the range
        of alternate alleles e.g., `hl.range(1, hl.len(mt.alleles))`
    :param treat_fields_as_allele_specific: Treat info fields as allele-specific.
        Defaults to False.
    :param retain_cdfs: If True, retains the cumulative distribution functions (CDFs)
        as an annotation for `median_agg_fields`. Keeping the CDFs is useful for
        annotations that require calculating the median across combined datasets at a
        later stage. Default is False.
    :param cdf_k: Parameter controlling the accuracy vs. memory usage tradeoff when
        retaining CDFs. A higher value of `cdf_k` results in a more accurate CDF
        approximation but increases memory usage and computation time. Default is 200.
    :return: Expression containing the AS info fields
    """
    if "DP" in list(sum_agg_fields) + list(int32_sum_agg_fields):
        logger.warning(
            "`DP` was included in allele-specific aggregation, however `DP` is"
            " typically not aggregated by allele; `VarDP` is.Note that the resulting"
            " `AS_DP` field will NOT include reference genotypes."
        )

    agg_expr = _get_info_agg_expr(
        mt=mt,
        sum_agg_fields=sum_agg_fields,
        int32_sum_agg_fields=int32_sum_agg_fields,
        median_agg_fields=median_agg_fields,
        array_sum_agg_fields=array_sum_agg_fields,
        prefix="" if treat_fields_as_allele_specific else "AS_",
        treat_fields_as_allele_specific=treat_fields_as_allele_specific,
        retain_cdfs=retain_cdfs,
        cdf_k=cdf_k,
    )

    if alt_alleles_range_array_field not in mt.row or mt[
        alt_alleles_range_array_field
    ].dtype != hl.dtype("array<int32>"):
        msg = (
            f"'get_as_info_expr' expected a row field '{alt_alleles_range_array_field}'"
            " of type array<int32>"
        )
        logger.error(msg)
        raise ValueError(msg)

    if not treat_fields_as_allele_specific:
        # Modify aggregations to aggregate per allele
        agg_expr = {
            f: hl.agg.array_agg(
                lambda ai: hl.agg.filter(mt.LA.contains(ai), expr),
                mt[alt_alleles_range_array_field],
            )
            for f, expr in agg_expr.items()
        }

    # Run aggregations
    info = hl.struct(**agg_expr)

    # Add FS and SOR if SB is present.
    if "AS_SB_TABLE" in info or "AS_SB" in info:
        drop = []
        # Rename AS_SB to AS_SB_TABLE if present and add SB Ax2 aggregation logic.
        if "AS_SB" in agg_expr:
            if "AS_SB_TABLE" in agg_expr:
                logger.warning(
                    "Both `AS_SB` and `AS_SB_TABLE` were specified for aggregation."
                    " `AS_SB` will be used for aggregation."
                )
            as_sb_table = hl.array(
                [
                    info.AS_SB.filter(lambda x: hl.is_defined(x)).fold(
                        lambda i, j: i[:2] + j[:2], [0, 0]
                    )  # ref
                ]
            ).extend(
                info.AS_SB.map(lambda x: x[2:])  # each alt
            )
            drop = ["AS_SB"]
        else:
            as_sb_table = info.AS_SB_TABLE
        info = info.annotate(
            AS_SB_TABLE=as_sb_table,
            AS_FS=hl.range(1, hl.len(mt.alleles)).map(
                lambda i: fs_from_sb(as_sb_table[0].extend(as_sb_table[i]))
            ),
            AS_SOR=hl.range(1, hl.len(mt.alleles)).map(
                lambda i: sor_from_sb(as_sb_table[0].extend(as_sb_table[i]))
            ),
        ).drop(*drop)

    return info


def get_site_info_expr(
    mt: hl.MatrixTable,
    sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["sum_agg_fields"],
    int32_sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["int32_sum_agg_fields"],
    median_agg_fields: Union[
        List[str], Dict[str, hl.expr.NumericExpression]
    ] = INFO_AGG_FIELDS["median_agg_fields"],
    array_sum_agg_fields: Union[
        List[str], Dict[str, hl.expr.ArrayNumericExpression]
    ] = INFO_AGG_FIELDS["array_sum_agg_fields"],
    retain_cdfs: bool = False,
    cdf_k: int = 200,
) -> hl.expr.StructExpression:
    """
    Create a site-level annotation Struct aggregating typical VCF INFO fields from GVCF INFO fields stored in the MT entries.

    .. note::

        - If `RAW_MQandDP` is specified in array_sum_agg_fields, it will be used for
          the `MQ` calculation and then dropped according to GATK recommendation.
        - If `RAW_MQ` and `MQ_DP` are given, they will be used for the `MQ` calculation
          and then dropped according to GATK recommendation.
        - If the fields to be aggregate (`sum_agg_fields`, `int32_sum_agg_fields`,
          `median_agg_fields`) are passed as list of str, then they should correspond
          to entry fields in `mt` or in `mt.gvcf_info`.
        - Priority is given to entry fields in `mt` over those in `mt.gvcf_info` in
          case of a name clash.

    :param mt: Input Matrix Table
    :param sum_agg_fields: Fields to aggregate using sum.
    :param int32_sum_agg_fields: Fields to aggregate using sum using int32.
    :param median_agg_fields: Fields to aggregate using (approximate) median.
    :param retain_cdfs: If True, retains the cumulative distribution functions (CDFs)
        as an annotation for `median_agg_fields`. Keeping the CDFs is useful for
        annotations that require calculating the median across combined datasets at a
        later stage. Default is False.
    :param cdf_k: Parameter controlling the accuracy vs. memory usage tradeoff when
        retaining CDFs. A higher value of `cdf_k` results in a more accurate CDF
        approximation but increases memory usage and computation time. Default is 200.
    :return: Expression containing the site-level info fields
    """
    if "DP" in list(sum_agg_fields) + list(int32_sum_agg_fields):
        logger.warning(
            "`DP` was included in site-level aggregation. This requires a densifying"
            " prior to running get_site_info_expr"
        )

    agg_expr = _get_info_agg_expr(
        mt=mt,
        sum_agg_fields=sum_agg_fields,
        int32_sum_agg_fields=int32_sum_agg_fields,
        median_agg_fields=median_agg_fields,
        array_sum_agg_fields=array_sum_agg_fields,
        retain_cdfs=retain_cdfs,
        cdf_k=cdf_k,
    )

    # Add FS and SOR if SB is present
    # This is done outside _get_info_agg_expr as the behavior is different
    # in site vs allele-specific versions
    if "SB" in agg_expr:
        agg_expr["FS"] = fs_from_sb(agg_expr["SB"])
        agg_expr["SOR"] = sor_from_sb(agg_expr["SB"])

    # Run aggregator on non-ref genotypes
    info = hl.agg.filter(
        mt.LGT.is_non_ref(),
        hl.struct(**{k: v for k, v in agg_expr.items() if k != "DP"}),
    )

    # Add DP, computed over both ref and non-ref genotypes, if present
    if "DP" in agg_expr:
        info = info.annotate(DP=agg_expr["DP"])

    return info


def default_compute_info(
    mt: hl.MatrixTable,
    site_annotations: bool = False,
    as_annotations: bool = False,
    # Set to True by default to prevent a breaking change.
    quasi_as_annotations: bool = True,
    n_partitions: Optional[int] = 5000,
    lowqual_indel_phred_het_prior: int = 40,
    ac_filter_groups: Optional[Dict[str, hl.Expression]] = None,
    retain_cdfs: bool = False,
    cdf_k: int = 200,
) -> hl.Table:
    """
    Compute a HT with the typical GATK allele-specific (AS) info fields as well as ACs and lowqual fields.

    .. note::

        - This table doesn't split multi-allelic sites.
        - At least one of `site_annotations`, `as_annotations` or `quasi_as_annotations`
          must be True.

    :param mt: Input MatrixTable. Note that this table should be filtered to nonref sites.
    :param site_annotations: Whether to generate site level info fields. Default is False.
    :param as_annotations: Whether to generate allele-specific info fields using
        allele-specific annotations in gvcf_info. Default is False.
    :param quasi_as_annotations: Whether to generate allele-specific info fields using
        non-allele-specific annotations in gvcf_info, but performing per allele
        aggregations. This method can be used in cases where genotype data doesn't
        contain allele-specific annotations to approximate allele-specific annotations.
        Default is True.
    :param n_partitions: Optional number of desired partitions for output Table. If
        specified, naive_coalesce is performed. Default is 5000.
    :param lowqual_indel_phred_het_prior: Phred-scaled prior for a het genotype at a
        site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to
        be more consistent with the filtering used by Broad's Data Sciences Platform
        for VQSR.
    :param ac_filter_groups: Optional dictionary of sample filter expressions to compute
        additional groupings of ACs. Default is None.
    :param retain_cdfs: If True, retains the cumulative distribution functions (CDFs)
        as an annotation for `median_agg_fields`. Keeping the CDFs is useful for
        annotations that require calculating the median across combined datasets at a
        later stage. Default is False.
    :param cdf_k: Parameter controlling the accuracy vs. memory usage tradeoff when
        retaining CDFs. A higher value of `cdf_k` results in a more accurate CDF
        approximation but increases memory usage and computation time. Default is 200.
    :return: Table with info fields
    :rtype: Table
    """
    if not site_annotations and not as_annotations and not quasi_as_annotations:
        raise ValueError(
            "At least one of `site_annotations`, `as_annotations`, or "
            "`quasi_as_annotations` must be True!"
        )

    # Add a temporary annotation for allele count groupings.
    ac_filter_groups = {"": True, **(ac_filter_groups or {})}
    mt = mt.annotate_cols(_ac_filter_groups=ac_filter_groups)

    # Move gvcf info entries out from nested struct.
    mt = mt.transmute_entries(**mt.gvcf_info)

    # Adding alt_alleles_range_array as a required annotation for
    # get_as_info_expr to reduce memory usage.
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    info_expr = None
    quasi_info_expr = None

    # Compute quasi-AS info expr.
    if quasi_as_annotations:
        info_expr = get_as_info_expr(mt, retain_cdfs=retain_cdfs, cdf_k=cdf_k)

    # Compute AS info expr using gvcf_info allele specific annotations.
    if as_annotations:
        if info_expr is not None:
            quasi_info_expr = info_expr
        info_expr = get_as_info_expr(
            mt,
            **AS_INFO_AGG_FIELDS,
            treat_fields_as_allele_specific=True,
            retain_cdfs=retain_cdfs,
            cdf_k=cdf_k,
        )

    if info_expr is not None:
        # Add allele specific pab_max
        info_expr = info_expr.annotate(
            AS_pab_max=pab_max_expr(mt.LGT, mt.LAD, mt.LA, hl.len(mt.alleles))
        )

    if site_annotations:
        site_expr = get_site_info_expr(mt, retain_cdfs=retain_cdfs, cdf_k=cdf_k)
        if info_expr is None:
            info_expr = site_expr
        else:
            info_expr = info_expr.annotate(**site_expr)

    # Add 'AC' and 'AC_raw' for each allele count filter group requested.
    # First compute ACs for each non-ref allele, grouped by adj.
    grp_ac_expr = {
        f: hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai) & mt._ac_filter_groups[f],
                hl.agg.group_by(
                    get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
                    hl.agg.sum(
                        mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[
                            mt.LA.index(ai)
                        ]
                    ),
                ),
            ),
            mt.alt_alleles_range_array,
        )
        for f in ac_filter_groups
    }

    # Then, for each non-ref allele, compute
    # 'AC' as the adj group
    # 'AC_raw' as the sum of adj and non-adj groups
    info_expr = info_expr.annotate(
        **{
            f"AC{'_' + f if f else f}_raw": grp.map(
                lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))
            )
            for f, grp in grp_ac_expr.items()
        },
        **{
            f"AC{'_' + f if f else f}": grp.map(lambda i: hl.int32(i.get(True, 0)))
            for f, grp in grp_ac_expr.items()
        },
    )

    ann_expr = {"info": info_expr}
    if quasi_info_expr is not None:
        ann_expr["quasi_info"] = quasi_info_expr

    info_ht = mt.select_rows(**ann_expr).rows()

    # Add AS lowqual flag
    info_ht = info_ht.annotate(
        AS_lowqual=get_lowqual_expr(
            info_ht.alleles,
            info_ht.info.AS_QUALapprox,
            indel_phred_het_prior=lowqual_indel_phred_het_prior,
        )
    )

    if site_annotations:
        # Add lowqual flag
        info_ht = info_ht.annotate(
            lowqual=get_lowqual_expr(
                info_ht.alleles,
                info_ht.info.QUALapprox,
                indel_phred_het_prior=lowqual_indel_phred_het_prior,
            )
        )

    if n_partitions is not None:
        info_ht = info_ht.naive_coalesce(n_partitions)

    info_ht = info_ht.annotate_globals(retain_cdfs=retain_cdfs)
    if retain_cdfs:
        info_ht = info_ht.annotate_globals(cdf_k=cdf_k)

    return info_ht


def split_info_annotation(
    info_expr: hl.expr.StructExpression, a_index: hl.expr.Int32Expression
) -> hl.expr.StructExpression:
    """
    Split multi-allelic allele-specific info fields.

    :param info_expr: Field containing info struct.
    :param a_index: Allele index. Output by hl.split_multi or hl.split_multi_hts.
    :return: Info struct with split annotations.
    """
    # Index AS annotations
    info_expr = info_expr.annotate(
        **{
            f: info_expr[f][a_index - 1]
            for f in info_expr
            if f.startswith("AC") or (f.startswith("AS_") and not f == "AS_SB_TABLE")
        }
    )
    if "AS_SB_TABLE" in info_expr:
        info_expr = info_expr.annotate(
            AS_SB_TABLE=info_expr.AS_SB_TABLE[0].extend(info_expr.AS_SB_TABLE[a_index])
        )

    return info_expr


def split_lowqual_annotation(
    lowqual_expr: hl.expr.ArrayExpression, a_index: hl.expr.Int32Expression
) -> hl.expr.BooleanExpression:
    """
    Split multi-allelic low QUAL annotation.

    :param lowqual_expr: Field containing low QUAL annotation.
    :param a_index: Allele index. Output by hl.split_multi or hl.split_multi_hts.
    :return: Low QUAL expression for particular allele.
    """
    return lowqual_expr[a_index - 1]


def impute_sex_ploidy(
    mt: hl.MatrixTable,
    excluded_calling_intervals: Optional[hl.Table] = None,
    included_calling_intervals: Optional[hl.Table] = None,
    normalization_contig: str = "chr20",
    chr_x: Optional[str] = None,
    chr_y: Optional[str] = None,
    use_only_variants: bool = False,
) -> hl.Table:
    """
    Impute sex ploidy from a sparse MatrixTable.

    Sex ploidy is imputed by normalizing the coverage of chromosomes X and Y using the coverage of an autosomal
    chromosome (by default chr20).

    Coverage is computed using the median block coverage (summed over the block size) and the non-ref coverage at
    non-ref genotypes unless the `use_only_variants` argument is set to True and then it will use the mean coverage
    defined by only the variants.

    :param mt: Input sparse Matrix Table
    :param excluded_calling_intervals: Optional table of intervals to exclude from the computation. Used only when
        determining contig size (not used when computing chromosome depth) when `use_only_variants` is False.
    :param included_calling_intervals: Optional table of intervals to use in the computation. Used only when
        determining contig size (not used when computing chromosome depth) when `use_only_variants` is False.
    :param normalization_contig: Which chromosome to normalize by
    :param chr_x: Optional X Chromosome contig name (by default uses the X contig in the reference)
    :param chr_y: Optional Y Chromosome contig name (by default uses the Y contig in the reference)
    :param use_only_variants: Whether to use depth of variant data within calling intervals instead of reference data.
        Default will only use reference data.

    :return: Table with mean coverage over chromosomes 20, X and Y and sex chromosomes ploidy based on normalized coverage.
    """
    ref = get_reference_genome(mt.locus, add_sequence=True)
    if chr_x is None:
        if len(ref.x_contigs) != 1:
            raise NotImplementedError(
                "Found {0} X chromosome contigs ({1}) in Genome reference."
                " sparse_impute_sex_ploidy currently only supports a single X"
                " chromosome contig. Please use the `chr_x` argument to  specify which"
                " X chromosome contig to use ".format(
                    len(ref.x_contigs), ",".join(ref.x_contigs)
                )
            )
        chr_x = ref.x_contigs[0]
    if chr_y is None:
        if len(ref.y_contigs) != 1:
            raise NotImplementedError(
                "Found {0} Y chromosome contigs ({1}) in Genome reference."
                " sparse_impute_sex_ploidy currently only supports a single Y"
                " chromosome contig. Please use the `chr_y` argument to  specify which"
                " Y chromosome contig to use ".format(
                    len(ref.y_contigs), ",".join(ref.y_contigs)
                )
            )
        chr_y = ref.y_contigs[0]

    def get_contig_size(contig: str) -> int:
        """
        Compute the size of the specified `contig` using the median block coverage (summed over the block size).

        The size of the contig will be determined using only non par regions if the contig is an X or Y reference contig
        and using the intervals specified by `included_calling_intervals` and excluding intervals specified by
        `excluded_calling_intervals` if either is defined in the outer function.

        :param contig: Contig to compute the size of
        :return: Integer of the contig size
        """
        logger.info("Working on %s", contig)
        contig_ht = hl.utils.range_table(
            ref.contig_length(contig),
            n_partitions=int(ref.contig_length(contig) / 500_000),
        )
        contig_ht = contig_ht.annotate(
            locus=hl.locus(contig=contig, pos=contig_ht.idx + 1, reference_genome=ref)
        )
        contig_ht = contig_ht.filter(contig_ht.locus.sequence_context().lower() != "n")

        if contig in ref.x_contigs:
            contig_ht = contig_ht.filter(contig_ht.locus.in_x_nonpar())
        if contig in ref.y_contigs:
            contig_ht = contig_ht.filter(contig_ht.locus.in_y_nonpar())

        contig_ht = contig_ht.key_by("locus")
        if included_calling_intervals is not None:
            contig_ht = contig_ht.filter(
                hl.is_defined(included_calling_intervals[contig_ht.key])
            )
        if excluded_calling_intervals is not None:
            contig_ht = contig_ht.filter(
                hl.is_missing(excluded_calling_intervals[contig_ht.key])
            )
        contig_size = contig_ht.count()
        logger.info("Contig %s has %d bases for coverage.", contig, contig_size)
        return contig_size

    def get_chr_dp_ann(chrom: str) -> hl.Table:
        """
        Compute the mean depth of the specified chromosome.

        The total depth will be determined using the sum DP of either reference and variant data or only variant data
        depending on the value of `use_only_variants` in the outer function.

        If `use_only_variants` is set to False then this value is computed using the median block coverage (summed over
        the block size). If `use_only_variants` is set to True, this value is computed using the sum of DP for  all
        variants divided by the total number of variants.

        The depth calculations will be determined using only non par regions if the contig is an X or Y reference contig
        and using the intervals specified by `included_calling_intervals` and excluding intervals specified by
        `excluded_calling_intervals` if either is defined in the outer function (when `use_only_variants` is not
        set this only applies to the contig size estimate and is not used when computing chromosome depth).

        :param chrom: Chromosome to compute the mean depth of
        :return: Table of a per sample mean depth of `chrom`
        """
        contig_size = get_contig_size(chrom)
        chr_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chrom)])

        if chrom in ref.x_contigs:
            chr_mt = chr_mt.filter_rows(chr_mt.locus.in_x_nonpar())
        if chrom in ref.y_contigs:
            chr_mt = chr_mt.filter_rows(chr_mt.locus.in_y_nonpar())

        if use_only_variants:
            if included_calling_intervals is not None:
                chr_mt = chr_mt.filter_rows(
                    hl.is_defined(included_calling_intervals[chr_mt.locus])
                )
            if excluded_calling_intervals is not None:
                chr_mt = chr_mt.filter_rows(
                    hl.is_missing(excluded_calling_intervals[chr_mt.locus])
                )
            return chr_mt.select_cols(
                **{
                    f"{chrom}_mean_dp": hl.agg.filter(
                        chr_mt.LGT.is_non_ref(),
                        hl.agg.sum(chr_mt.DP),
                    )
                    / hl.agg.filter(chr_mt.LGT.is_non_ref(), hl.agg.count())
                }
            ).cols()
        else:
            return chr_mt.select_cols(
                **{
                    f"{chrom}_mean_dp": (
                        hl.agg.sum(
                            hl.if_else(
                                chr_mt.LGT.is_hom_ref(),
                                chr_mt.DP * (1 + chr_mt.END - chr_mt.locus.position),
                                chr_mt.DP,
                            )
                        )
                        / contig_size
                    )
                }
            ).cols()

    normalization_chrom_dp = get_chr_dp_ann(normalization_contig)
    chrX_dp = get_chr_dp_ann(chr_x)
    chrY_dp = get_chr_dp_ann(chr_y)

    ht = normalization_chrom_dp.annotate(
        **chrX_dp[normalization_chrom_dp.key],
        **chrY_dp[normalization_chrom_dp.key],
    )

    return ht.annotate(
        **{
            f"{chr_x}_ploidy": ht[f"{chr_x}_mean_dp"]
            / (ht[f"{normalization_contig}_mean_dp"] / 2),
            f"{chr_y}_ploidy": ht[f"{chr_y}_mean_dp"]
            / (ht[f"{normalization_contig}_mean_dp"] / 2),
        }
    )


def densify_all_reference_sites(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    row_key_fields: Union[Tuple[str], List[str], Set[str]] = ("locus",),
    entry_keep_fields: Union[Tuple[str], List[str], Set[str]] = ("GT",),
) -> hl.MatrixTable:
    """
    Densify a VariantDataset or Sparse MatrixTable at all sites in a reference Table.

    :param mtds: Input sparse MatrixTable or VariantDataset.
    :param reference_ht: Table of reference sites.
    :param interval_ht: Optional Table of intervals to filter to.
    :param row_key_fields: Fields to use as row key. Defaults to locus.
    :param entry_keep_fields: Fields to keep in entries before performing the
        densification. Defaults to GT.
    :return: Densified MatrixTable.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)

    if interval_ht is not None and not is_vds:
        raise NotImplementedError(
            "Filtering to an interval list for a sparse Matrix Table is currently"
            " not supported."
        )

    # Filter datasets to interval list.
    if interval_ht is not None:
        reference_ht = reference_ht.filter(
            hl.is_defined(interval_ht[reference_ht.locus])
        )
        mtds = hl.vds.filter_intervals(
            vds=mtds, intervals=interval_ht, split_reference_blocks=False
        )

    entry_keep_fields = set(entry_keep_fields)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds
        entry_keep_fields.add("END")

    # Get the total number of samples.
    n_samples = mt.count_cols()
    mt_col_key_fields = list(mt.col_key)
    mt_row_key_fields = list(mt.row_key)
    ht = mt.select_entries(*entry_keep_fields).select_cols()

    # Localize entries and perform an outer join with the reference HT.
    ht = ht._localize_entries("__entries", "__cols")
    ht = ht.key_by(*row_key_fields)
    ht = ht.join(reference_ht.key_by(*row_key_fields).select(_in_ref=True), how="outer")
    ht = ht.key_by(*mt_row_key_fields)

    # Fill in missing entries with missing values for each entry field.
    ht = ht.annotate(
        __entries=hl.or_else(
            ht.__entries,
            hl.range(n_samples).map(
                lambda x: hl.missing(ht.__entries.dtype.element_type)
            ),
        )
    )

    # Unlocalize entries to turn the HT back to a MT.
    mt = ht._unlocalize_entries("__entries", "__cols", mt_col_key_fields)

    # Densify VDS/sparse MT at all sites.
    if is_vds:
        mt = hl.vds.to_dense_mt(
            hl.vds.VariantDataset(mtds.reference_data.select_cols().select_rows(), mt)
        )
    else:
        mt = hl.experimental.densify(mt)

    # Remove rows where the reference is missing.
    mt = mt.filter_rows(mt._in_ref)

    # Unfilter entries so that entries with no ref block overlap aren't null.
    mt = mt.unfilter_entries()

    # Rekey by requested row key field and drop unused keys.
    mt = mt.key_rows_by(*row_key_fields)
    mt = mt.drop(*[k for k in mt_row_key_fields if k not in row_key_fields])

    return mt


def compute_stats_per_ref_site(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    entry_agg_funcs: Dict[str, Tuple[Callable, Callable]],
    row_key_fields: Union[Tuple[str], List[str]] = ("locus",),
    interval_ht: Optional[hl.Table] = None,
    entry_keep_fields: Union[Tuple[str], List[str], Set[str]] = None,
    row_keep_fields: Union[Tuple[str], List[str], Set[str]] = None,
    entry_agg_group_membership: Optional[Dict[str, List[dict[str, str]]]] = None,
    strata_expr: Optional[List[Dict[str, hl.expr.StringExpression]]] = None,
    group_membership_ht: Optional[hl.Table] = None,
    sex_karyotype_field: Optional[str] = None,
    reduce_to_minimal_groups: bool = False,
    non_summable_strata: Optional[Set[str]] = None,
    reducible_aggs: Optional[Set[str]] = None,
) -> hl.Table:
    """
    Compute stats per site in a reference Table.

    .. rubric:: The `reduce_to_minimal_groups` parameter

    When True, the per-strata aggregation runs only on the "leaf"
    stratification groups (those that cannot be derived by summing other
    groups). The remaining groups are reconstructed by element-wise
    summation of leaves as a cheap post-processing step, so the returned
    `strata_meta` and per-site arrays are identical to the
    `reduce_to_minimal_groups=False` output. This is purely a cost
    optimization for large stratifications.

    The reduction assumes that the annotations subject to leaf expansion
    are summable (integers or struct-of-integers). By default every
    annotation in `entry_agg_funcs` is treated as reducible. Pass
    `reducible_aggs={...}` to opt only specific annotations into the
    expansion; the others must be narrowed via
    `entry_agg_group_membership` (otherwise their per-row arrays would
    have leaf shape, inconsistent with the full `strata_meta` global
    that the function emits). When leaf reduction is in effect, every
    entry in `entry_agg_funcs` must therefore appear in either
    `reducible_aggs` or as a key of `entry_agg_group_membership` —
    `compute_stats_per_ref_site` raises `ValueError` early in that case
    so no aggregation cost is wasted on a call whose output would be
    shape-inconsistent. This split lets a single call mix summable
    aggregations (e.g., AN) with non-summable ones (e.g., a coverage
    mean or a qual histogram restricted to a single global group via
    `entry_agg_group_membership`).

    Reduction is supported on both supplied paths:

        - When `strata_expr` is provided, the leaf reduction is performed
          inside the in-function call to
          `generate_freq_group_membership_array`.
        - When a pre-built `group_membership_ht` is supplied, this
          function honors any reduction already performed on it
          (detected by the `freq_reduced=True` global on
          `group_membership_ht`). Setting
          `reduce_to_minimal_groups=True` while passing a non-reduced
          `group_membership_ht` is a no-op.

    :param mtds: Input sparse Matrix Table or VariantDataset.
    :param reference_ht: Table of reference sites.
    :param entry_agg_funcs: Dict of entry aggregation functions to perform on the
        VariantDataset/MatrixTable. The keys of the dict are the names of the
        annotations and the values are tuples of functions. The first function is used
        to transform the `mt` entries in some way, and the second function is used to
        aggregate the output from the first function.
    :param row_key_fields: Fields to use as row key. Defaults to locus.
    :param interval_ht: Optional table of intervals to filter to.
    :param entry_keep_fields: Fields to keep in entries before performing the
        densification in `densify_all_reference_sites`. Should include any fields
        needed for the functions in `entry_agg_funcs`. By default, only GT or LGT is
        kept.
    :param row_keep_fields: Fields to keep in rows after performing the stats
        aggregation. By default, only the row key fields are kept.
    :param entry_agg_group_membership: Optional dict indicating the subset of group
        strata in 'freq_meta' to use the entry aggregation functions on. The keys of
        the dict can be any of the keys in `entry_agg_funcs` and the values are lists
        of dicts. Each dict in the list contains the strata in 'freq_meta' to use for
        the corresponding entry aggregation function. If provided, 'freq_meta' must be
        present in `group_membership_ht` and represent the same strata as those in
        'group_membership'. If not provided, all entries of the 'group_membership'
        annotation will have the entry aggregation functions applied to them.
    :param strata_expr: Optional list of dicts of expressions to stratify by.
    :param group_membership_ht: Optional Table of group membership annotations.
    :param sex_karyotype_field: Optional field to use to adjust genotypes for sex
        karyotype before stats aggregation. If provided, the field must be present in
        the columns of `mtds` (variant_data MT if `mtds` is a VDS) and use "XX" and
        "XY" as values. If not provided, no sex karyotype adjustment is performed.
        Default is None.
    :param reduce_to_minimal_groups: Whether to compute stats only on the
        leaf set of stratification groups and reconstruct the rest by
        element-wise summation. See the rubric above. Default is False.
    :param non_summable_strata: Strata names that should never be summed
        across their values when `reduce_to_minimal_groups` is True.
        Default is None, which resolves to `{"downsampling"}`.
    :param reducible_aggs: Optional set of annotation names from
        `entry_agg_funcs` that are summable and should be expanded from
        leaf shape to full shape via element-wise summation when leaf
        reduction is in effect. Defaults to all of `entry_agg_funcs`.
        Must be disjoint from the keys of `entry_agg_group_membership`.
        Ignored when no leaf reduction is in effect.
    :return: Table of stats per site.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds

    if sex_karyotype_field is not None and sex_karyotype_field not in mt.col:
        raise ValueError(
            f"The supplied 'sex_karyotype_field', {sex_karyotype_field}, is not present"
            " in the columns of the input!"
        )

    if group_membership_ht is not None and strata_expr is not None:
        raise ValueError(
            "Only one of 'group_membership_ht' or 'strata_expr' can be specified."
        )

    g = {} if group_membership_ht is None else group_membership_ht.globals
    if entry_agg_group_membership is not None and "freq_meta" not in g:
        raise ValueError(
            "The 'freq_meta' annotation must be present in 'group_membership_ht' if "
            "'entry_agg_group_membership' is specified."
        )

    # Determine reduction status from inputs so we can validate the
    # shape of `entry_agg_funcs` before any aggregation runs. Both
    # signals are inspectable without densify or aggregation:
    #   - When `group_membership_ht` is supplied, reduction is recorded
    #     by `generate_freq_group_membership_array` as the
    #     `freq_reduced` global on that table.
    #   - When `strata_expr` is supplied, the table is built later in
    #     this function with `reduce_to_minimal_groups` forwarded into
    #     `generate_freq_group_membership_array`, so the parameter
    #     directly determines reduction status. The pre-built case's
    #     "no-op when flag set but table not reduced" carve-out applies
    #     only to the supplied-table path.
    if group_membership_ht is not None:
        is_reduced = "freq_reduced" in group_membership_ht.globals and hl.eval(
            group_membership_ht.freq_reduced
        )
    else:
        is_reduced = reduce_to_minimal_groups

    if reducible_aggs is not None:
        if not reducible_aggs:
            raise ValueError(
                "`reducible_aggs` is empty; pass `None` to default to all of "
                "`entry_agg_funcs` or list at least one annotation."
            )
        unknown = set(reducible_aggs) - set(entry_agg_funcs)
        if unknown:
            raise ValueError(
                "`reducible_aggs` contains names not in `entry_agg_funcs`: "
                f"{sorted(unknown)}."
            )
        if entry_agg_group_membership is not None:
            overlap = set(reducible_aggs) & set(entry_agg_group_membership)
            if overlap:
                raise ValueError(
                    "`reducible_aggs` and `entry_agg_group_membership` must"
                    f" be disjoint; got overlap: {sorted(overlap)}. A"
                    " reducible annotation is computed on the leaf set and"
                    " reconstructed by summation, while"
                    " `entry_agg_group_membership` narrows an annotation to"
                    " a specific subset; the two are mutually exclusive for"
                    " the same annotation."
                )

    # When leaf reduction is in effect every annotation in
    # `entry_agg_funcs` must be accounted for: either it is reducible
    # (gets expanded leaf -> full post-aggregation) or it is narrowed
    # via `entry_agg_group_membership` to a specific subset of strata.
    # Anything in neither bucket would silently emit a leaf-shape
    # array against a full-shape `strata_meta` global, which is a
    # latent shape mismatch in downstream consumers; raise so the
    # caller fixes the call before paying the densify+aggregation cost.
    if is_reduced:
        reducible = (
            set(reducible_aggs) if reducible_aggs is not None else set(entry_agg_funcs)
        )
        grouped = (
            set(entry_agg_group_membership) if entry_agg_group_membership else set()
        )
        unaccounted = set(entry_agg_funcs) - (reducible | grouped)
        if unaccounted:
            raise ValueError(
                "Leaf reduction is in effect but the following entries of"
                f" `entry_agg_funcs` have no shape designation:"
                f" {sorted(unaccounted)}. Each annotation must be either"
                " (a) listed in `reducible_aggs` if its values are summable"
                " across strata, or (b) a key of"
                " `entry_agg_group_membership` to restrict it to a specific"
                " subset of groups; otherwise its per-row array would have"
                " leaf shape while `strata_meta` is emitted at full shape."
            )

    # Determine if the adj annotation is needed. It is only needed if "adj_groups" is
    # in the globals of the group_membership_ht and any entry is True, or "freq_meta"
    # is in the globals of the group_membership_ht and any entry has "group" == "adj".
    adj = hl.eval(
        hl.any(g.get("adj_groups", hl.empty_array("bool")))
        | hl.any(
            g.get("freq_meta", hl.empty_array("dict<str, str>")).map(
                lambda x: x.get("group", "NA") == "adj"
            )
        )
    )

    # Determine the entry fields on mt that should be densified.
    # "GT" or "LGT" is required for the genotype.
    # If the adj annotation is needed then "adj" must be present on mt, or AD/LAD, DP,
    # and GQ must be present.
    en = set(mt.entry)
    gt_field = en & {"GT"} or en & {"LGT"}
    ad_field = en & {"AD"} or en & {"LAD"}
    adj_fields = en & {"adj"} or ({"DP", "GQ"} | ad_field) if adj else set([])

    if not gt_field:
        raise ValueError("No genotype field found in entry fields!")

    if adj and not adj_fields.issubset(en):
        raise ValueError(
            "No 'adj' found in entry fields, and one of AD/LAD, DP, and GQ is missing "
            "so adj can't be computed!"
        )

    entry_keep_fields = set(entry_keep_fields or set([])) | gt_field | adj_fields

    # Write the sex karyotype field out to a temp HT so we can annotate the field back
    # onto the MT after 'densify_all_reference_sites' removes all column annotations.
    if sex_karyotype_field is not None:
        sex_karyotype_ht = (
            mt.cols()
            .select(sex_karyotype_field)
            .checkpoint(hl.utils.new_temp_file("sex_karyotype_ht", "ht"))
        )
    else:
        sex_karyotype_ht = None

    # Initialize no_strata and default strata_expr if neither group_membership_ht nor
    # strata_expr is provided.
    no_strata = group_membership_ht is None and strata_expr is None
    if no_strata:
        strata_expr = {}

    if group_membership_ht is None:
        logger.warning(
            "'group_membership_ht' is not specified, no stats are adj filtered."
        )

        # Annotate the MT cols with each of the expressions in strata_expr and redefine
        # strata_expr based on the column HT with added annotations.
        ht = mt.annotate_cols(
            **{k: v for d in strata_expr for k, v in d.items()}
        ).cols()
        strata_expr = [{k: ht[k] for k in d} for d in strata_expr]

        # Use 'generate_freq_group_membership_array' to create a group_membership Table
        # that gives stratification group membership info based on 'strata_expr'.
        # The per-ref-site stats path does no genotype-level filtering, so
        # label every freq_meta entry as "raw" rather than the default "adj".
        # `no_raw_group=True` skips the inserted raw-only entry; the
        # remaining entries cover all samples.
        group_membership_ht = generate_freq_group_membership_array(
            ht,
            strata_expr,
            no_raw_group=True,
            reduce_to_minimal_groups=reduce_to_minimal_groups,
            non_summable_strata=non_summable_strata,
            group_label="raw",
        )

    if is_vds:
        rmt = mtds.reference_data
        mtds = hl.vds.VariantDataset(
            rmt.select_entries(
                *((set(entry_keep_fields) & set(rmt.entry)) | {"END", "LEN"})
            ),
            mtds.variant_data,
        )

    mt = densify_all_reference_sites(
        mtds,
        reference_ht,
        interval_ht,
        row_key_fields,
        entry_keep_fields=entry_keep_fields,
    )

    if sex_karyotype_ht is not None:
        logger.info("Adjusting genotype ploidy based on sex karyotype.")
        gt_field = gt_field.pop()
        mt = mt.annotate_cols(
            sex_karyotype=sex_karyotype_ht[mt.col_key][sex_karyotype_field]
        )
        mt = mt.annotate_entries(
            **{
                gt_field: adjusted_sex_ploidy_expr(
                    mt.locus, mt[gt_field], mt.sex_karyotype
                )
            }
        )

    # Annotate with adj if needed.
    if adj and "adj" not in mt.entry:
        logger.info("Annotating the MT with adj.")
        mt = annotate_adj(mt)

    ht = agg_by_strata(
        mt,
        entry_agg_funcs,
        group_membership_ht=group_membership_ht,
        select_fields=row_keep_fields,
        entry_agg_group_membership=entry_agg_group_membership,
    )
    ht = ht.select_globals().checkpoint(hl.utils.new_temp_file("agg_stats", "ht"))

    group_globals = group_membership_ht.index_globals()

    # `is_reduced` was determined from the inputs above. If True, expand
    # each reducible per-strata aggregation back to the full set of
    # groups by element-wise summation of the leaf positions.
    if is_reduced:
        r = _read_reduction_globals(group_globals)
        anns_to_expand = (
            set(entry_agg_funcs) if reducible_aggs is None else set(reducible_aggs)
        )
        ht = ht.annotate(
            **{
                ann: expand_strata_array_from_leaves(
                    ht[ann], r["leaf_indices"], r["decomposition"], r["n_full"]
                )
                for ann in anns_to_expand
            }
        )

    global_expr = {}
    if no_strata:
        # If there was no stratification, move aggregated annotations to the top
        # level.
        ht = ht.select(**{ann: ht[ann][0] for ann in entry_agg_funcs})
        global_expr["sample_count"] = group_globals.freq_meta_sample_count[0]
    elif is_reduced:
        # Use the original (full) freq_meta/sample counts so the output is
        # indistinguishable from a non-reduced run.
        global_expr["strata_meta"] = r["freq_meta_full"]
        global_expr["strata_sample_count"] = r["freq_meta_sample_count_full"]
    else:
        global_expr["strata_meta"] = group_globals.freq_meta
        global_expr["strata_sample_count"] = group_globals.freq_meta_sample_count

    ht = ht.annotate_globals(**global_expr)

    return ht


def get_coverage_agg_func(
    dp_field: str = "DP", max_cov_bin: int = 100
) -> Tuple[Callable, Callable]:
    """
    Get a transformation and aggregation function for computing coverage.

    Can be used as an entry aggregation function in `compute_stats_per_ref_site`.

    :param dp_field: Depth field to use for computing coverage. Default is 'DP'.
    :param max_cov_bin: Maximum coverage bin (used when computing samples over X bin). Default is 100.
    :return: Tuple of functions to transform and aggregate coverage.
    """
    return (
        lambda t: hl.if_else(
            hl.is_missing(t[dp_field]) | hl.is_nan(t[dp_field]), 0, t[dp_field]
        ),
        lambda dp: hl.struct(
            # This expression creates a counter DP -> number of samples for DP
            # between 0 and max_cov_bin.
            coverage_counter=hl.agg.counter(hl.min(max_cov_bin, dp)),
            mean=hl.bind(
                lambda mean_dp: hl.if_else(hl.is_nan(mean_dp), 0, mean_dp),
                hl.agg.mean(dp),
            ),
            median_approx=hl.or_else(hl.agg.approx_median(dp), 0),
            total_DP=hl.agg.sum(dp),
        ),
    )


def get_coverage_agg_func_sparse(
    dp_field: str = "DP", max_cov_bin: int = 100
) -> Tuple[Callable, Callable]:
    """
    Sparse-path-compatible variant of `get_coverage_agg_func`.

    Returns a (transform, agg) pair with the same output schema as
    `get_coverage_agg_func`, but with `median_approx` set to a constant
    placeholder (``hl.int32(0)``). The exact median is recomputed at merge
    time by `merge_coverage_stats_array_expression` from the combined
    `coverage_counter`, so the placeholder value is never observed by callers.

    Why this exists: `hl.agg.approx_median` triggers an IR-level
    `NoSuchElementException: key not found: __cse_*` error when used inside
    `hl.experimental.densify`'s scan-densified context (the sparse path's
    densify), and there's no public Hail API for merging two
    `approx_median` sketches after aggregation. The standard
    `get_coverage_agg_func` continues to work on the dense path
    (`compute_stats_per_ref_site` with a VDS input) — use it there.

    :param dp_field: Depth field to use for computing coverage. Default 'DP'.
    :param max_cov_bin: Maximum coverage bin. Default 100.
    :return: ``(transform_fn, agg_fn)`` tuple. The output struct schema
        matches `get_coverage_agg_func`'s; use with
        `merge_coverage_stats_array_expression` to merge the ref-only and
        variant-only paths.
    """
    transform, _ = get_coverage_agg_func(dp_field, max_cov_bin)
    return (
        transform,
        lambda dp: hl.struct(
            coverage_counter=hl.agg.counter(hl.min(max_cov_bin, dp)),
            mean=hl.bind(
                lambda mean_dp: hl.if_else(hl.is_nan(mean_dp), 0, mean_dp),
                hl.agg.mean(dp),
            ),
            median_approx=hl.int32(0),
            total_DP=hl.agg.sum(dp),
        ),
    )


def compute_coverage_stats(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    row_key_fields: List[str] = ["locus"],
    strata_expr: Optional[List[Dict[str, hl.expr.StringExpression]]] = None,
    group_membership_ht: Optional[hl.Table] = None,
    dp_field: str = "DP",
) -> hl.Table:
    """
    Compute coverage statistics for every base of the `reference_ht` provided.

    The following coverage stats are calculated:
        - mean
        - median
        - total DP
        - fraction of samples with coverage above X, for each x in `coverage_over_x_bins`

    The `reference_ht` is a Table that contains a row for each locus coverage that should be
    computed on. It needs to be keyed by `locus`. The `reference_ht` can e.g. be
    created using `get_reference_ht`.

    :param mtds: Input sparse MT or VDS.
    :param reference_ht: Input reference HT.
    :param interval_ht: Optional Table containing intervals to filter to.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
    :param row_key_fields: List of row key fields to use for joining `mtds` with
        `reference_ht`.
    :param strata_expr: Optional list of dicts containing expressions to stratify the
        coverage stats by. Only one of `group_membership_ht` or `strata_expr` can be
        specified.
    :param group_membership_ht: Optional Table containing group membership annotations
        to stratify the coverage stats by. Only one of `group_membership_ht` or
        `strata_expr` can be specified.
    :param dp_field: Name of sample depth field. Default is DP.
    :return: Table with per-base coverage stats.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds

    # Determine the genotype field.
    en = set(mt.entry)
    gt_field = en & {"GT"} or en & {"LGT"}
    if not gt_field:
        raise ValueError("No genotype field found in entry fields!")

    gt_field = gt_field.pop()

    # Add function to compute coverage stats.
    cov_bins = sorted(coverage_over_x_bins)
    rev_cov_bins = list(reversed(cov_bins))
    max_cov_bin = cov_bins[-1]
    cov_bins = hl.array(cov_bins)
    entry_agg_funcs = {
        "coverage_stats": get_coverage_agg_func(
            dp_field=dp_field, max_cov_bin=max_cov_bin
        )
    }

    ht = compute_stats_per_ref_site(
        mtds,
        reference_ht,
        entry_agg_funcs,
        row_key_fields=row_key_fields,
        interval_ht=interval_ht,
        entry_keep_fields=[gt_field, dp_field],
        strata_expr=strata_expr,
        group_membership_ht=group_membership_ht,
    )

    # This expression aggregates the DP counter in reverse order of the cov_bins and
    # computes the cumulative sum over them. It needs to be in reverse order because we
    # want the sum over samples covered by > X.
    def _cov_stats(
        cov_stat: hl.expr.StructExpression, n: hl.expr.Int32Expression
    ) -> hl.expr.StructExpression:
        # The coverage was already floored to the max_coverage_bin, so no more
        # aggregation is needed for the max bin.
        count_expr = cov_stat.coverage_counter
        max_bin_expr = hl.int32(count_expr.get(max_cov_bin, 0))

        # For each of the other bins, coverage is summed between the boundaries.
        bin_expr = hl.range(hl.len(cov_bins) - 1, 0, step=-1)
        bin_expr = bin_expr.map(
            lambda i: hl.sum(
                hl.range(cov_bins[i - 1], cov_bins[i]).map(
                    lambda j: hl.int32(count_expr.get(j, 0))
                )
            )
        )
        bin_expr = hl.cumulative_sum(hl.array([max_bin_expr]).extend(bin_expr))

        bin_expr = {f"over_{x}": bin_expr[i] / n for i, x in enumerate(rev_cov_bins)}

        return cov_stat.annotate(**bin_expr).drop("coverage_counter")

    ht_globals = ht.index_globals()
    if isinstance(ht.coverage_stats, hl.expr.ArrayExpression):
        ht = ht.select_globals(
            coverage_stats_meta=ht_globals.strata_meta.map(
                lambda x: hl.dict(x.items().filter(lambda m: m[0] != "group"))
            ),
            coverage_stats_meta_sample_count=ht_globals.strata_sample_count,
        )
        cov_stats_expr = {
            "coverage_stats": hl.map(
                lambda c, n: _cov_stats(c, n),
                ht.coverage_stats,
                ht_globals.strata_sample_count,
            )
        }
    else:
        cov_stats_expr = _cov_stats(ht.coverage_stats, ht_globals.sample_count)

    ht = ht.transmute(**cov_stats_expr)

    return ht


def get_allele_number_agg_func(gt_field: str = "GT") -> Tuple[Callable, Callable]:
    """
    Get a transformation and aggregation function for computing the allele number.

    Can be used as an entry aggregation function in `compute_stats_per_ref_site`.

    :param gt_field: Genotype field to use for computing the allele number.
    :return: Tuple of functions to transform and aggregate the allele number.
    """
    return lambda t: t[gt_field].ploidy, hl.agg.sum


def compute_allele_number_per_ref_site(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    **kwargs,
) -> hl.Table:
    """
    Compute the allele number per reference site.

    .. note::

        This function supports the `reduce_to_minimal_groups` cost
        optimization (forwarded via `**kwargs` to
        `compute_stats_per_ref_site`). AN is summable across
        stratification groups (it is a sum of per-sample integer
        ploidies), so the optimization produces output identical to a
        non-reduced run. See `compute_stats_per_ref_site` for details.

    :param mtds: Input sparse Matrix Table or VariantDataset.
    :param reference_ht: Table of reference sites.
    :param kwargs: Keyword arguments to pass to `compute_stats_per_ref_site`.
    :return: Table of allele number per reference site.
    """
    if isinstance(mtds, hl.vds.VariantDataset):
        mt = mtds.variant_data
    else:
        mt = mtds

    # Determine the genotype field.
    en = set(mt.entry)
    gt_field = en & {"GT"} or en & {"LGT"}
    if not gt_field:
        raise ValueError(
            "No genotype field found in entry fields, needed for ploidy calculation!"
        )

    # Use ploidy to determine the number of alleles for each sample at each site.
    entry_agg_funcs = {"AN": get_allele_number_agg_func(gt_field.pop())}

    return compute_stats_per_ref_site(mtds, reference_ht, entry_agg_funcs, **kwargs)


def merge_sum_array_expression(
    ref_arr: hl.expr.ArrayExpression,
    var_arr: hl.expr.ArrayExpression,
    sample_counts: hl.expr.ArrayExpression,
) -> hl.expr.ArrayExpression:
    """
    Merge two per-stratum arrays via element-wise sum.

    Suitable for any annotation whose dense-path aggregation is a sum (e.g.,
    AN, raw counts). `var_arr` may be missing (no variant_data row at the
    locus); treated as a zeros array. `sample_counts` is accepted for API
    uniformity with the other merge helpers and is not used here — element-
    wise sum needs no per-stratum sample-count context.

    :param ref_arr: Per-stratum array from the reference-data aggregation.
    :param var_arr: Per-stratum array from the variant-data aggregation; may be
        missing if the locus has no variant_data row.
    :param sample_counts: Per-stratum sample counts; unused.
    :return: Element-wise sum array.
    """
    del sample_counts  # API uniformity.
    # Zeros array of the same per-element type as `ref_arr` — `x - x` is the
    # idiomatic way to produce a zero of arbitrary numeric type in Hail.
    zeros = ref_arr.map(lambda x: x - x)
    var_arr = hl.coalesce(var_arr, zeros)
    return hl.zip(ref_arr, var_arr).map(lambda x: x[0] + x[1])


def _median_from_counter(
    counter: hl.expr.DictExpression, n_total: hl.expr.NumericExpression
) -> hl.expr.Int32Expression:
    """
    Return the integer median DP value from a binned coverage_counter histogram.

    The counter maps DP value to sample count at that DP. The median is the
    smallest DP whose cumulative count from below reaches ``ceil(n_total / 2)``
    (the lower median for even n_total). Returns 0 when the counter is empty
    or n_total is 0.
    """
    half = (hl.int64(n_total) + 1) // 2
    sorted_items = hl.sorted(counter.items(), key=lambda kv: kv[0])
    result = hl.fold(
        lambda acc, item: hl.struct(
            cum=acc.cum + item[1],
            median=hl.if_else(
                hl.is_missing(acc.median) & ((acc.cum + item[1]) >= half),
                item[0],
                acc.median,
            ),
        ),
        hl.struct(cum=hl.int64(0), median=hl.missing(hl.tint32)),
        sorted_items,
    )
    return hl.or_else(result.median, hl.int32(0))


def _merge_one_coverage_stats(
    ref_cs: hl.expr.StructExpression,
    var_cs: hl.expr.StructExpression,
    n_total: hl.expr.NumericExpression,
) -> hl.expr.StructExpression:
    """
    Merge a single pair of ``coverage_stats`` structs into the dense-equivalent.

    Implements the per-stratum merge for `merge_coverage_stats_array_expression`;
    see that function's docstring for the full semantic model and trade-offs.
    """
    ref_counter = ref_cs.coverage_counter
    var_counter = var_cs.coverage_counter
    n_data = hl.sum(ref_counter.values()) + hl.sum(var_counter.values())
    n_no_data = hl.int64(n_total) - n_data

    # Sum non-zero keys directly; rebuild the (0, count) entry with the
    # contribution from samples that had no entry in either path (they
    # contribute DP=0 in the dense path via the missing→0 transform).
    # `set.map` returns a set, but we need an array to append the (0, count)
    # entry — explicitly convert via `hl.array(...)`.
    non_zero_keys = hl.array(ref_counter.key_set().union(var_counter.key_set())).filter(
        lambda k: k != 0
    )
    non_zero_kv = non_zero_keys.map(
        lambda k: (
            k,
            ref_counter.get(k, hl.int64(0)) + var_counter.get(k, hl.int64(0)),
        )
    )
    zero_count = (
        ref_counter.get(0, hl.int64(0)) + var_counter.get(0, hl.int64(0)) + n_no_data
    )
    # Only include the (0, zero_count) entry when it's actually populated —
    # the dense path's counter doesn't carry a 0 key when no sample has DP=0.
    merged_counter = hl.dict(
        hl.if_else(
            zero_count > 0,
            non_zero_kv.append((hl.int32(0), zero_count)),
            non_zero_kv,
        )
    )

    total_dp = ref_cs.total_DP + var_cs.total_DP
    mean = hl.if_else(
        hl.int64(n_total) > 0,
        hl.float64(total_dp) / hl.float64(n_total),
        hl.float64(0),
    )
    median = _median_from_counter(merged_counter, n_total)

    return hl.struct(
        coverage_counter=merged_counter,
        mean=mean,
        median_approx=median,
        total_DP=total_dp,
    )


def merge_coverage_stats_array_expression(
    ref_arr: hl.expr.ArrayExpression,
    var_arr: hl.expr.ArrayExpression,
    sample_counts: hl.expr.ArrayExpression,
) -> hl.expr.ArrayExpression:
    """
    Merge two per-stratum arrays of ``coverage_stats`` structs.

    The dense path's coverage_stats struct has fields ``coverage_counter``,
    ``mean``, ``median_approx``, and ``total_DP``. Merging the ref-only and
    variant-only aggregations into the dense-equivalent requires care because
    the dense aggregation iterates over *every* sample at each locus (no-data
    samples have missing entries → DP transformed to 0), while the sparse
    paths each only iterate over samples present in that path.

    Per-stratum merge logic:
        - ``coverage_counter``: element-wise sum of the two counter dicts, then
          add ``n_no_data = n_total − n_ref − n_var`` to the count at key 0,
          where ``n_ref`` and ``n_var`` are recovered as
          ``sum(counter_ref.values())`` and ``sum(counter_var.values())``.
          This restores the no-data samples' DP=0 contributions that the dense
          path's missing→0 transform would have produced.
        - ``total_DP``: simple sum (no-data samples contribute 0).
        - ``mean``: recomputed as ``total_DP / n_total`` (matches the dense
          semantic; no-data samples contribute 0 to total_DP).
        - ``median_approx``: exact median over the merged counter (which now
          fully represents the dense per-sample DP distribution). This is a
          semantic shift from the dense path's ``hl.agg.approx_median``
          (a non-deterministic streaming estimator), but on the same binned
          input the two are equivalent up to approximation noise.

    `var_arr` may be missing (no variant_data row at the locus); treated as an
    array of zero-coverage_stats structs.

    :param ref_arr: Per-stratum coverage_stats array from the ref-data path.
    :param var_arr: Per-stratum coverage_stats array from the variant-data
        path; may be missing.
    :param sample_counts: Per-stratum total sample counts. Must align with the
        index space of `ref_arr` / `var_arr`. Used to derive ``n_no_data`` per
        stratum and to recompute ``mean``.
    :return: Per-stratum merged coverage_stats array.
    """
    elem_type = ref_arr.dtype.element_type
    zero_cs = hl.struct(
        coverage_counter=hl.empty_dict(
            elem_type["coverage_counter"].key_type,
            elem_type["coverage_counter"].value_type,
        ),
        mean=hl.float64(0),
        median_approx=hl.int32(0),
        total_DP=hl.int64(0),
    )
    n_groups = hl.len(ref_arr)
    var_arr = hl.coalesce(var_arr, hl.range(n_groups).map(lambda _: zero_cs))
    return hl.zip(ref_arr, var_arr, sample_counts).map(
        lambda x: _merge_one_coverage_stats(x[0], x[1], x[2])
    )


def _merge_qual_hists_struct(
    ref_struct: hl.expr.StructExpression,
    var_struct: hl.expr.StructExpression,
) -> hl.expr.StructExpression:
    """
    Recursively merge two nested histogram structs.

    A "histogram leaf" is any struct with fields {bin_edges, bin_freq,
    n_smaller, n_larger} — the shape produced by `hl.agg.hist`. Leaves merge
    via `merge_histograms` with `operation="sum"`. Other structs recurse into
    their fields. This shape mirrors the nested layout produced by
    `qual_hist_expr` (with `split_adj_and_raw=True`, which wraps adj/raw under
    `qual_hists` and `raw_qual_hists` sub-structs).
    """
    fields = list(ref_struct.dtype.fields)
    if {"bin_edges", "bin_freq", "n_smaller", "n_larger"}.issubset(set(fields)):
        return merge_histograms([ref_struct, var_struct], operation="sum")
    return hl.struct(
        **{f: _merge_qual_hists_struct(ref_struct[f], var_struct[f]) for f in fields}
    )


def _zero_qual_hists_struct(struct_type: hl.expr.types.HailType):
    """
    Build a zero-value struct mirroring `struct_type`'s nested histogram shape.

    For histogram leaves, both array-valued fields (`bin_edges`, `bin_freq`)
    are emitted as missing so `merge_histograms` falls back to the non-missing
    side via its built-in `or_else` handling. Counters (`n_smaller`,
    `n_larger`) are zeroed.
    """
    fields = list(struct_type.fields)
    if {"bin_edges", "bin_freq", "n_smaller", "n_larger"}.issubset(set(fields)):
        return hl.struct(
            bin_edges=hl.missing(struct_type["bin_edges"]),
            bin_freq=hl.missing(struct_type["bin_freq"]),
            n_smaller=hl.int64(0),
            n_larger=hl.int64(0),
        )
    return hl.struct(**{f: _zero_qual_hists_struct(struct_type[f]) for f in fields})


def merge_qual_hists_array_expression(
    ref_arr: hl.expr.ArrayExpression,
    var_arr: hl.expr.ArrayExpression,
    sample_counts: hl.expr.ArrayExpression,
) -> hl.expr.ArrayExpression:
    """
    Merge two per-stratum arrays of nested histogram structs.

    Used for `qual_hists`-style annotations (e.g., the struct returned by
    `qual_hist_expr`). Histograms naturally exclude missing values, so no-data
    samples don't appear in either path and the merge is a simple element-wise
    sum of bin counts (delegated to `merge_histograms`). `sample_counts` is
    accepted for API uniformity and unused.

    :param ref_arr: Per-stratum nested-histogram array from the ref-data path.
    :param var_arr: Per-stratum nested-histogram array from the variant-data
        path; may be missing.
    :param sample_counts: Unused; present for API uniformity.
    :return: Per-stratum merged nested-histogram array.
    """
    del sample_counts
    elem_type = ref_arr.dtype.element_type
    zero_struct = _zero_qual_hists_struct(elem_type)
    n_groups = hl.len(ref_arr)
    var_arr = hl.coalesce(var_arr, hl.range(n_groups).map(lambda _: zero_struct))
    return hl.zip(ref_arr, var_arr).map(lambda x: _merge_qual_hists_struct(x[0], x[1]))


def _wrap_transforms_for_sparse_path(
    entry_agg_funcs: Dict[str, Tuple[Callable, Callable]],
    gt_field: str,
) -> Dict[str, Tuple[Callable, Callable]]:
    """
    Wrap each (transform, agg) so the aggregation runs only over present samples.

    The sparse paths each iterate over only the samples present in that path
    (variant_data or reference_data). Each (transform, agg) pair gets a
    Boolean-filtered wrapper applied at the aggregation level:

        - The transform is wrapped to return missing when the entry's genotype
          field is missing (marking a sample not in this path).
        - The agg is wrapped in `hl.agg.filter(hl.is_defined(expr), agg(expr))`
          so the aggregator only consumes values from samples in this path.

    Wrapping at the agg level — not just the transform — is required because
    user-supplied transforms (e.g., `get_coverage_agg_func`'s missing→0 rule)
    can map a missing input to a non-missing constant, and some downstream
    Hail builders (e.g., `hl.min`) treat one missing argument as the other's
    value, which would let no-data samples leak into the counter at the
    `max_cov_bin` key.

    For transforms that already return missing for no-data samples (e.g., the
    AN transform ``t.GT.ploidy``), the transform wrap is a no-op, and the
    agg-level `hl.agg.filter` retains the same set of samples that hail
    aggregators would skip on their own.
    """

    def _wrap(transform, agg):
        def transform_wrapped(t):
            return hl.or_missing(hl.is_defined(t[gt_field]), transform(t))

        def agg_wrapped(expr):
            return hl.agg.filter(hl.is_defined(expr), agg(expr))

        return (transform_wrapped, agg_wrapped)

    return {
        ann: _wrap(transform, agg) for ann, (transform, agg) in entry_agg_funcs.items()
    }


def _aggregate_stats_on_variant_data(
    vmt: hl.MatrixTable,
    reference_ht: hl.Table,
    entry_agg_funcs: Dict[str, Tuple[Callable, Callable]],
    group_membership_ht: hl.Table,
    sex_karyotype_field: Optional[str] = None,
    entry_keep_fields: Optional[Set[str]] = None,
    entry_agg_group_membership: Optional[Dict[str, List[dict]]] = None,
) -> hl.Table:
    """
    Aggregate `entry_agg_funcs` on `variant_data` only — no densification.

    variant_data carries rows only at variant-call sites, so the function
    simply filters to `reference_ht` loci, optionally adjusts ploidy on
    chrX/chrY non-PAR for XY samples, optionally annotates `adj`, and then
    runs `agg_by_strata` per (locus, alleles). The transforms in
    `entry_agg_funcs` are wrapped to return missing for undefined-entry
    samples (see `_wrap_transforms_for_sparse_path`).

    Caller must ensure the VDS is not split-multiallelic. In a split VDS a
    heterozygous multi-allelic call appears in two rows at the same locus and
    would be double-counted at the per-locus rekeying step below.

    :param vmt: `variant_data` MatrixTable of a VDS.
    :param reference_ht: Reference sites Table; the output is filtered to
        these loci.
    :param entry_agg_funcs: As in `compute_stats_per_ref_site`.
    :param group_membership_ht: Pre-built group membership Table — must be
        the same one used for the ref-data path so strata indices align.
    :param sex_karyotype_field: Optional sex-karyotype col-field on `vmt`.
    :param entry_keep_fields: Entry fields to retain on `vmt` for the
        aggregation. The genotype field (GT/LGT) and `adj` (if needed) are
        added automatically.
    :param entry_agg_group_membership: Optional per-annotation narrowing —
        forwarded to `agg_by_strata`.
    :return: Table keyed by `locus` with one column per annotation in
        `entry_agg_funcs`. The (locus, alleles) row key from `vmt` is
        collapsed to locus assuming one row per locus (unsplit VDS).
    """
    en = set(vmt.entry)
    gt_field_set = en & {"GT"} or en & {"LGT"}
    if not gt_field_set:
        raise ValueError(
            "`variant_data` has no GT or LGT entry field; cannot aggregate."
        )
    gt_field = next(iter(gt_field_set))

    if sex_karyotype_field is not None:
        vmt = vmt.annotate_entries(
            **{
                gt_field: adjusted_sex_ploidy_expr(
                    vmt.locus, vmt[gt_field], vmt[sex_karyotype_field]
                )
            }
        )

    # Determine if `adj` is needed (mirrors `compute_stats_per_ref_site`'s
    # detection logic).
    g = group_membership_ht.index_globals()
    adj_needed = hl.eval(
        hl.any(g.get("adj_groups", hl.empty_array("bool")))
        | hl.any(
            g.get("freq_meta", hl.empty_array("dict<str, str>")).map(
                lambda x: x.get("group", "NA") == "adj"
            )
        )
    )
    if adj_needed and "adj" not in vmt.entry:
        vmt = annotate_adj(vmt)

    entry_keep = set(entry_keep_fields or set()) | {gt_field}
    if adj_needed:
        entry_keep |= {"adj"}
    entry_keep &= set(vmt.entry)
    vmt = vmt.select_entries(*entry_keep)
    vmt = vmt.filter_rows(hl.is_defined(reference_ht[vmt.locus]))

    wrapped_funcs = _wrap_transforms_for_sparse_path(entry_agg_funcs, gt_field)
    per_row_ht = agg_by_strata(
        vmt,
        wrapped_funcs,
        group_membership_ht=group_membership_ht,
        entry_agg_group_membership=entry_agg_group_membership,
    )

    # Re-key by locus only. For an unsplit VDS there is at most one row per
    # locus in variant_data, so dropping `alleles` collapses cleanly.
    per_row_ht = per_row_ht.key_by("locus")
    if "alleles" in per_row_ht.row:
        per_row_ht = per_row_ht.drop("alleles")
    return per_row_ht


def _resolve_sample_counts_for_annotation(
    ann: str,
    entry_agg_group_membership: Optional[Dict[str, List[dict]]],
    strata_meta_full: List[Dict[str, str]],
    strata_sample_count_full: List[int],
    reduction_globals: Optional[Dict[str, object]],
) -> List[int]:
    """
    Compute the per-stratum sample counts that align with an annotation's array.

    For an annotation NOT narrowed via `entry_agg_group_membership`, the
    returned list is the full `strata_sample_count`. For a narrowed annotation
    each entry in `entry_agg_group_membership[ann]` is resolved to one or more
    indices in the (full) `strata_meta`; the sample count for that target is
    the leaf entry's count, or the SUM of its leaf children's counts when the
    target is a non-leaf parent under reduction.
    """
    if not entry_agg_group_membership or ann not in entry_agg_group_membership:
        return list(strata_sample_count_full)

    decomposition = reduction_globals["decomposition"] if reduction_globals else {}
    counts = []
    for target in entry_agg_group_membership[ann]:
        target_dict = dict(target)
        found = False
        for i, m in enumerate(strata_meta_full):
            if dict(m) == target_dict:
                if i in decomposition and decomposition[i]:
                    counts.append(
                        sum(strata_sample_count_full[c] for c in decomposition[i])
                    )
                else:
                    counts.append(strata_sample_count_full[i])
                found = True
                break
        if not found:
            raise ValueError(
                f"`entry_agg_group_membership[{ann!r}]` target {target} is"
                " not in `strata_meta` (or `freq_meta_full` under reduction)."
            )
    return counts


def compute_stats_per_ref_site_sparse(
    vds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    entry_agg_funcs: Dict[str, Tuple[Callable, Callable]],
    include_variant_data: bool = False,
    merge_funcs: Optional[Dict[str, Callable]] = None,
    row_key_fields: Union[Tuple[str], List[str]] = ("locus",),
    interval_ht: Optional[hl.Table] = None,
    entry_keep_fields: Union[Tuple[str], List[str], Set[str]] = None,
    row_keep_fields: Union[Tuple[str], List[str], Set[str]] = None,
    entry_agg_group_membership: Optional[Dict[str, List[dict]]] = None,
    strata_expr: Optional[List[Dict[str, hl.expr.StringExpression]]] = None,
    group_membership_ht: Optional[hl.Table] = None,
    sex_karyotype_field: Optional[str] = None,
    reduce_to_minimal_groups: bool = False,
    non_summable_strata: Optional[Set[str]] = None,
    reducible_aggs: Optional[Set[str]] = None,
) -> hl.Table:
    """
    Compute per-reference-site stats from a VDS without a full densify.

    Sparse-aware counterpart to `compute_stats_per_ref_site` for VDS inputs.
    The standard function densifies the VDS via `hl.vds.to_dense_mt`, which
    performs a row union of `variant_data`, `reference_data`, and
    `reference_ht` and is the dominant cost when aggregating across a large
    reference panel. This function instead aggregates on `reference_data`
    alone via `hl.experimental.densify` (which uses `hl.scan._densify`
    internally) and — optionally, via `include_variant_data=True` — also
    aggregates directly on `variant_data` (no densify needed there because
    `variant_data` only has rows at variant-call sites) and merges the two
    results per annotation.

    .. rubric:: Annotation transforms run in "sparse mode"

    Each transform in `entry_agg_funcs` is wrapped to return missing for
    samples whose entry has an undefined genotype (GT/LGT) — i.e., samples
    that are NOT in the current path. Hail aggregators skip missing inputs by
    standard semantics, so each path's aggregation naturally restricts to the
    samples present in that path. This is required for correctness when the
    user's transform converts missing values to a non-missing constant (e.g.,
    `get_coverage_agg_func`'s `missing→0` rule, which would otherwise double-
    count "other path" samples as DP=0 contributions). For transforms that
    already return missing for no-data samples (e.g., AN's `t.GT.ploidy`),
    the wrap is a no-op.

    .. rubric:: Merging the two paths

    When `include_variant_data=True`, `merge_funcs` must supply a callable
    `(ref_arr, var_arr, sample_counts) -> merged_arr` for every annotation
    in `entry_agg_funcs`. Three helpers are provided for the gnomAD use
    cases:

        - `merge_sum_array_expression`: element-wise sum (AN-style).
        - `merge_coverage_stats_array_expression`: per-stratum merge of the
          `coverage_stats` struct produced by `get_coverage_agg_func`.
        - `merge_qual_hists_array_expression`: recursive element-wise sum of
          nested histogram structs (e.g., from `qual_hist_expr`).

    `sample_counts` passed to each merge function is the per-stratum
    sample-count array that aligns with that annotation's per-row array
    shape — full `strata_sample_count` for non-narrowed annotations, and
    the resolved counts for annotations narrowed via
    `entry_agg_group_membership`.

    .. rubric:: Constraints

    - `vds.variant_data` must be unsplit (one row per locus). In a split VDS
      a multi-allelic call lands on multiple rows at the same locus and would
      be double-counted by the per-locus rekeying step.
    - When `include_variant_data=True`, either `strata_expr` or
      `group_membership_ht` must be supplied so the two intermediate arrays
      share aligned strata indices.

    :param vds: Input `hl.vds.VariantDataset`.
    :param reference_ht: Table of reference sites.
    :param entry_agg_funcs: As in `compute_stats_per_ref_site` — dict of
        ``ann_name -> (transform_fn, agg_fn)``.
    :param include_variant_data: When True, also aggregate on `variant_data`
        and merge the two paths per annotation using `merge_funcs`. Default
        False (ref-only path; cheapest but undercounts at variant-called
        samples).
    :param merge_funcs: Required when `include_variant_data=True`; dict of
        ``ann_name -> Callable[[ref_arr, var_arr, sample_counts], merged_arr]``.
        See the three helpers exported alongside this function.
    :param row_key_fields: Forwarded to the inner aggregation calls.
    :param interval_ht: Optional interval Table. When supplied, the VDS is
        filtered with `hl.vds.filter_intervals(split_reference_blocks=False)`
        and `reference_ht` is filtered to the same intervals up front.
    :param entry_keep_fields: Forwarded to the inner aggregation calls.
    :param row_keep_fields: Forwarded to the inner aggregation calls.
    :param entry_agg_group_membership: Per-annotation narrowing dict
        (forwarded to both paths' inner aggregations).
    :param strata_expr: Mutually exclusive with `group_membership_ht`. When
        supplied, a single `group_membership_ht` is built up front and used
        for both paths.
    :param group_membership_ht: Pre-built group membership Table. Mutually
        exclusive with `strata_expr`.
    :param sex_karyotype_field: Optional sex-karyotype col-field on
        `variant_data`. Propagated to `reference_data.cols()` automatically.
    :param reduce_to_minimal_groups: Forwarded to the ref-data path. When
        True under `include_variant_data=True`, the variant-data path's
        leaf-shape per-row arrays are expanded back to full shape via
        `expand_strata_array_from_leaves` so the merge sees aligned full-
        shape arrays on both sides.
    :param non_summable_strata: Forwarded.
    :param reducible_aggs: Forwarded to the ref-data path. Controls which
        annotations are leaf-expanded on the variant-data path as well.
    :return: Table keyed by `row_key_fields` with one column per annotation
        in `entry_agg_funcs` and `strata_meta` / `strata_sample_count`
        globals matching `compute_stats_per_ref_site`'s output.
    """
    if not isinstance(vds, hl.vds.VariantDataset):
        raise ValueError("`vds` must be a `hl.vds.VariantDataset`.")
    if strata_expr is not None and group_membership_ht is not None:
        raise ValueError(
            "Only one of `strata_expr` or `group_membership_ht` can be" " specified."
        )
    if include_variant_data:
        if strata_expr is None and group_membership_ht is None:
            raise ValueError(
                "When `include_variant_data=True`, one of `strata_expr` or"
                " `group_membership_ht` must be supplied so the two"
                " intermediate per-stratum arrays share aligned indices for"
                " the merge."
            )
        if merge_funcs is None:
            raise ValueError(
                "When `include_variant_data=True`, `merge_funcs` must be"
                " supplied with one merge callable per annotation in"
                " `entry_agg_funcs`."
            )
        missing = set(entry_agg_funcs) - set(merge_funcs)
        if missing:
            raise ValueError(
                "`merge_funcs` is missing entries for:"
                f" {sorted(missing)}. Each annotation in `entry_agg_funcs`"
                " needs a merge callable."
            )

    if interval_ht is not None:
        vds = hl.vds.filter_intervals(vds, interval_ht, split_reference_blocks=False)
        reference_ht = reference_ht.filter(
            hl.is_defined(interval_ht[reference_ht.locus])
        )

    rmt = vds.reference_data
    vmt = vds.variant_data

    if sex_karyotype_field is not None:
        if sex_karyotype_field not in vmt.col:
            raise ValueError(
                f"`sex_karyotype_field` {sex_karyotype_field!r} not found"
                " in `variant_data` cols."
            )
        if sex_karyotype_field not in rmt.col:
            sex_ht = vmt.cols().select(sex_karyotype_field)
            rmt = rmt.annotate_cols(
                **{sex_karyotype_field: sex_ht[rmt.col_key][sex_karyotype_field]}
            )

    # Build a single group_membership_ht from strata_expr up front so both
    # paths share the same group set. Resolve expressions against vmt.cols
    # (rmt and vmt share the same sample set by construction).
    if strata_expr is not None:
        cols_ht = vmt.annotate_cols(
            **{k: v for d in strata_expr for k, v in d.items()}
        ).cols()
        strata_resolved = [{k: cols_ht[k] for k in d} for d in strata_expr]
        group_membership_ht = generate_freq_group_membership_array(
            cols_ht,
            strata_resolved,
            no_raw_group=True,
            reduce_to_minimal_groups=reduce_to_minimal_groups,
            non_summable_strata=non_summable_strata,
            group_label="raw",
        )
        strata_expr = None

    # Determine the genotype field on the reference_data MT for the auto-wrap.
    rmt_en = set(rmt.entry)
    rmt_gt_field_set = rmt_en & {"GT"} or rmt_en & {"LGT"}
    if not rmt_gt_field_set:
        raise ValueError(
            "`vds.reference_data` has no GT or LGT entry field; cannot" " aggregate."
        )
    rmt_gt_field = next(iter(rmt_gt_field_set))

    # Step 1: Aggregate on reference_data via the sparse-MT path.
    ref_ht_out = compute_stats_per_ref_site(
        rmt,
        reference_ht,
        _wrap_transforms_for_sparse_path(entry_agg_funcs, rmt_gt_field),
        row_key_fields=row_key_fields,
        entry_keep_fields=entry_keep_fields,
        row_keep_fields=row_keep_fields,
        entry_agg_group_membership=entry_agg_group_membership,
        strata_expr=None,
        group_membership_ht=group_membership_ht,
        sex_karyotype_field=sex_karyotype_field,
        reducible_aggs=reducible_aggs,
        # Reduction is already applied (or not) in the supplied
        # group_membership_ht; forward False and rely on the supplied table's
        # `freq_reduced` global for detection.
        reduce_to_minimal_groups=False,
        non_summable_strata=non_summable_strata,
    )

    if not include_variant_data:
        return ref_ht_out

    # Step 2: Aggregate on variant_data directly (no densify).
    var_ht_out = _aggregate_stats_on_variant_data(
        vmt,
        reference_ht,
        entry_agg_funcs,
        group_membership_ht=group_membership_ht,
        sex_karyotype_field=sex_karyotype_field,
        entry_keep_fields=set(entry_keep_fields or set()),
        entry_agg_group_membership=entry_agg_group_membership,
    )

    # Under reduction, the variant-data path emits leaf-shape arrays for the
    # reducible annotations; expand them to full shape so the merge sees
    # aligned arrays on both sides. Narrowed annotations (those listed in
    # `entry_agg_group_membership`) are already at their final length-k shape
    # on both sides and need no expansion.
    g = group_membership_ht.index_globals()
    is_reduced = "freq_reduced" in group_membership_ht.globals and hl.eval(
        g.freq_reduced
    )
    if is_reduced:
        r = _read_reduction_globals(g)
        grouped = (
            set(entry_agg_group_membership) if entry_agg_group_membership else set()
        )
        if reducible_aggs is None:
            anns_to_expand = set(entry_agg_funcs) - grouped
        else:
            anns_to_expand = set(reducible_aggs)
        var_ht_out = var_ht_out.annotate(
            **{
                ann: expand_strata_array_from_leaves(
                    var_ht_out[ann],
                    r["leaf_indices"],
                    r["decomposition"],
                    r["n_full"],
                )
                for ann in anns_to_expand
            }
        )

    # Step 3: Merge per annotation.
    strata_meta_full = hl.eval(ref_ht_out.strata_meta)
    strata_sample_count_full = hl.eval(ref_ht_out.strata_sample_count)
    reduction_info = _read_reduction_globals(g) if is_reduced else None

    var_indexed = var_ht_out[ref_ht_out[row_key_fields[0]]]
    merged_annotations = {}
    for ann in entry_agg_funcs:
        per_ann_counts = _resolve_sample_counts_for_annotation(
            ann,
            entry_agg_group_membership,
            strata_meta_full,
            strata_sample_count_full,
            reduction_info,
        )
        sample_counts_expr = hl.array([hl.int64(c) for c in per_ann_counts])
        merged_annotations[ann] = merge_funcs[ann](
            ref_ht_out[ann],
            var_indexed[ann],
            sample_counts_expr,
        )
    return ref_ht_out.annotate(**merged_annotations)


def compute_allele_number_per_ref_site_sparse(
    vds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    include_variant_data: bool = False,
    strata_expr: Optional[List[Dict[str, hl.expr.StringExpression]]] = None,
    group_membership_ht: Optional[hl.Table] = None,
    sex_karyotype_field: Optional[str] = None,
    interval_ht: Optional[hl.Table] = None,
    row_key_fields: Union[Tuple[str], List[str]] = ("locus",),
    row_keep_fields: Union[Tuple[str], List[str], Set[str]] = None,
    entry_keep_fields: Union[Tuple[str], List[str], Set[str]] = None,
    reduce_to_minimal_groups: bool = False,
    non_summable_strata: Optional[Set[str]] = None,
) -> hl.Table:
    """
    Compute AN per reference site without a full VDS densify.

    Thin wrapper around `compute_stats_per_ref_site_sparse` for the common
    AN-only case. Uses `merge_sum_array_expression` to merge the ref-only
    and variant-only paths when `include_variant_data=True`. See
    `compute_stats_per_ref_site_sparse` for the full semantic model.

    :param vds: Input VDS.
    :param reference_ht: Table of reference sites.
    :param include_variant_data: If True, also aggregate AN from
        `variant_data` (no densify) and merge via element-wise sum to
        recover the full AN matching the dense-path output.
    :param strata_expr: Mutually exclusive with `group_membership_ht`.
    :param group_membership_ht: Pre-built group membership Table. Mutually
        exclusive with `strata_expr`.
    :param sex_karyotype_field: Optional sex-karyotype col-field; applied to
        both paths.
    :param interval_ht: Optional interval filter Table.
    :param row_key_fields: Forwarded.
    :param row_keep_fields: Forwarded.
    :param entry_keep_fields: Forwarded.
    :param reduce_to_minimal_groups: Forwarded.
    :param non_summable_strata: Forwarded.
    :return: Table keyed by `row_key_fields` with an `AN` array per stratum
        and `strata_meta` / `strata_sample_count` globals.
    """
    if not isinstance(vds, hl.vds.VariantDataset):
        raise ValueError("`vds` must be a `hl.vds.VariantDataset`.")

    # The generic function dispatches the same `entry_agg_funcs` to both the
    # ref-data and variant-data paths. In VDSs where the two halves use
    # different genotype field names (e.g., test fixtures with LGT on
    # reference_data and GT on variant_data), a hardcoded field name in the
    # transform would fail on one side. Inspect the MT at expression-build
    # time and pick whichever of GT/LGT is present.
    def _an_transform(mt):
        gt = "GT" if "GT" in mt.entry else "LGT"
        return mt[gt].ploidy

    entry_agg_funcs = {"AN": (_an_transform, hl.agg.sum)}
    merge_funcs = {"AN": merge_sum_array_expression} if include_variant_data else None
    return compute_stats_per_ref_site_sparse(
        vds,
        reference_ht,
        entry_agg_funcs,
        include_variant_data=include_variant_data,
        merge_funcs=merge_funcs,
        row_key_fields=row_key_fields,
        interval_ht=interval_ht,
        entry_keep_fields=entry_keep_fields,
        row_keep_fields=row_keep_fields,
        strata_expr=strata_expr,
        group_membership_ht=group_membership_ht,
        sex_karyotype_field=sex_karyotype_field,
        reduce_to_minimal_groups=reduce_to_minimal_groups,
        non_summable_strata=non_summable_strata,
        reducible_aggs={"AN"} if include_variant_data else None,
    )


def filter_ref_blocks(
    t: Union[hl.MatrixTable, hl.Table],
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter ref blocks out of the Table or MatrixTable.

    :param t: Input MT/HT
    :return: MT/HT with ref blocks removed
    """
    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows((hl.len(t.alleles) > 1))
    else:
        t = t.filter((hl.len(t.alleles) > 1))

    return t
