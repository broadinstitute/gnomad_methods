# noqa: D100

import logging
from typing import Callable, Dict, List, Optional, Set, Tuple, Union

import hail as hl

from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    agg_by_strata,
    annotate_adj,
    fs_from_sb,
    generate_freq_group_membership_array,
    get_adj_expr,
    get_lowqual_expr,
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
    :param treat_fields_as_allele_specific: Treat info fields as allele-specific. Defaults to False.
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
                if k.startswith("AS_RAW_") and k.endswith("RankSum"):
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
        info_expr = get_as_info_expr(mt)

    # Compute AS info expr using gvcf_info allele specific annotations.
    if as_annotations:
        if info_expr is not None:
            quasi_info_expr = info_expr
        info_expr = get_as_info_expr(
            mt,
            **AS_INFO_AGG_FIELDS,
            treat_fields_as_allele_specific=True,
        )

    if info_expr is not None:
        # Add allele specific pab_max
        info_expr = info_expr.annotate(
            AS_pab_max=pab_max_expr(mt.LGT, mt.LAD, mt.LA, hl.len(mt.alleles))
        )

    if site_annotations:
        site_expr = get_site_info_expr(mt)
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
) -> hl.Table:
    """
    Compute stats per site in a reference Table.

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
        # that gives stratification group membership info based on 'strata_expr'. The
        # returned Table has the following annotations: 'freq_meta',
        # 'freq_meta_sample_count', and 'group_membership'. By default, this
        # function returns annotations where the second element is a placeholder for the
        # "raw" frequency of all samples, where the first 2 elements are the same sample
        # set, but 'freq_meta' starts with [{"group": "adj", "group": "raw", ...]. Use
        # `no_raw_group` to exclude the "raw" group so there is a single annotation
        # representing the full samples set. Update all 'freq_meta' entries' "group"
        # to "raw" because `generate_freq_group_membership_array` will return them all
        # as "adj" since it was built for frequency computation, but for the coverage
        # computation we don't want to do any filtering.
        group_membership_ht = generate_freq_group_membership_array(
            ht, strata_expr, no_raw_group=True
        )
        group_membership_ht = group_membership_ht.annotate_globals(
            freq_meta=group_membership_ht.freq_meta.map(
                lambda x: hl.dict(
                    x.items().map(
                        lambda m: hl.if_else(m[0] == "group", ("group", "raw"), m)
                    )
                )
            )
        )

    if is_vds:
        rmt = mtds.reference_data
        mtds = hl.vds.VariantDataset(
            rmt.select_entries(*((set(entry_keep_fields) & set(rmt.entry)) | {"END"})),
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
    global_expr = {}
    if no_strata:
        # If there was no stratification, move aggregated annotations to the top
        # level.
        ht = ht.select(**{ann: ht[ann][0] for ann in entry_agg_funcs})
        global_expr["sample_count"] = group_globals.freq_meta_sample_count[0]
    else:
        # If there was stratification, add the metadata and sample count info for the
        # stratification to the globals.
        global_expr["strata_meta"] = group_globals.freq_meta
        global_expr["strata_sample_count"] = group_globals.freq_meta_sample_count

    ht = ht.annotate_globals(**global_expr)

    return ht


def compute_coverage_stats(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    row_key_fields: List[str] = ["locus"],
    strata_expr: Optional[List[Dict[str, hl.expr.StringExpression]]] = None,
    group_membership_ht: Optional[hl.Table] = None,
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

    :param mtds: Input sparse MT or VDS
    :param reference_ht: Input reference HT
    :param interval_ht: Optional Table containing intervals to filter to
    :param coverage_over_x_bins: List of boundaries for computing samples over X
    :param row_key_fields: List of row key fields to use for joining `mtds` with
        `reference_ht`
    :param strata_expr: Optional list of dicts containing expressions to stratify the
        coverage stats by. Only one of `group_membership_ht` or `strata_expr` can be
        specified.
    :param group_membership_ht: Optional Table containing group membership annotations
        to stratify the coverage stats by. Only one of `group_membership_ht` or
        `strata_expr` can be specified.
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
        "coverage_stats": (
            lambda t: hl.if_else(hl.is_missing(t.DP) | hl.is_nan(t.DP), 0, t.DP),
            lambda dp: hl.struct(
                # This expression creates a counter DP -> number of samples for DP
                # between 0 and max_cov_bin.
                coverage_counter=hl.agg.counter(hl.min(max_cov_bin, dp)),
                mean=hl.if_else(hl.is_nan(hl.agg.mean(dp)), 0, hl.agg.mean(dp)),
                median_approx=hl.or_else(hl.agg.approx_median(dp), 0),
                total_DP=hl.agg.sum(dp),
            ),
        )
    }

    ht = compute_stats_per_ref_site(
        mtds,
        reference_ht,
        entry_agg_funcs,
        row_key_fields=row_key_fields,
        interval_ht=interval_ht,
        entry_keep_fields=[gt_field, "DP"],
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
