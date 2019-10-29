from .generic import *

INFO_SUM_AGG_FIELDS= ['RAW_MQ', 'QUALapprox']
INFO_INT32_SUM_AGG_FIELDS = ['DP', 'MQ_DP', 'VarDP']
INFO_MEDIAN_AGG_FIELDS = ['ReadPosRankSum', 'MQRankSum']

def compute_last_ref_block_end(mt: hl.MatrixTable) -> hl.Table:
    """
    This function takes a sparse MT and computes for each row the genomic position of the
    most upstream reference block overlapping that row.
    Note that since reference blocks do not extend beyond contig boundaries, only the position is kept.

    This function returns a Table with that annotation.  (`last_END_position`).

    :param MatrixTable mt: Input MatrixTable
    :return: Output Table with `last_END_position` annotation
    :rtype: Table
    """
    mt = mt.select_entries('END')
    ht = mt._localize_entries('__entries', '__cols')
    ht = ht.select(
        last_END_position=hl.or_else(
            hl.min(
                hl.scan.array_agg(
                    lambda entry: hl.scan._prev_nonnull(
                        hl.or_missing(
                            hl.is_defined(entry.END),
                            hl.tuple([
                                ht.locus,
                                entry.END
                            ])
                        )
                    ),
                    ht.__entries
                ).map(
                    lambda x: hl.or_missing(
                        (x[1] >= ht.locus.position) & (x[0].contig == ht.locus.contig),
                        x[0].position
                    )
                )
            ),
            ht.locus.position
        )
    )
    return ht.select_globals()


def densify_sites(
        mt: hl.MatrixTable,
        sites_ht: hl.Table,
        last_END_positions_ht: hl.Table,
        semi_join_rows: bool = True
) -> hl.MatrixTable:
    """
    Densifies the input sparse MT at the sites in `sites_ht` reading the minimal amount of data required.
    Note that only rows that appear both in `mt` and `sites_ht` are returned.

    :param MatrixTable mt: Input sparse MT
    :param Table sites_ht: Desired sites to densify
    :param Table last_END_positions_ht: Table storing positions of the furthest ref block (END tag)
    :param bool semi_join_rows: Whether to filter the MT rows based on semi-join (default, better if sites_ht is large) or based on filter_intervals (better if sites_ht only contains a few sites)
    :return: Dense MT filtered to the sites in `sites_ht`
    :rtype: MatrixTable
    """
    logger.info("Computing intervals to densify from sites Table.")
    sites_ht = sites_ht.key_by('locus')
    sites_ht = sites_ht.annotate(
        interval=hl.locus_interval(
            sites_ht.locus.contig,
            last_END_positions_ht[sites_ht.key].last_END_position,
            end=sites_ht.locus.position,
            includes_end=True,
            reference_genome=sites_ht.locus.dtype.reference_genome
        )
    )
    sites_ht = sites_ht.filter(hl.is_defined(sites_ht.interval))

    if semi_join_rows:
        mt = mt.filter_rows(hl.is_defined(sites_ht.key_by('interval')[mt.locus]))
    else:
        logger.info("Collecting intervals to densify.")
        intervals = sites_ht.interval.collect()

        print("Found {0} intervals, totalling {1} bp in the dense Matrix.".format(
            len(intervals),
            sum([interval_length(interval) for interval in union_intervals(intervals)])
        ))

        mt = hl.filter_intervals(mt, intervals)

    mt = hl.experimental.densify(mt)

    return mt.filter_rows(
        hl.is_defined(sites_ht[mt.locus])
    )


def _get_info_agg_expr(
        mt: hl.MatrixTable,
        sum_agg_fields: Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_SUM_AGG_FIELDS,
        int32_sum_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_INT32_SUM_AGG_FIELDS,
        median_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_MEDIAN_AGG_FIELDS,
        prefix: str = ''
) -> Dict[str, hl.expr.Aggregation]:
    """
    Helper function containing code to create Aggregators for both site or AS info expression aggregations.

    :param MatrixTable mt: Input MT
    :param list of str or dict of str -> NumericExpression sum_agg_fields: Fields to aggregate using sum.
    :param list of str or dict of str -> NumericExpression int32_sum_agg_fields: Fields to aggregate using sum using int32.
    :param list of str or dict of str -> NumericExpression median_agg_fields: Fields to aggregate using median.
    :param str prefix: Optional prefix for the fields. Used for adding 'AS_' in the AS case.

    :return: Dictionary of expression names and their corresponding aggregation Expression
    :rtype: dict of str -> Aggregation
    """

    def _agg_list_to_dict(mt: hl.MatrixTable, fields: List[str]) -> Dict[str, hl.expr.NumericExpression]:
        out_fields = {}
        if 'gvcf_info' in mt.entry:
            out_fields = {f: mt.gvcf_info[f] for f in fields if f in mt.gvcf_info}

        out_fields.update(
            {f: mt[f] for f in fields if f in mt.entry}
        )

        #Check that all fields were found
        missing_fields = [f for f in fields if f not in out_fields]
        if missing_fields:
            raise  ValueError("Could not find the following field(s)in the MT entry schema (or nested under mt.gvcf_info: {}".format(
                ",".join(missing_fields)
            ))

        return out_fields

    # Map str to expressions where needed
    if isinstance(sum_agg_fields, list):
        sum_agg_fields = _agg_list_to_dict(mt, sum_agg_fields)

    if isinstance(int32_sum_agg_fields, list):
        int32_sum_agg_fields = _agg_list_to_dict(mt, int32_sum_agg_fields)

    if isinstance(median_agg_fields, list):
        median_agg_fields = _agg_list_to_dict(mt, median_agg_fields)

    # Create aggregators
    agg_expr = {}

    agg_expr.update({
        f'{prefix}{k}': hl.agg.approx_quantiles(expr, 0.5)
        for k, expr in median_agg_fields.items()
    })
    agg_expr.update({
        f'{prefix}{k}': hl.agg.sum(expr)
        for k, expr in sum_agg_fields.items()
    })
    agg_expr.update({
        f'{prefix}{k}': hl.int32(hl.agg.sum(expr))
        for k, expr in int32_sum_agg_fields.items()
    })

    return agg_expr


def get_as_info_expr(
        mt: hl.MatrixTable,
        sum_agg_fields: Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_SUM_AGG_FIELDS,
        int32_sum_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_INT32_SUM_AGG_FIELDS,
        median_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_MEDIAN_AGG_FIELDS,
) -> hl.expr.StructExpression:
    """
    Returns an allele-specific annotation Struct containing typical VCF INFO fields from GVCF INFO fields stored in the MT entries.
    Note that while not a parameter, if `SB` is found, it will be aggregated in `AS_SB_TABLE`

    :param MatrixTable mt: Input Matrix Table
    :param list of str or dict of str -> NumericExpression sum_agg_fields: Fields to aggregate using sum.
    :param list of str or dict of str -> NumericExpression int32_sum_agg_fields: Fields to aggregate using sum using int32.
    :param list of str or dict of str -> NumericExpression median_agg_fields: Fields to aggregate using median.
    :return: Expression containing the AS info fields
    :rtype: StructExpression
    """

    agg_expr = _get_info_agg_expr(
        mt=mt,
        sum_agg_fields=sum_agg_fields,
        int32_sum_agg_fields=int32_sum_agg_fields,
        median_agg_fields=median_agg_fields,
        prefix='AS_'
    )

    # Modify aggregations to aggregate per allele
    agg_expr = {
        f: hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai),
                expr
            ),
            hl.range(1, hl.len(mt.alleles))
        )
        for f, expr in agg_expr.items()
    }

    # Add SB aggregation
    sb_expr = mt['SB'] if 'SB' in mt.entry else mt['gvcf_info']['SB'] if 'gvcf_info' in mt and 'SB' in 'gvcf_info' else None
    if sb_expr is not None:
        agg_expr['AS_SB_TABLE'] = hl.array([
            hl.agg.array_agg(lambda x: hl.agg.sum(x), sb_expr[:2]).map(lambda x: hl.int32(x)) # ref
        ]).extend(
            hl.agg.array_agg(
                lambda ai: hl.agg.filter(
                    mt.LA.contains(ai),
                    hl.agg.array_agg(lambda x: hl.agg.sum(x), sb_expr[2:]).map(lambda x: hl.int32(x)) # each alt
                ),
                hl.range(1, hl.len(mt.alleles))
            )
        )

    # Run aggregations
    info = hl.struct(
        **agg_expr
    )

    # Add metrics that combine aggregations
    additional_metrics = {}
    if 'AS_RAW_MQ' in info and 'AS_MQ_DP' in info:
        additional_metrics['AS_MQ'] = (info.AS_RAW_MQ / info.AS_MQ_DP).map(lambda x: hl.sqrt(x))
    if 'AS_QUALapprox' in info and 'AS_VarDP' in info:
        additional_metrics['AS_QD'] = info.AS_QUALapprox/ info.AS_VarDP
    if 'AS_SB_TABLE' in info:
        additional_metrics['AS_FS'] = hl.range(1, hl.len(mt.alleles)).map(
            lambda i: fs_from_sb(info.AS_SB_TABLE[0].extend(info.AS_SB_TABLE[i]))
        )

    return info.annotate(**additional_metrics)


def get_site_info_expr(
        mt: hl.MatrixTable,
        sum_agg_fields: Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_SUM_AGG_FIELDS,
        int32_sum_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_INT32_SUM_AGG_FIELDS,
        median_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_MEDIAN_AGG_FIELDS
) -> hl.expr.StructExpression:
    """
    Creates a site-level annotation Struct aggregating typical VCF INFO fields from GVCF INFO fields stored in the MT entries.

    If the fields to be aggregate (`sum_agg_fields`, `int32_sum_agg_fields`, `median_agg_fields`) are passed as list of str,
    then they should correspond to entry fields in `mt` or in `mt.gvcf_info`.
    Priority is given to entry fields in `mt` to those in `mt.gvcf_info` in case of a name clash.
    Note that `SB` will also be aggregated if present.

    :param MatrixTable mt: Input Matrix Table
    :param list of str or dict of str -> NumericExpression sum_agg_fields: Fields to aggregate using sum.
    :param list of str or dict of str -> NumericExpression int32_sum_agg_fields: Fields to aggregate using sum using int32.
    :param list of str or dict of str -> NumericExpression median_agg_fields: Fields to aggregate using median.
    :return: Expression containing the site-level info fields
    :rtype: StructExpression
    """

    agg_expr = _get_info_agg_expr(
        mt=mt,
        sum_agg_fields=sum_agg_fields,
        int32_sum_agg_fields=int32_sum_agg_fields,
        median_agg_fields=median_agg_fields
    )

    # Run aggregator on non-ref genotypes
    info = hl.agg.filter(
        mt.LGT.is_non_ref(),
        hl.struct(
            **agg_expr
        )
    )

    # Add metrics that combine aggregations
    additional_metrics = {}
    if 'RAW_MQ' in info and 'MQ_DP' in info:
        additional_metrics['MQ'] = hl.sqrt(info.RAW_MQ / info.MQ_DP)
    if 'QUALapprox' in info and 'VarDP'in info:
        additional_metrics['QD'] = info.QUALapprox / info.VarDP
    if 'SB' in info:
        additional_metrics['FS'] = fs_from_sb(hl.array(info.SB))

    return info.annotate(**additional_metrics)
