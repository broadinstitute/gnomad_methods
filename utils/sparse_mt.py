from .generic import *

INFO_SUM_AGG_FIELDS= ['RAW_MQ', 'QUALapprox']
INFO_INT32_SUM_AGG_FIELDS = ['DP', 'MQ_DP', 'VarDP']
INFO_MEDIAN_AGG_FIELDS = ['ReadPosRankSum', 'MQRankSum']
INFO_VCF_AS_PIPE_DELIMITED_FIELDS = ['AS_QUALapprox', 'AS_VarDP', 'AS_MQ_DP', 'AS_RAW_MQ', 'AS_SB_TABLE']

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
    t = mt._localize_entries('__entries', '__cols')
    return t.select(
        last_END_position=hl.or_else(
            hl.min(
                hl.scan.array_agg(
                    lambda entry: hl.scan._prev_nonnull(
                        hl.or_missing(
                            hl.is_defined(entry.END),
                            hl.tuple([
                                t.locus,
                                entry.END
                            ])
                        )
                    ),
                    t.__entries
                ).map(
                    lambda x: hl.or_missing(
                        (x[1] >= t.locus.position) & (x[0].contig == t.locus.contig),
                        x[0].position
                    )
                )
            ),
            t.locus.position
        )
    )


def _explode_info_array(
        info_array_expr: hl.expr.ArrayExpression,
        a_index: Optional[hl.expr.Int32Expression] = None
) -> hl.expr.StructExpression:
    """
    Explodes the info array created by annotate_info with option `explode_info` set to False.
    It expects `info_array_expr` to be an array with one entry per allele (ref included) and will return
    a struct expression with each info metric as an array with one entry per non-reference allele,
    except for `SB` (Strand bias table) for which the array has one entry per allele (ref included).

    If `a_index` is set, then the resulting struct expression will have one field per info metric with a single value
    corresponding to that of the allele at `a_index`,  except for `SB` that will be an array with two vales for [ref, a_index] resp.

    Note: The struct expression returned is what is typically found in VCF format.

    :param ArrayExpression info_array_expr: The info array expression to explode
    :param a_index: Optional non-ref allele to return info for
    :return: Exploded info Struct
    :rtype: StructExpression
    """
    fields = get_array_element_type(info_array_expr).fields

    if a_index is not None:
        if a_index < 1:
            raise ValueError("a_index in explode_info_array should be >1 (index of the desired non-reference allele in the MT alleles)")

        fields = {field: info_array_expr[a_index][field] for field in fields if field != 'SB' and not field.startswith('site_')}
        fields['SB'] = hl.array([info_array_expr[0].SB, info_array_expr[a_index].SB])

    else:
        fields = {field: info_array_expr[1:].map(lambda x: x[field]) for field in fields if field != 'SB' and not field.startswith('site_')}
        fields['SB'] = info_array_expr.map(lambda x: x.SB)
    return  hl.struct(**fields)


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


def get_as_info_expr(
        mt: hl.MatrixTable,
        sum_agg_fields: List[str] = INFO_SUM_AGG_FIELDS,
        int32_sum_agg_fields: List[str] = INFO_INT32_SUM_AGG_FIELDS,
        median_agg_fields: List[str] = INFO_MEDIAN_AGG_FIELDS
) -> hl.expr.StructExpression:
    """
    Returns an allele-specific annotation Struct containing typical VCF INFO fields from GVCF INFO fields stored in the MT entries.
    Note that while not a parameter, if `SB` is found, it will be aggregated in `AS_SB_TABLE`


    :param MatrixTable mt: Input Matrix Table
    :param list of str sum_agg_fields: Fields to aggregate that should be summed.
    :param list of str int32_sum_agg_fields: Fields to aggregate that should be summed as int32.
    :param list of str median_agg_fields: Fields to aggregate that should be aggregated using median.
    :return:
    """

    # Map str to expressions where needed
    if isinstance(sum_agg_fields, list):
        sum_agg_fields = _agg_list_to_dict(mt, sum_agg_fields)

    if isinstance(sum_agg_fields, list):
        int32_sum_agg_fields = _agg_list_to_dict(mt, int32_sum_agg_fields)

    if isinstance(median_agg_fields, list):
        median_agg_fields = _agg_list_to_dict(mt, median_agg_fields)

    # Create aggregators
    agg_expr = {}
    sb_expr = mt['SB'] if 'SB' in mt.entry else mt['gvcf_info']['SB'] if 'gvcf_info' in mt and 'SB' in 'gvcf_info' else None
    if sb_expr is not None:
        agg_expr = {'AS_SB_TABLE': hl.agg.array_agg(lambda x: hl.agg.sum(x), sb_expr[2:]).map(lambda x: hl.int32(x))}

    agg_expr.update({
                    f'AS_{k}': hl.agg.approx_quantiles(expr, 0.5)
                    for k, expr in median_agg_fields.items()
                })
    agg_expr.update({
                    f'AS_{k}': hl.agg.sum(expr)
                    for k, expr in sum_agg_fields.items()
                })
    agg_expr.update({
                    f'AS_{k}': hl.int32(hl.agg.sum(expr))
                    for k, expr in int32_sum_agg_fields.items()
                })


    # Aggregate values per allele
    info_expr = hl.agg.array_agg(
        lambda ai: hl.agg.filter(
            mt.LA.contains(ai),
            hl.struct(
                **agg_expr
            )
        ),
        hl.range(1, hl.len(mt.alleles))
    )

    # TODO: check whether bind can  be replaced by annotate. Also, separate 3 expressions so that they can be input-dependent
    info_expr = hl.bind(
            lambda info_arr: hl.range(0, hl.len(info_arr)).map(
                lambda ai: info_arr[ai].annotate(
                    AS_MQ=hl.or_missing(ai > 0, hl.sqrt(info_arr[ai].RAW_MQ / info_arr[ai].MQ_DP)),
                    AS_QD=hl.or_missing(ai > 0, info_arr[ai].QUALapprox / info_arr[ai].VarDP),
                    AS_FS=hl.or_missing(ai > 0, fs_from_sb(hl.array([info_arr[0].SB, info_arr[ai].SB])))
                )
            ),
            info_expr
        )


    #  Add reference allele annotations (really just `SB`)
    info_element_type = get_array_element_type(info_expr)
    info_expr = hl.array([
        hl.struct(
            AS_SB_TABLE=hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.SB[:2]).map(lambda x: hl.int32(x)),
            **{f: hl.null(info_element_type.get(f)) for f in info_element_type.fields if f != 'AS_SB_TABLE'}
        )
    ]).extend(info_expr)



    info_expr = hl.bind(
                _explode_info_array,
                info_expr
            )

    return info_expr


def get_site_info_expr(
        mt: hl.MatrixTable,
        sum_agg_fields: Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_SUM_AGG_FIELDS,
        int32_sum_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_INT32_SUM_AGG_FIELDS,
        median_agg_fields:  Union[List[str], Dict[str, hl.expr.NumericExpression]] = INFO_MEDIAN_AGG_FIELDS,

) -> hl.expr.StructExpression:
    """
    Creates a site-level annotation Struct aggregating typical VCF INFO fields from GVCF INFO fields stored in the MT entries.

    If the fields to be aggregate (`sum_agg_fields`, `int32_sum_agg_fields`, `median_agg_fields`) are passed as list of str,
    then they should correspond to entry fields in `mt` or in `mt.gvcf_info`.
    Priority is given to entry fields in `mt` to those in `mt.gvcf_info` in case of a name clash.
    Note that `SB` will also be aggregated if present.

    :param MatrixTable mt: Input Matrix Table
    :param list of str or dict of str ->  NumericExpression sum_agg_fields: Fields to aggregate that should be summed.
    :param list of str or dict of str ->  NumericExpression int32_sum_agg_fields: Fields to aggregate that should be summed as int32.
    :param list of str or dict of str ->  NumericExpression median_agg_fields: Fields to aggregate that should be aggregated using median.
    :return:
    """

    # Map str to expressions if needed
    if isinstance(sum_agg_fields, list):
        sum_agg_fields = _agg_list_to_dict(mt, sum_agg_fields)

    if isinstance(int32_sum_agg_fields, list):
        int32_sum_agg_fields = _agg_list_to_dict(mt, int32_sum_agg_fields)

    if isinstance(median_agg_fields, list):
        median_agg_fields = _agg_list_to_dict(mt, median_agg_fields)


    # Create aggregators
    agg_expr = {}
    sb_expr = mt['SB'] if 'SB' in mt.entry else mt['gvcf_info']['SB'] if 'gvcf_info' in mt and 'SB' in 'gvcf_info' else None
    if sb_expr is not None:
        sb_agg = hl.agg.array_agg(lambda x: hl.agg.sum(x), sb_expr).map(lambda x: hl.int32(x))
        agg_expr = {'SB': [sb_agg[:2], sb_agg[2:]]}

    agg_expr.update({
                    k: hl.agg.approx_quantiles(expr, 0.5)
                    for k, expr in median_agg_fields.items()
                })
    agg_expr.update({
                    k: hl.agg.sum(expr)
                    for k, expr in sum_agg_fields.items()
                })
    agg_expr.update({
                    k: hl.int32(hl.agg.sum(expr))
                    for k, expr in int32_sum_agg_fields.items()
                })

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


def export_info_vcf(
        info_ht: hl.Table,
        info_vcf_path: str,
        pipe_delimited_annotations  : List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS
):

    info_expr = {}

    # Creates R-based, pipe-delimited AS strings
    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f in [f in info_ht.info]:

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, but it is overwritten below
    for f in pipe_delimited_annotations:
        if f in info_ht.info:
            info_expr[f] = "|" + get_pipe_expr(info_ht.info[f])

    # Take care of strand balance field
    if 'SB' in info_ht.info:
        info_expr['SB'] = info_ht.info.SB[0].extend(info_ht.info.SB[1])

    if 'AS_SB_TABLE' in info_ht.info:
        info_expr['AS_SB_TABLE'] = get_pipe_expr(info_ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ","))),

    # Annotate with new expression and add 's' empty string field required to cast HT to MT
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(**info_expr),
        s=hl.null(hl.tstr)
    )

    logger.info("Exporting VCF from HT with following schema:")
    info_ht.describe()

    # Create an MT with no cols so that we acn export to VCF
    info_mt = info_ht.to_matrix_table_row_major(columns=['s'], entry_field_name='s')
    info_mt = info_mt.filter_cols(False)
    hl.export_vcf(info_mt, info_vcf_path)