import base64
import gzip
import logging
import os
import subprocess
import uuid
from typing import List, Optional, Tuple, Union

import hail as hl

INFO_VCF_AS_PIPE_DELIMITED_FIELDS = ['AS_QUALapprox', 'AS_VarDP', 'AS_MQ_DP', 'AS_RAW_MQ', 'AS_SB_TABLE']

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def file_exists(fname: str) -> bool:
    """
    Check whether a file exists.
    Supports either local or Google cloud (gs://) paths.
    If the file is a Hail file (.ht, .mt extensions), it checks that _SUCCESS is present.

    :param fname: File name
    :return: Whether the file exists
    """
    fext = os.path.splitext(fname)[1]
    if fext in ['.ht', '.mt']:
        fname += '/_SUCCESS'
    if fname.startswith('gs://'):
        return hl.hadoop_exists(fname)
    else:
        return os.path.isfile(fname)


def write_temp_gcs(t: Union[hl.MatrixTable, hl.Table], gcs_path: str,
                   overwrite: bool = False, temp_path: Optional[str] = None) -> None:
    if not temp_path:
        temp_path = f'/tmp_{uuid.uuid4()}.h'
    t.write(temp_path, overwrite=True)
    t = hl.read_matrix_table(temp_path) if isinstance(t, hl.MatrixTable) else hl.read_table(temp_path)
    t.write(gcs_path, overwrite=overwrite)


def select_primitives_from_ht(ht: hl.Table) -> hl.Table:
    """
    Select only primitive types (string, int, float, bool) from a Table.
    Particularly useful for exporting a Table.

    :param ht: Input Table
    :return: Table with only primitive types selected
    """
    return ht.select(**{x: v for x, v in ht.row_value.items() if
                        v.dtype in {hl.tstr, hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64, hl.tbool}})


def ht_to_vcf_mt(
        info_ht: hl.Table,
        pipe_delimited_annotations  : List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS
 ) -> hl.MatrixTable:
    """
    Creates a MT ready for vcf export from a HT. In particular, the following conversions are done:
    - All int64 are coerced to int32
    - Fields specified by `pipe_delimited_annotations` will be converted from arrays to pipe-delimited strings

    .. note::

        The MT returned has no cols.

    :param info_ht: Input HT
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct)
    :return: MatrixTable ready for VCF export
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    info_ht = info_ht.key_by('locus', 'alleles')

    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in info_ht.info.dtype.items():
        if ft == hl.dtype('int64'):
            logger.warning(f"Coercing field info.{f} from int64 to int32 for VCF output. Value will be capped at int32 max value.")
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(**{f: hl.int32(hl.min(2**31 - 1, info_ht.info[f]))})
            )
        elif ft == hl.dtype('array<int64>'):
            logger.warning(f"Coercing field info.{f} from array<int64> to array<int32> for VCF output. Array values will be capped at int32 max value.")
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(**{f: info_ht.info[f].map(lambda x: hl.int32(hl.min(2**31 - 1, x)))})
            )

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, but it is overwritten below
    for f in pipe_delimited_annotations:
        if f in info_ht.info:
            info_expr[f] = "|" + get_pipe_expr(info_ht.info[f])

    # Flatten SB if it is an array of arrays
    if 'SB' in info_ht.info and not isinstance(info_ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr['SB'] = info_ht.info.SB[0].extend(info_ht.info.SB[1])

    if 'AS_SB_TABLE' in info_ht.info:
        info_expr['AS_SB_TABLE'] = get_pipe_expr(info_ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ",")))

    # Annotate with new expression and add 's' empty string field required to cast HT to MT
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(**info_expr),
        s=hl.null(hl.tstr)
    )

    # Create an MT with no cols so that we acn export to VCF
    info_mt = info_ht.to_matrix_table_row_major(columns=['s'], entry_field_name='s')
    return info_mt.filter_cols(False)


def rep_on_read(path: str, n_partitions: int) -> hl.MatrixTable:
    """
    Repartitions a MatrixTable on read. Currently the best way to increase the number of partitions in a MatrixTable.

    :param path: Path to input MatrixTable
    :param n_partitions: Number of desired partitions
    :return: MatrixTable with the number of desired partitions
    """
    mt = hl.read_matrix_table(path)
    intervals = mt._calculate_new_partitions(n_partitions)
    return hl.read_matrix_table(path, _intervals=intervals)



def get_file_stats(url: str) -> Tuple[int, str, str]:
    """
    Gets size (as both int and str) and md5 for file at specified URL.
    Typically used to get stats on VCFs.

    :param url: Path to file of interest.
    :return: Tuple of file size and md5.
    """
    one_gibibyte = 2 ** 30
    one_mebibyte = 2 ** 20

    output = subprocess.check_output(["gsutil", "stat", url]).decode("utf8")
    lines = output.split("\n")

    info = {}
    for line in lines:
        if not line:
            continue

        label, value = [s.strip() for s in line.split(":", 1)]
        if label == "Content-Length":
            size = int(value)
            if size >= one_gibibyte:
                info["size"] = f"{round(size / one_gibibyte, 2)} GiB"
            else:
                info["size"] = f"{round(size / one_mebibyte, 2)} MiB"

        if label == "Hash (md5)":
            info["md5"] = base64.b64decode(value).hex()

    return (size, info["size"], info["md5"])


def read_list_data(input_file_path: str) -> List[str]:
    """
    Reads a file input into a python list (each line will be an element).
    Supports Google storage paths and .gz compression.

    :param input_file_path: File path
    :return: List of lines
    """
    if input_file_path.startswith('gs://'):
        hl.hadoop_copy(input_file_path, 'file:///' + input_file_path.split("/")[-1])
        f = gzip.open("/" + os.path.basename(input_file_path)) if input_file_path.endswith('gz') else open("/" + os.path.basename(input_file_path))
    else:
        f = gzip.open(input_file_path) if input_file_path.endswith('gz') else open(input_file_path)
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def unphase_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Generate unphased version of MatrixTable (assumes call is in mt.GT and is diploid or haploid only)
    """
    return mt.annotate_entries(GT=hl.case()
                               .when(mt.GT.is_diploid(), hl.call(mt.GT[0], mt.GT[1], phased=False))
                               .when(mt.GT.is_haploid(), hl.call(mt.GT[0], phased=False))
                               .default(hl.null(hl.tcall))
                               )
