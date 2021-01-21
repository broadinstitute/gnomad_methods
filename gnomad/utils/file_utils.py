import base64
import gzip
import logging
import os
import subprocess
import uuid
from typing import List, Optional, Tuple, Union

import hail as hl

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
    if fext in [".ht", ".mt"]:
        fname += "/_SUCCESS"
    if fname.startswith("gs://"):
        return hl.hadoop_exists(fname)
    else:
        return os.path.isfile(fname)


def write_temp_gcs(
    t: Union[hl.MatrixTable, hl.Table],
    gcs_path: str,
    overwrite: bool = False,
    temp_path: Optional[str] = None,
) -> None:
    if not temp_path:
        temp_path = f"/tmp_{uuid.uuid4()}.h"
    t.write(temp_path, overwrite=True)
    t = (
        hl.read_matrix_table(temp_path)
        if isinstance(t, hl.MatrixTable)
        else hl.read_table(temp_path)
    )
    t.write(gcs_path, overwrite=overwrite)


def select_primitives_from_ht(ht: hl.Table) -> hl.Table:
    """
    Select only primitive types (string, int, float, bool) from a Table.
    Particularly useful for exporting a Table.

    :param ht: Input Table
    :return: Table with only primitive types selected
    """
    return ht.select(
        **{
            x: v
            for x, v in ht.row_value.items()
            if v.dtype
            in {hl.tstr, hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64, hl.tbool}
        }
    )


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
    if input_file_path.startswith("gs://"):
        hl.hadoop_copy(input_file_path, "file:///" + input_file_path.split("/")[-1])
        f = (
            gzip.open("/" + os.path.basename(input_file_path))
            if input_file_path.endswith("gz")
            else open("/" + os.path.basename(input_file_path))
        )
    else:
        f = (
            gzip.open(input_file_path)
            if input_file_path.endswith("gz")
            else open(input_file_path)
        )
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output
