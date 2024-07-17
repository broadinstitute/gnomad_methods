# noqa: D100

import base64
import gzip
import logging
import os
import subprocess
import uuid
from typing import List, Optional, Tuple, Union

import hail as hl

from gnomad.resources.resource_utils import DataException

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def file_exists(fname: str) -> bool:
    """
    Check whether a file exists.

    Supports either local or Google cloud (gs://) paths.
    If the file is a Hail file (.ht, .mt, .bm, .parquet, .he, and .vds extensions), it
    checks that _SUCCESS is present.

    :param fname: File name.
    :return: Whether the file exists.
    """
    fext = os.path.splitext(fname)[1]
    if fext in {".ht", ".mt", ".bm", ".parquet", ".he"}:
        paths = [f"{fname}/_SUCCESS"]
    elif fext == ".vds":
        paths = [f"{fname}/reference_data/_SUCCESS", f"{fname}/variant_data/_SUCCESS"]
    else:
        paths = [fname]

    if fname.startswith("gs://"):
        exists_func = hl.hadoop_exists
    else:
        exists_func = os.path.isfile

    exists = all([exists_func(p) for p in paths])

    return exists


def check_file_exists_raise_error(
    fname: Union[str, List[str]],
    error_if_exists: bool = False,
    error_if_not_exists: bool = False,
    error_if_exists_msg: str = "The following files already exist: ",
    error_if_not_exists_msg: str = "The following files do not exist: ",
) -> bool:
    """
    Check whether the file or all files in a list of files exist and optionally raise an exception.

    This can be useful when writing out to files at the end of a pipeline to first check if the file already
    exists and therefore requires the file to be removed or overwrite specified so the pipeline doesn't fail.

    :param fname: File path, or list of file paths to check the existence of.
    :param error_if_exists: Whether to raise an exception if any of the files exist. Default is True.
    :param error_if_not_exists: Whether to raise an exception if any of the files do not exist. Default is False.
    :param error_if_exists_msg: String of the error message to print if any of the files exist.
    :param error_if_not_exists_msg: String of the error message to print if any of the files do not exist.
    :return: Boolean indicating if `fname` or all files in `fname` exist.
    """
    if isinstance(fname, str):
        fname = [fname]

    all_exist = True
    exist = []
    not_exist = []
    for f in fname:
        exists = file_exists(f)
        all_exist &= exists
        if exists and error_if_exists:
            exist.append(f)
        if not exists and error_if_not_exists:
            not_exist.append(f)

    error_msg = ""
    if exist:
        error_msg = error_if_exists_msg + ", ".join(exist)
    if not_exist:
        error_msg = error_msg + "\n" + error_if_not_exists_msg + ", ".join(not_exist)
    if error_msg:
        raise DataException(error_msg)

    return all_exist


def write_temp_gcs(
    t: Union[hl.MatrixTable, hl.Table],
    gcs_path: str,
    overwrite: bool = False,
    temp_path: Optional[str] = None,
) -> None:  # noqa: D103
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


def get_file_stats(url: str, project_id: Optional[str] = None) -> Tuple[int, str, str]:
    """
    Get size (as both int and str) and md5 for file at specified URL.

    Typically used to get stats on VCFs.

    :param url: Path to file of interest.
    :param project_id: Google project ID. Specify if URL points to a requester-pays bucket.
    :return: Tuple of file size and md5.
    """
    one_gibibyte = 2**30
    one_mebibyte = 2**20

    if project_id:
        output = subprocess.check_output(
            ["gsutil", "-u", project_id, "stat", url]
        ).decode("utf8")
    else:
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
    Read a file input into a python list (each line will be an element).

    Supports Google storage paths and .gz compression.

    :param input_file_path: File path
    :return: List of lines
    """
    if input_file_path.startswith("gs://"):
        hl.hadoop_copy(input_file_path, "file:///" + input_file_path.split("/")[-1])
        f = (
            gzip.open("/" + os.path.basename(input_file_path), encoding="utf-8")
            if input_file_path.endswith("gz")
            else open("/" + os.path.basename(input_file_path), encoding="utf-8")
        )
    else:
        f = (
            gzip.open(input_file_path, encoding="utf-8")
            if input_file_path.endswith("gz")
            else open(input_file_path, encoding="utf-8")
        )
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def repartition_for_join(
    ht_path: str,
    new_partition_percent: float = 1.1,
) -> List[hl.expr.IntervalExpression]:
    """
    Calculate new partition intervals using input Table.

    Reading in all Tables using the same partition intervals (via
    `_intervals`) before they are joined makes the joins much more efficient.
    For more information, see:
    https://discuss.hail.is/t/room-for-improvement-when-joining-multiple-hts/2278/8

    :param ht_path: Path to Table to use for interval partition calculation.
    :param new_partition_percent: Percent of initial dataset partitions to use.
        Value should be greater than 1 so that input Table will have more
        partitions for the join. Defaults to 1.1.
    :return: List of IntervalExpressions calculated over new set of partitions
        (number of partitions in HT * desired percent increase).
    """
    ht = hl.read_table(ht_path)
    if new_partition_percent < 1:
        logger.warning(
            "new_partition_percent value is less than 1! The new HT will have fewer"
            " partitions than the original HT!"
        )
    return ht._calculate_new_partitions(ht.n_partitions() * new_partition_percent)


def create_vds(
    gvcfs: str,
    output_path: str,
    temp_path: str,
    save_path: Optional[str] = None,
    use_genome_default_intervals: bool = False,
    use_exome_default_intervals: bool = False,
    intervals: Optional[str] = None,
    gvcf_batch_size: Optional[int] = None,
    reference_genome: str = "GRCh38",
) -> hl.vds.VariantDataset:
    """
    Combine GVCFs into a single VDS.

    :param gvcfs: Path to file containing GVCF paths with no header.
    :param output_path: Path to write output VDS.
    :param temp_path: Directory path to write temporary files. A bucket with a life-cycle
        policy is recommended.
    :param save_path: Path to write combiner to on failure. Can be used to restart
        combiner from a failed state. If not specified, defaults to temp_path +
        combiner_plan.json.
    :param use_genome_default_intervals: Use the default genome intervals.
    :param use_exome_default_intervals: Use the default exome intervals.
    :param intervals: Path to text file with intervals to use for VDS creation.
    :param gvcf_batch_size: Number of GVCFs to combine into a Variant Dataset at once.
    :param reference_genome: Reference genome to use. Default is GRCh38.
    :return: Combined VDS.
    """
    if not save_path and temp_path:
        save_path = temp_path + "combiner_plan.json"

    gvcfs = read_list_data(gvcfs)
    intervals = (
        hl.import_locus_intervals(
            intervals, reference_genome=reference_genome
        ).interval.collect()
        if intervals
        else None
    )

    if not len(gvcfs) > 0:
        raise DataException("No GVCFs provided in file")

    if intervals and not len(intervals) > 0:
        raise DataException("No intervals provided in passed intervals file")

    logger.info("Combining %s GVCFs into a single VDS", len(gvcfs))
    combiner = hl.vds.new_combiner(
        output_path=output_path,
        temp_path=temp_path,
        save_path=save_path,
        gvcf_paths=gvcfs,
        use_genome_default_intervals=use_genome_default_intervals,
        use_exome_default_intervals=use_exome_default_intervals,
        intervals=intervals,
        gvcf_batch_size=gvcf_batch_size,
    )
    combiner.run()
    vds = hl.vds.read_vds(output_path)
    return vds
