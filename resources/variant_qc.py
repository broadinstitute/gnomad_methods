from resources.basics import *


def rf_run_hash_path(data_type: str, hail_version: str = CURRENT_HAIL_VERSION):
    """
    Returns the path to the json file containing the RF runs list.

    :param str data_type: One of 'exomes' or 'genomes'
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to json file
    :rtype: str
    """

    return 'gs://gnomad/annotations/hail-{0}/mt/{1}/rf/gnomad.{1}.runs.json'.format(hail_version, data_type)


def rf_annotated_path(
        data_type: str,
        adj: bool = False,
        hail_version: str = CURRENT_HAIL_VERSION
) -> str:
    """
    Returns the path to the RF-ready annotated HT

    :param str data_type: One of 'exomes' or 'genomes'
    :param bool adj: Whether to load 'adj' or 'raw'
    :param hail_version: One of the HAIL_VERSIONs
    :return:
    """
    return 'gs://gnomad/annotations/hail-{0}/mt/{1}/rf/gnomad.{1}.{2}.ht'.format(
        hail_version,
        data_type,
        "adj" if adj else 'raw')


def rf_path(data_type: str,
            data: str = 'rf_result',
            run_hash: str = None,
            hail_version: str = CURRENT_HAIL_VERSION
            ) -> str:
    """

    Gets the path to the desired RF data.
    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering

    :param str data_type: One of 'exomes' or 'genomes'
    :param str data: One of 'pre_rf', 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :param str hail_version: One of the HAIL_VERSIONs
    :return:
    """

    extension = 'model' if data == 'model' else 'ht'
    return 'gs://gnomad/annotations/hail-{0}/mt/{1}/rf/{2}/gnomad.{1}.{2}.{3}.{4}'.format(hail_version, data_type, run_hash, data, extension)
