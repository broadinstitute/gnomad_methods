from gnomad_hail.resources.resource_utils import (
    DataException,
    PedigreeResource,
    TableResource,
    VersionedTableResource,
)

DATA_TYPES = ["exomes", "genomes"]

CURRENT_EXOME_RELEASE = "2.1.1"
CURRENT_GENOME_RELEASE = "2.1.1"

CURRENT_EXOME_META = "2018-10-11"
CURRENT_GENOME_META = "2018-10-11"

CURRENT_EXOME_FAM = "2018-04-12"
CURRENT_GENOME_FAM = "2018-04-12"

EXOME_RELEASES = ["2.1", "2.1.1"]
GENOME_RELEASES = ["2.1", "2.1.1"]

EXOME_META_DATES = [
    "2018-03-29",
    "2018-03-30",
    "2018-04-03",
    "2018-04-10",
    "2018-04-11",
    "2018-04-12",
    "2018-04-13",
    "2018-04-18",
    "2018-04-19",
    "2018-04-25",
    "2018-05-11",
    "2018-06-08",
    "2018-06-10",
    "2018-08-04",
    "2018-08-16",
    "2018-09-12",
    "2018-10-10",
    "2018-10-11",
]
GENOME_META_DATES = [
    "2018-03-29",
    "2018-03-30",
    "2018-04-03",
    "2018-04-10",
    "2018-04-11",
    "2018-04-12",
    "2018-04-13",
    "2018-04-18",
    "2018-04-19",
    "2018-04-25",
    "2018-05-11",
    "2018-06-08",
    "2018-06-10",
    "2018-08-04",
    "2018-08-16",
    "2018-09-12",
    "2018-10-10",
    "2018-10-11",
]

EXOME_FAM_DATES = [
    "2017-10-03",
    "2018-04-09",
    "2018-04-12",
    "2018-04-12.true_trios",
]
GEOME_FAM_DATES = [
    "2017-10-03",
    "2018-04-09",
    "2018-04-12",
    "2018-04-12.true_trios",
]

SUBPOPS = {
    "NFE": ["BGR", "EST", "NWE", "SEU", "SWE", "ONF"],
    "EAS": ["KOR", "JPN", "OEA"],
}
GENOME_POPS = ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH"]
EXOME_POPS = ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "SAS"]
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

POP_NAMES = {
    "oth": "Other",
    "afr": "African-American/African",
    "ami": "Amish",
    "amr": "Latino",
    "eas": "East Asian",
    "fin": "Finnish",
    "eur": "European",
    "nfe": "Non-Finnish European",
    "sas": "South Asian",
    "mde": "Middle Eastern",
    "asj": "Ashkenazi Jewish",
    "uniform": "Uniform",
    "sas_non_consang": "South Asian (F < 0.05)",
    "consanguineous": "South Asian (F > 0.05)",
    "exac": "ExAC",
    "bgr": "Bulgarian (Eastern European)",
    "deu": "German",
    "est": "Estonian",
    "esp": "Spanish",
    "gbr": "British",
    "nwe": "North-Western European",
    "seu": "Southern European",
    "ita": "Italian",
    "swe": "Swedish",
    "chn": "Chinese",
    "kor": "Korean",
    "hkg": "Hong Kong",
    "sgp": "Singaporean",
    "twn": "Taiwanese",
    "jpn": "Japanese",
    "oea": "Other East Asian",
    "oeu": "Other European",
    "onf": "Other Non-Finnish European",
    "unk": "Unknown",
}


def _public_release_ht_path(data_type: str, version: str) -> str:
    """
    Get public release table path

    :param str data_type: One of "exomes" or "genomes"
    :param str version: One of the release versions of gnomAD on GRCh37
    :return: Path to release Table
    :rtype: str
    """
    return f"gs://gnomad-public/release/{version}/ht/{data_type}/gnomad.{data_type}.r{version}.sites.ht"


def _public_coverage_ht_path(data_type: str, version: str) -> str:
    """
    Get public coverage hail table

    :param str data_type: One of "exomes" or "genomes"
    :param str version: One of the release versions of gnomAD on GRCh37
    :return: path to coverage Table
    :rtype: str
    """
    return f"gs://gnomad-public/release/{version}/coverage/{data_type}/gnomad.{data_type}.r{version}.coverage.ht"


def _public_pca_ht_path(subpop: str) -> str:
    """
    Get public pca loadings path

    :param subpop: Can be empty ("") -> global, "eas" or "nfe"
    :return: Path to release Table
    :rtype: str
    """
    subpop = f".{subpop}" if subpop else ""
    return f"gs://gnomad-public/release/2.1/pca/gnomad.r2.1.pca_loadings{subpop}.ht"


def _metadata_ht_path(data_type: str, date: str) -> str:
    """
    Metadata ht path for supplied date

    :param str data_type: One of "exomes" or "genomes"
    :param str date: One of the dates in data_types meta date array
    :return: Path to metadata
    :rtype: str
    """
    return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.metadata.{date}.ht"


def _fam_path(data_type: str, date: str, true_trios=False) -> str:
    """
    Pedigree path for supplied date

    :param str data_type: One of "exomes" or "genomes"
    :param str date: Date from fam file generation
    :param bool true_trios: Whether to include only true trios
    :return:
    """
    true_trios = ".true_trios" if true_trios else ""
    return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.metadata.{date}{true_trios}.ht"


def _liftover_data_path(data_type: str, version: str) -> str:
    """
    Paths to liftover gnomAD Table

    :param str data_type: One of `exomes` or `genomes`
    :param str version: One of the release versions of gnomAD on GRCh37
    :return: Path to chosen Table
    :rtype: str
    """
    return f"gs://gnomad-public/release/{version}/liftover_grch38/ht/{data_type}/gnomad.{data_type}.r{version}.sites.liftover_grch38.ht"


def public_release(data_type: str) -> VersionedTableResource:
    """
    Retrieves publicly released versioned table resource

    :param str data_type: One of "exomes" or "genomes"
    :return: Release Table
    :rtype: VersionedTableResource
    """

    if data_type not in DATA_TYPES:
        raise DataException(f'{data_type} not in {DATA_TYPES}')

    if data_type == "exomes":
        current_release = CURRENT_EXOME_RELEASE
        releases = EXOME_RELEASES
    else:
        current_release = CURRENT_GENOME_RELEASE
        releases = GENOME_RELEASES

    return VersionedTableResource(
        current_release,
        {release: TableResource(path=_public_release_ht_path(data_type, release)) for release in releases},
    )


def coverage(data_type: str) -> VersionedTableResource:
    """
    Retrieves gnomAD's coverage table by data_type

    :param str data_type: One of "exomes" or "genomes"
    :return: Coverage Table
    :rtype: VersionedTableResource
    """
    if data_type not in DATA_TYPES:
        raise DataException(f'{data_type} not in {DATA_TYPES}')

    if data_type == "exomes":
        current_release = "2.1"
        releases = EXOME_RELEASES
        releases.remove("2.1.1")
    else:
        current_release = CURRENT_GENOME_RELEASE
        releases = GENOME_RELEASES

    return VersionedTableResource(
        current_release,
        {release: TableResource(path=_public_coverage_ht_path(data_type, release)) for release in releases},
    )


def liftover(data_type: str) -> VersionedTableResource:
    """
    Get the 38 liftover of gnomad v2.1.1

    :param str data_type: One of "exomes" or "genomes"
    :return: Release Table
    :rtype: VersionedTableResource
    """
    if data_type not in DATA_TYPES:
        raise DataException(f'{data_type} not in {DATA_TYPES}')

    if data_type == "exomes":
        current_release = CURRENT_EXOME_RELEASE
        releases = EXOME_RELEASES
        releases.remove("2.1")
    else:
        current_release = CURRENT_GENOME_RELEASE
        releases = GENOME_RELEASES

    return VersionedTableResource(
        current_release,
        {release: TableResource(path=_public_coverage_ht_path(data_type, release)) for release in releases},
    )


def public_pca_loadings(subpop: str = "") -> TableResource:
    """
    Returns the TableResource containing sites and loadings from population PCA

    :param str subpop: Can be empty ("") -> global, "eas" or "nfe"
    :return: gnomAD public PCA loadings TableResource
    :rtype: TableResource
    """
    if subpop not in ["", "eas", "nfe"]:
        raise DataException(
            'Available subpops are "eas" or "nfe", default value "" for global'
        )

    return TableResource(path=_public_pca_ht_path(subpop))


def release_vcf_path(data_type: str, version: str, contig: str) -> str:
    """
    Publically released VCF. Provide specific contig, i.e. "20", to retrieve contig
    specific VCF
    :param str data_type: One of "exomes" or "genomes"
    :param str version: One of the release versions of gnomAD on GRCh37
    :param str contig: Single contig "1" to "Y"
    :return: Path to VCF
    :rtype: str
    """
    contig = f".{contig}" if contig else ""
    return f"gs://gnomad-public/release/{version}/vcf/{data_type}/gnomad.{data_type}.r{version}.sites{contig}.vcf.bgz"


def metadata(data_type: str) -> VersionedTableResource:
    """
    Sample metadata associated with gnomAD release

    :param str data_type: One of "exomes" or "genomes"
    :return: Versioned metadata Tables
    :rtype: VersionedTableResource
    """
    if data_type not in DATA_TYPES:
        raise DataException(f'{data_type} not in {DATA_TYPES}')

    if data_type == "exomes":
        current_meta = CURRENT_EXOME_META
        dates = EXOME_META_DATES
    else:
        current_meta = CURRENT_GENOME_META
        dates = GENOME_META_DATES

    return VersionedTableResource(
        current_meta,
        {date: TableResource(path=_metadata_ht_path(data_type, date)) for date in dates},
    )


def fam(data_type: str, true_trios: bool = False) -> PedigreeResource:
    """
    Returns the path to gnomAD pedigree file.

    :param str data_type: One of 'exomes' or 'genomes'
    :param bool true_trios: If set, removes all families with more than one offspring
    :return: Path to fam file
    :rtype: str
    """
    if data_type not in DATA_TYPES:
        raise DataException(f'{data_type} not in {DATA_TYPES}')

    if data_type == "exomes":
        current_fam = CURRENT_EXOME_FAM
    else:
        current_fam = CURRENT_GENOME_FAM

    return PedigreeResource(path=_fam_path(data_type, current_fam, true_trios))


def vep_config_path():
    return "gs://hail-common/vep/vep/vep85-loftee-gcloud.json"
