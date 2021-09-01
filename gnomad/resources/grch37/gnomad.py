# noqa: D100

from gnomad.resources.resource_utils import (
    DataException,
    GnomadPublicTableResource,
    VersionedTableResource,
)

DATA_TYPES = ["exomes", "genomes"]

CURRENT_EXOME_RELEASE = "2.1.1"
CURRENT_GENOME_RELEASE = "2.1.1"

EXOME_RELEASES = ["2.1", "2.1.1"]
GENOME_RELEASES = ["2.1", "2.1.1"]

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
    Get public release table path.

    :param data_type: One of "exomes" or "genomes"
    :param version: One of the release versions of gnomAD on GRCh37
    :return: Path to release Table
    """
    return f"gs://gnomad-public-requester-pays/release/{version}/ht/{data_type}/gnomad.{data_type}.r{version}.sites.ht"


def _public_coverage_ht_path(data_type: str, version: str) -> str:
    """
    Get public coverage hail table.

    :param data_type: One of "exomes" or "genomes"
    :param version: One of the release versions of gnomAD on GRCh37
    :return: path to coverage Table
    """
    return f"gs://gnomad-public-requester-pays/release/{version}/coverage/{data_type}/gnomad.{data_type}.r{version}.coverage.ht"


def _public_pca_ht_path(subpop: str) -> str:
    """
    Get public pca loadings path.

    :param subpop: Can be empty ("") -> global, "eas" or "nfe"
    :return: Path to release Table
    """
    subpop = f".{subpop}" if subpop else ""
    return f"gs://gnomad-public-requester-pays/release/2.1/pca/gnomad.r2.1.pca_loadings{subpop}.ht"


def _liftover_data_path(data_type: str, version: str) -> str:
    """
    Paths to liftover gnomAD Table.

    :param data_type: One of `exomes` or `genomes`
    :param version: One of the release versions of gnomAD on GRCh37
    :return: Path to chosen Table
    """
    return f"gs://gnomad-public-requester-pays/release/{version}/liftover_grch38/ht/{data_type}/gnomad.{data_type}.r{version}.sites.liftover_grch38.ht"


def public_release(data_type: str) -> VersionedTableResource:
    """
    Retrieve publicly released versioned table resource.

    :param data_type: One of "exomes" or "genomes"
    :return: Release Table
    """
    if data_type not in DATA_TYPES:
        raise DataException(f"{data_type} not in {DATA_TYPES}")

    if data_type == "exomes":
        current_release = CURRENT_EXOME_RELEASE
        releases = EXOME_RELEASES
    else:
        current_release = CURRENT_GENOME_RELEASE
        releases = GENOME_RELEASES

    return VersionedTableResource(
        current_release,
        {
            release: GnomadPublicTableResource(
                path=_public_release_ht_path(data_type, release)
            )
            for release in releases
        },
    )


def coverage(data_type: str) -> VersionedTableResource:
    """
    Retrieve gnomAD's coverage table by data_type.

    :param data_type: One of "exomes" or "genomes"
    :return: Coverage Table
    """
    if data_type not in DATA_TYPES:
        raise DataException(f"{data_type} not in {DATA_TYPES}")

    if data_type == "exomes":
        current_release = "2.1"
        releases = [r for r in EXOME_RELEASES if r != "2.1.1"]
    else:
        current_release = "2.1"
        releases = [r for r in GENOME_RELEASES if r != "2.1.1"]

    return VersionedTableResource(
        current_release,
        {
            release: GnomadPublicTableResource(
                path=_public_coverage_ht_path(data_type, release)
            )
            for release in releases
        },
    )


def liftover(data_type: str) -> VersionedTableResource:
    """
    Get the 38 liftover of gnomad v2.1.1.

    :param data_type: One of "exomes" or "genomes"
    :return: Release Table
    """
    if data_type not in DATA_TYPES:
        raise DataException(f"{data_type} not in {DATA_TYPES}")

    if data_type == "exomes":
        current_release = CURRENT_EXOME_RELEASE
        releases = [r for r in EXOME_RELEASES if r != "2.1"]
    else:
        current_release = CURRENT_GENOME_RELEASE
        releases = [r for r in GENOME_RELEASES if r != "2.1"]

    return VersionedTableResource(
        current_release,
        {
            release: GnomadPublicTableResource(
                path=_liftover_data_path(data_type, release)
            )
            for release in releases
        },
    )


def public_pca_loadings(subpop: str = "") -> GnomadPublicTableResource:
    """
    Return the TableResource containing sites and loadings from population PCA.

    :param subpop: Can be empty ("") -> global, "eas" or "nfe"
    :return: gnomAD public PCA loadings TableResource
    """
    if subpop not in ["", "eas", "nfe"]:
        raise DataException(
            'Available subpops are "eas" or "nfe", default value "" for global'
        )

    return GnomadPublicTableResource(path=_public_pca_ht_path(subpop))


def release_vcf_path(data_type: str, version: str, contig: str) -> str:
    """
    Publically released VCF. Provide specific contig, i.e. "20", to retrieve contig specific VCF.

    :param data_type: One of "exomes" or "genomes"
    :param version: One of the release versions of gnomAD on GRCh37
    :param contig: Single contig "1" to "Y"
    :return: Path to VCF
    """
    if not version.startswith("2"):
        raise DataException(
            f"gnomAD version {version} is not available on reference genome GRCh37"
        )

    contig = f".{contig}" if contig else ""
    return f"gs://gcp-public-data--gnomad/release/{version}/vcf/{data_type}/gnomad.{data_type}.r{version}.sites{contig}.vcf.bgz"
