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

EXOME_COVERAGE_RELEASES = ["2.1"]
GENOME_COVERAGE_RELEASES = ["2.1"]

CURRENT_EXOME_COVERAGE_RELEASE = "2.1"
CURRENT_GENOME_COVERAGE_RELEASE = "2.1"

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


def _public_constraint_ht_path() -> str:
    """
    Get public gene constraint Table path.

    :return: Path to constraint Table.
    """
    return "gs://gnomad-public-requester-pays/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.ht"


def _public_pext_path(pext_type: str = "base_level") -> str:
    """
    Get public proportion expressed across transcripts (pext) data.

    :param pext_type: One of "annotation_level" or "base_level". Default is "base_level".
    :return: Path to pext data.
    :raises DataException: If the provided pext_type is invalid.
    """
    pext_paths = {
        "annotation_level": "gs://gnomad-public-requester-pays/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht",
        "base_level": "gs://gnomad-public-requester-pays/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.ht",
    }

    if pext_type not in pext_paths:
        valid_types = list(pext_paths.keys())
        raise DataException(
            f"Invalid pext_type: '{pext_type}'. Valid options are {valid_types}."
        )

    return pext_paths[pext_type]


def _public_browser_gene_ht_path() -> str:
    """
    Get public browser gene table path.

    .. note::

       This table has smaller number of partitions (n=100) for faster computation and
       contains pext data compared to gnomad.genes.GRCh37.GENCODEv19.ht (which was
       used by the browser for ES export) under the same path.

    :return: Path to browser gene Table.
    """
    return "gs://gnomad-public-requester-pays/resources/grch37/browser/gnomad.genes.GRCh37.GENCODEv19.pext.ht"


def _public_mnv_path(distance: str = "1") -> str:
    """
    Get path to public multinucleotide variant (MNV) data.

    :param distance: Distance between two SNVs in MNV. Default is '1'.
    :return: Path to MNV data.
    :raises DataException: If the provided distance is invalid (outside 1-10).
    """
    distances = set(map(str, range(1, 11)))
    if distance not in distances:
        raise DataException(
            f"Invalid distance: '{distance}'. Valid options are {distances}."
        )
    return f"gs://gnomad-public-requester-pays/release/2.1/mnv/gnomad_mnv_genome_d{distance}.ht"


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
        current_release = CURRENT_EXOME_COVERAGE_RELEASE
        releases = EXOME_COVERAGE_RELEASES
    else:
        current_release = CURRENT_GENOME_COVERAGE_RELEASE
        releases = GENOME_COVERAGE_RELEASES

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


def pext(pext_type: str = "base_level") -> GnomadPublicTableResource:
    """
    Retrieve proportion expressed across transcripts (pext) data.

    :param pext_type: One of "annotation_level" or "base_level". Default is "base_level".
    :return: Pext Table.
    """
    return GnomadPublicTableResource(path=_public_pext_path(pext_type))


def constraint() -> GnomadPublicTableResource:
    """
    Retrieve gene constraint table.

    :return: Gene constraint Table.
    """
    return GnomadPublicTableResource(path=_public_constraint_ht_path())


def mnv(distance: str = "1") -> GnomadPublicTableResource:
    """
    Retrieve multinucleotide variant table.

    :param distance: Distance between two SNVs in MNV. Default is '1'.
    :return: MNV Table.
    """
    return GnomadPublicTableResource(path=_public_mnv_path(distance))


def browser_gene() -> GnomadPublicTableResource:
    """
    Retrieve browser gene table.

    :return: Browser gene Table.
    """
    return GnomadPublicTableResource(path=_public_browser_gene_ht_path())
