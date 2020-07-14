from gnomad.resources.resource_utils import TableResource, BlockMatrixResource
from gnomad.resources.grch37.gnomad import CURRENT_EXOME_RELEASE, CURRENT_GENOME_RELEASE
from typing import Optional


def _ld_matrix_path(
    data_type: str,
    pop: str,
    common_only: bool = True,
    adj: bool = True,
    version: Optional[str] = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    subdir = "sv/" if data_type == "genomes_snv_sv" else ""
    return f'gs://gnomad-public-requester-pays/release/{version}/ld/{subdir}gnomad.{data_type}.r{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.bm'


def _ld_index_path(
    data_type: str,
    pop: str,
    common_only: bool = True,
    adj: bool = True,
    version: Optional[str] = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    subdir = "sv/" if data_type == "genomes_snv_sv" else ""
    return f'gs://gnomad-public-requester-pays/release/{version}/ld/{subdir}gnomad.{data_type}.r{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.variant_indices.ht'


def _ld_snv_sv_path(pop):
    return f"gs://gnomad-public-requester-pays/release/2.1.1/ld/sv/gnomad.genomes_snv_sv.r2.1.1.{pop}.snv_sv.ld.ht"


def _ld_snv_sv_index_path(pop, type):
    return f"gs://gnomad-public-requester-pays/release/2.1.1/ld/sv/gnomad.genomes_snv_sv.r2.1.1.{pop}.snv_sv.ld.{type}.txt.bgz"


def _cross_pop_ld_scores_path(
    data_type: str,
    pop1: str,
    pop2: str,
    adj: bool = True,
    version: Optional[str] = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-public-requester-pays/release/{version}/ld/scores/gnomad.{data_type}.r{version}.{pop1}.{pop2}.{"adj." if adj else ""}ld_scores.ht'


def _ld_scores_path(
    data_type: str, pop: str, adj: bool = True, version: Optional[str] = None
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-public-requester-pays/release/{version}/ld/scores/gnomad.{data_type}.r{version}.{pop}.{"adj." if adj else ""}ld_scores.ht'


def ld_matrix(pop: str) -> BlockMatrixResource:
    return BlockMatrixResource(path=_ld_matrix_path("genomes", pop))


def ld_index(pop: str) -> TableResource:
    return TableResource(path=_ld_index_path("genomes", pop))


def ld_scores(pop: str) -> TableResource:
    return TableResource(path=_ld_scores_path("genomes", pop))
