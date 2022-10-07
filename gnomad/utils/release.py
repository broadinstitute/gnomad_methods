# noqa: D100

from typing import Dict, List, Optional

from gnomad.resources.grch38.gnomad import GROUPS, POPS, SEXES, SUBSETS
from gnomad.utils.vcf import index_globals


def make_faf_index_dict(
    faf_meta: List[Dict[str, str]],
    groups: List[str] = ["adj"],
    pops: List[str] = POPS,
    sexes: List[str] = SEXES,
    label_delimiter: str = "_",
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array.

    :param faf_meta: Global annotation containing the set of groupings for each element of the faf array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :param groups: List of sample groups [adj, raw]. Default is GROUPS
    :param pops: List of sample global population names for gnomAD genomes. Default is POPS
    :param sexes: List of sample sexes used in VCF export. Default is SEXES
    :param label_delimiter: String used as delimiter when making group label combinations
    :return: Dictionary of faf annotation population groupings, where values are the corresponding 0-based indices for the
        groupings in the faf_meta array
    """

    def _get_index(label_groups):
        return index_globals(faf_meta, label_groups, label_delimiter)

    index_dict = {
        **_get_index(dict(group=groups)),
        **_get_index(dict(group=groups, pop=pops)),
        **_get_index(dict(group=groups, sex=sexes)),
        **_get_index(dict(group=groups, pop=pops, sex=sexes)),
    }
    return index_dict


def make_freq_index_dict(
    freq_meta: List[Dict[str, str]],
    groups: List[str] = GROUPS,
    pops: List[str] = POPS,
    sexes: List[str] = SEXES,
    subsets: List[str] = SUBSETS,
    downsamplings: Optional[List[int]] = None,
    label_delimiter: str = "_",
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array.

    .. note:

        Downsampling groupings are only computed on 'adj'-filtered genotypes

    :param freq_meta: List containing the set of groupings for each element of the freq array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :param groups: List of sample groups [adj, raw]. Default is GROUPS
    :param pops: List of sample global population names for gnomAD genomes. Default is POPS
    :param sexes: List of sample sexes used in VCF export. Default is SEXES
    :param subsets: List of sample subsets in dataset. Default is SUBSETS
    :param downsamplings: List of downsampling cohort sizes present in global frequency array
    :param label_delimiter: String used as delimiter when making group label combinations
    :return: Dictionary keyed by the grouping combinations found in the frequency array, where values are the corresponding
        0-based indices for the groupings in the freq_meta array
    """

    def _get_index(label_groups):
        return index_globals(freq_meta, label_groups, label_delimiter)

    index_dict = {
        **_get_index(dict(group=groups)),
        **_get_index(dict(group=groups, pop=pops)),
        **_get_index(dict(group=groups, sex=sexes)),
        **_get_index(dict(group=groups, pop=pops, sex=sexes)),
        **_get_index(dict(group=groups, subset=subsets)),
        **_get_index(dict(group=groups, subset=subsets, pop=pops)),
        **_get_index(dict(group=groups, subset=subsets, sex=sexes)),
        **_get_index(dict(group=groups, subset=subsets, pop=pops, sex=sexes)),
    }

    if downsamplings:
        index_dict.update(
            {**_get_index(dict(downsampling=downsamplings, group=["adj"], pop=pops))}
        )

    return index_dict
