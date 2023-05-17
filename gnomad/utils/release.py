# noqa: D100
import copy
from typing import Dict, List, Optional

from gnomad.resources.grch38.gnomad import (
    CURRENT_MAJOR_RELEASE,
    GROUPS,
    POPS,
    SEXES,
    SUBSETS,
)
from gnomad.utils.vcf import SORT_ORDER, index_globals


def make_faf_index_dict(
    faf_meta: List[Dict[str, str]],
    groups: List[str] = ["adj"],
    pops: List[str] = POPS[CURRENT_MAJOR_RELEASE],
    sexes: List[str] = SEXES,
    label_delimiter: str = "_",
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array.

    :param faf_meta: Global annotation containing the set of groupings for each element of the faf array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :param groups: List of sample groups [adj, raw]. Default is GROUPS
    :param pops: List of sample global population names for gnomAD genomes. Default is POPS[CURRENT_MAJOR_RELEASE]
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
    pops: List[str] = POPS[CURRENT_MAJOR_RELEASE],
    sexes: List[str] = SEXES,
    subsets: List[str] = SUBSETS[CURRENT_MAJOR_RELEASE],
    downsamplings: Optional[List[int]] = None,
    label_delimiter: str = "_",
    sort_order: List[str] = SORT_ORDER,
    additional_strata: Optional[List[Dict[str, List[str]]]] = None,
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array.

    .. note:

        Downsampling groupings are only computed on 'adj'-filtered genotypes

    :param freq_meta: List containing the set of groupings for each element of the freq array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :param groups: List of sample groups [adj, raw]. Default is GROUPS
    :param pops: List of sample global population names for gnomAD genomes. Default is POPS[CURRENT_MAJOR_RELEASE]
    :param sexes: List of sample sexes used in VCF export. Default is SEXES
    :param subsets: List of sample subsets in dataset. Default is SUBSETS[CURRENT_MAJOR_RELEASE]
    :param downsamplings: List of downsampling cohort sizes present in global frequency array
    :param label_delimiter: String used as delimiter when making group label combinations
    :param sort_order: List of strings specifying the order to sort subgroupings in frequency dictionary.
    :param additional_strata: Optional List of additional strata as dictionaries to include in the index dictionary.
        Key is strata, value is list of strata values.
    :return: Dictionary keyed by the grouping combinations found in the frequency array, where values are the corresponding
        0-based indices for the groupings in the freq_meta array
    """
    sort_order = copy.deepcopy(sort_order)

    def _get_index(label_groups):
        return index_globals(freq_meta, label_groups, label_delimiter, sort_order)

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

    if additional_strata:
        for strata in additional_strata:
            for k, v in strata.items():
                # Note: Tack the new strata onto the end of the sort order so the labels
                # can be made
                if k not in sort_order:
                    sort_order.append(
                        k
                    )  # TODO: Should we add the additional strata before group, aka adj? sort_order.insert(-1,k)
            index_dict.update({**_get_index(dict(group=groups, **strata))})

    return index_dict
