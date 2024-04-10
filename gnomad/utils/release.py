# noqa: D100

import logging
from typing import Dict, List, Optional

import hail as hl

from gnomad.resources.grch38.gnomad import (
    CURRENT_MAJOR_RELEASE,
    GROUPS,
    POPS,
    SEXES,
    SUBSETS,
)
from gnomad.utils.vcf import SORT_ORDER, index_globals

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def make_faf_index_dict(
    faf_meta: List[Dict[str, str]],
    groups: List[str] = ["adj"],
    pops: List[str] = POPS[CURRENT_MAJOR_RELEASE]["exomes"],
    sexes: List[str] = SEXES,
    label_delimiter: str = "_",
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array.

    :param faf_meta: Global annotation containing the set of groupings for each element of the faf array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :param groups: List of sample groups [adj, raw]. Default is GROUPS
    :param pops: List of sample global population names for gnomAD data type. Default is POPS[CURRENT_MAJOR_RELEASE]["exomes"].
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
    pops: List[str] = POPS[CURRENT_MAJOR_RELEASE]["exomes"],
    sexes: List[str] = SEXES,
    subsets: List[str] = SUBSETS[CURRENT_MAJOR_RELEASE],
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
    :param pops: List of sample global population names for gnomAD data type. Default is POPS[CURRENT_MAJOR_RELEASE]["exomes"].
    :param sexes: List of sample sexes used in VCF export. Default is SEXES
    :param subsets: List of sample subsets in dataset. Default is SUBSETS[CURRENT_MAJOR_RELEASE]
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


def make_freq_index_dict_from_meta(
    freq_meta: List[Dict[str, str]],
    label_delimiter: str = "_",
    sort_order: Optional[List[str]] = SORT_ORDER,
) -> Dict[str, int]:
    """
    Create a dictionary for accessing frequency array.

    The dictionary is keyed by the grouping combinations found in the frequency metadata
    array, where values are the corresponding 0-based indices for the groupings in the
    frequency array. For example, if the `freq_meta` entry [{'pop': 'nfe'}, {'sex': 'XX'}]
    corresponds to the 5th entry in the frequency array, the returned dictionary entry
    would be {'nfe_XX': 4}.

    :param freq_meta: List of dictionaries containing frequency metadata.
    :param label_delimiter: Delimiter to use when joining frequency metadata labels.
    :param sort_order: List of frequency metadata labels to use when sorting the dictionary.
    :return: Dictionary of frequency metadata.
    """
    # Confirm all groups in freq_meta are in sort_order. Warn user if not.
    if sort_order is not None:
        diff = hl.eval(hl.set(freq_meta.flatmap(lambda i: i.keys()))) - set(sort_order)
        if diff:
            logger.warning(
                "Found unexpected frequency metadata groupings: %s. These groupings"
                " are not present in the provided sort_order: %s. These groupings"
                " will not be included in the returned dictionary.",
                diff,
                sort_order,
            )

    index_dict = {}
    for i, f in enumerate(hl.eval(freq_meta)):
        if sort_order is None or len(set(f.keys()) - set(sort_order)) < 1:
            index_dict[
                label_delimiter.join(
                    [
                        f[g]
                        for g in sorted(
                            f.keys(),
                            key=(lambda x: sort_order.index(x)) if sort_order else None,
                        )
                    ]
                )
            ] = i

    return index_dict
