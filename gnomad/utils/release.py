from typing import Dict, List, Optional, Union

import hail as hl

from gnomad.resources.grch38.gnomad import (
    GROUPS,
    POPS,
    SEXES,
    SUBSETS,
)
from gnomad.utils.vcf import index_globals, make_label_combos


def make_faf_index_dict(faf_meta: List[Dict[str, str]]) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array.

    .. note: Filtering allele frequencies are only computed on 'adj'-filtered genotypes

    :param faf_meta: Global annotation containing the set of groupings for each element of the faf array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :return: Dictionary of faf annotation population groupings, where values are the corresponding 0-based indices for the
        groupings in the faf_meta array
    """

    index_dict = {
        **index_globals(faf_meta, dict(group=["adj"])),
        **index_globals(faf_meta, dict(group=["adj"], pop=POPS)),
        **index_globals(faf_meta, dict(group=["adj"], sex=SEXES)),
        **index_globals(faf_meta, dict(group=["adj"], pop=POPS, sex=SEXES)),
    }

    return index_dict


def make_freq_index_dict(
    freq_meta: List[Dict[str, str]], downsamplings: Optional[List[int]]
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array.

    .. note: Downsampling groupings are only computed on 'adj'-filtered genotypes

    :param freq_meta: Global annotation continaing the set of groupings for each element of the freq array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :param downsamplings: List of downsampling cohort sizes present in global frequency array
    :return: Dictionary keyed by the grouping combinations found in the frequency array, where values are the corresponding
        0-based indices for the groupings in the freq_meta array
    """

    index_dict = {
        **index_globals(freq_meta, dict(group=GROUPS)),
        **index_globals(freq_meta, dict(group=GROUPS, pop=POPS)),
        **index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)),
        **index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, pop=POPS)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, sex=SEXES)),
        **index_globals(
            freq_meta, dict(group=GROUPS, subset=SUBSETS, pop=POPS, sex=SEXES)
        ),
    }
    if downsamplings:
        index_dict.update(
            {
                **index_globals(
                    freq_meta, dict(downsampling=downsamplings, group=["adj"], pop=POPS)
                )
            }
        )

    return index_dict


