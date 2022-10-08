# noqa: D100

import logging
import random
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

import hail as hl

from gnomad.utils.annotations import annotate_adj

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

UNRELATED = "unrelated"
"""
String representation for a pair of unrelated individuals in this module.
Typically >2nd degree relatives, but the threshold is user-dependant.
"""

SECOND_DEGREE_RELATIVES = "second degree relatives"
"""
String representation for a pair of 2nd degree relatives in this module.
"""

PARENT_CHILD = "parent-child"
"""
String representation for a parent-child pair in this module.
"""

SIBLINGS = "siblings"
"""
String representation for a sibling pair in this module.
"""

DUPLICATE_OR_TWINS = "duplicate/twins"
"""
String representation for a pair of samples who are identical (either MZ twins of duplicate) in this module.
"""

AMBIGUOUS_RELATIONSHIP = "ambiguous"
"""
String representation for a pair of samples whose relationship is ambiguous.
This is used in the case of a pair of samples which kinship/IBD values do not correspond to any biological relationship between two individuals.
"""


def get_duplicated_samples(
    relationship_ht: hl.Table,
    i_col: str = "i",
    j_col: str = "j",
    rel_col: str = "relationship",
) -> List[Set[str]]:
    """
    Extract the list of duplicate samples using a Table ouput from pc_relate.

    :param relationship_ht: Table with relationships between pairs of samples
    :param i_col: Column containing the 1st sample
    :param j_col: Column containing the 2nd sample
    :param rel_col: Column containing the sample pair relationship annotated with get_relationship_expr
    :return: List of sets of samples that are duplicates
    """

    def get_all_dups(
        s: str, dups: Set[str], samples_duplicates: Dict[str, Set[str]]
    ) -> Tuple[Set[str], Dict[str, Set[str]]]:
        """
        Create the set of all duplicated samples corresponding to `s` that are found in `sample_duplicates`.

        Also return the remaining sample duplicates after removing all duplicated samples corresponding to `s`.

        Works by recursively adding duplicated samples to the set.

        :param s: sample to identify duplicates for
        :param dups: set of corresponding samples already identified
        :param samples_duplicates: dict of sample -> duplicate-pair left to assign
        :return: (set of duplicates corresponding to s found in samples_duplicates, remaining samples_duplicates)
        """
        if s in samples_duplicates:
            dups.add(s)
            s_dups = samples_duplicates.pop(s)
            for s_dup in s_dups:
                if s_dup not in dups:
                    dups, samples_duplicates = get_all_dups(
                        s_dup, dups, samples_duplicates
                    )
        return dups, samples_duplicates

    logger.info("Computing duplicate sets")
    dup_pairs = relationship_ht.aggregate(
        hl.agg.filter(
            relationship_ht[rel_col] == DUPLICATE_OR_TWINS,
            hl.agg.collect(hl.tuple([relationship_ht[i_col], relationship_ht[j_col]])),
        )
    )

    samples_duplicates = defaultdict(set)
    for i, j in dup_pairs:
        samples_duplicates[i].add(j)
        samples_duplicates[j].add(i)

    duplicated_samples = []
    while len(samples_duplicates) > 0:
        dup_set, samples_duplicates = get_all_dups(
            list(samples_duplicates)[0], set(), samples_duplicates
        )
        duplicated_samples.append(dup_set)

    return duplicated_samples


def get_duplicated_samples_ht(
    duplicated_samples: List[Set[str]],
    samples_rankings_ht: hl.Table,
    rank_ann: str = "rank",
):
    """
    Create a HT with duplicated samples sets.

    Each row is indexed by the sample that is kept and also contains the set of duplicate samples that should be filtered.

    `samples_rankings_ht` is a HT containing a global rank for each of the samples (smaller is better).

    :param duplicated_samples: List of sets of duplicated samples
    :param samples_rankings_ht: HT with global rank for each sample
    :param rank_ann: Annotation in `samples_ranking_ht` containing each sample global rank (smaller is better).
    :return: HT with duplicate sample sets, including which to keep/filter
    """
    dups_ht = hl.Table.parallelize(
        [
            hl.struct(dup_set=i, dups=duplicated_samples[i])
            for i in range(0, len(duplicated_samples))
        ]
    )
    dups_ht = dups_ht.explode(dups_ht.dups, name="_dup")
    dups_ht = dups_ht.key_by("_dup")
    dups_ht = dups_ht.annotate(rank=samples_rankings_ht[dups_ht.key][rank_ann])
    dups_cols = hl.bind(
        lambda x: hl.struct(kept=x[0], filtered=x[1:]),
        hl.sorted(
            hl.agg.collect(hl.tuple([dups_ht._dup, dups_ht.rank])), key=lambda x: x[1]
        ).map(lambda x: x[0]),
    )
    dups_ht = dups_ht.group_by(dups_ht.dup_set).aggregate(**dups_cols)

    if isinstance(dups_ht.kept, hl.expr.StructExpression):
        dups_ht = dups_ht.key_by(**dups_ht.kept).drop("kept")
    else:
        dups_ht = dups_ht.key_by(
            s=dups_ht.kept
        )  # Since there is no defined name in the case of a non-struct type, use `s`
    return dups_ht


def explode_duplicate_samples_ht(dups_ht: hl.Table) -> hl.Table:
    """
    Explode the result of `get_duplicated_samples_ht`, so that each line contains a single sample.

    An additional annotation is added: `dup_filtered` indicating which of the duplicated samples was kept.
    Requires a field `filtered` which type should be the same as the input duplicated samples Table key.

    :param dups_ht: Input HT
    :return: Flattened HT
    """

    def get_dups_to_keep_expr():
        if dups_ht.filtered.dtype.element_type == dups_ht.key.dtype:
            return (dups_ht.key, False)
        elif (len(dups_ht.key) == 1) & (
            dups_ht.filtered.dtype.element_type == dups_ht.key[0].dtype
        ):
            return (dups_ht.key[0], False)
        else:
            raise TypeError(
                "Cannot explode table as types of the filtered field"
                f" ({dups_ht.filtered.dtype}) and the key ({dups_ht.key.dtype}) are"
                " incompatible."
            )

    dups_ht = dups_ht.annotate(
        dups=hl.array([get_dups_to_keep_expr()]).extend(
            dups_ht.filtered.map(lambda x: (x, True))
        )
    )
    dups_ht = dups_ht.explode("dups")
    dups_ht = dups_ht.key_by()
    return dups_ht.select(s=dups_ht.dups[0], dup_filtered=dups_ht.dups[1]).key_by("s")


def get_relationship_expr(  # TODO: The threshold detection could be easily automated by fitting distributions over the data.
    kin_expr: hl.expr.NumericExpression,
    ibd0_expr: hl.expr.NumericExpression,
    ibd1_expr: hl.expr.NumericExpression,
    ibd2_expr: hl.expr.NumericExpression,
    first_degree_kin_thresholds: Tuple[float, float] = (0.19, 0.4),
    second_degree_min_kin: float = 0.1,
    ibd0_0_max: float = 0.025,
    ibd0_25_thresholds: Tuple[float, float] = (0.1, 0.425),
    # ibd0_50_thresholds = [0.37, 0.625], Not useful for relationship inference
    # ibd0_100_threshold = 0.625  , Not useful for relationship inference
    ibd1_0_thresholds: Tuple[float, float] = (-0.15, 0.1),
    # ibd1_25_thresholds: Tuple[float, float] = (0.1, 0.37), Not useful for
    # relationship inference
    ibd1_50_thresholds: Tuple[float, float] = (0.275, 0.75),
    ibd1_100_min: float = 0.75,
    ibd2_0_max: float = 0.125,
    ibd2_25_thresholds: Tuple[float, float] = (0.1, 0.5),
    ibd2_100_thresholds: Tuple[float, float] = (0.75, 1.25),
) -> hl.expr.StringExpression:
    """
    Return an expression indicating the relationship between a pair of samples given their kin coefficient and IBDO, IBD1, IBD2 values.

    The kinship coefficient values in the defaults are in line with those output from
    `hail.methods.pc_relate <https://hail.is/docs/0.2/methods/genetics.html?highlight=pc_relate#hail.methods.pc_relate>`.

    :param kin_expr: Kin coefficient expression
    :param ibd0_expr: IBDO expression
    :param ibd1_expr: IBD1 expression
    :param ibd2_expr: IDB2 expression
    :param first_degree_kin_thresholds: (min, max) kinship threshold for 1st degree relatives
    :param second_degree_min_kin: min kinship threshold for 2nd degree relatives
    :param ibd0_0_max: max IBD0 threshold for 0 IBD0 sharing
    :param ibd0_25_thresholds: (min, max) thresholds for 0.25 IBD0 sharing
    :param ibd1_0_thresholds: (min, max) thresholds for 0 IBD1 sharing. Note that the min is there because pc_relate can output large negative values in some corner cases.
    :param ibd1_50_thresholds: (min, max) thresholds for 0.5 IBD1 sharing
    :param ibd1_100_min: min IBD1 threshold for 1.0 IBD1 sharing
    :param ibd2_0_max: max IBD2 threshold for 0 IBD2 sharing
    :param ibd2_25_thresholds: (min, max) thresholds for 0.25 IBD2 sharing
    :param ibd2_100_thresholds: (min, max) thresholds for 1.00 IBD2 sharing. Note that the min is there because pc_relate can output much larger IBD2 values in some corner cases.
    :return: The relationship annotation using the constants defined in this module.
    """
    return (
        hl.case()
        .when(kin_expr < second_degree_min_kin, UNRELATED)
        .when((kin_expr < first_degree_kin_thresholds[0]), SECOND_DEGREE_RELATIVES)
        .when(
            (kin_expr < first_degree_kin_thresholds[1])
            & (ibd0_expr <= ibd0_0_max)
            & (ibd1_expr >= ibd1_100_min)
            & (ibd2_expr <= ibd2_0_max),
            PARENT_CHILD,
        )
        .when(
            (kin_expr < first_degree_kin_thresholds[1])
            & (ibd0_expr >= ibd0_25_thresholds[0])
            & (ibd0_expr <= ibd0_25_thresholds[1])
            & (ibd1_expr >= ibd1_50_thresholds[0])
            & (ibd1_expr <= ibd1_50_thresholds[1])
            & (ibd2_expr >= ibd2_25_thresholds[0])
            & (ibd2_expr <= ibd2_25_thresholds[1]),
            SIBLINGS,
        )
        .when(
            (kin_expr > first_degree_kin_thresholds[1])
            & (ibd0_expr < ibd0_0_max)
            & (ibd1_expr >= ibd1_0_thresholds[0])
            & (ibd1_expr <= ibd1_0_thresholds[1])
            & (ibd2_expr >= ibd2_100_thresholds[0])
            & (ibd2_expr <= ibd2_100_thresholds[1]),
            DUPLICATE_OR_TWINS,
        )
        .default(AMBIGUOUS_RELATIONSHIP)
    )


def infer_families(
    relationship_ht: hl.Table,
    sex: Union[hl.Table, Dict[str, bool]],
    duplicate_samples_ht: hl.Table,
    i_col: str = "i",
    j_col: str = "j",
    relationship_col: str = "relationship",
) -> hl.Pedigree:
    """
    Generate a pedigree containing trios inferred from the `relationship_ht`.

    This function takes a hail Table with a row for each pair of related individuals i, j in the data (it's OK to have
    unrelated samples too).

    The `relationship_col` should be a column specifying the relationship between each two samples as defined in this
     module's constants.

    This function returns a pedigree containing trios inferred from the data. Family ID can be the same for multiple
    trios if one or more members of the trios are related (e.g. sibs, multi-generational family). Trios are ordered by family ID.

    .. note::

        This function only returns complete trios defined as: one child, one father and one mother (sex is required for both parents).

    :param relationship_ht: Input relationship table
    :param sex: A Table or dict giving the sex for each sample (`TRUE`=female, `FALSE`=male). If a Table is given, it should have a field `is_female`.
    :param duplicated_samples: All duplicated samples TO REMOVE (If not provided, this function won't work as it assumes that each child has exactly two parents)
    :param i_col: Column containing the 1st sample of the pair in the relationship table
    :param j_col: Column containing the 2nd sample of the pair in the relationship table
    :param relationship_col: Column contatining the relationship for the sample pair as defined in this module constants.
    :return: Pedigree of complete trios
    """

    def group_parent_child_pairs_by_fam(
        parent_child_pairs: Iterable[Tuple[str, str]]
    ) -> List[List[Tuple[str, str]]]:
        """
        Group parent-child pairs into a list of families.

        A family here is defined as a list of sample-pairs which all share at least one sample with at least one other
        sample-pair in the list.

        :param parent_child_pairs: All the parent-children pairs
        :return: A list of families, where each element of the list is a list of the parent-children pairs
        """
        fam_id = 1  # stores the current family id
        s_fam = dict()  # stores the family id for each sample
        fams = defaultdict(list)  # stores fam_id -> sample-pairs
        for pair in parent_child_pairs:
            if pair[0] in s_fam:
                if pair[1] in s_fam:
                    if (
                        s_fam[pair[0]] != s_fam[pair[1]]
                    ):  # If both samples are in different families, merge the families
                        new_fam_id = s_fam[pair[0]]
                        fam_id_to_merge = s_fam[pair[1]]
                        for s in s_fam:
                            if s_fam[s] == fam_id_to_merge:
                                s_fam[s] = new_fam_id
                        fams[new_fam_id].extend(fams.pop(fam_id_to_merge))
                else:  # If only the 1st sample in the pair is already in a family, assign the 2nd sample in the pair to the same family
                    s_fam[pair[1]] = s_fam[pair[0]]
                fams[s_fam[pair[0]]].append(pair)
            elif (
                pair[1] in s_fam
            ):  # If only the 2nd sample in the pair is already in a family, assign the 1st sample in the pair to the same family
                s_fam[pair[0]] = s_fam[pair[1]]
                fams[s_fam[pair[1]]].append(pair)
            else:  # If none of the samples in the pair is already in a family, create a new family
                s_fam[pair[0]] = fam_id
                s_fam[pair[1]] = fam_id
                fams[fam_id].append(pair)
                fam_id += 1

        return list(fams.values())

    def get_trios(
        fam_id: str,
        parent_child_pairs: List[Tuple[str, str]],
        related_pairs: Dict[Tuple[str, str], str],
    ) -> List[hl.Trio]:
        """
        Generate trios based on the list of parent-child pairs in the family and all related pairs in the data.

        Only complete parent/offspring trios are included in the results.

        The trios are assembled as follows:
        1. All pairs of unrelated samples with different sexes within the family are extracted as possible parent pairs
        2. For each possible parent pair, a list of all children is constructed (each child in the list has a parent-offspring pair with each parent)
        3. If there are multiple children for a given parent pair, all children should be siblings with each other
        4. Check that each child was only assigned a single pair of parents. If a child is found to have multiple parent pairs, they are ALL discarded.

        :param fam_id: The family ID
        :param parent_child_pairs: The parent-child pairs for this family
        :param related_pairs: All related sample pairs in the data
        :return: List of trios in the family
        """

        def get_possible_parents(samples: List[str]) -> List[Tuple[str, str]]:
            """
            Return all pairs of unrelated samples with different sexes within the family are extracted as possible parent pairs.

            :param samples: All samples in the family
            :return: Possible parent pairs
            """
            possible_parents = []
            for i in range(len(samples)):
                for j in range(i + 1, len(samples)):
                    if (
                        related_pairs.get(tuple(sorted([samples[i], samples[j]])))
                        is None
                    ):
                        if sex.get(samples[i]) is False and sex.get(samples[j]) is True:
                            possible_parents.append((samples[i], samples[j]))
                        elif (
                            sex.get(samples[i]) is True and sex.get(samples[j]) is False
                        ):
                            possible_parents.append((samples[j], samples[i]))
            return possible_parents

        def get_children(possible_parents: Tuple[str, str]) -> List[str]:
            """
            Construct a list of all children for a given possible parent pair.

            Each child in the list has a parent-offspring pair with each parent.

            :param possible_parents: A pair of possible parents
            :return: The list of all children (if any) corresponding to the possible parents
            """
            possible_offsprings = defaultdict(
                set
            )  # stores sample -> set of parents in the possible_parents where (sample, parent) is found in possible_child_pairs
            for pair in parent_child_pairs:
                if possible_parents[0] == pair[0]:
                    possible_offsprings[pair[1]].add(possible_parents[0])
                elif possible_parents[0] == pair[1]:
                    possible_offsprings[pair[0]].add(possible_parents[0])
                elif possible_parents[1] == pair[0]:
                    possible_offsprings[pair[1]].add(possible_parents[1])
                elif possible_parents[1] == pair[1]:
                    possible_offsprings[pair[0]].add(possible_parents[1])

            return [
                s for s, parents in possible_offsprings.items() if len(parents) == 2
            ]

        def check_sibs(children: List[str]) -> bool:
            """
            Confirm that all children of a parent pair are siblings with each other.

            If there are multiple children for a given parent pair, all children should be siblings with each other.

            :param children: List of all children for a given parent pair
            :return: Whether all children in the list are siblings
            """
            for i in range(len(children)):
                for j in range(i + 1, len(children)):
                    if (
                        related_pairs[tuple(sorted([children[i], children[j]]))]
                        != SIBLINGS
                    ):
                        return False
            return True

        def discard_multi_parents_children(trios: List[hl.Trio]):
            """
            Check that each child was only assigned a single pair of parents.

            If a child is found to have multiple parent pairs, they are ALL discarded.

            :param trios: All trios formed for this family
            :return: The list of trios for which each child has a single parents pair.
            """
            children_trios = defaultdict(list)
            for trio in trios:
                children_trios[trio.s].append(trio)

            for s, s_trios in children_trios.items():
                if len(s_trios) > 1:
                    logger.warning(
                        "Discarded duplicated child %s found multiple in trios: %s",
                        s,
                        ", ".join([str(trio) for trio in s_trios]),
                    )

            return [trios[0] for trios in children_trios.values() if len(trios) == 1]

        # Get all possible pairs of parents in (father, mother) order
        all_possible_parents = get_possible_parents(
            list({s for pair in parent_child_pairs for s in pair})
        )

        trios = []
        for possible_parents in all_possible_parents:
            children = get_children(possible_parents)
            if check_sibs(children):
                trios.extend(
                    [
                        hl.Trio(
                            s=s,
                            fam_id=fam_id,
                            pat_id=possible_parents[0],
                            mat_id=possible_parents[1],
                            is_female=sex.get(s),
                        )
                        for s in children
                    ]
                )
            else:
                logger.warning(
                    "Discarded family with same parents, and multiple offspring that"
                    " weren't siblings:\nMother: %s\nFather:%s\nChildren:%s",
                    possible_parents[0],
                    possible_parents[1],
                    ", ".join(children),
                )

        return discard_multi_parents_children(trios)

    # Get all the relations we care about:
    # => Remove unrelateds and duplicates
    dups = duplicate_samples_ht.aggregate(
        hl.agg.explode(
            lambda dup: hl.agg.collect_as_set(dup), duplicate_samples_ht.filtered
        ),
        _localize=False,
    )
    relationship_ht = relationship_ht.filter(
        ~dups.contains(relationship_ht[i_col])
        & ~dups.contains(relationship_ht[j_col])
        & (relationship_ht[relationship_col] != UNRELATED)
    )

    # Check relatedness table format
    if not relationship_ht[i_col].dtype == relationship_ht[j_col].dtype:
        logger.error(
            "i_col and j_col of the relatedness table need to be of the same type."
        )

    # If i_col and j_col aren't str, then convert them
    if not isinstance(relationship_ht[i_col], hl.expr.StringExpression):
        logger.warning(
            "Pedigrees can only be constructed from string IDs, but your relatedness_ht"
            " ID column is of type: %s. Expression will be converted to string in"
            " Pedigrees.",
            relationship_ht[i_col].dtype,
        )
        if isinstance(relationship_ht[i_col], hl.expr.StructExpression):
            logger.warning(
                "Struct fields %s will be joined by underscores to use as sample names"
                " in Pedigree.",
                list(relationship_ht[i_col]),
            )
            relationship_ht = relationship_ht.key_by(
                **{
                    i_col: hl.delimit(
                        hl.array(
                            [
                                hl.str(relationship_ht[i_col][x])
                                for x in relationship_ht[i_col]
                            ]
                        ),
                        "_",
                    ),
                    j_col: hl.delimit(
                        hl.array(
                            [
                                hl.str(relationship_ht[j_col][x])
                                for x in relationship_ht[j_col]
                            ]
                        ),
                        "_",
                    ),
                }
            )
        else:
            raise NotImplementedError(
                "The `i_col` and `j_col` columns of the `relationship_ht` argument"
                " passed to infer_families are not of type StringExpression or Struct."
            )

    # If sex is a Table, extract sex information as a Dict
    if isinstance(sex, hl.Table):
        sex = dict(hl.tuple([sex.s, sex.is_female]).collect())

    # Collect all related sample pairs and
    # create a dictionnary with pairs as keys and relationships as values
    # Sample-pairs are tuples ordered by sample name
    related_pairs = {
        tuple(sorted([i, j])): rel
        for i, j, rel in hl.tuple(
            [relationship_ht.i, relationship_ht.j, relationship_ht.relationship]
        ).collect()
    }

    parent_child_pairs_by_fam = group_parent_child_pairs_by_fam(
        [pair for pair, rel in related_pairs.items() if rel == PARENT_CHILD]
    )
    return hl.Pedigree(
        [
            trio
            for fam_index, parent_child_pairs in enumerate(parent_child_pairs_by_fam)
            for trio in get_trios(str(fam_index), parent_child_pairs, related_pairs)
        ]
    )


def create_fake_pedigree(
    n: int,
    sample_list: List[str],
    exclude_real_probands: bool = False,
    max_tries: int = 10,
    real_pedigree: Optional[hl.Pedigree] = None,
) -> hl.Pedigree:
    """
    Generate a pedigree made of trios created by sampling 3 random samples in the sample list.

    - If `real_pedigree` is given, then children in the resulting fake trios will not include any trio with proband - parents
      that are in the real ones.
    - Each sample can be used only once as a proband in the resulting trios.
    - Sex of probands in fake trios is random.

    :param n: Number of fake trios desired in the pedigree
    :param sample_list: List of samples
    :param exclude_real_probands: If set, then fake trios probands cannot be in the real trios probands.
    :param max_tries: Maximum number of sampling to try before bailing out (preventing infinite loop if `n` is too large w.r.t. the number of samples)
    :param real_pedigree: Optional pedigree to exclude children from
    :return: Fake pedigree
    """
    real_trios = (
        {trio.s: trio for trio in real_pedigree.trios}
        if real_pedigree is not None
        else dict()
    )

    if exclude_real_probands and len(real_trios) == len(set(sample_list)):
        logger.warning(
            "All samples are in the real probands list; cannot create any fake"
            " pedigrees with exclude_real_probands=True. Returning an empty Pedigree."
        )
        return hl.Pedigree([])

    fake_trios = {}
    tries = 0
    while len(fake_trios) < n and tries < max_tries:
        s, mat_id, pat_id = random.sample(sample_list, 3)
        if (
            s in real_trios
            and (
                exclude_real_probands
                or {mat_id, pat_id} == {real_trios[s].mat_id, real_trios[s].pat_id}
            )
        ) or s in fake_trios:
            tries += 1
        else:
            tries = 0
            fake_trios[s] = hl.Trio(
                s=s,
                pat_id=pat_id,
                mat_id=mat_id,
                fam_id=f"fake_{str(len(fake_trios))}",
                is_female=bool(random.getrandbits(1)),
            )

    if tries == max_tries:
        logger.warning(
            "Only returning %d fake trios; random trio sampling stopped after reaching"
            " the maximum %d iterations",
            len(fake_trios),
            max_tries,
        )

    return hl.Pedigree(list(fake_trios.values()))


def compute_related_samples_to_drop(
    relatedness_ht: hl.Table,
    rank_ht: hl.Table,
    kin_threshold: float,
    filtered_samples: Optional[hl.expr.SetExpression] = None,
    min_related_hard_filter: Optional[int] = None,
) -> hl.Table:
    """
    Compute a Table with the list of samples to drop (and their global rank) to get the maximal independent set of unrelated samples.

    .. note::

        - `relatedness_ht` should be keyed by exactly two fields of the same type, identifying the pair of samples for each row.
        - `rank_ht` should be keyed by a single key of the same type as a single sample identifier in `relatedness_ht`.

    :param relatedness_ht: relatedness HT, as produced by e.g. pc-relate
    :param kin_threshold: Kinship threshold to consider two samples as related
    :param rank_ht: Table with a global rank for each sample (smaller is preferred)
    :param filtered_samples: An optional set of samples to exclude (e.g. these samples were hard-filtered)  These samples will then appear in the resulting samples to drop.
    :param min_related_hard_filter: If provided, any sample that is related to more samples than this parameter will be filtered prior to computing the maximal independent set and appear in the results.
    :return: A Table with the list of the samples to drop along with their rank.
    """
    # Make sure that the key types are valid
    assert len(list(relatedness_ht.key)) == 2
    assert relatedness_ht.key[0].dtype == relatedness_ht.key[1].dtype
    assert len(list(rank_ht.key)) == 1
    assert relatedness_ht.key[0].dtype == rank_ht.key[0].dtype

    logger.info("Filtering related samples using a kin threshold of %f", kin_threshold)
    relatedness_ht = relatedness_ht.filter(relatedness_ht.kin > kin_threshold)

    filtered_samples_rel = set()
    if min_related_hard_filter is not None:
        logger.info(
            "Computing samples related to too many individuals (>%d) for exclusion",
            min_related_hard_filter,
        )
        gbi = relatedness_ht.annotate(s=list(relatedness_ht.key))
        gbi = gbi.explode(gbi.s)
        gbi = gbi.group_by(gbi.s).aggregate(n=hl.agg.count())
        filtered_samples_rel = gbi.aggregate(
            hl.agg.filter(gbi.n > min_related_hard_filter, hl.agg.collect_as_set(gbi.s))
        )
        logger.info(
            "Found %d samples with too many 1st/2nd degree relatives. These samples"
            " will be excluded.",
            len(filtered_samples_rel),
        )

    if filtered_samples is not None:
        filtered_samples_rel = filtered_samples_rel.union(
            relatedness_ht.aggregate(
                hl.agg.explode(
                    lambda s: hl.agg.collect_as_set(s),
                    hl.array(list(relatedness_ht.key)).filter(
                        lambda s: filtered_samples.contains(s)
                    ),
                )
            )
        )

    if len(filtered_samples_rel) > 0:
        filtered_samples_lit = hl.literal(filtered_samples_rel)
        relatedness_ht = relatedness_ht.filter(
            filtered_samples_lit.contains(relatedness_ht.key[0])
            | filtered_samples_lit.contains(relatedness_ht.key[1]),
            keep=False,
        )

    logger.info("Annotating related sample pairs with rank.")
    i, j = list(relatedness_ht.key)
    relatedness_ht = relatedness_ht.key_by(s=relatedness_ht[i])
    relatedness_ht = relatedness_ht.annotate(
        **{i: hl.struct(s=relatedness_ht.s, rank=rank_ht[relatedness_ht.key].rank)}
    )
    relatedness_ht = relatedness_ht.key_by(s=relatedness_ht[j])
    relatedness_ht = relatedness_ht.annotate(
        **{j: hl.struct(s=relatedness_ht.s, rank=rank_ht[relatedness_ht.key].rank)}
    )
    relatedness_ht = relatedness_ht.key_by(i, j)
    relatedness_ht = relatedness_ht.drop("s")
    relatedness_ht = relatedness_ht.persist()

    related_samples_to_drop_ht = hl.maximal_independent_set(
        relatedness_ht[i],
        relatedness_ht[j],
        keep=False,
        tie_breaker=lambda l, r: l.rank - r.rank,
    )
    related_samples_to_drop_ht = related_samples_to_drop_ht.key_by()
    related_samples_to_drop_ht = related_samples_to_drop_ht.select(
        **related_samples_to_drop_ht.node
    )
    related_samples_to_drop_ht = related_samples_to_drop_ht.key_by("s")

    if len(filtered_samples_rel) > 0:
        related_samples_to_drop_ht = related_samples_to_drop_ht.union(
            hl.Table.parallelize(
                [hl.struct(s=s, rank=hl.null(hl.tint64)) for s in filtered_samples_rel],
                key="s",
            )
        )

    return related_samples_to_drop_ht


def filter_mt_to_trios(mt: hl.MatrixTable, fam_ht: hl.Table) -> hl.MatrixTable:
    """
    Filter a MatrixTable to a set of trios in `fam_ht` and annotates with adj.

    :param mt: A Matrix Table to filter to only trios
    :param fam_ht: A Table of trios to filter to, loaded using `hl.import_fam`
    :return: A MT filtered to trios and adj annotated
    """
    # Filter MT to samples present in any of the trios
    fam_ht = fam_ht.annotate(fam_members=[fam_ht.id, fam_ht.pat_id, fam_ht.mat_id])
    fam_ht = fam_ht.explode("fam_members", name="s")
    fam_ht = fam_ht.key_by("s").select().distinct()

    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.col_key]))
    if "adj" not in mt.entry:
        mt = annotate_adj(mt)

    return mt


def generate_trio_stats_expr(
    trio_mt: hl.MatrixTable,
    transmitted_strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    de_novo_strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    ac_strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    proband_is_female_expr: Optional[hl.expr.BooleanExpression] = None,
) -> hl.expr.StructExpression:
    """
    Generate a row-wise expression containing trio transmission stats.

    The expression will generate the following counts:
        - Number of alleles in het parents transmitted to the proband
        - Number of alleles in het parents not transmitted to the proband
        - Number of de novo mutations
        - Parent allele count
        - Proband allele count

    Transmission and de novo mutation metrics and allele counts can be stratified using additional filters.
    `transmitted_strata`, `de_novo_strata`, and `ac_strata` all expect a dictionary of filtering expressions keyed
    by their desired suffix to append for labeling. The default will perform counts using all genotypes and append
    'raw' to the label.

    .. note::

        Expects that `mt` is dense if dealing with a sparse MT `hl.experimental.densify` must be run first.

    :param trio_mt: A trio standard trio MT (with the format as produced by hail.methods.trio_matrix)
    :param transmitted_strata: Strata for the transmission counts
    :param de_novo_strata: Strata for the de novo counts
    :param ac_strata: Strata for the parent and child allele counts
    :param proband_is_female_expr: An optional expression giving the sex the proband. If not given, DNMs are only computed for autosomes.
    :return: An expression with the counts
    """
    # Create map for transmitted, untransmitted and DNM
    hom_ref = 0
    het = 1
    hom_var = 2

    auto_or_par = 2
    hemi_x = 1
    hemi_y = 0

    trans_config_counts = {
        # kid, dad, mom, copy -> t, u
        (hom_ref, het, het, auto_or_par): (0, 2),
        (hom_ref, hom_ref, het, auto_or_par): (0, 1),
        (hom_ref, het, hom_ref, auto_or_par): (0, 1),
        (het, het, het, auto_or_par): (1, 1),
        (het, hom_ref, het, auto_or_par): (1, 0),
        (het, het, hom_ref, auto_or_par): (1, 0),
        (het, hom_var, het, auto_or_par): (0, 1),
        (het, het, hom_var, auto_or_par): (0, 1),
        (hom_var, het, het, auto_or_par): (2, 0),
        (hom_var, het, hom_var, auto_or_par): (1, 0),
        (hom_var, hom_var, het, auto_or_par): (1, 0),
        (hom_ref, hom_ref, het, hemi_x): (0, 1),
        (hom_ref, hom_var, het, hemi_x): (0, 1),
        (hom_var, hom_ref, het, hemi_x): (1, 0),
        (hom_var, hom_var, het, hemi_x): (1, 0),
    }

    trans_count_map = hl.literal(trans_config_counts)

    def _get_copy_state(locus: hl.expr.LocusExpression) -> hl.expr.Int32Expression:
        """Get copy-state int from LocusExpression for indexing into trans_count_map."""
        return (
            hl.case()
            .when(locus.in_autosome_or_par(), auto_or_par)
            .when(locus.in_x_nonpar(), hemi_x)
            .when(locus.in_y_nonpar(), hemi_y)
            .or_missing()
        )

    def _is_dnm(
        proband_gt: hl.expr.CallExpression,
        father_gt: hl.expr.CallExpression,
        mother_gt: hl.expr.CallExpression,
        locus: hl.expr.LocusExpression,
        proband_is_female: Optional[hl.expr.BooleanExpression],
    ) -> hl.expr.BooleanExpression:
        """Determine whether a trio genotype combination is a DNM."""
        if proband_is_female is None:
            logger.warning(
                "Since no proband sex expression was given to generate_trio_stats_expr,"
                " only DNMs in autosomes will be counted."
            )
            return hl.or_missing(
                locus.in_autosome(),
                proband_gt.is_het() & father_gt.is_hom_ref() & mother_gt.is_hom_ref(),
            )
        return hl.cond(
            locus.in_autosome_or_par() | (proband_is_female & locus.in_x_nonpar()),
            proband_gt.is_het() & father_gt.is_hom_ref() & mother_gt.is_hom_ref(),
            hl.or_missing(
                ~proband_is_female, proband_gt.is_hom_var() & father_gt.is_hom_ref()
            ),
        )

    def _ac_an_parent_child_count(
        proband_gt: hl.expr.CallExpression,
        father_gt: hl.expr.CallExpression,
        mother_gt: hl.expr.CallExpression,
    ) -> Dict[str, hl.expr.Int64Expression]:
        """Get AC and AN for parents and children."""
        ac_parent_expr = hl.agg.sum(
            father_gt.n_alt_alleles() + mother_gt.n_alt_alleles()
        )
        an_parent_expr = hl.agg.sum(
            (hl.is_defined(father_gt) + hl.is_defined(mother_gt)) * 2
        )
        ac_child_expr = hl.agg.sum(proband_gt.n_alt_alleles())
        an_child_expr = hl.agg.sum(hl.is_defined(proband_gt) * 2)

        return {
            "ac_parents": ac_parent_expr,
            "an_parents": an_parent_expr,
            "ac_children": ac_child_expr,
            "an_children": an_child_expr,
        }

    # Create transmission counters
    trio_stats = hl.struct(
        **{
            f"{name2}_{name}": hl.agg.filter(
                (
                    trio_mt.proband_entry.GT.is_non_ref()
                    | trio_mt.father_entry.GT.is_non_ref()
                    | trio_mt.mother_entry.GT.is_non_ref()
                )
                & expr,
                hl.agg.sum(
                    trans_count_map.get(
                        (
                            trio_mt.proband_entry.GT.n_alt_alleles(),
                            trio_mt.father_entry.GT.n_alt_alleles(),
                            trio_mt.mother_entry.GT.n_alt_alleles(),
                            _get_copy_state(trio_mt.locus),
                        ),
                        default=(0, 0),
                    )[i]
                ),
            )
            for name, expr in transmitted_strata.items()
            for i, name2 in enumerate(["n_transmitted", "n_untransmitted"])
        }
    )

    # Create de novo counters
    trio_stats = trio_stats.annotate(
        **{
            f"n_de_novos_{name}": hl.agg.filter(
                _is_dnm(
                    trio_mt.proband_entry.GT,
                    trio_mt.father_entry.GT,
                    trio_mt.mother_entry.GT,
                    trio_mt.locus,
                    proband_is_female_expr,
                )
                & expr,
                hl.agg.count(),
            )
            for name, expr in de_novo_strata.items()
        }
    )

    trio_stats = trio_stats.annotate(
        **{
            f"{name2}_{name}": hl.agg.filter(
                expr,
                _ac_an_parent_child_count(
                    trio_mt.proband_entry.GT,
                    trio_mt.father_entry.GT,
                    trio_mt.mother_entry.GT,
                )[name2],
            )
            for name, expr in ac_strata.items()
            for name2 in ["ac_parents", "an_parents", "ac_children", "an_children"]
        }
    )

    return trio_stats


def generate_sib_stats_expr(
    mt: hl.MatrixTable,
    sib_ht: hl.Table,
    i_col: str = "i",
    j_col: str = "j",
    strata: Dict[str, hl.expr.BooleanExpression] = {"raw": True},
    is_female: Optional[hl.expr.BooleanExpression] = None,
) -> hl.expr.StructExpression:
    """
    Generate a row-wise expression containing the number of alternate alleles in common between sibling pairs.

    The sibling sharing counts can be stratified using additional filters using `stata`.

    .. note::

        This function expects that the `mt` has either been split or filtered to only bi-allelics
        If a sample has multiple sibling pairs, only one pair will be counted

    :param mt: Input matrix table
    :param sib_ht: Table defining sibling pairs with one sample in a col (`i_col`) and the second in another col (`j_col`)
    :param i_col: Column containing the 1st sample of the pair in the relationship table
    :param j_col: Column containing the 2nd sample of the pair in the relationship table
    :param strata: Dict with additional strata to use when computing shared sibling variant counts
    :param is_female: An optional column in mt giving the sample sex. If not given, counts are only computed for autosomes.
    :return: A Table with the sibling shared variant counts
    """

    def _get_alt_count(locus, gt, is_female):
        """Calculate alt allele count with sex info if present."""
        if is_female is None:
            return hl.or_missing(locus.in_autosome(), gt.n_alt_alleles())
        return (
            hl.case()
            .when(locus.in_autosome_or_par(), gt.n_alt_alleles())
            .when(
                ~is_female & (locus.in_x_nonpar() | locus.in_y_nonpar()),
                hl.min(1, gt.n_alt_alleles()),
            )
            .when(is_female & locus.in_y_nonpar(), 0)
            .default(0)
        )

    if is_female is None:
        logger.warning(
            "Since no sex expression was given to generate_sib_stats_expr, only"
            " variants in autosomes will be counted."
        )

    # If a sample is in sib_ht more than one time, keep only one of the sibling pairs
    # First filter to only samples found in mt to keep as many pairs as possible
    s_to_keep = mt.aggregate_cols(hl.agg.collect_as_set(mt.s), _localize=False)
    sib_ht = sib_ht.filter(
        s_to_keep.contains(sib_ht[i_col].s) & s_to_keep.contains(sib_ht[j_col].s)
    )
    sib_ht = sib_ht.add_index("sib_idx")
    sib_ht = sib_ht.annotate(sibs=[sib_ht[i_col].s, sib_ht[j_col].s])
    sib_ht = sib_ht.explode("sibs")
    sib_ht = sib_ht.group_by("sibs").aggregate(
        sib_idx=(hl.agg.take(sib_ht.sib_idx, 1, ordering=sib_ht.sib_idx)[0])
    )
    sib_ht = sib_ht.group_by(sib_ht.sib_idx).aggregate(sibs=hl.agg.collect(sib_ht.sibs))
    sib_ht = sib_ht.filter(hl.len(sib_ht.sibs) == 2).persist()

    logger.info(
        "Generating sibling variant sharing counts using %d pairs.", sib_ht.count()
    )
    sib_ht = sib_ht.explode("sibs").key_by("sibs")[mt.s]

    # Create sibling sharing counters
    sib_stats = hl.struct(
        **{
            f"n_sib_shared_variants_{name}": hl.sum(
                hl.agg.filter(
                    expr,
                    hl.agg.group_by(
                        sib_ht.sib_idx,
                        hl.or_missing(
                            hl.agg.sum(hl.is_defined(mt.GT)) == 2,
                            hl.agg.min(_get_alt_count(mt.locus, mt.GT, is_female)),
                        ),
                    ),
                ).values()
            )
            for name, expr in strata.items()
        }
    )

    sib_stats = sib_stats.annotate(
        **{
            f"ac_sibs_{name}": hl.agg.filter(
                expr & hl.is_defined(sib_ht.sib_idx), hl.agg.sum(mt.GT.n_alt_alleles())
            )
            for name, expr in strata.items()
        }
    )

    return sib_stats
