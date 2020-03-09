import hail as hl
import logging
from typing import Dict, List, Tuple, Set, Union, Iterable, Optional
from collections import defaultdict
import random


logger = logging.getLogger("gnomad.utils")


UNRELATED = 'Unrelated'
"""
String representation for a pair of unrelated individuals in this module.
Typically >2nd degree relatives, but the threshold is user-dependant.
"""

SECOND_DEGREE_RELATIVES = '2nd degree relatives'
"""
String representation for a pair of 2nd degree relatives in this module.
"""

PARENT_CHILD = 'Parent-child'
"""
String representation for a parent-child pair in this module.
"""

SIBLINGS = 'Siblings'
"""
String representation for a sibling pair in this module.
"""

DUPLICATE_OR_TWINS = 'Duplicate/twins'
"""
String representation for a pair of samples who are identical (either MZ twins of duplicate) in this module.
"""

AMBIGUOUS_RELATIONSHIP = 'Ambiguous'
"""
String representation for a pair of samples whose relationship is ambiguous.
This is used in the case of a pair of samples which kinship/IBD values do not correspond to any biological relationship between two individuals.
"""


def get_duplicated_samples(
        relationship_ht: hl.Table,
        i_col: str = 'i',
        j_col: str = 'j',
        rel_col: str = 'relationship'
) -> List[Set[str]]:
    """
    Given a pc_relate output Table, extract the list of duplicate samples. Returns a list of set of samples that are duplicates.
    :param relationship_ht: Table with relationships between pairs of samples
    :param i_col: Column containing the 1st sample
    :param j_col: Column containing the 2nd sample
    :param rel_col: Column containing the sample pair relationship annotated with get_relationship_expr
    :return: List of sets of samples that are duplicates
    """

    def get_all_dups(s: str, dups: Set[str], samples_duplicates: Dict[str, Set[str]]) -> Tuple[Set[str], Dict[str, Set[str]]]:
        """
        Creates the set of all duplicated samples corresponding to `s` that are found in `sample_duplicates`.
        Also returns the remaining sample duplicates after removing all duplicated samples corresponding to `s`.

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
                    dups, samples_duplicates = get_all_dups(s_dup, dups, samples_duplicates)
        return dups, samples_duplicates

    logger.info("Computing duplicate sets")
    dup_pairs = relationship_ht.aggregate(
        hl.agg.filter(
            relationship_ht[rel_col] == DUPLICATE_OR_TWINS,
            hl.agg.collect(hl.tuple([relationship_ht[i_col], relationship_ht[j_col]]))
        )
    )

    samples_duplicates = defaultdict(set)
    for i, j in dup_pairs:
        samples_duplicates[i].add(j)
        samples_duplicates[j].add(i)

    duplicated_samples = []
    while len(samples_duplicates) > 0:
        dup_set, samples_duplicates = get_all_dups(list(samples_duplicates)[0], set(), samples_duplicates)
        duplicated_samples.append(dup_set)

    return duplicated_samples


def get_duplicated_samples_ht(
        duplicated_samples: List[Set[str]],
        samples_rankings_ht: hl.Table,
        rank_ann: str = 'rank'
):
    """
    Creates a HT with duplicated samples sets.
    Each row is indexed by the sample that is kept and also contains the set of duplicate samples that should be filtered.

    `samples_rankings_ht` is a HT containing a global rank for each of the samples (smaller is better).

    :param duplicated_samples: List of sets of duplicated samples
    :param samples_rankings_ht: HT with global rank for each sample
    :param rank_ann: Annotation in `samples_ranking_ht` containing each sample global rank (smaller is better).
    :return: HT with duplicate sample sets, including which to keep/filter
    """
    dups_ht = hl.Table.parallelize([hl.struct(dup_set=i, dups=duplicated_samples[i]) for i in range(0, len(duplicated_samples))])
    dups_ht = dups_ht.explode(dups_ht.dups, name='_dup')
    dups_ht = dups_ht.key_by('_dup')
    dups_ht = dups_ht.annotate(rank=samples_rankings_ht[dups_ht.key][rank_ann])
    dups_cols = hl.bind(
        lambda x: hl.struct(
            kept=x[0],
            filtered=x[1:]
        ),
        hl.sorted(hl.agg.collect(hl.tuple([dups_ht._dup, dups_ht.rank])), key=lambda x: x[1]).map(lambda x: x[0])
    )
    dups_ht = dups_ht.group_by(dups_ht.dup_set).aggregate(
        **dups_cols
    )

    if isinstance(dups_ht.kept, hl.expr.StructExpression):
        dups_ht = dups_ht.key_by(**dups_ht.kept).drop('kept')
    else:
        dups_ht = dups_ht.key_by(s=dups_ht.kept)  # Since there is no defined name in the case of a non-struct type, use `s`
    return dups_ht


def explode_duplicate_samples_ht(dups_ht: hl.Table) -> hl.Table:
    """
    Flattens the result of `filter_duplicate_samples`, so that each line contains a single sample.
    An additional annotation is added: `dup_filtered` indicating which of the duplicated samples was kept.

    :param dups_ht: Input HT
    :return: Flattened HT
    """
    dups_ht = dups_ht.annotate(
        dups=hl.array([(dups_ht.key, False)]).extend(
            dups_ht.filtered.map(lambda x: (x, True))
        )
    )
    dups_ht = dups_ht.explode('dups')
    dups_ht = dups_ht.key_by()
    return dups_ht.select(s=dups_ht.dups[0], dup_filtered=dups_ht.dups[1]).key_by('s')


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
        # ibd1_25_thresholds: Tuple[float, float] = (0.1, 0.37), Not useful for relationship inference
        ibd1_50_thresholds: Tuple[float, float] = (0.275, 0.75),
        ibd1_100_min: float = 0.75,
        ibd2_0_max: float = 0.125,
        ibd2_25_thresholds: Tuple[float, float] = (0.1, 0.5),
        ibd2_100_thresholds: Tuple[float, float] = (0.75, 1.25)
) -> hl.expr.StringExpression:
    """
    Returns an expression that gives the relationship between a pair of samples given their kin coefficient and IBDO, IBD1, IBD2 values.
    The kinship coefficient values in the defaults are in line with those output from `hail.methods.pc_relate <https://hail.is/docs/0.2/methods/genetics.html?highlight=pc_relate#hail.methods.pc_relate>`.
    
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
            .when(
            (kin_expr < first_degree_kin_thresholds[0]),
            SECOND_DEGREE_RELATIVES
        )
            .when(
            (kin_expr < first_degree_kin_thresholds[1]) &
            (ibd0_expr <= ibd0_0_max) &
            (ibd1_expr >= ibd1_100_min) &
            (ibd2_expr <= ibd2_0_max),
            PARENT_CHILD
        )
            .when(
            (kin_expr < first_degree_kin_thresholds[1]) &
            (ibd0_expr >= ibd0_25_thresholds[0]) &
            (ibd0_expr <= ibd0_25_thresholds[1]) &
            (ibd1_expr >= ibd1_50_thresholds[0]) &
            (ibd1_expr <= ibd1_50_thresholds[1]) &
            (ibd2_expr >= ibd2_25_thresholds[0]) &
            (ibd2_expr <= ibd2_25_thresholds[1]),
            SIBLINGS
        )
            .when(
            (kin_expr > first_degree_kin_thresholds[1]) &
            (ibd0_expr < ibd0_0_max) &
            (ibd1_expr >= ibd1_0_thresholds[0]) &
            (ibd1_expr <= ibd1_0_thresholds[1]) &
            (ibd2_expr >= ibd2_100_thresholds[0]) &
            (ibd2_expr <= ibd2_100_thresholds[1]),
            DUPLICATE_OR_TWINS
        )
            .default(AMBIGUOUS_RELATIONSHIP)
    )


def infer_families(relationship_ht: hl.Table,
                   sex: Union[hl.Table, Dict[str, bool]],
                   duplicate_samples_ht: hl.Table,
                   i_col: str = 'i',
                   j_col: str = 'j',
                   relationship_col: str = 'relationship'
                   ) -> hl.Pedigree:
    """
    This function takes a hail Table with a row for each pair of individuals i,j in the data that are related (it's OK to have unrelated samples too).
    The `relationship_col` should be a column specifying the relationship between each two samples as defined in this module's constants.

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

    def group_parent_child_pairs_by_fam(parent_child_pairs: Iterable[Tuple[str, str]]) -> List[List[Tuple[str, str]]]:
        """
        Takes all parent-children pairs and groups them by family.
        A family here is defined as a list of sample-pairs which all share at least one sample with at least one other sample-pair in the list.

        :param parent_child_pairs: All the parent-children pairs
        :return: A list of families, where each element of the list is a list of the parent-children pairs
        """
        fam_id = 1  # stores the current family id
        s_fam = dict()  # stores the family id for each sample
        fams = defaultdict(list)  # stores fam_id -> sample-pairs
        for pair in parent_child_pairs:
            if pair[0] in s_fam:
                if pair[1] in s_fam:
                    if s_fam[pair[0]] != s_fam[pair[1]]:  # If both samples are in different families, merge the families
                        new_fam_id = s_fam[pair[0]]
                        fam_id_to_merge = s_fam[pair[1]]
                        for s in s_fam:
                            if s_fam[s] == fam_id_to_merge:
                                s_fam[s] = new_fam_id
                        fams[new_fam_id].extend(fams.pop(fam_id_to_merge))
                else:  # If only the 1st sample in the pair is already in a family, assign the 2nd sample in the pair to the same family
                    s_fam[pair[1]] = s_fam[pair[0]]
                fams[s_fam[pair[0]]].append(pair)
            elif pair[1] in s_fam:   # If only the 2nd sample in the pair is already in a family, assign the 1st sample in the pair to the same family
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
            related_pairs: Dict[Tuple[str, str], str]
    ) -> List[hl.Trio]:
        """
        Generates trios based from the list of parent-child pairs in the family and
        all related pairs in the data. Only complete parent/offspring trios are included in the results.

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
            1. All pairs of unrelated samples with different sexes within the family are extracted as possible parent pairs

            :param samples: All samples in the family
            :return: Possible parent pairs
            """
            possible_parents = []
            for i in range(len(samples)):
                for j in range(i + 1, len(samples)):
                    if related_pairs.get(tuple(sorted([samples[i], samples[j]]))) is None:
                        if sex.get(samples[i]) is True and sex.get(samples[j]) is False:
                            possible_parents.append((samples[i], samples[j]))
                        elif sex.get(samples[i]) is False and sex.get(samples[j]) is True:
                            possible_parents.append((samples[j], samples[i]))
            return possible_parents

        def get_children(possible_parents: Tuple[str, str]) -> List[str]:
            """
            2. For a given possible parent pair, a list of all children is constructed (each child in the list has a parent-offspring pair with each parent)

            :param possible_parents: A pair of possible parents
            :return: The list of all children (if any) corresponding to the possible parents
            """
            possible_offsprings = defaultdict(set)  # stores sample -> set of parents in the possible_parents where (sample, parent) is found in possible_child_pairs
            for pair in parent_child_pairs:
                if possible_parents[0] == pair[0]:
                    possible_offsprings[pair[1]].add(possible_parents[0])
                elif possible_parents[0] == pair[1]:
                    possible_offsprings[pair[0]].add(possible_parents[0])
                elif possible_parents[1] == pair[0]:
                    possible_offsprings[pair[1]].add(possible_parents[1])
                elif possible_parents[1] == pair[1]:
                    possible_offsprings[pair[0]].add(possible_parents[1])

            return [s for s, parents in possible_offsprings.items() if len(parents) == 2]

        def check_sibs(children: List[str]) -> bool:
            """
            3. If there are multiple children for a given parent pair, all children should be siblings with each other

            :param children: List of all children for a given parent pair
            :return: Whether all children in the list are siblings
            """
            for i in range(len(children)):
                for j in range(i + 1, len(children)):
                    if related_pairs[tuple(sorted([children[i], children[j]]))] != SIBLINGS:
                        return False
            return True

        def discard_multi_parents_children(trios: List[hl.Trio]):
            """
            4. Check that each child was only assigned a single pair of parents. If a child is found to have multiple parent pairs, they are ALL discarded.

            :param trios: All trios formed for this family
            :return: The list of trios for which each child has a single parents pair.
            """
            children_trios = defaultdict(list)
            for trio in trios:
                children_trios[trio.s].append(trio)

            for s, s_trios in children_trios.items():
                if len(s_trios) > 1:
                    logger.warning("Discarded duplicated child {0} found multiple in trios: {1}".format(
                        s,
                        ", ".join([str(trio) for trio in s_trios])
                    ))

            return [trios[0] for trios in children_trios.values() if len(trios) == 1]

        # Get all possible pairs of parents in (father, mother) order
        all_possible_parents = get_possible_parents(list({s for pair in parent_child_pairs for s in pair}))

        trios = []
        for possible_parents in all_possible_parents:
            children = get_children(possible_parents)
            if check_sibs(children):
                trios.extend([
                    hl.Trio(
                        s=s,
                        fam_id=fam_id,
                        pat_id=possible_parents[0],
                        mat_id=possible_parents[1],
                        is_female=sex.get(s)
                    )
                    for s in children
                ])
            else:
                logger.warn(
                    "Discarded family with same parents, and multiple offspring that weren't siblings:"
                    "\nMother: {}\nFather:{}\nChildren:{}".format(
                        possible_parents[0],
                        possible_parents[1],
                        ", ".join(children)
                    )
                )

        return discard_multi_parents_children(trios)

    # Get all the relations we care about:
    # => Remove unrelateds and duplicates
    dups = duplicate_samples_ht.aggregate(
        hl.agg.explode(
            lambda dup: hl.agg.collect_as_set(dup),
            duplicate_samples_ht.filtered
        ),
        _localize=False
    )
    relationship_ht = relationship_ht.filter(
        ~dups.contains(relationship_ht[i_col]) &
        ~dups.contains(relationship_ht[j_col]) &
        (relationship_ht[relationship_col] != UNRELATED)
    )

    # Check relatedness table format
    if not relationship_ht[i_col].dtype == relationship_ht[j_col].dtype:
        logger.error("i_col and j_col of the relatedness table need to be of the same type.")

    # If i_col and j_col aren't str, then convert them
    if not isinstance(relationship_ht[i_col], hl.expr.StringExpression):
        logger.warning(f"Pedigrees can only be constructed from string IDs, but your relatedness_ht ID column is of type: {relationship_ht[i_col].dtype}. Expression will be converted to string in Pedigrees.")
        if isinstance(relationship_ht[i_col], hl.expr.StructExpression):
            logger.warning(f"Struct fields {list(relationship_ht[i_col])} will be joined by underscores to use as sample names in Pedigree.")
            relationship_ht = relationship_ht.key_by(
                **{
                    i_col: hl.delimit(hl.array([hl.str(relationship_ht[i_col][x]) for x in relationship_ht[i_col]]), "_"),
                    j_col: hl.delimit(hl.array([hl.str(relationship_ht[j_col][x]) for x in relationship_ht[j_col]]), "_")
                }
            )
        else:
            raise NotImplementedError("The `i_col` and `j_col` columns of the `relationship_ht` argument passed to infer_families are not of type StringExpression or Struct.")

    # If sex is a Table, extract sex information as a Dict
    if isinstance(sex, hl.Table):
        sex = dict(hl.tuple([sex.s, sex.is_female]).collect())

    # Collect all related sample pairs and
    # create a dictionnary with pairs as keys and relationships as values
    # Sample-pairs are tuples ordered by sample name
    related_pairs = {
        tuple(sorted([i, j])): rel for i, j, rel in
        hl.tuple([relationship_ht.i, relationship_ht.j, relationship_ht.relationship]).collect()
    }

    parent_child_pairs_by_fam = group_parent_child_pairs_by_fam([pair for pair, rel in related_pairs.items() if rel == PARENT_CHILD])
    return hl.Pedigree([
        trio for fam_index, parent_child_pairs in enumerate(parent_child_pairs_by_fam)
        for trio in get_trios(str(fam_index), parent_child_pairs, related_pairs)
    ])


def create_fake_pedigree(
        n: int,
        sample_list: List[str],
        exclude_real_probands: bool = False,
        max_tries: int = 10,
        real_pedigree: Optional[hl.Pedigree] = None
) -> hl.Pedigree:
    """
    Generates a pedigree made of trios created by sampling 3 random samples in the sample list.
    If `real_pedigree` is given, then children the resulting fake trios will not include any trio with proband - parents that are in the real ones.
    Each sample can be used only once as a proband in the resulting trios.
    Sex of probands in fake trios is random.

    :param n: Number of fake trios desired in the pedigree
    :param sample_list: List of samples
    :param exclude_real_probands: If set, then fake trios probands cannot be in the real trios probands.
    :param max_tries: Maximum number of sampling to try before bailing out (preventing infinite loop if `n` is too large w.r.t. the number of samples)
    :param real_pedigree: Optional pedigree to exclude children from
    :return: Fake pedigree
    """
    real_trios = {trio.s: trio for trio in real_pedigree.trios} if real_pedigree is not None else dict()

    if exclude_real_probands and len(real_trios) == len(set(sample_list)):
        logger.warning("All samples are in the real probands list; cannot create any fake pedigrees with exclude_real_probands=True. Returning an empty Pedigree.")
        return hl.Pedigree([])

    fake_trios = {}
    tries = 0
    while len(fake_trios) < n and tries < max_tries:
        s, mat_id, pat_id = random.sample(sample_list, 3)
        if (
                (s in real_trios and (
                        exclude_real_probands or
                        {mat_id, pat_id} == {real_trios[s].mat_id, real_trios[s].pat_id}
                )) or
                s in fake_trios
        ):
            tries += 1
        else:
            tries = 0
            fake_trios[s] = hl.Trio(
                s=s,
                pat_id=pat_id,
                mat_id=mat_id,
                fam_id=f"fake_{str(len(fake_trios))}",
                is_female=bool(random.getrandbits(1))
            )

    if tries == max_tries:
        logger.warning(f"Only returning {len(fake_trios)} fake trios; random trio sampling stopped after reaching the maximum {max_tries} iterations")

    return hl.Pedigree(list(fake_trios.values()))
