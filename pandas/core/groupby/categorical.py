from __future__ import annotations

import dataclasses
from typing import final

import numpy as np

from pandas._typing import npt

from pandas.core.algorithms import unique1d
from pandas.core.arrays.categorical import (
    Categorical,
    CategoricalDtype,
    recode_for_categories,
)
from pandas.core.indexes.api import (
    CategoricalIndex,
    Index,
)


@final
@dataclasses.dataclass
class CategoricalGrouper:
    original_grouping_vector: Categorical | None
    new_grouping_vector: Categorical
    observed: bool
    sort: bool

    def result_index(self, group_index: Index) -> Index:
        if self.original_grouping_vector is None:
            return group_index
        assert isinstance(group_index, CategoricalIndex)
        return recode_from_groupby(
            self.original_grouping_vector, self.sort, group_index
        )

    @classmethod
    def make(cls, grouping_vector, sort: bool, observed: bool) -> CategoricalGrouper:
        new_grouping_vector, original_grouping_vector = recode_for_groupby(
            grouping_vector, sort, observed
        )
        return cls(
            original_grouping_vector=original_grouping_vector,
            new_grouping_vector=new_grouping_vector,
            observed=observed,
            sort=sort,
        )

    def codes_and_uniques(
        self, cat: Categorical
    ) -> tuple[npt.NDArray[np.signedinteger], Categorical]:
        """
        This does a version of `algorithms.factorize()` that works on categoricals.
        """
        categories = cat.categories

        if self.observed:
            ucodes = unique1d(cat.codes)
            ucodes = ucodes[ucodes != -1]
            if self.sort or cat.ordered:
                ucodes = np.sort(ucodes)
        else:
            ucodes = np.arange(len(categories))

        uniques = Categorical.from_codes(
            codes=ucodes, categories=categories, ordered=cat.ordered
        )
        return cat.codes, uniques


def recode_for_groupby(
    c: Categorical, sort: bool, observed: bool
) -> tuple[Categorical, Categorical | None]:
    """
    Code the categories to ensure we can groupby for categoricals.

    If observed=True, we return a new Categorical with the observed
    categories only.

    If sort=False, return a copy of self, coded with categories as
    returned by .unique(), followed by any categories not appearing in
    the data. If sort=True, return self.

    This method is needed solely to ensure the categorical index of the
    GroupBy result has categories in the order of appearance in the data
    (GH-8868).

    Parameters
    ----------
    c : Categorical
    sort : bool
        The value of the sort parameter groupby was called with.
    observed : bool
        Account only for the observed values

    Returns
    -------
    Categorical
        If sort=False, the new categories are set to the order of
        appearance in codes (unless ordered=True, in which case the
        original order is preserved), followed by any unrepresented
        categories in the original order.
    Categorical or None
        If we are observed, return the original categorical, otherwise None
    """
    # we only care about observed values
    if observed:
        # In cases with c.ordered, this is equivalent to
        #  return c.remove_unused_categories(), c

        unique_codes = unique1d(c.codes)

        take_codes = unique_codes[unique_codes != -1]
        if c.ordered:
            take_codes = np.sort(take_codes)

        # we recode according to the uniques
        categories = c.categories.take(take_codes)
        codes = recode_for_categories(c.codes, c.categories, categories)

        # return a new categorical that maps our new codes
        # and categories
        dtype = CategoricalDtype(categories, ordered=c.ordered)
        return Categorical(codes, dtype=dtype, fastpath=True), c

    # Already sorted according to c.categories; all is fine
    if sort:
        return c, None

    # sort=False should order groups in as-encountered order (GH-8868)
    cat = c.unique()

    # See GH-38140 for block below
    # exclude nan from indexer for categories
    take_codes = cat.codes[cat.codes != -1]
    if cat.ordered:
        take_codes = np.sort(take_codes)
    cat = cat.set_categories(cat.categories.take(take_codes))

    # But for groupby to work, all categories should be present,
    # including those missing from the data (GH-13179), which .unique()
    # above dropped
    cat = cat.add_categories(c.categories[~c.categories.isin(cat.categories)])

    return c.reorder_categories(cat.categories), None


def recode_from_groupby(
    c: Categorical, sort: bool, ci: CategoricalIndex
) -> CategoricalIndex:
    """
    Reverse the codes_to_groupby to account for sort / observed.

    Parameters
    ----------
    c : Categorical
    sort : bool
        The value of the sort parameter groupby was called with.
    ci : CategoricalIndex
        The codes / categories to recode

    Returns
    -------
    CategoricalIndex
    """
    # we re-order to the original category orderings
    if sort:
        # error: "CategoricalIndex" has no attribute "set_categories"
        return ci.set_categories(c.categories)  # type: ignore[attr-defined]

    # we are not sorting, so add unobserved to the end
    new_cats = c.categories[~c.categories.isin(ci.categories)]
    # error: "CategoricalIndex" has no attribute "add_categories"
    return ci.add_categories(new_cats)  # type: ignore[attr-defined]
