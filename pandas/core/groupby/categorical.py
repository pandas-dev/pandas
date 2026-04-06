from __future__ import annotations

import numpy as np

from pandas.core.algorithms import unique1d
from pandas.core.arrays.categorical import (
    Categorical,
    CategoricalDtype,
)


def recode_for_groupby(c: Categorical, sort: bool, observed: bool) -> Categorical:
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
    """
    # we only care about observed values
    if observed:
        # In cases with c.ordered, this is equivalent to
        #  return c.remove_unused_categories(), c

        take_codes = unique1d(c.codes[c.codes != -1])

        if sort:
            take_codes = np.sort(take_codes)

    elif sort:
        # Already sorted according to c.categories; all is fine
        return c

    else:
        # sort=False should order groups in as-encountered order (GH-8868)

        # GH:46909: Re-ordering codes faster than using (set|add|reorder)_categories
        # GH 38140: exclude nan from indexer for categories
        unique_notnan_codes = unique1d(c.codes[c.codes != -1])
        if (num_cat := len(c.categories)) > len(unique_notnan_codes):
            # GH 13179: All categories need to be present, even if missing
            # from the data
            missing_codes = np.setdiff1d(
                np.arange(num_cat), unique_notnan_codes, assume_unique=True
            )
            take_codes = np.concatenate((unique_notnan_codes, missing_codes))
        else:
            take_codes = unique_notnan_codes

    # GH#48749: Remap codes using direct array indexing instead of
    # hash-based recode_for_categories / Categorical constructor,
    # which use get_indexer_for (hash table lookup).
    reverse_indexer = np.empty(len(c.categories), dtype=np.intp)
    reverse_indexer[take_codes] = np.arange(len(take_codes))

    mask = c.codes >= 0
    new_codes = np.full_like(c.codes, fill_value=-1)
    if mask.any():
        new_codes[mask] = reverse_indexer[c.codes[mask]]

    new_cats = c.categories.take(take_codes)
    dtype = CategoricalDtype(new_cats, ordered=c.ordered)
    return Categorical._simple_new(new_codes, dtype=dtype)
