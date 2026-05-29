from __future__ import annotations

import numpy as np

from pandas.core.dtypes.cast import coerce_indexer_dtype

from pandas.core.algorithms import (
    take_nd,
    unique1d,
)
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
            if len(take_codes) == len(c.categories):
                # Every category is observed; sorted codes are already
                # 0..n-1, so c is unchanged.
                return c

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

    # GH#48749: Build the old-code -> new-code mapping by scatter-assignment
    # (cheap integer indexing) rather than recode_for_categories / the
    # Categorical constructor, which call get_indexer_for -- a hash-table
    # lookup over the categories that is expensive when there are many
    # (especially string) categories. The mapping is then applied with
    # take_nd, which maps the NaN sentinel (-1) to -1 via fill_value.
    reverse_indexer = np.empty(len(c.categories), dtype=np.intp)
    reverse_indexer[take_codes] = np.arange(len(take_codes))
    reverse_indexer = coerce_indexer_dtype(reverse_indexer, c.categories)

    new_codes = take_nd(reverse_indexer, c.codes, fill_value=-1)

    new_cats = c.categories.take(take_codes)
    dtype = CategoricalDtype(new_cats, ordered=c.ordered)
    return Categorical._simple_new(new_codes, dtype=dtype)
