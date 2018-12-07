'''
The tests in this package are to ensure the proper resultant dtypes of
set operations.
'''
import itertools as it

import numpy as np
import pytest

from pandas.core.dtypes.common import is_dtype_equal

import pandas as pd
from pandas import Int64Index, RangeIndex
from pandas.tests.indexes.conftest import indices
import pandas.util.testing as tm

COMPATIBLE_INCONSISTENT_PAIRS = {
    (Int64Index, RangeIndex): (tm.makeIntIndex, tm.makeRangeIndex)
}


def test_union_same_types(indices):
    # Union with a non-unique, non-monotonic index raises error
    # Only needed for bool index factory
    idx1 = indices.sort_values()
    idx2 = indices.sort_values()
    assert idx1.union(idx2).dtype == idx1.dtype

    # Note: catIndex reflects only left dtype, should it reflect both?


@pytest.mark.parametrize(
    'idx1,idx2',
    list(it.combinations(indices._pytestfixturefunction.params, 2))
)
def test_union_different_types(idx1, idx2):
    # GH 23525
    pair = tuple(sorted([type(idx1), type(idx2)], key=lambda x: str(x)))
    if pair in COMPATIBLE_INCONSISTENT_PAIRS:
        return

    if any(isinstance(idx, pd.MultiIndex) for idx in [idx1, idx2]):
        return

    if is_dtype_equal(idx1.dtype, idx2.dtype):
        return

    # A union with a CategoricalIndex (even as dtype('O')) and a
    # non-CategoricalIndex can only be made if both indices are monotonic.
    # This is true before this PR as well.

    # Union with a non-unique, non-monotonic index raises error
    # This applies to the boolean index
    idx1 = idx1.sort_values()
    idx2 = idx2.sort_values()

    assert idx1.union(idx2).dtype == np.dtype('O')
    assert idx2.union(idx1).dtype == np.dtype('O')


@pytest.mark.parametrize('idx_fact1,idx_fact2',
                         COMPATIBLE_INCONSISTENT_PAIRS.values())
def test_compatible_inconsistent_pairs(idx_fact1, idx_fact2):
    # GH 23525
    idx1 = idx_fact1(10)
    idx2 = idx_fact2(20)

    res1 = idx1.union(idx2)
    res2 = idx2.union(idx1)

    assert res1.dtype in (idx1.dtype, idx2.dtype)
    assert res2.dtype in (idx1.dtype, idx2.dtype)
