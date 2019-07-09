"""
The tests in this package are to ensure the proper resultant dtypes of
set operations.
"""
from collections import OrderedDict
import itertools as it

import numpy as np
import pytest

from pandas.core.dtypes.common import is_dtype_equal

import pandas as pd
from pandas import Float64Index, Int64Index, RangeIndex, UInt64Index
from pandas.api.types import pandas_dtype
from pandas.tests.indexes.conftest import indices_list
import pandas.util.testing as tm

COMPATIBLE_INCONSISTENT_PAIRS = OrderedDict(
    [
        ((Int64Index, RangeIndex), (tm.makeIntIndex, tm.makeRangeIndex)),
        ((Float64Index, Int64Index), (tm.makeFloatIndex, tm.makeIntIndex)),
        ((Float64Index, RangeIndex), (tm.makeFloatIndex, tm.makeIntIndex)),
        ((Float64Index, UInt64Index), (tm.makeFloatIndex, tm.makeUIntIndex)),
    ]
)


@pytest.fixture(
    params=list(it.combinations(indices_list, 2)),
    ids=lambda x: type(x[0]).__name__ + type(x[1]).__name__,
)
def index_pair(request):
    """
    Create all combinations of 2 index types.
    """
    return request.param


def test_union_same_types(indices):
    # Union with a non-unique, non-monotonic index raises error
    # Only needed for bool index factory
    idx1 = indices.sort_values()
    idx2 = indices.sort_values()
    assert idx1.union(idx2).dtype == idx1.dtype


def test_union_different_types(index_pair):
    # GH 23525
    idx1, idx2 = index_pair
    type_pair = tuple(sorted([type(idx1), type(idx2)], key=lambda x: str(x)))
    if type_pair in COMPATIBLE_INCONSISTENT_PAIRS:
        pytest.xfail("This test only considers non compatible indexes.")

    if any(isinstance(idx, pd.MultiIndex) for idx in index_pair):
        pytest.xfail("This test doesn't consider multiindixes.")

    if is_dtype_equal(idx1.dtype, idx2.dtype):
        pytest.xfail("This test only considers non matching dtypes.")

    # A union with a CategoricalIndex (even as dtype('O')) and a
    # non-CategoricalIndex can only be made if both indices are monotonic.
    # This is true before this PR as well.

    # Union with a non-unique, non-monotonic index raises error
    # This applies to the boolean index
    idx1 = idx1.sort_values()
    idx2 = idx2.sort_values()

    assert idx1.union(idx2).dtype == np.dtype("O")
    assert idx2.union(idx1).dtype == np.dtype("O")


@pytest.mark.parametrize("idx_fact1,idx_fact2", COMPATIBLE_INCONSISTENT_PAIRS.values())
def test_compatible_inconsistent_pairs(idx_fact1, idx_fact2):
    # GH 23525
    idx1 = idx_fact1(10)
    idx2 = idx_fact2(20)

    res1 = idx1.union(idx2)
    res2 = idx2.union(idx1)

    assert res1.dtype in (idx1.dtype, idx2.dtype)
    assert res2.dtype in (idx1.dtype, idx2.dtype)


@pytest.mark.parametrize(
    "left, right, expected",
    [
        ("int64", "int64", "int64"),
        ("int64", "uint64", "object"),
        ("int64", "float64", "float64"),
        ("uint64", "float64", "float64"),
        ("uint64", "uint64", "uint64"),
        ("float64", "float64", "float64"),
        ("datetime64[ns]", "int64", "object"),
        ("datetime64[ns]", "uint64", "object"),
        ("datetime64[ns]", "float64", "object"),
        ("datetime64[ns, CET]", "int64", "object"),
        ("datetime64[ns, CET]", "uint64", "object"),
        ("datetime64[ns, CET]", "float64", "object"),
        ("Period[D]", "int64", "object"),
        ("Period[D]", "uint64", "object"),
        ("Period[D]", "float64", "object"),
    ],
)
def test_union_dtypes(left, right, expected):
    left = pandas_dtype(left)
    right = pandas_dtype(right)
    a = pd.Index([], dtype=left)
    b = pd.Index([], dtype=right)
    result = (a | b).dtype
    assert result == expected
