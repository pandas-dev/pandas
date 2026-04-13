import numpy as np
import pytest

from pandas import (
    CategoricalIndex,
    Index,
)
import pandas._testing as tm


@pytest.mark.parametrize("method", ["union", "intersection"])
def test_setop_mismatched_category_order(method):
    # GH#55335 unordered CategoricalIndex with same categories in different
    # order should give correct results for union/intersection
    cats1 = ["a", "b", "c", "d"]
    cats2 = ["d", "c", "b", "a"]
    id1 = CategoricalIndex(["a", "c"], categories=cats1)
    id2 = CategoricalIndex(["d", "b"], categories=cats2)

    result = getattr(id1, method)(id2)
    if method == "union":
        # sort=None (default) sorts the result
        expected = CategoricalIndex(["a", "b", "c", "d"], categories=cats1)
    else:
        expected = CategoricalIndex([], categories=cats1)
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("na_value", [None, np.nan])
def test_difference_with_na(na_value):
    # GH 57318
    ci = CategoricalIndex(["a", "b", "c", None])
    other = Index(["c", na_value])
    result = ci.difference(other)
    expected = CategoricalIndex(["a", "b"], categories=["a", "b", "c"])
    tm.assert_index_equal(result, expected)
