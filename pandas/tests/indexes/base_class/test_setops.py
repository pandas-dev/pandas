import numpy as np
import pytest

import pandas as pd
from pandas import Index, Series
import pandas._testing as tm
from pandas.core.algorithms import safe_sort


class TestIndexSetOps:
    @pytest.mark.parametrize(
        "method", ["union", "intersection", "difference", "symmetric_difference"]
    )
    def test_setops_disallow_true(self, method):
        idx1 = pd.Index(["a", "b"])
        idx2 = pd.Index(["b", "c"])

        with pytest.raises(ValueError, match="The 'sort' keyword only takes"):
            getattr(idx1, method)(idx2, sort=True)

    def test_setops_preserve_object_dtype(self):
        idx = pd.Index([1, 2, 3], dtype=object)
        result = idx.intersection(idx[1:])
        expected = idx[1:]
        tm.assert_index_equal(result, expected)

        # if other is not monotonic increasing, intersection goes through
        #  a different route
        result = idx.intersection(idx[1:][::-1])
        tm.assert_index_equal(result, expected)

        result = idx._union(idx[1:], sort=None)
        expected = idx
        tm.assert_index_equal(result, expected)

        result = idx.union(idx[1:], sort=None)
        tm.assert_index_equal(result, expected)

        # if other is not monotonic increasing, _union goes through
        #  a different route
        result = idx._union(idx[1:][::-1], sort=None)
        tm.assert_index_equal(result, expected)

        result = idx.union(idx[1:][::-1], sort=None)
        tm.assert_index_equal(result, expected)

    def test_union_base(self):
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[3:]
        second = index[:5]

        result = first.union(second)

        expected = Index([0, 1, 2, "a", "b", "c"])
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("klass", [np.array, Series, list])
    def test_union_different_type_base(self, klass):
        # GH 10149
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[3:]
        second = index[:5]

        result = first.union(klass(second.values))

        assert tm.equalContents(result, index)

    def test_union_sort_other_incomparable(self):
        # https://github.com/pandas-dev/pandas/issues/24959
        idx = pd.Index([1, pd.Timestamp("2000")])
        # default (sort=None)
        with tm.assert_produces_warning(RuntimeWarning):
            result = idx.union(idx[:1])

        tm.assert_index_equal(result, idx)

        # sort=None
        with tm.assert_produces_warning(RuntimeWarning):
            result = idx.union(idx[:1], sort=None)
        tm.assert_index_equal(result, idx)

        # sort=False
        result = idx.union(idx[:1], sort=False)
        tm.assert_index_equal(result, idx)

    @pytest.mark.xfail(reason="Not implemented")
    def test_union_sort_other_incomparable_true(self):
        # TODO decide on True behaviour
        # sort=True
        idx = pd.Index([1, pd.Timestamp("2000")])
        with pytest.raises(TypeError, match=".*"):
            idx.union(idx[:1], sort=True)

    def test_intersection_base(self, sort):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[:5]
        second = index[:3]

        expected = Index([0, 1, "a"]) if sort is None else Index([0, "a", 1])
        result = first.intersection(second, sort=sort)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("klass", [np.array, Series, list])
    def test_intersection_different_type_base(self, klass, sort):
        # GH 10149
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[:5]
        second = index[:3]

        result = first.intersection(klass(second.values), sort=sort)
        assert tm.equalContents(result, second)

    def test_intersect_nosort(self):
        result = pd.Index(["c", "b", "a"]).intersection(["b", "a"])
        expected = pd.Index(["b", "a"])
        tm.assert_index_equal(result, expected)

    def test_intersection_equal_sort(self):
        idx = pd.Index(["c", "a", "b"])
        tm.assert_index_equal(idx.intersection(idx, sort=False), idx)
        tm.assert_index_equal(idx.intersection(idx, sort=None), idx)

    def test_difference_base(self, sort):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[:4]
        second = index[3:]

        result = first.difference(second, sort)
        expected = Index([0, "a", 1])
        if sort is None:
            expected = Index(safe_sort(expected))
        tm.assert_index_equal(result, expected)

    def test_symmetric_difference(self):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[:4]
        second = index[3:]

        result = first.symmetric_difference(second)
        expected = Index([0, 1, 2, "a", "c"])
        tm.assert_index_equal(result, expected)
