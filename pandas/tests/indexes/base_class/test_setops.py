import numpy as np
import pytest

from pandas import Index, Series
import pandas._testing as tm
from pandas.core.algorithms import safe_sort


class TestIndexSetOps:
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

    @pytest.mark.parametrize("sort", [None, False])
    def test_intersection_base(self, sort):
        # (same results for py2 and py3 but sortedness not tested elsewhere)
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[:5]
        second = index[:3]

        expected = Index([0, 1, "a"]) if sort is None else Index([0, "a", 1])
        result = first.intersection(second, sort=sort)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("klass", [np.array, Series, list])
    @pytest.mark.parametrize("sort", [None, False])
    def test_intersection_different_type_base(self, klass, sort):
        # GH 10149
        index = Index([0, "a", 1, "b", 2, "c"])
        first = index[:5]
        second = index[:3]

        result = first.intersection(klass(second.values), sort=sort)
        assert tm.equalContents(result, second)

    @pytest.mark.parametrize("sort", [None, False])
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
