"""
The tests in this package are to ensure the proper resultant dtypes of
set operations.
"""
import numpy as np
import pytest

from pandas.core.dtypes.common import is_dtype_equal

import pandas as pd
from pandas import (
    CategoricalIndex,
    DatetimeIndex,
    Float64Index,
    Int64Index,
    MultiIndex,
    RangeIndex,
    Series,
    TimedeltaIndex,
    UInt64Index,
)
import pandas._testing as tm
from pandas.api.types import is_datetime64tz_dtype, pandas_dtype

COMPATIBLE_INCONSISTENT_PAIRS = {
    (Int64Index, RangeIndex): (tm.makeIntIndex, tm.makeRangeIndex),
    (Float64Index, Int64Index): (tm.makeFloatIndex, tm.makeIntIndex),
    (Float64Index, RangeIndex): (tm.makeFloatIndex, tm.makeIntIndex),
    (Float64Index, UInt64Index): (tm.makeFloatIndex, tm.makeUIntIndex),
}


def test_union_same_types(index):
    # Union with a non-unique, non-monotonic index raises error
    # Only needed for bool index factory
    idx1 = index.sort_values()
    idx2 = index.sort_values()
    assert idx1.union(idx2).dtype == idx1.dtype


def test_union_different_types(index, index_fixture2):
    # This test only considers combinations of indices
    # GH 23525
    idx1, idx2 = index, index_fixture2
    type_pair = tuple(sorted([type(idx1), type(idx2)], key=lambda x: str(x)))
    if type_pair in COMPATIBLE_INCONSISTENT_PAIRS:
        pytest.xfail("This test only considers non compatible indexes.")

    if any(isinstance(idx, pd.MultiIndex) for idx in (idx1, idx2)):
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
@pytest.mark.parametrize("names", [("foo", "foo", "foo"), ("foo", "bar", None)])
def test_union_dtypes(left, right, expected, names):
    left = pandas_dtype(left)
    right = pandas_dtype(right)
    a = pd.Index([], dtype=left, name=names[0])
    b = pd.Index([], dtype=right, name=names[1])
    result = a.union(b)
    assert result.dtype == expected
    assert result.name == names[2]

    # Testing name retention
    # TODO: pin down desired dtype; do we want it to be commutative?
    result = a.intersection(b)
    assert result.name == names[2]


def test_dunder_inplace_setops_deprecated(index):
    # GH#37374 these will become logical ops, not setops

    with tm.assert_produces_warning(FutureWarning):
        index |= index

    with tm.assert_produces_warning(FutureWarning):
        index &= index

    with tm.assert_produces_warning(FutureWarning):
        index ^= index


@pytest.mark.parametrize("values", [[1, 2, 2, 3], [3, 3]])
def test_intersection_duplicates(values):
    # GH#31326
    a = pd.Index(values)
    b = pd.Index([3, 3])
    result = a.intersection(b)
    expected = pd.Index([3])
    tm.assert_index_equal(result, expected)


class TestSetOps:
    # Set operation tests shared by all indexes in the `index` fixture
    @pytest.mark.parametrize("case", [0.5, "xxx"])
    @pytest.mark.parametrize(
        "method", ["intersection", "union", "difference", "symmetric_difference"]
    )
    def test_set_ops_error_cases(self, case, method, index):
        # non-iterable input
        msg = "Input must be Index or array-like"
        with pytest.raises(TypeError, match=msg):
            getattr(index, method)(case)

    def test_intersection_base(self, index):
        if isinstance(index, CategoricalIndex):
            return

        first = index[:5]
        second = index[:3]
        intersect = first.intersection(second)
        assert tm.equalContents(intersect, second)

        if is_datetime64tz_dtype(index.dtype):
            # The second.values below will drop tz, so the rest of this test
            #  is not applicable.
            return

        # GH#10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.intersection(case)
            assert tm.equalContents(result, second)

        if isinstance(index, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.intersection([1, 2, 3])

    def test_union_base(self, index):
        first = index[3:]
        second = index[:5]
        everything = index
        union = first.union(second)
        assert tm.equalContents(union, everything)

        if is_datetime64tz_dtype(index.dtype):
            # The second.values below will drop tz, so the rest of this test
            #  is not applicable.
            return

        # GH#10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            if not isinstance(index, CategoricalIndex):
                result = first.union(case)
                assert tm.equalContents(result, everything), (
                    result,
                    everything,
                    type(case),
                )

        if isinstance(index, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.union([1, 2, 3])

    def test_difference_base(self, sort, index):
        first = index[2:]
        second = index[:4]
        if isinstance(index, CategoricalIndex) or index.is_boolean():
            answer = []
        else:
            answer = index[4:]
        result = first.difference(second, sort)
        assert tm.equalContents(result, answer)

        # GH#10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(index, (DatetimeIndex, TimedeltaIndex)):
                assert type(result) == type(answer)
                tm.assert_numpy_array_equal(
                    result.sort_values().asi8, answer.sort_values().asi8
                )
            else:
                result = first.difference(case, sort)
                assert tm.equalContents(result, answer)

        if isinstance(index, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.difference([1, 2, 3], sort)

    def test_symmetric_difference(self, index):
        if isinstance(index, CategoricalIndex):
            return
        if len(index) < 2:
            return
        if index[0] in index[1:] or index[-1] in index[:-1]:
            # index fixture has e.g. an index of bools that does not satisfy this,
            #  another with [0, 0, 1, 1, 2, 2]
            return

        first = index[1:]
        second = index[:-1]
        answer = index[[0, -1]]
        result = first.symmetric_difference(second)
        assert tm.equalContents(result, answer)

        # GH#10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            if is_datetime64tz_dtype(first):
                with pytest.raises(ValueError, match="Tz-aware"):
                    # `second.values` casts to tznaive
                    # TODO: should the symmetric_difference then be the union?
                    first.symmetric_difference(case)
                continue
            result = first.symmetric_difference(case)
            assert tm.equalContents(result, answer)

        if isinstance(index, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with pytest.raises(TypeError, match=msg):
                first.symmetric_difference([1, 2, 3])

    @pytest.mark.parametrize(
        "fname, sname, expected_name",
        [
            ("A", "A", "A"),
            ("A", "B", None),
            ("A", None, None),
            (None, "B", None),
            (None, None, None),
        ],
    )
    def test_corner_union(self, index, fname, sname, expected_name):
        # GH#9943, GH#9862
        # Test unions with various name combinations
        # Do not test MultiIndex or repeats

        if isinstance(index, MultiIndex) or not index.is_unique:
            pytest.skip("Not for MultiIndex or repeated indices")

        # Test copy.union(copy)
        first = index.copy().set_names(fname)
        second = index.copy().set_names(sname)
        union = first.union(second)
        expected = index.copy().set_names(expected_name)
        tm.assert_index_equal(union, expected)

        # Test copy.union(empty)
        first = index.copy().set_names(fname)
        second = index.drop(index).set_names(sname)
        union = first.union(second)
        expected = index.copy().set_names(expected_name)
        tm.assert_index_equal(union, expected)

        # Test empty.union(copy)
        first = index.drop(index).set_names(fname)
        second = index.copy().set_names(sname)
        union = first.union(second)
        expected = index.copy().set_names(expected_name)
        tm.assert_index_equal(union, expected)

        # Test empty.union(empty)
        first = index.drop(index).set_names(fname)
        second = index.drop(index).set_names(sname)
        union = first.union(second)
        expected = index.drop(index).set_names(expected_name)
        tm.assert_index_equal(union, expected)

    @pytest.mark.parametrize(
        "fname, sname, expected_name",
        [
            ("A", "A", "A"),
            ("A", "B", None),
            ("A", None, None),
            (None, "B", None),
            (None, None, None),
        ],
    )
    def test_union_unequal(self, index, fname, sname, expected_name):
        if isinstance(index, MultiIndex) or not index.is_unique:
            pytest.skip("Not for MultiIndex or repeated indices")

        # test copy.union(subset) - need sort for unicode and string
        first = index.copy().set_names(fname)
        second = index[1:].set_names(sname)
        union = first.union(second).sort_values()
        expected = index.set_names(expected_name).sort_values()
        tm.assert_index_equal(union, expected)

    @pytest.mark.parametrize(
        "fname, sname, expected_name",
        [
            ("A", "A", "A"),
            ("A", "B", None),
            ("A", None, None),
            (None, "B", None),
            (None, None, None),
        ],
    )
    def test_corner_intersect(self, index, fname, sname, expected_name):
        # GH#35847
        # Test intersections with various name combinations

        if isinstance(index, MultiIndex) or not index.is_unique:
            pytest.skip("Not for MultiIndex or repeated indices")

        # Test copy.intersection(copy)
        first = index.copy().set_names(fname)
        second = index.copy().set_names(sname)
        intersect = first.intersection(second)
        expected = index.copy().set_names(expected_name)
        tm.assert_index_equal(intersect, expected)

        # Test copy.intersection(empty)
        first = index.copy().set_names(fname)
        second = index.drop(index).set_names(sname)
        intersect = first.intersection(second)
        expected = index.drop(index).set_names(expected_name)
        tm.assert_index_equal(intersect, expected)

        # Test empty.intersection(copy)
        first = index.drop(index).set_names(fname)
        second = index.copy().set_names(sname)
        intersect = first.intersection(second)
        expected = index.drop(index).set_names(expected_name)
        tm.assert_index_equal(intersect, expected)

        # Test empty.intersection(empty)
        first = index.drop(index).set_names(fname)
        second = index.drop(index).set_names(sname)
        intersect = first.intersection(second)
        expected = index.drop(index).set_names(expected_name)
        tm.assert_index_equal(intersect, expected)

    @pytest.mark.parametrize(
        "fname, sname, expected_name",
        [
            ("A", "A", "A"),
            ("A", "B", None),
            ("A", None, None),
            (None, "B", None),
            (None, None, None),
        ],
    )
    def test_intersect_unequal(self, index, fname, sname, expected_name):
        if isinstance(index, MultiIndex) or not index.is_unique:
            pytest.skip("Not for MultiIndex or repeated indices")

        # test copy.intersection(subset) - need sort for unicode and string
        first = index.copy().set_names(fname)
        second = index[1:].set_names(sname)
        intersect = first.intersection(second).sort_values()
        expected = index[1:].set_names(expected_name).sort_values()
        tm.assert_index_equal(intersect, expected)

    def test_intersection_name_retention_with_nameless(self, index):
        if isinstance(index, MultiIndex):
            index = index.rename(list(range(index.nlevels)))
        else:
            index = index.rename("foo")

        other = np.asarray(index)

        result = index.intersection(other)
        assert result.name == index.name

        # empty other, same dtype
        result = index.intersection(other[:0])
        assert result.name == index.name

        # empty `self`
        result = index[:0].intersection(other)
        assert result.name == index.name

    def test_difference_preserves_type_empty(self, index, sort):
        # GH#20040
        # If taking difference of a set and itself, it
        # needs to preserve the type of the index
        if not index.is_unique:
            return
        result = index.difference(index, sort=sort)
        expected = index[:0]
        tm.assert_index_equal(result, expected, exact=True)

    def test_difference_name_retention_equals(self, index, sort, names):
        if isinstance(index, MultiIndex):
            names = [[x] * index.nlevels for x in names]
        index = index.rename(names[0])
        other = index.rename(names[1])

        assert index.equals(other)

        result = index.difference(other)
        expected = index[:0].rename(names[2])
        tm.assert_index_equal(result, expected)

    def test_intersection_difference_match_empty(self, index, sort):
        # GH#20040
        # Test that the intersection of an index with an
        # empty index produces the same index as the difference
        # of an index with itself.  Test for all types
        if not index.is_unique:
            return
        inter = index.intersection(index[:0])
        diff = index.difference(index, sort=sort)
        tm.assert_index_equal(inter, diff, exact=True)
