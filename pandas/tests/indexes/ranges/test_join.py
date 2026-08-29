import numpy as np
import pytest

from pandas import (
    Index,
    RangeIndex,
)
import pandas._testing as tm


class TestJoin:
    def test_join_outer(self):
        # join with Index[int64]
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index(np.arange(25, 14, -1, dtype=np.int64))

        res, lidx, ridx = index.join(other, how="outer", return_indexers=True)
        noidx_res = index.join(other, how="outer")
        tm.assert_index_equal(res, noidx_res)

        eres = Index(
            [0, 2, 4, 6, 8, 10, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
        )
        elidx = np.array(
            [0, 1, 2, 3, 4, 5, 6, 7, -1, 8, -1, 9, -1, -1, -1, -1, -1, -1, -1],
            dtype=np.intp,
        )
        eridx = np.array(
            [-1, -1, -1, -1, -1, -1, -1, -1, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
            dtype=np.intp,
        )

        assert isinstance(res, Index) and res.dtype == np.dtype(np.int64)
        assert not isinstance(res, RangeIndex)
        tm.assert_index_equal(res, eres, exact=True)
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # join with RangeIndex
        other = RangeIndex(25, 14, -1)

        res, lidx, ridx = index.join(other, how="outer", return_indexers=True)
        noidx_res = index.join(other, how="outer")
        tm.assert_index_equal(res, noidx_res)

        assert isinstance(res, Index) and res.dtype == np.int64
        assert not isinstance(res, RangeIndex)
        tm.assert_index_equal(res, eres)
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_inner(self):
        # Join with non-RangeIndex
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index(np.arange(25, 14, -1, dtype=np.int64))

        res, lidx, ridx = index.join(other, how="inner", return_indexers=True)

        # no guarantee of sortedness, so sort for comparison purposes
        ind = res.argsort()
        res = res.take(ind)
        lidx = lidx.take(ind)
        ridx = ridx.take(ind)

        eres = Index([16, 18])
        elidx = np.array([8, 9], dtype=np.intp)
        eridx = np.array([9, 7], dtype=np.intp)

        assert isinstance(res, Index) and res.dtype == np.int64
        tm.assert_index_equal(res, eres, exact="equiv")
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # Join two RangeIndex
        other = RangeIndex(25, 14, -1)

        res, lidx, ridx = index.join(other, how="inner", return_indexers=True)

        assert isinstance(res, RangeIndex)
        tm.assert_index_equal(res, eres, exact="equiv")
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_left(self):
        # Join with Index[int64]
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index(np.arange(25, 14, -1, dtype=np.int64))

        res, lidx, ridx = index.join(other, how="left", return_indexers=True)
        eres = index
        eridx = np.array([-1, -1, -1, -1, -1, -1, -1, -1, 9, 7], dtype=np.intp)

        assert isinstance(res, RangeIndex)
        tm.assert_index_equal(res, eres)
        assert lidx is None
        tm.assert_numpy_array_equal(ridx, eridx)

        # Join withRangeIndex
        other = Index(np.arange(25, 14, -1, dtype=np.int64))

        res, lidx, ridx = index.join(other, how="left", return_indexers=True)

        assert isinstance(res, RangeIndex)
        tm.assert_index_equal(res, eres)
        assert lidx is None
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_right(self):
        # Join with Index[int64]
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index(np.arange(25, 14, -1, dtype=np.int64))

        res, lidx, ridx = index.join(other, how="right", return_indexers=True)
        eres = other
        elidx = np.array([-1, -1, -1, -1, -1, -1, -1, 9, -1, 8, -1], dtype=np.intp)

        assert isinstance(other, Index) and other.dtype == np.int64
        tm.assert_index_equal(res, eres)
        tm.assert_numpy_array_equal(lidx, elidx)
        assert ridx is None

        # Join withRangeIndex
        other = RangeIndex(25, 14, -1)

        res, lidx, ridx = index.join(other, how="right", return_indexers=True)
        eres = other

        assert isinstance(other, RangeIndex)
        tm.assert_index_equal(res, eres)
        tm.assert_numpy_array_equal(lidx, elidx)
        assert ridx is None

    def test_join_non_int_index(self):
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index([3, 6, 7, 8, 10], dtype=object)

        outer = index.join(other, how="outer")
        outer2 = other.join(index, how="outer")
        expected = Index([0, 2, 3, 4, 6, 7, 8, 10, 12, 14, 16, 18])
        tm.assert_index_equal(outer, outer2)
        tm.assert_index_equal(outer, expected)

        inner = index.join(other, how="inner")
        inner2 = other.join(index, how="inner")
        expected = Index([6, 8, 10])
        tm.assert_index_equal(inner, inner2)
        tm.assert_index_equal(inner, expected)

        left = index.join(other, how="left")
        tm.assert_index_equal(left, index.astype(object))

        left2 = other.join(index, how="left")
        tm.assert_index_equal(left2, other)

        right = index.join(other, how="right")
        tm.assert_index_equal(right, other)

        right2 = other.join(index, how="right")
        tm.assert_index_equal(right2, index.astype(object))

    def test_join_non_unique(self):
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index([4, 4, 3, 3])

        res, lidx, ridx = index.join(other, return_indexers=True)

        eres = Index([0, 2, 4, 4, 6, 8, 10, 12, 14, 16, 18])
        elidx = np.array([0, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.intp)
        eridx = np.array([-1, -1, 0, 1, -1, -1, -1, -1, -1, -1, -1], dtype=np.intp)

        tm.assert_index_equal(res, eres)
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_self(self, join_type):
        index = RangeIndex(start=0, stop=20, step=2)
        joined = index.join(index, how=join_type)
        assert index is joined


@pytest.mark.parametrize(
    "left, right, expected, expected_lidx, expected_ridx, how",
    [
        [RangeIndex(2), RangeIndex(3), RangeIndex(2), None, [0, 1], "left"],
        [RangeIndex(2), RangeIndex(2), RangeIndex(2), None, None, "left"],
        [RangeIndex(2), RangeIndex(20, 22), RangeIndex(2), None, [-1, -1], "left"],
        [RangeIndex(2), RangeIndex(3), RangeIndex(3), [0, 1, -1], None, "right"],
        [RangeIndex(2), RangeIndex(2), RangeIndex(2), None, None, "right"],
        [
            RangeIndex(2),
            RangeIndex(20, 22),
            RangeIndex(20, 22),
            [-1, -1],
            None,
            "right",
        ],
        [RangeIndex(2), RangeIndex(3), RangeIndex(2), [0, 1], [0, 1], "inner"],
        [RangeIndex(2), RangeIndex(2), RangeIndex(2), None, None, "inner"],
        [RangeIndex(2), RangeIndex(1, 3), RangeIndex(1, 2), [1], [0], "inner"],
        [RangeIndex(2), RangeIndex(3), RangeIndex(3), [0, 1, -1], [0, 1, 2], "outer"],
        [RangeIndex(2), RangeIndex(2), RangeIndex(2), None, None, "outer"],
        [
            RangeIndex(2),
            RangeIndex(2, 4),
            RangeIndex(4),
            [0, 1, -1, -1],
            [-1, -1, 0, 1],
            "outer",
        ],
        [RangeIndex(2), RangeIndex(0), RangeIndex(2), None, [-1, -1], "left"],
        [RangeIndex(2), RangeIndex(0), RangeIndex(0), [], None, "right"],
        [RangeIndex(2), RangeIndex(0), RangeIndex(0), [], None, "inner"],
        [RangeIndex(2), RangeIndex(0), RangeIndex(2), None, [-1, -1], "outer"],
    ],
)
@pytest.mark.parametrize(
    "right_type", [RangeIndex, lambda x: Index(list(x), dtype=x.dtype)]
)
def test_join_preserves_rangeindex(
    left, right, expected, expected_lidx, expected_ridx, how, right_type
):
    result, lidx, ridx = left.join(right_type(right), how=how, return_indexers=True)
    tm.assert_index_equal(result, expected, exact=True)

    if expected_lidx is None:
        assert lidx is expected_lidx
    else:
        exp_lidx = np.array(expected_lidx, dtype=np.intp)
        tm.assert_numpy_array_equal(lidx, exp_lidx)

    if expected_ridx is None:
        assert ridx is expected_ridx
    else:
        exp_ridx = np.array(expected_ridx, dtype=np.intp)
        tm.assert_numpy_array_equal(ridx, exp_ridx)


def test_rangeindex_can_use_libjoin_is_false():
    # GH#54646 - RangeIndex should not use libjoin fastpath
    # because _get_join_target allocates memory via np.arange
    index = RangeIndex(5)
    assert not index._can_use_libjoin


def test_join_disjoint_rangeindexes_outer():
    # GH#54646 - outer join of disjoint RangeIndexes should be sorted
    idx1 = RangeIndex(3, 6)  # [3, 4, 5]
    idx2 = RangeIndex(0, 3)  # [0, 1, 2]

    result, lidx, ridx = idx1.join(idx2, how="outer", return_indexers=True)

    expected = RangeIndex(0, 6)
    tm.assert_index_equal(result, expected, exact=True)
    exp_lidx = np.array([-1, -1, -1, 0, 1, 2], dtype=np.intp)
    exp_ridx = np.array([0, 1, 2, -1, -1, -1], dtype=np.intp)
    tm.assert_numpy_array_equal(lidx, exp_lidx)
    tm.assert_numpy_array_equal(ridx, exp_ridx)


def test_join_overlapping_rangeindexes_outer():
    # GH#54646
    idx1 = RangeIndex(0, 5)  # [0, 1, 2, 3, 4]
    idx2 = RangeIndex(3, 8)  # [3, 4, 5, 6, 7]

    result, lidx, ridx = idx1.join(idx2, how="outer", return_indexers=True)

    expected = RangeIndex(0, 8)
    tm.assert_index_equal(result, expected, exact=True)
    exp_lidx = np.array([0, 1, 2, 3, 4, -1, -1, -1], dtype=np.intp)
    exp_ridx = np.array([-1, -1, -1, 0, 1, 2, 3, 4], dtype=np.intp)
    tm.assert_numpy_array_equal(lidx, exp_lidx)
    tm.assert_numpy_array_equal(ridx, exp_ridx)


def test_join_rangeindex_with_int_index():
    # GH#54646 - RangeIndex join with non-RangeIndex should still work
    idx1 = RangeIndex(5)
    idx2 = Index(np.arange(3, 8, dtype=np.int64))

    result, lidx, ridx = idx1.join(idx2, how="inner", return_indexers=True)

    expected = RangeIndex(3, 5)
    tm.assert_index_equal(result, expected, exact=True)


def test_join_rangeindex_self():
    # GH#54646 - joining a RangeIndex with itself should return self
    index = RangeIndex(10)
    for how in ["left", "right", "inner", "outer"]:
        result = index.join(index, how=how)
        assert result is index


def test_join_rangeindex_unconvertible_fallback():
    # GH#54646 - when other cannot be converted to RangeIndex,
    # RangeIndex._join_monotonic raises NotImplementedError and
    # falls back to _join_via_get_indexer
    idx1 = RangeIndex(5)
    idx2 = Index([0, 1, 3, 7, 9], dtype=np.int64)  # not a valid range

    result, lidx, ridx = idx1.join(idx2, how="inner", return_indexers=True)

    expected = Index([0, 1, 3], dtype=np.int64)
    tm.assert_index_equal(result, expected)
