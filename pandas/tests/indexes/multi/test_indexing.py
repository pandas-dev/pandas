from collections import namedtuple
from datetime import timedelta
import re

import numpy as np
import pytest

from pandas._libs import index as libindex
from pandas.compat import PY312
from pandas.errors import InvalidIndexError

import pandas as pd
import pandas._testing as tm


class TestSliceLocs:
    def test_slice_locs_partial(self, idx):
        sorted_idx, _ = idx.sortlevel(0)

        result = sorted_idx.slice_locs(("foo", "two"), ("qux", "one"))
        assert result == (1, 5)

        result = sorted_idx.slice_locs(None, ("qux", "one"))
        assert result == (0, 5)

        result = sorted_idx.slice_locs(("foo", "two"), None)
        assert result == (1, len(sorted_idx))

        result = sorted_idx.slice_locs("bar", "baz")
        assert result == (2, 4)

    def test_slice_locs(self):
        df = pd.DataFrame(
            np.random.default_rng(2).standard_normal((50, 4)),
            columns=pd.Index(list("ABCD"), dtype=object),
            index=pd.date_range("2000-01-01", periods=50, freq="B"),
        )
        stacked = df.stack()
        idx = stacked.index

        slob = slice(*idx.slice_locs(df.index[5], df.index[15]))
        sliced = stacked[slob]
        expected = df[5:16].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

        slob = slice(
            *idx.slice_locs(
                df.index[5] + timedelta(seconds=30),
                df.index[15] - timedelta(seconds=30),
            )
        )
        sliced = stacked[slob]
        expected = df[6:15].stack()
        tm.assert_almost_equal(sliced.values, expected.values)

    def test_slice_locs_with_type_mismatch(self):
        df = pd.DataFrame(
            np.random.default_rng(2).standard_normal((10, 4)),
            columns=pd.Index(list("ABCD"), dtype=object),
            index=pd.date_range("2000-01-01", periods=10, freq="B"),
        )
        stacked = df.stack()
        idx = stacked.index
        with pytest.raises(TypeError, match="^Level type mismatch"):
            idx.slice_locs((1, 3))
        with pytest.raises(TypeError, match="^Level type mismatch"):
            idx.slice_locs(df.index[5] + timedelta(seconds=30), (5, 2))
        df = pd.DataFrame(
            np.ones((5, 5)),
            index=pd.Index([f"i-{i}" for i in range(5)], name="a"),
            columns=pd.Index([f"i-{i}" for i in range(5)], name="a"),
        )
        stacked = df.stack()
        idx = stacked.index
        with pytest.raises(TypeError, match="^Level type mismatch"):
            idx.slice_locs(timedelta(seconds=30))
        # TODO: Try creating a UnicodeDecodeError in exception message
        with pytest.raises(TypeError, match="^Level type mismatch"):
            idx.slice_locs(df.index[1], (16, "a"))

    def test_slice_locs_not_sorted(self):
        index = pd.MultiIndex(
            levels=[
                pd.Index(np.arange(4)),
                pd.Index(np.arange(4)),
                pd.Index(np.arange(4)),
            ],
            codes=[
                np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                np.array([1, 0, 1, 1, 0, 0, 1, 0]),
            ],
        )
        msg = "Key length.*was greater than MultiIndex lexsort depth"
        with pytest.raises(KeyError, match=msg):
            index.slice_locs((1, 0, 1), (2, 1, 0))

        # works
        sorted_index, _ = index.sortlevel(0)
        # should there be a test case here???
        sorted_index.slice_locs((1, 0, 1), (2, 1, 0))

    def test_slice_locs_not_contained(self):
        # some searchsorted action

        index = pd.MultiIndex(
            levels=[[0, 2, 4, 6], [0, 2, 4]],
            codes=[[0, 0, 0, 1, 1, 2, 3, 3, 3], [0, 1, 2, 1, 2, 2, 0, 1, 2]],
        )

        result = index.slice_locs((1, 0), (5, 2))
        assert result == (3, 6)

        result = index.slice_locs(1, 5)
        assert result == (3, 6)

        result = index.slice_locs((2, 2), (5, 2))
        assert result == (3, 6)

        result = index.slice_locs(2, 5)
        assert result == (3, 6)

        result = index.slice_locs((1, 0), (6, 3))
        assert result == (3, 8)

        result = index.slice_locs(-1, 10)
        assert result == (0, len(index))

    @pytest.mark.parametrize(
        "index_arr,expected,start_idx,end_idx",
        [
            ([[np.nan, "a", "b"], ["c", "d", "e"]], (0, 3), np.nan, None),
            ([[np.nan, "a", "b"], ["c", "d", "e"]], (0, 3), np.nan, "b"),
            ([[np.nan, "a", "b"], ["c", "d", "e"]], (0, 3), np.nan, ("b", "e")),
            ([["a", "b", "c"], ["d", np.nan, "e"]], (1, 3), ("b", np.nan), None),
            ([["a", "b", "c"], ["d", np.nan, "e"]], (1, 3), ("b", np.nan), "c"),
            ([["a", "b", "c"], ["d", np.nan, "e"]], (1, 3), ("b", np.nan), ("c", "e")),
        ],
    )
    def test_slice_locs_with_missing_value(
        self, index_arr, expected, start_idx, end_idx
    ):
        # issue 19132
        idx = pd.MultiIndex.from_arrays(index_arr)
        result = idx.slice_locs(start=start_idx, end=end_idx)
        assert result == expected


class TestPutmask:
    def test_putmask_with_wrong_mask(self, idx):
        # GH18368

        msg = "putmask: mask and data must be the same size"
        with pytest.raises(ValueError, match=msg):
            idx.putmask(np.ones(len(idx) + 1, np.bool_), 1)

        with pytest.raises(ValueError, match=msg):
            idx.putmask(np.ones(len(idx) - 1, np.bool_), 1)

        with pytest.raises(ValueError, match=msg):
            idx.putmask("foo", 1)

    def test_putmask_multiindex_other(self):
        # GH#43212 `value` is also a MultiIndex

        left = pd.MultiIndex.from_tuples([(np.nan, 6), (np.nan, 6), ("a", 4)])
        right = pd.MultiIndex.from_tuples([("a", 1), ("a", 1), ("d", 1)])
        mask = np.array([True, True, False])

        result = left.putmask(mask, right)

        expected = pd.MultiIndex.from_tuples([right[0], right[1], left[2]])
        tm.assert_index_equal(result, expected)

    def test_putmask_keep_dtype(self, any_numeric_ea_dtype):
        # GH#49830
        midx = pd.MultiIndex.from_arrays(
            [pd.Series([1, 2, 3], dtype=any_numeric_ea_dtype), [10, 11, 12]]
        )
        midx2 = pd.MultiIndex.from_arrays(
            [pd.Series([5, 6, 7], dtype=any_numeric_ea_dtype), [-1, -2, -3]]
        )
        result = midx.putmask([True, False, False], midx2)
        expected = pd.MultiIndex.from_arrays(
            [pd.Series([5, 2, 3], dtype=any_numeric_ea_dtype), [-1, 11, 12]]
        )
        tm.assert_index_equal(result, expected)

    def test_putmask_keep_dtype_shorter_value(self, any_numeric_ea_dtype):
        # GH#49830
        midx = pd.MultiIndex.from_arrays(
            [pd.Series([1, 2, 3], dtype=any_numeric_ea_dtype), [10, 11, 12]]
        )
        midx2 = pd.MultiIndex.from_arrays(
            [pd.Series([5], dtype=any_numeric_ea_dtype), [-1]]
        )
        result = midx.putmask([True, False, False], midx2)
        expected = pd.MultiIndex.from_arrays(
            [pd.Series([5, 2, 3], dtype=any_numeric_ea_dtype), [-1, 11, 12]]
        )
        tm.assert_index_equal(result, expected)


class TestGetIndexer:
    def test_get_indexer(self):
        major_axis = pd.Index(np.arange(4))
        minor_axis = pd.Index(np.arange(2))

        major_codes = np.array([0, 0, 1, 2, 2, 3, 3], dtype=np.intp)
        minor_codes = np.array([0, 1, 0, 0, 1, 0, 1], dtype=np.intp)

        index = pd.MultiIndex(
            levels=[major_axis, minor_axis], codes=[major_codes, minor_codes]
        )
        idx1 = index[:5]
        idx2 = index[[1, 3, 5]]

        r1 = idx1.get_indexer(idx2)
        tm.assert_almost_equal(r1, np.array([1, 3, -1], dtype=np.intp))

        r1 = idx2.get_indexer(idx1, method="pad")
        e1 = np.array([-1, 0, 0, 1, 1], dtype=np.intp)
        tm.assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method="pad")
        tm.assert_almost_equal(r2, e1[::-1])

        rffill1 = idx2.get_indexer(idx1, method="ffill")
        tm.assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method="backfill")
        e1 = np.array([0, 0, 1, 1, 2], dtype=np.intp)
        tm.assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method="backfill")
        tm.assert_almost_equal(r2, e1[::-1])

        rbfill1 = idx2.get_indexer(idx1, method="bfill")
        tm.assert_almost_equal(r1, rbfill1)

        # pass non-MultiIndex
        r1 = idx1.get_indexer(idx2.values)
        rexp1 = idx1.get_indexer(idx2)
        tm.assert_almost_equal(r1, rexp1)

        r1 = idx1.get_indexer([1, 2, 3])
        assert (r1 == [-1, -1, -1]).all()

        # create index with duplicates
        idx1 = pd.Index(list(range(10)) + list(range(10)))
        idx2 = pd.Index(list(range(20)))

        msg = "Reindexing only valid with uniquely valued Index objects"
        with pytest.raises(InvalidIndexError, match=msg):
            idx1.get_indexer(idx2)

    def test_get_indexer_nearest(self):
        midx = pd.MultiIndex.from_tuples([("a", 1), ("b", 2)])
        msg = (
            "method='nearest' not implemented yet for MultiIndex; see GitHub issue 9365"
        )
        with pytest.raises(NotImplementedError, match=msg):
            midx.get_indexer(["a"], method="nearest")
        msg = "tolerance not implemented yet for MultiIndex"
        with pytest.raises(NotImplementedError, match=msg):
            midx.get_indexer(["a"], method="pad", tolerance=2)

    def test_get_indexer_categorical_time(self):
        # https://github.com/pandas-dev/pandas/issues/21390
        midx = pd.MultiIndex.from_product(
            [
                pd.Categorical(["a", "b", "c"]),
                pd.Categorical(pd.date_range("2012-01-01", periods=3, freq="h")),
            ]
        )
        result = midx.get_indexer(midx)
        tm.assert_numpy_array_equal(result, np.arange(9, dtype=np.intp))

    @pytest.mark.parametrize(
        "index_arr,labels,expected",
        [
            (
                [[1, np.nan, 2], [3, 4, 5]],
                [1, np.nan, 2],
                np.array([-1, -1, -1], dtype=np.intp),
            ),
            ([[1, np.nan, 2], [3, 4, 5]], [(np.nan, 4)], np.array([1], dtype=np.intp)),
            ([[1, 2, 3], [np.nan, 4, 5]], [(1, np.nan)], np.array([0], dtype=np.intp)),
            (
                [[1, 2, 3], [np.nan, 4, 5]],
                [np.nan, 4, 5],
                np.array([-1, -1, -1], dtype=np.intp),
            ),
        ],
    )
    def test_get_indexer_with_missing_value(self, index_arr, labels, expected):
        # issue 19132
        idx = pd.MultiIndex.from_arrays(index_arr)
        result = idx.get_indexer(labels)
        tm.assert_numpy_array_equal(result, expected)

    def test_get_indexer_methods(self):
        # https://github.com/pandas-dev/pandas/issues/29896
        # test getting an indexer for another index with different methods
        # confirms that getting an indexer without a filling method, getting an
        # indexer and backfilling, and getting an indexer and padding all behave
        # correctly in the case where all of the target values fall in between
        # several levels in the MultiIndex into which they are getting an indexer
        #
        # visually, the MultiIndexes used in this test are:
        # mult_idx_1:
        #  0: -1 0
        #  1:    2
        #  2:    3
        #  3:    4
        #  4:  0 0
        #  5:    2
        #  6:    3
        #  7:    4
        #  8:  1 0
        #  9:    2
        # 10:    3
        # 11:    4
        #
        # mult_idx_2:
        #  0: 0 1
        #  1:   3
        #  2:   4
        mult_idx_1 = pd.MultiIndex.from_product([[-1, 0, 1], [0, 2, 3, 4]])
        mult_idx_2 = pd.MultiIndex.from_product([[0], [1, 3, 4]])

        indexer = mult_idx_1.get_indexer(mult_idx_2)
        expected = np.array([-1, 6, 7], dtype=indexer.dtype)
        tm.assert_almost_equal(expected, indexer)

        backfill_indexer = mult_idx_1.get_indexer(mult_idx_2, method="backfill")
        expected = np.array([5, 6, 7], dtype=backfill_indexer.dtype)
        tm.assert_almost_equal(expected, backfill_indexer)

        # ensure the legacy "bfill" option functions identically to "backfill"
        backfill_indexer = mult_idx_1.get_indexer(mult_idx_2, method="bfill")
        expected = np.array([5, 6, 7], dtype=backfill_indexer.dtype)
        tm.assert_almost_equal(expected, backfill_indexer)

        pad_indexer = mult_idx_1.get_indexer(mult_idx_2, method="pad")
        expected = np.array([4, 6, 7], dtype=pad_indexer.dtype)
        tm.assert_almost_equal(expected, pad_indexer)

        # ensure the legacy "ffill" option functions identically to "pad"
        pad_indexer = mult_idx_1.get_indexer(mult_idx_2, method="ffill")
        expected = np.array([4, 6, 7], dtype=pad_indexer.dtype)
        tm.assert_almost_equal(expected, pad_indexer)

    @pytest.mark.parametrize("method", ["pad", "ffill", "backfill", "bfill", "nearest"])
    def test_get_indexer_methods_raise_for_non_monotonic(self, method):
        # 53452
        mi = pd.MultiIndex.from_arrays([[0, 4, 2], [0, 4, 2]])
        if method == "nearest":
            err = NotImplementedError
            msg = "not implemented yet for MultiIndex"
        else:
            err = ValueError
            msg = "index must be monotonic increasing or decreasing"
        with pytest.raises(err, match=msg):
            mi.get_indexer([(1, 1)], method=method)

    def test_get_indexer_three_or_more_levels(self):
        # https://github.com/pandas-dev/pandas/issues/29896
        # tests get_indexer() on MultiIndexes with 3+ levels
        # visually, these are
        # mult_idx_1:
        #  0: 1 2 5
        #  1:     7
        #  2:   4 5
        #  3:     7
        #  4:   6 5
        #  5:     7
        #  6: 3 2 5
        #  7:     7
        #  8:   4 5
        #  9:     7
        # 10:   6 5
        # 11:     7
        #
        # mult_idx_2:
        #  0: 1 1 8
        #  1: 1 5 9
        #  2: 1 6 7
        #  3: 2 1 6
        #  4: 2 7 6
        #  5: 2 7 8
        #  6: 3 6 8
        mult_idx_1 = pd.MultiIndex.from_product([[1, 3], [2, 4, 6], [5, 7]])
        mult_idx_2 = pd.MultiIndex.from_tuples(
            [
                (1, 1, 8),
                (1, 5, 9),
                (1, 6, 7),
                (2, 1, 6),
                (2, 7, 7),
                (2, 7, 8),
                (3, 6, 8),
            ]
        )
        # sanity check
        assert mult_idx_1.is_monotonic_increasing
        assert mult_idx_1.is_unique
        assert mult_idx_2.is_monotonic_increasing
        assert mult_idx_2.is_unique

        # show the relationships between the two
        assert mult_idx_2[0] < mult_idx_1[0]
        assert mult_idx_1[3] < mult_idx_2[1] < mult_idx_1[4]
        assert mult_idx_1[5] == mult_idx_2[2]
        assert mult_idx_1[5] < mult_idx_2[3] < mult_idx_1[6]
        assert mult_idx_1[5] < mult_idx_2[4] < mult_idx_1[6]
        assert mult_idx_1[5] < mult_idx_2[5] < mult_idx_1[6]
        assert mult_idx_1[-1] < mult_idx_2[6]

        indexer_no_fill = mult_idx_1.get_indexer(mult_idx_2)
        expected = np.array([-1, -1, 5, -1, -1, -1, -1], dtype=indexer_no_fill.dtype)
        tm.assert_almost_equal(expected, indexer_no_fill)

        # test with backfilling
        indexer_backfilled = mult_idx_1.get_indexer(mult_idx_2, method="backfill")
        expected = np.array([0, 4, 5, 6, 6, 6, -1], dtype=indexer_backfilled.dtype)
        tm.assert_almost_equal(expected, indexer_backfilled)

        # now, the same thing, but forward-filled (aka "padded")
        indexer_padded = mult_idx_1.get_indexer(mult_idx_2, method="pad")
        expected = np.array([-1, 3, 5, 5, 5, 5, 11], dtype=indexer_padded.dtype)
        tm.assert_almost_equal(expected, indexer_padded)

        # now, do the indexing in the other direction
        assert mult_idx_2[0] < mult_idx_1[0] < mult_idx_2[1]
        assert mult_idx_2[0] < mult_idx_1[1] < mult_idx_2[1]
        assert mult_idx_2[0] < mult_idx_1[2] < mult_idx_2[1]
        assert mult_idx_2[0] < mult_idx_1[3] < mult_idx_2[1]
        assert mult_idx_2[1] < mult_idx_1[4] < mult_idx_2[2]
        assert mult_idx_2[2] == mult_idx_1[5]
        assert mult_idx_2[5] < mult_idx_1[6] < mult_idx_2[6]
        assert mult_idx_2[5] < mult_idx_1[7] < mult_idx_2[6]
        assert mult_idx_2[5] < mult_idx_1[8] < mult_idx_2[6]
        assert mult_idx_2[5] < mult_idx_1[9] < mult_idx_2[6]
        assert mult_idx_2[5] < mult_idx_1[10] < mult_idx_2[6]
        assert mult_idx_2[5] < mult_idx_1[11] < mult_idx_2[6]

        indexer = mult_idx_2.get_indexer(mult_idx_1)
        expected = np.array(
            [-1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1], dtype=indexer.dtype
        )
        tm.assert_almost_equal(expected, indexer)

        backfill_indexer = mult_idx_2.get_indexer(mult_idx_1, method="bfill")
        expected = np.array(
            [1, 1, 1, 1, 2, 2, 6, 6, 6, 6, 6, 6], dtype=backfill_indexer.dtype
        )
        tm.assert_almost_equal(expected, backfill_indexer)

        pad_indexer = mult_idx_2.get_indexer(mult_idx_1, method="pad")
        expected = np.array(
            [0, 0, 0, 0, 1, 2, 5, 5, 5, 5, 5, 5], dtype=pad_indexer.dtype
        )
        tm.assert_almost_equal(expected, pad_indexer)

    def test_get_indexer_crossing_levels(self):
        # https://github.com/pandas-dev/pandas/issues/29896
        # tests a corner case with get_indexer() with MultiIndexes where, when we
        # need to "carry" across levels, proper tuple ordering is respected
        #
        # the MultiIndexes used in this test, visually, are:
        # mult_idx_1:
        #  0: 1 1 1 1
        #  1:       2
        #  2:     2 1
        #  3:       2
        #  4: 1 2 1 1
        #  5:       2
        #  6:     2 1
        #  7:       2
        #  8: 2 1 1 1
        #  9:       2
        # 10:     2 1
        # 11:       2
        # 12: 2 2 1 1
        # 13:       2
        # 14:     2 1
        # 15:       2
        #
        # mult_idx_2:
        #  0: 1 3 2 2
        #  1: 2 3 2 2
        mult_idx_1 = pd.MultiIndex.from_product([[1, 2]] * 4)
        mult_idx_2 = pd.MultiIndex.from_tuples([(1, 3, 2, 2), (2, 3, 2, 2)])

        # show the tuple orderings, which get_indexer() should respect
        assert mult_idx_1[7] < mult_idx_2[0] < mult_idx_1[8]
        assert mult_idx_1[-1] < mult_idx_2[1]

        indexer = mult_idx_1.get_indexer(mult_idx_2)
        expected = np.array([-1, -1], dtype=indexer.dtype)
        tm.assert_almost_equal(expected, indexer)

        backfill_indexer = mult_idx_1.get_indexer(mult_idx_2, method="bfill")
        expected = np.array([8, -1], dtype=backfill_indexer.dtype)
        tm.assert_almost_equal(expected, backfill_indexer)

        pad_indexer = mult_idx_1.get_indexer(mult_idx_2, method="ffill")
        expected = np.array([7, 15], dtype=pad_indexer.dtype)
        tm.assert_almost_equal(expected, pad_indexer)

    def test_get_indexer_kwarg_validation(self):
        # GH#41918
        mi = pd.MultiIndex.from_product([range(3), ["A", "B"]])

        msg = "limit argument only valid if doing pad, backfill or nearest"
        with pytest.raises(ValueError, match=msg):
            mi.get_indexer(mi[:-1], limit=4)

        msg = "tolerance argument only valid if doing pad, backfill or nearest"
        with pytest.raises(ValueError, match=msg):
            mi.get_indexer(mi[:-1], tolerance="piano")

    def test_get_indexer_nan(self):
        # GH#37222
        idx1 = pd.MultiIndex.from_product([["A"], [1.0, 2.0]], names=["id1", "id2"])
        idx2 = pd.MultiIndex.from_product([["A"], [np.nan, 2.0]], names=["id1", "id2"])
        expected = np.array([-1, 1])
        result = idx2.get_indexer(idx1)
        tm.assert_numpy_array_equal(result, expected, check_dtype=False)
        result = idx1.get_indexer(idx2)
        tm.assert_numpy_array_equal(result, expected, check_dtype=False)


def test_getitem(idx):
    # scalar
    assert idx[2] == ("bar", "one")

    # slice
    result = idx[2:5]
    expected = idx[[2, 3, 4]]
    assert result.equals(expected)

    # boolean
    result = idx[[True, False, True, False, True, True]]
    result2 = idx[np.array([True, False, True, False, True, True])]
    expected = idx[[0, 2, 4, 5]]
    assert result.equals(expected)
    assert result2.equals(expected)


def test_getitem_group_select(idx):
    sorted_idx, _ = idx.sortlevel(0)
    assert sorted_idx.get_loc("baz") == slice(3, 4)
    assert sorted_idx.get_loc("foo") == slice(0, 2)


@pytest.mark.parametrize("box", [list, pd.Index])
def test_getitem_bool_index_all(box):
    # GH#22533
    ind1 = box([True] * 5)
    idx = pd.MultiIndex.from_tuples([(10, 1), (20, 2), (30, 3), (40, 4), (50, 5)])
    tm.assert_index_equal(idx[ind1], idx)

    ind2 = box([True, False, True, False, False])
    expected = pd.MultiIndex.from_tuples([(10, 1), (30, 3)])
    tm.assert_index_equal(idx[ind2], expected)


@pytest.mark.parametrize("box", [list, pd.Index])
def test_getitem_bool_index_single(box):
    # GH#22533
    ind1 = box([True])
    idx = pd.MultiIndex.from_tuples([(10, 1)])
    tm.assert_index_equal(idx[ind1], idx)

    ind2 = box([False])
    expected = pd.MultiIndex(
        levels=[np.array([], dtype=np.int64), np.array([], dtype=np.int64)],
        codes=[[], []],
    )
    tm.assert_index_equal(idx[ind2], expected)


class TestGetLoc:
    def test_get_loc(self, idx):
        assert idx.get_loc(("foo", "two")) == 1
        assert idx.get_loc(("baz", "two")) == 3
        with pytest.raises(KeyError, match=r"^\('bar', 'two'\)$"):
            idx.get_loc(("bar", "two"))
        with pytest.raises(KeyError, match=r"^'quux'$"):
            idx.get_loc("quux")

        # 3 levels
        index = pd.MultiIndex(
            levels=[
                pd.Index(np.arange(4)),
                pd.Index(np.arange(4)),
                pd.Index(np.arange(4)),
            ],
            codes=[
                np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                np.array([1, 0, 1, 1, 0, 0, 1, 0]),
            ],
        )
        with pytest.raises(KeyError, match=r"^\(1, 1\)$"):
            index.get_loc((1, 1))
        assert index.get_loc((2, 0)) == slice(3, 5)

    def test_get_loc_duplicates(self):
        index = pd.Index([2, 2, 2, 2])
        result = index.get_loc(2)
        expected = slice(0, 4)
        assert result == expected

        index = pd.Index(["c", "a", "a", "b", "b"])
        rs = index.get_loc("c")
        xp = 0
        assert rs == xp

        with pytest.raises(KeyError, match="2"):
            index.get_loc(2)

    def test_get_loc_level(self):
        index = pd.MultiIndex(
            levels=[
                pd.Index(np.arange(4)),
                pd.Index(np.arange(4)),
                pd.Index(np.arange(4)),
            ],
            codes=[
                np.array([0, 0, 1, 2, 2, 2, 3, 3]),
                np.array([0, 1, 0, 0, 0, 1, 0, 1]),
                np.array([1, 0, 1, 1, 0, 0, 1, 0]),
            ],
        )
        loc, new_index = index.get_loc_level((0, 1))
        expected = slice(1, 2)
        exp_index = index[expected].droplevel(0).droplevel(0)
        assert loc == expected
        assert new_index.equals(exp_index)

        loc, new_index = index.get_loc_level((0, 1, 0))
        expected = 1
        assert loc == expected
        assert new_index is None

        with pytest.raises(KeyError, match=r"^\(2, 2\)$"):
            index.get_loc_level((2, 2))
        # GH 22221: unused label
        with pytest.raises(KeyError, match=r"^2$"):
            index.drop(2).get_loc_level(2)
        # Unused label on unsorted level:
        with pytest.raises(KeyError, match=r"^2$"):
            index.drop(1, level=2).get_loc_level(2, level=2)

        index = pd.MultiIndex(
            levels=[[2000], list(range(4))],
            codes=[np.array([0, 0, 0, 0]), np.array([0, 1, 2, 3])],
        )
        result, new_index = index.get_loc_level((2000, slice(None, None)))
        expected = slice(None, None)
        assert result == expected
        assert new_index.equals(index.droplevel(0))

    @pytest.mark.parametrize("dtype1", [int, float, bool, str])
    @pytest.mark.parametrize("dtype2", [int, float, bool, str])
    def test_get_loc_multiple_dtypes(self, dtype1, dtype2):
        # GH 18520
        levels = [np.array([0, 1]).astype(dtype1), np.array([0, 1]).astype(dtype2)]
        idx = pd.MultiIndex.from_product(levels)
        assert idx.get_loc(idx[2]) == 2

    @pytest.mark.parametrize("level", [0, 1])
    @pytest.mark.parametrize("dtypes", [[int, float], [float, int]])
    def test_get_loc_implicit_cast(self, level, dtypes):
        # GH 18818, GH 15994 : as flat index, cast int to float and vice-versa
        levels = [["a", "b"], ["c", "d"]]
        key = ["b", "d"]
        lev_dtype, key_dtype = dtypes
        levels[level] = np.array([0, 1], dtype=lev_dtype)
        key[level] = key_dtype(1)
        idx = pd.MultiIndex.from_product(levels)
        assert idx.get_loc(tuple(key)) == 3

    @pytest.mark.parametrize("dtype", [bool, object])
    def test_get_loc_cast_bool(self, dtype):
        # GH 19086 : int is casted to bool, but not vice-versa (for object dtype)
        #  With bool dtype, we don't cast in either direction.
        levels = [pd.Index([False, True], dtype=dtype), np.arange(2, dtype="int64")]
        idx = pd.MultiIndex.from_product(levels)

        if dtype is bool:
            with pytest.raises(KeyError, match=r"^\(0, 1\)$"):
                assert idx.get_loc((0, 1)) == 1
            with pytest.raises(KeyError, match=r"^\(1, 0\)$"):
                assert idx.get_loc((1, 0)) == 2
        else:
            # We use python object comparisons, which treat 0 == False and 1 == True
            assert idx.get_loc((0, 1)) == 1
            assert idx.get_loc((1, 0)) == 2

        with pytest.raises(KeyError, match=r"^\(False, True\)$"):
            idx.get_loc((False, True))
        with pytest.raises(KeyError, match=r"^\(True, False\)$"):
            idx.get_loc((True, False))

    @pytest.mark.parametrize("level", [0, 1])
    def test_get_loc_nan(self, level, nulls_fixture):
        # GH 18485 : NaN in MultiIndex
        levels = [["a", "b"], ["c", "d"]]
        key = ["b", "d"]
        levels[level] = np.array([0, nulls_fixture], dtype=type(nulls_fixture))
        key[level] = nulls_fixture
        idx = pd.MultiIndex.from_product(levels)
        assert idx.get_loc(tuple(key)) == 3

    def test_get_loc_missing_nan(self):
        # GH 8569
        idx = pd.MultiIndex.from_arrays([[1.0, 2.0], [3.0, 4.0]])
        assert isinstance(idx.get_loc(1), slice)
        with pytest.raises(KeyError, match=r"^3$"):
            idx.get_loc(3)
        with pytest.raises(KeyError, match=r"^nan$"):
            idx.get_loc(np.nan)
        with pytest.raises(InvalidIndexError, match=r"\[nan\]"):
            # listlike/non-hashable raises TypeError
            idx.get_loc([np.nan])

    def test_get_loc_with_values_including_missing_values(self):
        # issue 19132
        idx = pd.MultiIndex.from_product([[np.nan, 1]] * 2)
        expected = slice(0, 2, None)
        assert idx.get_loc(np.nan) == expected

        idx = pd.MultiIndex.from_arrays([[np.nan, 1, 2, np.nan]])
        expected = np.array([True, False, False, True])
        tm.assert_numpy_array_equal(idx.get_loc(np.nan), expected)

        idx = pd.MultiIndex.from_product([[np.nan, 1]] * 3)
        expected = slice(2, 4, None)
        assert idx.get_loc((np.nan, 1)) == expected

    def test_get_loc_duplicates2(self):
        # TODO: de-duplicate with test_get_loc_duplicates above?
        index = pd.MultiIndex(
            levels=[["D", "B", "C"], [0, 26, 27, 37, 57, 67, 75, 82]],
            codes=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2], [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
            names=["tag", "day"],
        )

        assert index.get_loc("D") == slice(0, 3)

    def test_get_loc_past_lexsort_depth(self, performance_warning):
        # GH#30053
        idx = pd.MultiIndex(
            levels=[["a"], [0, 7], [1]],
            codes=[[0, 0], [1, 0], [0, 0]],
            names=["x", "y", "z"],
            sortorder=0,
        )
        key = ("a", 7)

        msg = "indexing past lexsort depth may impact performance"
        with tm.assert_produces_warning(performance_warning, match=msg):
            # PerformanceWarning: indexing past lexsort depth may impact performance
            result = idx.get_loc(key)

        assert result == slice(0, 1, None)

    def test_multiindex_get_loc_list_raises(self):
        # GH#35878
        idx = pd.MultiIndex.from_tuples([("a", 1), ("b", 2)])
        msg = r"\[\]"
        with pytest.raises(InvalidIndexError, match=msg):
            idx.get_loc([])

    def test_get_loc_nested_tuple_raises_keyerror(self):
        # raise KeyError, not TypeError
        mi = pd.MultiIndex.from_product([range(3), range(4), range(5), range(6)])
        key = ((2, 3, 4), "foo")

        with pytest.raises(KeyError, match=re.escape(str(key))):
            mi.get_loc(key)

    def test_get_loc_string_key_on_numeric_level(self):
        # GH#60104 - string key on integer level should raise KeyError,
        # not return wrong results due to numpy searchsorted cross-dtype bug
        mi = pd.MultiIndex.from_product([[2000, 2001], ["a", "b"], ["x", "y"]])
        df = pd.DataFrame({"value": range(8)}, index=mi)

        with pytest.raises(KeyError, match="2000"):
            df.loc[("2000",)]
        with pytest.raises(KeyError, match="a"):
            df.loc[("2000", "a")]
        with pytest.raises(KeyError, match="GIBBERISH"):
            df.loc[("2000", "GIBBERISH")]

        # Valid lookups should still work
        result = df.loc[(2000, "a")]
        assert len(result) == 2
        result = df.loc[(2000,)]
        assert len(result) == 4


class TestWhere:
    def test_where(self):
        i = pd.MultiIndex.from_tuples([("A", 1), ("A", 2)])

        msg = r"\.where is not supported for MultiIndex operations"
        with pytest.raises(NotImplementedError, match=msg):
            i.where(True)

    def test_where_array_like(self, listlike_box):
        mi = pd.MultiIndex.from_tuples([("A", 1), ("A", 2)])
        cond = [False, True]
        msg = r"\.where is not supported for MultiIndex operations"
        with pytest.raises(NotImplementedError, match=msg):
            mi.where(listlike_box(cond))


class TestContains:
    def test_contains_top_level(self):
        midx = pd.MultiIndex.from_product([["A", "B"], [1, 2]])
        assert "A" in midx
        assert "A" not in midx._engine

    def test_contains_with_nat(self):
        # MI with a NaT
        mi = pd.MultiIndex(
            levels=[["C"], pd.date_range("2012-01-01", periods=5)],
            codes=[[0, 0, 0, 0, 0, 0], [-1, 0, 1, 2, 3, 4]],
            names=[None, "B"],
        )
        assert ("C", pd.Timestamp("2012-01-01")) in mi
        for val in mi.values:
            assert val in mi

    def test_contains(self, idx):
        assert ("foo", "two") in idx
        assert ("bar", "two") not in idx
        assert None not in idx

    def test_contains_with_missing_value(self):
        # GH#19132
        idx = pd.MultiIndex.from_arrays([[1, np.nan, 2]])
        assert np.nan in idx

        idx = pd.MultiIndex.from_arrays([[1, 2], [np.nan, 3]])
        assert np.nan not in idx
        assert (1, np.nan) in idx

    def test_multiindex_contains_dropped(self):
        # GH#19027
        # test that dropped MultiIndex levels are not in the MultiIndex
        # despite continuing to be in the MultiIndex's levels
        idx = pd.MultiIndex.from_product([[1, 2], [3, 4]])
        assert 2 in idx
        idx = idx.drop(2)

        # drop implementation keeps 2 in the levels
        assert 2 in idx.levels[0]
        # but it should no longer be in the index itself
        assert 2 not in idx

        # also applies to strings
        idx = pd.MultiIndex.from_product([["a", "b"], ["c", "d"]])
        assert "a" in idx
        idx = idx.drop("a")
        assert "a" in idx.levels[0]
        assert "a" not in idx

    def test_contains_td64_level(self):
        # GH#24570
        tx = pd.timedelta_range("09:30:00", "16:00:00", freq="30 min")
        idx = pd.MultiIndex.from_arrays([tx, np.arange(len(tx))])
        assert tx[0] in idx
        assert "element_not_exit" not in idx
        assert "0 day 09:30:00" in idx

    def test_large_mi_contains(self, monkeypatch):
        # GH#10645
        with monkeypatch.context():
            monkeypatch.setattr(libindex, "_SIZE_CUTOFF", 10)
            result = pd.MultiIndex.from_arrays([range(10), range(10)])
            assert (10, 0) not in result


def test_timestamp_multiindex_indexer():
    # https://github.com/pandas-dev/pandas/issues/26944
    idx = pd.MultiIndex.from_product(
        [
            pd.date_range("2019-01-01T00:15:33", periods=100, freq="h", name="date"),
            ["x"],
            [3],
        ]
    )
    df = pd.DataFrame({"foo": np.arange(len(idx))}, idx)
    result = df.loc[pd.IndexSlice["2019-1-2":, "x", :], "foo"]
    qidx = pd.MultiIndex.from_product(
        [
            pd.date_range(
                start="2019-01-02T00:15:33",
                end="2019-01-05T03:15:33",
                freq="h",
                name="date",
            ),
            ["x"],
            [3],
        ]
    )
    should_be = pd.Series(data=np.arange(24, len(qidx) + 24), index=qidx, name="foo")
    tm.assert_series_equal(result, should_be)


@pytest.mark.parametrize(
    "index_arr,expected,target,algo",
    [
        ([[np.nan, "a", "b"], ["c", "d", "e"]], 0, np.nan, "left"),
        ([[np.nan, "a", "b"], ["c", "d", "e"]], 1, (np.nan, "c"), "right"),
        ([["a", "b", "c"], ["d", np.nan, "d"]], 1, ("b", np.nan), "left"),
    ],
)
def test_get_slice_bound_with_missing_value(index_arr, expected, target, algo):
    # issue 19132
    idx = pd.MultiIndex.from_arrays(index_arr)
    result = idx.get_slice_bound(target, side=algo)
    assert result == expected


@pytest.mark.parametrize(
    "index_arr,expected,start_idx,end_idx",
    [
        ([[np.nan, 1, 2], [3, 4, 5]], slice(0, 2, None), np.nan, 1),
        ([[np.nan, 1, 2], [3, 4, 5]], slice(0, 3, None), np.nan, (2, 5)),
        ([[1, 2, 3], [4, np.nan, 5]], slice(1, 3, None), (2, np.nan), 3),
        ([[1, 2, 3], [4, np.nan, 5]], slice(1, 3, None), (2, np.nan), (3, 5)),
    ],
)
def test_slice_indexer_with_missing_value(index_arr, expected, start_idx, end_idx):
    # issue 19132
    idx = pd.MultiIndex.from_arrays(index_arr)
    result = idx.slice_indexer(start=start_idx, end=end_idx)
    assert result == expected


@pytest.mark.parametrize(
    "N, expected_dtype",
    [
        (1, "uint8"),  # 2*4*N = 8
        (2, "uint16"),  # 2*4*N = 16
        (4, "uint32"),  # 2*4*N = 32
        (8, "uint64"),  # 2*4*N = 64
        (10, "object"),  # 2*4*N = 80
    ],
)
def test_pyint_engine(N, expected_dtype):
    # GH#18519 : when combinations of codes cannot be represented in 64
    # bits, the index underlying the MultiIndex engine works with Python
    # integers, rather than uint64.
    keys = [
        tuple(arr)
        for arr in [
            [0] * 4 * N,
            [1] * 4 * N,
            [np.nan] * N + [0] * 3 * N,
            [0] * N + [1] * 3 * N,
            [np.nan] * N + [1] * 2 * N + [0] * N,
        ]
    ]
    # Each level contains 3 elements (NaN, 0, 1), and it's represented
    # in 2 bits to store 4 possible values (0=notfound, 1=NaN, 2=0, 3=1), for
    # a total of 2*N*4 = 80 > 64 bits where N=10 and the number of levels is N*4.
    # If we were using a 64 bit engine and truncating the first levels, the
    # fourth and fifth keys would collide; if truncating the last levels, the
    # fifth and sixth; if rotating bits rather than shifting, the third and fifth.

    index = pd.MultiIndex.from_tuples(keys)
    assert index._engine.values.dtype == expected_dtype

    for idx, key_value in enumerate(keys):
        assert index.get_loc(key_value) == idx

        expected = np.arange(idx + 1, dtype=np.intp)
        result = index.get_indexer([keys[i] for i in expected])
        tm.assert_numpy_array_equal(result, expected)

    # With missing key:
    idces = range(len(keys))
    expected = np.array([-1, *list(idces)], dtype=np.intp)
    missing = tuple([0, 1, 0, 1] * N)
    result = index.get_indexer([missing] + [keys[i] for i in idces])
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize(
    "keys,expected",
    [
        ((slice(None), [5, 4]), [1, 0]),
        ((slice(None), [4, 5]), [0, 1]),
        (([True, False, True], [4, 6]), [0, 2]),
        (([True, False, True], [6, 4]), [0, 2]),
        ((2, [4, 5]), [0, 1]),
        ((2, [5, 4]), [1, 0]),
        (([2], [4, 5]), [0, 1]),
        (([2], [5, 4]), [1, 0]),
    ],
)
def test_get_locs_reordering(keys, expected):
    # GH48384
    idx = pd.MultiIndex.from_arrays(
        [
            [2, 2, 1],
            [4, 5, 6],
        ]
    )
    result = idx.get_locs(keys)
    expected = np.array(expected, dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_basic():
    # GH#55786 - vectorized path for list-like keys
    idx = pd.MultiIndex.from_product([["a", "b"], [1, 2, 3]])
    result = idx.get_locs((["a"], [1, 3]))
    expected = np.array([0, 2], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_with_nan():
    # GH#55786 - NaN labels are stored as code -1
    idx = pd.MultiIndex.from_arrays([["a", "a", "b"], [1.0, np.nan, 2.0]])
    result = idx.get_locs((["a"], [np.nan]))
    expected = np.array([1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_with_nan_and_valid():
    # GH#55786 - mix of NaN and valid labels
    idx = pd.MultiIndex.from_arrays([["a", "a", "a"], [1.0, np.nan, 2.0]])
    result = idx.get_locs((slice(None), [1.0, np.nan]))
    expected = np.array([0, 1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_missing_raises():
    # GH#55786 - missing labels raise KeyError
    idx = pd.MultiIndex.from_product([["a", "b"], [1, 2]])
    with pytest.raises(KeyError, match="99"):
        idx.get_locs((["a"], [1, 99]))


def test_get_locs_list_like_nan_not_in_level():
    # GH#55786 - NaN in query but not in index raises KeyError
    idx = pd.MultiIndex.from_product([["a"], [1, 2]])
    with pytest.raises(KeyError, match="nan"):
        idx.get_locs((["a"], [np.nan]))


@pytest.mark.parametrize(
    "level",
    [
        pd.date_range("2020-01-01", "2020-03-31", freq="D"),
        pd.period_range("2020-01-01", "2020-03-31", freq="D"),
    ],
)
@pytest.mark.parametrize("keys", [["2020-02"], ["2020-02", "2020-03"]])
def test_get_locs_list_like_partial_string(level, keys):
    # GH#64807 - a list-like of partial-string keys on a level that supports
    #  partial indexing (Datetime/Period) must match the same rows as the
    #  scalar key, not just the exact-coerced timestamp.
    idx = pd.MultiIndex.from_product([level, ["a", "b"]])

    result = idx.get_locs([keys])
    expected = np.concatenate([idx.get_locs([key]) for key in keys])
    expected.sort()
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_partial_string_missing_raises():
    # GH#64807 - a partial-string key absent from the level still raises KeyError
    level = pd.date_range("2020-01-01", "2020-03-31", freq="D")
    idx = pd.MultiIndex.from_product([level, ["a", "b"]])
    with pytest.raises(KeyError, match="2021-01"):
        idx.get_locs([["2021-01"]])


@pytest.mark.parametrize(
    "make_key",
    [
        lambda: range(1, 3),
        lambda: (val for val in [1, 2]),
        lambda: {1, 2},
        lambda: {1: None, 2: None},
    ],
)
def test_get_locs_list_like_not_a_list(make_key):
    # GH#64807 - a level key that is list-like but neither a list nor an array
    #  (range/generator/set/dict) must select the same rows as the equivalent
    #  list. Ordering is a separate matter: reordering an unordered container
    #  raises here, both before and after GH#64807.
    idx = pd.MultiIndex.from_product([["a", "b"], [1, 2, 3]])
    result = idx.get_locs((["a"], make_key()))
    expected = np.array([0, 1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_not_a_list_with_na():
    # GH#64807 - the NA inside a non-list list-like has to survive into the code
    #  bookkeeping; isna() on the container itself is a scalar, not elementwise
    idx = pd.MultiIndex.from_arrays([["a", "a", "a"], [1.0, np.nan, 3.0]])
    result = idx.get_locs((slice(None), {1.0: 0, np.nan: 0}))
    tm.assert_numpy_array_equal(result, np.array([0, 1], dtype=np.intp))


@pytest.mark.parametrize(
    "labels", [[frozenset({1}), frozenset({2})], [range(2), range(3)]]
)
def test_get_locs_hashable_list_like_label(labels):
    # GH#64807 - a hashable list-like may be a level label rather than a
    #  sequence of labels, so it must not be materialized before the lookup
    idx = pd.MultiIndex.from_arrays([labels, ["x", "y"]])
    result = idx.get_locs((labels[0], slice(None)))
    tm.assert_numpy_array_equal(result, np.array([0], dtype=np.intp))


def test_get_locs_generator_reordering():
    # GH#64807 - the generator must survive long enough to reorder the result
    idx = pd.MultiIndex.from_product([["a", "b"], [1, 2, 3]])
    result = idx.get_locs((["a"], (val for val in [3, 1])))
    tm.assert_numpy_array_equal(result, np.array([2, 0], dtype=np.intp))


def test_get_locs_list_like_loc():
    # GH#64807 - .loc with a non-list list-like level key
    idx = pd.MultiIndex.from_product([["a", "b"], [0, 1, 2]])
    df = pd.DataFrame({"x": range(6)}, index=idx)
    result = df.loc[(["a"], range(2)), :]
    expected = df.iloc[[0, 1]]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("keys", [([], [1, 2]), (["a"], []), ([], [])])
def test_get_locs_list_like_empty(keys):
    # GH#64807 - an empty level key matches nothing rather than raising
    idx = pd.MultiIndex.from_product([["a", "b"], [1, 2, 3]])
    result = idx.get_locs(keys)
    expected = np.array([], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_locs_list_like_empty_typed():
    # GH#64807 - an empty key that is not object dtype reaches the vectorized
    #  branch, where it has to short-circuit the levels after it rather than let
    #  them raise. An empty *list* is object dtype and takes the strict path.
    idx = pd.MultiIndex.from_product([[1, 2], [1, 2, 3]])
    result = idx.get_locs((np.array([], dtype=np.int64), [1]))
    tm.assert_numpy_array_equal(result, np.array([], dtype=np.intp))


@pytest.mark.skipif(not PY312, reason="a slice is unhashable before 3.12")
@pytest.mark.parametrize("dtype", [None, object])
def test_get_locs_list_like_slice_label(dtype):
    # GH#64807 - a slice is hashable from Python 3.12 on, so it arrives as a
    #  level label; resolving it must not depend on the level's dtype
    idx = pd.MultiIndex.from_product([["a"], pd.Index([1, 2, 3], dtype=dtype)])
    result = idx.get_locs((["a"], [slice(1, 2)]))
    tm.assert_numpy_array_equal(result, np.array([0, 1], dtype=np.intp))


@pytest.mark.skipif(not PY312, reason="a slice is unhashable before 3.12")
def test_get_locs_list_like_slice_label_mixed_level():
    # GH#64807 - a level that infers as "mixed", just as the slice key does,
    #  gets past the dtype gate, so the slice needs recognizing on its own
    level = pd.Index([(1,), (2,), (3,)], dtype=object)
    assert level.inferred_type == "mixed"
    idx = pd.MultiIndex.from_product([["a"], level])
    result = idx.get_locs((["a"], [slice((1,), (2,))]))
    tm.assert_numpy_array_equal(result, np.array([0, 1], dtype=np.intp))


@pytest.mark.parametrize("dtype", [None, object])
def test_get_locs_list_like_unhashable_after_missing(dtype):
    # GH#64807 - an unhashable entry is not a label, but it must not preempt an
    #  earlier label that is merely missing. The tuple level infers as "mixed",
    #  like the key, so nothing else routes it off the vectorized path.
    labels = [(1,), (2,), (3,)] if dtype is object else [1, 2, 3]
    missing, unhashable = ((9,), [1, 2]) if dtype is object else ("missing", [1, 2])
    idx = pd.MultiIndex.from_product([["a"], pd.Index(labels, dtype=dtype)])

    with pytest.raises(KeyError, match=re.escape(str(missing))):
        idx.get_locs((["a"], [missing, unhashable]))
    with pytest.raises(InvalidIndexError, match=re.escape(str([unhashable, missing]))):
        idx.get_locs((["a"], [unhashable, missing]))


def test_get_locs_list_like_unsupported_timedelta_unit():
    # GH#64807 - an Index cannot hold a month-unit timedelta, so building one
    #  from the key must not turn a lookup into ValueError
    level = pd.Index(
        np.array([np.timedelta64(1, "M"), np.timedelta64(2, "M")], dtype=object),
        dtype=object,
    )
    idx = pd.MultiIndex.from_arrays([["x", "y"], level])
    result = idx.get_locs((slice(None), [np.timedelta64(1, "M")]))
    tm.assert_numpy_array_equal(result, np.array([0], dtype=np.intp))

    idx = pd.MultiIndex.from_product([["x"], pd.date_range("2020", periods=3)])
    with pytest.raises(KeyError, match="1"):
        idx.get_locs((slice(None), [np.timedelta64(1, "M")]))


def test_get_locs_list_like_interval_level():
    # GH#64807 - an IntervalIndex level matches a numeric key by containment,
    #  the same as a flat IntervalIndex does for a list-like key
    level = pd.IntervalIndex.from_breaks([0.0, 1.0, 2.0, 3.0])
    idx = pd.MultiIndex.from_product([["a"], level])
    result = idx.get_locs((["a"], [0.5, 2.5]))
    tm.assert_numpy_array_equal(result, np.array([0, 2], dtype=np.intp))


def test_get_locs_list_like_categorical_overlapping_interval_level():
    # GH#64807 - a CategoricalIndex reports itself unique whatever its
    #  categories do, so overlapping interval categories still need the
    #  per-label path
    categories = pd.IntervalIndex.from_tuples([(0, 2), (1, 3)])
    level = pd.CategoricalIndex(
        pd.Categorical.from_codes([0, 1], dtype=pd.CategoricalDtype(categories))
    )
    assert level._index_as_unique and not level.categories._index_as_unique
    idx = pd.MultiIndex.from_arrays([["a", "b"], level])
    # the lookups the per-label path makes work
    tm.assert_numpy_array_equal(
        idx.get_locs((slice(None), [pd.Interval(0, 2)])), np.array([0], dtype=np.intp)
    )
    tm.assert_numpy_array_equal(
        idx.get_locs((slice(None), [0.5])), np.array([0], dtype=np.intp)
    )


def test_get_locs_list_like_float16_key():
    # GH#64807 - an Index cannot hold float16 at all, so building one from the
    #  key must not turn a lookup into NotImplementedError
    idx = pd.MultiIndex.from_product([["a"], [1.0, 2.0, 3.0]])
    result = idx.get_locs((["a"], [np.float16(1.0)]))
    tm.assert_numpy_array_equal(result, np.array([0], dtype=np.intp))

    idx = pd.MultiIndex.from_product([["a"], pd.date_range("2020", periods=3)])
    with pytest.raises(KeyError, match="1.0"):
        idx.get_locs((["a"], [np.float16(1.0)]))


def test_get_locs_list_like_nan_key_without_nan_rows():
    # GH#64807 - a NaN key must raise when the data carries no NA rows, even
    #  though NaN is looked up without needing a level entry
    idx = pd.MultiIndex.from_arrays([["a", "a"], [1.0, 2.0]])
    with pytest.raises(KeyError, match="nan"):
        idx.get_locs((slice(None), [np.nan]))


def test_get_locs_list_like_categorical_level():
    # GH#64807 - smoke test for a CategoricalIndex level, which no other case
    #  in this file covers
    idx = pd.MultiIndex.from_product([pd.CategoricalIndex(list("abc")), [1, 2]])
    result = idx.get_locs((["a", "c"], slice(None)))
    tm.assert_numpy_array_equal(result, np.array([0, 1, 4, 5], dtype=np.intp))


@pytest.mark.parametrize(
    "keys, match",
    [
        ((["a", "b"], slice(None)), "['a', 'b']"),
        ((["b"], slice(None)), "['b']"),
    ],
)
def test_get_locs_list_like_unused_level_raises(keys, match):
    # GH#64807 - a label present in the level but unused by any code is missing,
    #  just as it is when passed as a scalar
    idx = pd.MultiIndex.from_product([["a", "b"], [1, 2, 3]])[:3]
    assert "b" in idx.levels[0]
    with pytest.raises(KeyError, match=re.escape(match)):
        idx.get_locs(keys)


def test_get_locs_list_like_nan_code_vs_last_level_label():
    # GH#64807 - the -1 code standing for NaN must not wrap around onto the
    #  last level entry, which would both leak the NaN row in...
    idx = pd.MultiIndex.from_arrays([["a"] * 3, [1.0, np.nan, 2.0]])
    result = idx.get_locs((["a"], [2.0]))
    tm.assert_numpy_array_equal(result, np.array([2], dtype=np.intp))

    # ...and make an unused last label look used
    idx = pd.MultiIndex.from_arrays([["a"] * 4, [1.0, np.nan, 1.0, 2.0]])[:3]
    assert 2.0 in idx.levels[1]
    with pytest.raises(KeyError, match=re.escape("[2.0]")):
        idx.get_locs((["a"], [2.0]))


def test_get_locs_list_like_na_label_against_na_in_level():
    # GH#64807 - a level can carry an NA entry of its own that no code uses,
    #  while the NA rows carry code -1 (e.g. concatenating a
    #  groupby(dropna=False) result with an ordinary NaN-indexed frame, then
    #  dropping rows). get_indexer matches an NA label to that entry, so it must
    #  not be mistaken for a label the level does not have.
    idx = pd.MultiIndex(
        levels=[pd.Index(["a"]), pd.Index([1.0, np.nan, 2.0])], codes=[[0, 0], [-1, 0]]
    )
    expected = idx.get_locs((slice(None), np.nan))
    result = idx.get_locs((slice(None), [np.nan]))
    tm.assert_numpy_array_equal(result, expected)

    # and a genuinely missing label alongside it must still raise
    with pytest.raises(KeyError, match=re.escape("[nan, 9.9]")):
        idx.get_locs((slice(None), [np.nan, 9.9]))


def test_get_locs_list_like_na_label_not_matched_to_level_entry():
    # GH#64807 - get_indexer can match an NA label to an NA entry of the level,
    #  but the scalar path always resolves NA to code -1, so the list path has
    #  to select the -1 rows only and agree with it
    idx = pd.MultiIndex.from_arrays([["a"] * 3, [1.0, np.nan, 3.0]])
    idx = idx.set_levels(pd.Index([1.0, np.nan]), level=1, verify_integrity=False)
    assert idx.levels[1].hasnans

    expected = np.array([1], dtype=np.intp)
    tm.assert_numpy_array_equal(idx.get_locs((slice(None), [np.nan])), expected)
    tm.assert_numpy_array_equal(idx.get_locs((slice(None), np.nan)), expected)


def test_get_locs_list_like_multiindex_key():
    # GH#64807 - tuple-valued level labels, so the key itself is a MultiIndex.
    #  isna() answers per tuple *element*, which for a single label broadcasts
    #  against the codes silently and for several of them does not broadcast at
    #  all, so both lengths are worth covering.
    idx = pd.MultiIndex.from_arrays([[("x", 1), ("y", 2), ("z", 3)], list("abc")])
    result = idx.get_locs((pd.Index([("x", 1), ("y", 2), ("z", 3)]), slice(None)))
    tm.assert_numpy_array_equal(result, np.array([0, 1, 2], dtype=np.intp))

    result = idx.get_locs((pd.Index([("x", 1)]), slice(None)))
    tm.assert_numpy_array_equal(result, np.array([0], dtype=np.intp))


def test_get_locs_list_like_mixed_type_key():
    # GH#64807 - a list-like key resolves through get_indexer, which cannot
    #  convert an object-dtype target mixing strings and Timedeltas, so it
    #  reports the string as missing. A flat Index rejects the same key, even
    #  though the scalar path parses "1 days".
    idx = pd.MultiIndex.from_product([["a"], pd.timedelta_range("1 day", periods=3)])
    key = ["1 days", pd.Timedelta("2 days")]
    with pytest.raises(KeyError, match="1 days"):
        idx.get_locs((["a"], key))
    with pytest.raises(KeyError, match="1 days"):
        pd.Series(range(3), index=idx.levels[1]).loc[key]

    tm.assert_numpy_array_equal(
        idx.get_locs((["a"], "1 days")), np.array([0], dtype=np.intp)
    )


def test_get_locs_list_like_float32_level():
    # GH#64807 - a float64 key never equals a float32 label once get_indexer
    #  widens the level to hold it, matching what a flat Index does. The scalar
    #  path is not held to that; the two are not expected to agree here.
    idx = pd.MultiIndex.from_arrays(
        [["a", "a", "b"], np.array([1.1, 2.2, 1.1], dtype=np.float32)]
    )
    with pytest.raises(KeyError, match="1.1"):
        idx.get_locs((slice(None), [1.1]))
    with pytest.raises(KeyError, match="1.1"):
        pd.Series(range(2), index=idx.levels[1]).loc[[1.1]]

    tm.assert_numpy_array_equal(
        idx.get_locs((slice(None), 1.1)), np.array([0, 2], dtype=np.intp)
    )


def test_get_locs_list_like_na_in_object_container():
    # GH#64807 - pd.NA alongside a real label leaves the key container object
    #  dtype; the NA still has to be recognized and matched to the -1 codes
    idx = pd.MultiIndex(
        levels=[pd.Index(["a", "b"]), pd.Index([1.1, 2.2])],
        codes=[[0, 0, 1], [0, -1, 0]],
    )
    # the key order puts the 1.1 rows ahead of the NA row
    result = idx.get_locs((slice(None), [1.1, pd.NA]))
    tm.assert_numpy_array_equal(result, np.array([0, 2, 1], dtype=np.intp))


def test_get_locs_list_like_signed_unsigned_level():
    # GH#64807 - get_indexer reconciles a signed level and an unsigned key
    #  through float64, colliding two labels above 2**53 and selecting the
    #  wrong row; a flat Index resolves the same key correctly
    idx = pd.MultiIndex.from_arrays([["a", "b"], [2**53, 2**53 + 1]])
    key = pd.array([2**53], dtype="UInt64")
    tm.assert_numpy_array_equal(
        idx.get_locs((slice(None), key)), np.array([0], dtype=np.intp)
    )
    assert pd.Series(range(2), index=idx.levels[1]).loc[key].tolist() == [0]


def test_get_locs_list_like_categorical_key():
    # GH#64807 - a Categorical key against an object level of tuple labels:
    #  a label the level lacks still raises rather than matching
    level = pd.Index([(1,), (2,), (3,)], dtype=object)
    idx = pd.MultiIndex(
        levels=[pd.Index(list("abc")), level],
        codes=[[0, 1, 2], [0, 1, 2]],
        verify_integrity=False,
    )
    with pytest.raises(KeyError, match="9"):
        idx.get_locs((slice(None), pd.Categorical([(1,), 9])))


def test_get_locs_list_like_narrow_int_level():
    # GH#64807 - a wider integer holds every value of a narrower one, so an
    #  int64 key still matches an int32 level exactly (no float widening)
    idx = pd.MultiIndex.from_product([["a"], pd.Index(np.arange(5, dtype="int32"))])
    result = idx.get_locs((["a"], [1, 3]))
    tm.assert_numpy_array_equal(result, np.array([1, 3], dtype=np.intp))


def test_get_locs_list_like_ragged_tuple_labels():
    # GH#64807 - a tuple key of tuple-valued labels must not be tupleized into
    #  a single label by Index()
    idx = pd.MultiIndex.from_arrays([[("x", 1), ("y", 2, 3)], ["a", "b"]])
    expected = np.array([0, 1], dtype=np.intp)
    tm.assert_numpy_array_equal(
        idx.get_locs(((("x", 1), ("y", 2, 3)), slice(None))), expected
    )
    tm.assert_numpy_array_equal(
        idx.get_locs(([("x", 1), ("y", 2, 3)], slice(None))), expected
    )


def test_get_locs_list_like_wide_level():
    # GH#64807 - slicing a MultiIndex leaves its levels untrimmed, so the level
    #  can be bigger than the codes; unused labels must still raise there
    idx = pd.MultiIndex.from_product([list("abcde"), [1]])[:2]
    assert len(idx.levels[0]) > len(idx.codes[0])

    result = idx.get_locs((["a", "b"], slice(None)))
    tm.assert_numpy_array_equal(result, np.array([0, 1], dtype=np.intp))

    with pytest.raises(KeyError, match=re.escape("['a', 'c']")):
        idx.get_locs((["a", "c"], slice(None)))


@pytest.mark.parametrize(
    "keys, expected", [([1.5], [0, 1]), ([pd.Interval(0, 2)], [0])]
)
def test_get_locs_list_like_overlapping_interval_level(keys, expected):
    # GH#64807 - an overlapping IntervalIndex level cannot use get_indexer; a
    #  label there can match several level entries, e.g. 1.5 is in both
    level = pd.IntervalIndex.from_tuples([(0, 2), (1, 3)])
    idx = pd.MultiIndex.from_product([["a", "b"], level])
    result = idx.get_locs((["a"], keys))
    tm.assert_numpy_array_equal(result, np.array(expected, dtype=np.intp))


def test_get_indexer_for_multiindex_with_nans(nulls_fixture):
    # GH37222
    idx1 = pd.MultiIndex.from_product([["A"], [1.0, 2.0]], names=["id1", "id2"])
    idx2 = pd.MultiIndex.from_product(
        [["A"], [nulls_fixture, 2.0]], names=["id1", "id2"]
    )

    result = idx2.get_indexer(idx1)
    expected = np.array([-1, 1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)

    result = idx1.get_indexer(idx2)
    expected = np.array([-1, 1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)


def test_get_loc_namedtuple_behaves_like_tuple():
    # GH57922
    NamedIndex = namedtuple("NamedIndex", ("a", "b"))
    multi_idx = pd.MultiIndex.from_tuples(
        [NamedIndex("i1", "i2"), NamedIndex("i3", "i4"), NamedIndex("i5", "i6")]
    )
    for idx in (multi_idx, multi_idx.to_flat_index()):
        assert idx.get_loc(NamedIndex("i1", "i2")) == 0
        assert idx.get_loc(NamedIndex("i3", "i4")) == 1
        assert idx.get_loc(NamedIndex("i5", "i6")) == 2
        assert idx.get_loc(("i1", "i2")) == 0
        assert idx.get_loc(("i3", "i4")) == 1
        assert idx.get_loc(("i5", "i6")) == 2
    multi_idx = pd.MultiIndex.from_tuples([("i1", "i2"), ("i3", "i4"), ("i5", "i6")])
    for idx in (multi_idx, multi_idx.to_flat_index()):
        assert idx.get_loc(NamedIndex("i1", "i2")) == 0
        assert idx.get_loc(NamedIndex("i3", "i4")) == 1
        assert idx.get_loc(NamedIndex("i5", "i6")) == 2
        assert idx.get_loc(("i1", "i2")) == 0
        assert idx.get_loc(("i3", "i4")) == 1
        assert idx.get_loc(("i5", "i6")) == 2
