from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
from pandas import (
    Categorical,
    CategoricalIndex,
    Index,
    IntervalIndex,
    MultiIndex,
    date_range,
)
import pandas._testing as tm
from pandas.core.indexes.base import InvalidIndexError


def test_slice_locs_partial(idx):
    sorted_idx, _ = idx.sortlevel(0)

    result = sorted_idx.slice_locs(("foo", "two"), ("qux", "one"))
    assert result == (1, 5)

    result = sorted_idx.slice_locs(None, ("qux", "one"))
    assert result == (0, 5)

    result = sorted_idx.slice_locs(("foo", "two"), None)
    assert result == (1, len(sorted_idx))

    result = sorted_idx.slice_locs("bar", "baz")
    assert result == (2, 4)


def test_slice_locs():
    df = tm.makeTimeDataFrame()
    stacked = df.stack()
    idx = stacked.index

    slob = slice(*idx.slice_locs(df.index[5], df.index[15]))
    sliced = stacked[slob]
    expected = df[5:16].stack()
    tm.assert_almost_equal(sliced.values, expected.values)

    slob = slice(
        *idx.slice_locs(
            df.index[5] + timedelta(seconds=30), df.index[15] - timedelta(seconds=30)
        )
    )
    sliced = stacked[slob]
    expected = df[6:15].stack()
    tm.assert_almost_equal(sliced.values, expected.values)


def test_slice_locs_with_type_mismatch():
    df = tm.makeTimeDataFrame()
    stacked = df.stack()
    idx = stacked.index
    with pytest.raises(TypeError, match="^Level type mismatch"):
        idx.slice_locs((1, 3))
    with pytest.raises(TypeError, match="^Level type mismatch"):
        idx.slice_locs(df.index[5] + timedelta(seconds=30), (5, 2))
    df = tm.makeCustomDataframe(5, 5)
    stacked = df.stack()
    idx = stacked.index
    with pytest.raises(TypeError, match="^Level type mismatch"):
        idx.slice_locs(timedelta(seconds=30))
    # TODO: Try creating a UnicodeDecodeError in exception message
    with pytest.raises(TypeError, match="^Level type mismatch"):
        idx.slice_locs(df.index[1], (16, "a"))


def test_slice_locs_not_sorted():
    index = MultiIndex(
        levels=[Index(np.arange(4)), Index(np.arange(4)), Index(np.arange(4))],
        codes=[
            np.array([0, 0, 1, 2, 2, 2, 3, 3]),
            np.array([0, 1, 0, 0, 0, 1, 0, 1]),
            np.array([1, 0, 1, 1, 0, 0, 1, 0]),
        ],
    )
    msg = "[Kk]ey length.*greater than MultiIndex lexsort depth"
    with pytest.raises(KeyError, match=msg):
        index.slice_locs((1, 0, 1), (2, 1, 0))

    # works
    sorted_index, _ = index.sortlevel(0)
    # should there be a test case here???
    sorted_index.slice_locs((1, 0, 1), (2, 1, 0))


def test_slice_locs_not_contained():
    # some searchsorted action

    index = MultiIndex(
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


def test_putmask_with_wrong_mask(idx):
    # GH18368

    msg = "putmask: mask and data must be the same size"
    with pytest.raises(ValueError, match=msg):
        idx.putmask(np.ones(len(idx) + 1, np.bool), 1)

    with pytest.raises(ValueError, match=msg):
        idx.putmask(np.ones(len(idx) - 1, np.bool), 1)

    with pytest.raises(ValueError, match=msg):
        idx.putmask("foo", 1)


def test_get_indexer():
    major_axis = Index(np.arange(4))
    minor_axis = Index(np.arange(2))

    major_codes = np.array([0, 0, 1, 2, 2, 3, 3], dtype=np.intp)
    minor_codes = np.array([0, 1, 0, 0, 1, 0, 1], dtype=np.intp)

    index = MultiIndex(
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
    idx1 = Index(list(range(10)) + list(range(10)))
    idx2 = Index(list(range(20)))

    msg = "Reindexing only valid with uniquely valued Index objects"
    with pytest.raises(InvalidIndexError, match=msg):
        idx1.get_indexer(idx2)


def test_get_indexer_nearest():
    midx = MultiIndex.from_tuples([("a", 1), ("b", 2)])
    msg = "method='nearest' not implemented yet for MultiIndex; see GitHub issue 9365"
    with pytest.raises(NotImplementedError, match=msg):
        midx.get_indexer(["a"], method="nearest")
    msg = "tolerance not implemented yet for MultiIndex"
    with pytest.raises(NotImplementedError, match=msg):
        midx.get_indexer(["a"], method="pad", tolerance=2)


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


def test_get_indexer_consistency(idx):
    # See GH 16819
    if isinstance(idx, IntervalIndex):
        pass

    if idx.is_unique or isinstance(idx, CategoricalIndex):
        indexer = idx.get_indexer(idx[0:2])
        assert isinstance(indexer, np.ndarray)
        assert indexer.dtype == np.intp
    else:
        e = "Reindexing only valid with uniquely valued Index objects"
        with pytest.raises(InvalidIndexError, match=e):
            idx.get_indexer(idx[0:2])

    indexer, _ = idx.get_indexer_non_unique(idx[0:2])
    assert isinstance(indexer, np.ndarray)
    assert indexer.dtype == np.intp


@pytest.mark.parametrize("ind1", [[True] * 5, pd.Index([True] * 5)])
@pytest.mark.parametrize(
    "ind2",
    [[True, False, True, False, False], pd.Index([True, False, True, False, False])],
)
def test_getitem_bool_index_all(ind1, ind2):
    # GH#22533
    idx = MultiIndex.from_tuples([(10, 1), (20, 2), (30, 3), (40, 4), (50, 5)])
    tm.assert_index_equal(idx[ind1], idx)

    expected = MultiIndex.from_tuples([(10, 1), (30, 3)])
    tm.assert_index_equal(idx[ind2], expected)


@pytest.mark.parametrize("ind1", [[True], pd.Index([True])])
@pytest.mark.parametrize("ind2", [[False], pd.Index([False])])
def test_getitem_bool_index_single(ind1, ind2):
    # GH#22533
    idx = MultiIndex.from_tuples([(10, 1)])
    tm.assert_index_equal(idx[ind1], idx)

    expected = pd.MultiIndex(
        levels=[np.array([], dtype=np.int64), np.array([], dtype=np.int64)],
        codes=[[], []],
    )
    tm.assert_index_equal(idx[ind2], expected)


def test_get_loc(idx):
    assert idx.get_loc(("foo", "two")) == 1
    assert idx.get_loc(("baz", "two")) == 3
    with pytest.raises(KeyError, match=r"^10$"):
        idx.get_loc(("bar", "two"))
    with pytest.raises(KeyError, match=r"^'quux'$"):
        idx.get_loc("quux")

    msg = "only the default get_loc method is currently supported for MultiIndex"
    with pytest.raises(NotImplementedError, match=msg):
        idx.get_loc("foo", method="nearest")

    # 3 levels
    index = MultiIndex(
        levels=[Index(np.arange(4)), Index(np.arange(4)), Index(np.arange(4))],
        codes=[
            np.array([0, 0, 1, 2, 2, 2, 3, 3]),
            np.array([0, 1, 0, 0, 0, 1, 0, 1]),
            np.array([1, 0, 1, 1, 0, 0, 1, 0]),
        ],
    )
    with pytest.raises(KeyError, match=r"^\(1, 1\)$"):
        index.get_loc((1, 1))
    assert index.get_loc((2, 0)) == slice(3, 5)


def test_get_loc_duplicates():
    index = Index([2, 2, 2, 2])
    result = index.get_loc(2)
    expected = slice(0, 4)
    assert result == expected
    # pytest.raises(Exception, index.get_loc, 2)

    index = Index(["c", "a", "a", "b", "b"])
    rs = index.get_loc("c")
    xp = 0
    assert rs == xp


def test_get_loc_level():
    index = MultiIndex(
        levels=[Index(np.arange(4)), Index(np.arange(4)), Index(np.arange(4))],
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

    index = MultiIndex(
        levels=[[2000], list(range(4))],
        codes=[np.array([0, 0, 0, 0]), np.array([0, 1, 2, 3])],
    )
    result, new_index = index.get_loc_level((2000, slice(None, None)))
    expected = slice(None, None)
    assert result == expected
    assert new_index.equals(index.droplevel(0))


@pytest.mark.parametrize("dtype1", [int, float, bool, str])
@pytest.mark.parametrize("dtype2", [int, float, bool, str])
def test_get_loc_multiple_dtypes(dtype1, dtype2):
    # GH 18520
    levels = [np.array([0, 1]).astype(dtype1), np.array([0, 1]).astype(dtype2)]
    idx = pd.MultiIndex.from_product(levels)
    assert idx.get_loc(idx[2]) == 2


@pytest.mark.parametrize("level", [0, 1])
@pytest.mark.parametrize("dtypes", [[int, float], [float, int]])
def test_get_loc_implicit_cast(level, dtypes):
    # GH 18818, GH 15994 : as flat index, cast int to float and vice-versa
    levels = [["a", "b"], ["c", "d"]]
    key = ["b", "d"]
    lev_dtype, key_dtype = dtypes
    levels[level] = np.array([0, 1], dtype=lev_dtype)
    key[level] = key_dtype(1)
    idx = MultiIndex.from_product(levels)
    assert idx.get_loc(tuple(key)) == 3


def test_get_loc_cast_bool():
    # GH 19086 : int is casted to bool, but not vice-versa
    levels = [[False, True], np.arange(2, dtype="int64")]
    idx = MultiIndex.from_product(levels)

    assert idx.get_loc((0, 1)) == 1
    assert idx.get_loc((1, 0)) == 2

    with pytest.raises(KeyError, match=r"^\(False, True\)$"):
        idx.get_loc((False, True))
    with pytest.raises(KeyError, match=r"^\(True, False\)$"):
        idx.get_loc((True, False))


@pytest.mark.parametrize("level", [0, 1])
def test_get_loc_nan(level, nulls_fixture):
    # GH 18485 : NaN in MultiIndex
    levels = [["a", "b"], ["c", "d"]]
    key = ["b", "d"]
    levels[level] = np.array([0, nulls_fixture], dtype=type(nulls_fixture))
    key[level] = nulls_fixture
    idx = MultiIndex.from_product(levels)
    assert idx.get_loc(tuple(key)) == 3


def test_get_loc_missing_nan():
    # GH 8569
    idx = MultiIndex.from_arrays([[1.0, 2.0], [3.0, 4.0]])
    assert isinstance(idx.get_loc(1), slice)
    with pytest.raises(KeyError, match=r"^3\.0$"):
        idx.get_loc(3)
    with pytest.raises(KeyError, match=r"^nan$"):
        idx.get_loc(np.nan)
    with pytest.raises(TypeError, match=r"'\[nan\]' is an invalid key"):
        # listlike/non-hashable raises TypeError
        idx.get_loc([np.nan])


def test_get_indexer_categorical_time():
    # https://github.com/pandas-dev/pandas/issues/21390
    midx = MultiIndex.from_product(
        [
            Categorical(["a", "b", "c"]),
            Categorical(date_range("2012-01-01", periods=3, freq="H")),
        ]
    )
    result = midx.get_indexer(midx)
    tm.assert_numpy_array_equal(result, np.arange(9, dtype=np.intp))


def test_timestamp_multiindex_indexer():
    # https://github.com/pandas-dev/pandas/issues/26944
    idx = pd.MultiIndex.from_product(
        [
            pd.date_range("2019-01-01T00:15:33", periods=100, freq="H", name="date"),
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
                end="2019-01-05T02:15:33",
                freq="H",
                name="date",
            ),
            ["x"],
            [3],
        ]
    )
    should_be = pd.Series(data=np.arange(24, len(qidx) + 24), index=qidx, name="foo")
    tm.assert_series_equal(result, should_be)


def test_get_loc_with_values_including_missing_values():
    # issue 19132
    idx = MultiIndex.from_product([[np.nan, 1]] * 2)
    expected = slice(0, 2, None)
    assert idx.get_loc(np.nan) == expected

    idx = MultiIndex.from_arrays([[np.nan, 1, 2, np.nan]])
    expected = np.array([True, False, False, True])
    tm.assert_numpy_array_equal(idx.get_loc(np.nan), expected)

    idx = MultiIndex.from_product([[np.nan, 1]] * 3)
    expected = slice(2, 4, None)
    assert idx.get_loc((np.nan, 1)) == expected


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
def test_get_indexer_with_missing_value(index_arr, labels, expected):
    # issue 19132
    idx = MultiIndex.from_arrays(index_arr)
    result = idx.get_indexer(labels)
    tm.assert_numpy_array_equal(result, expected)


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
    idx = MultiIndex.from_arrays(index_arr)
    result = idx.get_slice_bound(target, side=algo, kind="loc")
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
    idx = MultiIndex.from_arrays(index_arr)
    result = idx.slice_indexer(start=start_idx, end=end_idx)
    assert result == expected


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
def test_slice_locs_with_missing_value(index_arr, expected, start_idx, end_idx):
    # issue 19132
    idx = MultiIndex.from_arrays(index_arr)
    result = idx.slice_locs(start=start_idx, end=end_idx)
    assert result == expected
