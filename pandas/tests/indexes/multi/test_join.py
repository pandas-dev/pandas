import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("other", [["three", "one", "two"], ["one"], ["one", "three"]])
def test_join_level(idx, other, join_type):
    other = pd.Index(other)
    join_index, lidx, ridx = other.join(
        idx, how=join_type, level="second", return_indexers=True
    )

    exp_level = other.join(idx.levels[1], how=join_type)
    assert join_index.levels[0].equals(idx.levels[0])
    assert join_index.levels[1].equals(exp_level)

    # pare down levels
    mask = np.array([x[1] in exp_level for x in idx], dtype=bool)
    exp_values = idx.values[mask]
    tm.assert_numpy_array_equal(join_index.values, exp_values)

    if join_type in ("outer", "inner"):
        join_index2, ridx2, lidx2 = idx.join(
            other, how=join_type, level="second", return_indexers=True
        )

        assert join_index.equals(join_index2)
        tm.assert_numpy_array_equal(lidx, lidx2)
        if ridx is None:
            assert ridx == ridx2
        else:
            tm.assert_numpy_array_equal(ridx, ridx2)
        tm.assert_numpy_array_equal(join_index2.values, exp_values)


def test_join_level_corner_case(idx):
    # some corner cases
    index = pd.Index(["three", "one", "two"])
    result = index.join(idx, level="second")
    assert isinstance(result, pd.MultiIndex)

    with pytest.raises(TypeError, match="Join.*MultiIndex.*ambiguous"):
        idx.join(idx, level=1)


def test_join_self(idx, join_type):
    result = idx.join(idx, how=join_type)
    expected = idx
    if join_type == "outer":
        expected = expected.sort_values()
    tm.assert_index_equal(result, expected)


def test_join_multi():
    # GH 10665
    midx = pd.MultiIndex.from_product([np.arange(4), np.arange(4)], names=["a", "b"])
    idx = pd.Index([1, 2, 5], name="b")

    # inner
    jidx, lidx, ridx = midx.join(idx, how="inner", return_indexers=True)
    exp_idx = pd.MultiIndex.from_product([np.arange(4), [1, 2]], names=["a", "b"])
    exp_lidx = np.array([1, 2, 5, 6, 9, 10, 13, 14], dtype=np.intp)
    exp_ridx = np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=np.intp)
    tm.assert_index_equal(jidx, exp_idx)
    tm.assert_numpy_array_equal(lidx, exp_lidx)
    tm.assert_numpy_array_equal(ridx, exp_ridx)
    # flip
    jidx, ridx, lidx = idx.join(midx, how="inner", return_indexers=True)
    tm.assert_index_equal(jidx, exp_idx)
    tm.assert_numpy_array_equal(lidx, exp_lidx)
    tm.assert_numpy_array_equal(ridx, exp_ridx)

    # keep MultiIndex
    jidx, lidx, ridx = midx.join(idx, how="left", return_indexers=True)
    exp_ridx = np.array(
        [-1, 0, 1, -1, -1, 0, 1, -1, -1, 0, 1, -1, -1, 0, 1, -1], dtype=np.intp
    )
    tm.assert_index_equal(jidx, midx)
    assert lidx is None
    tm.assert_numpy_array_equal(ridx, exp_ridx)
    # flip
    jidx, ridx, lidx = idx.join(midx, how="right", return_indexers=True)
    tm.assert_index_equal(jidx, midx)
    assert lidx is None
    tm.assert_numpy_array_equal(ridx, exp_ridx)


@pytest.mark.parametrize("flat_first", [True, False])
@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
def test_join_flat_with_multi_unmatched_flat_keys(how, flat_first):
    # GH#34292 a key of the flat side with no MultiIndex row is kept when `how`
    # keeps that side, with missing values in the other levels
    midx = pd.MultiIndex.from_arrays([[1, 3, 5], [2, 4, 6]], names=["a", "b"])
    idx = pd.Index([2, 5], name="a")
    left, right = (idx, midx) if flat_first else (midx, idx)

    result, lidx, ridx = left.join(right, how=how, return_indexers=True)

    keep_flat = how in ("outer", "left" if flat_first else "right")
    keep_multi = how in ("outer", "right" if flat_first else "left")
    multi_rows = [0, 1, 2] if keep_multi else [2]
    flat_rows = [-1, -1, 1] if keep_multi else [1]
    codes = [multi_rows, multi_rows]
    if keep_flat:
        codes = [[*multi_rows, 3], [*multi_rows, -1]]
        multi_rows, flat_rows = [*multi_rows, -1], [*flat_rows, 0]
    expected = pd.MultiIndex(
        levels=[[1, 3, 5, 2], [2, 4, 6]], codes=codes, names=["a", "b"]
    )
    tm.assert_index_equal(result, expected)

    flat_idx, multi_idx = (lidx, ridx) if flat_first else (ridx, lidx)
    if multi_idx is None:
        # None means the result is the MultiIndex itself
        assert result.equals(midx)
        multi_idx = np.arange(len(midx), dtype=np.intp)
    tm.assert_numpy_array_equal(flat_idx, np.array(flat_rows, dtype=np.intp))
    tm.assert_numpy_array_equal(multi_idx, np.array(multi_rows, dtype=np.intp))


def test_join_flat_with_multi_unmatched_missing_and_object_keys():
    # GH#34292 appended keys keep the level's dtype; a missing key stays missing
    midx = pd.MultiIndex.from_arrays(
        [pd.Index(["a", "b"], dtype=object), [1, 2]], names=["k", "o"]
    )
    idx = pd.Index(["b", np.nan, "d"], dtype=object, name="k")

    result = idx.join(midx, how="left")

    expected = pd.MultiIndex(
        levels=[pd.Index(["b", np.nan, "d"], dtype=object), [1, 2]],
        codes=[[0, -1, 2], [1, -1, -1]],
        names=["k", "o"],
    )
    tm.assert_index_equal(result, expected, exact=True)
    tm.assert_index_equal(result.levels[0], expected.levels[0], exact=True)
    tm.assert_numpy_array_equal(result.codes[0], expected.codes[0])


def test_join_flat_with_multi_unmatched_middle_level():
    # GH#34292 shared level not first; the sliced MultiIndex keeps the unused
    # level value 5, so the level join returns no MultiIndex indexer
    midx = pd.MultiIndex.from_arrays(
        [["x", "y", "z"], [1, 3, 5], [7, 8, 9]], names=["x", "a", "z"]
    )[:2]
    idx = pd.Index([1, 3, 5], name="a")

    result, lidx, ridx = idx.join(midx, how="left", return_indexers=True)

    expected = pd.MultiIndex(
        levels=[["x", "y", "z"], [1, 3, 5], [7, 8, 9]],
        codes=[[0, 1, -1], [0, 1, 2], [0, 1, -1]],
        names=["x", "a", "z"],
    )
    tm.assert_index_equal(result, expected)
    tm.assert_numpy_array_equal(lidx, np.array([0, 1, 2], dtype=np.intp))
    tm.assert_numpy_array_equal(ridx, np.array([0, 1, -1], dtype=np.intp))


def test_join_flat_with_multi_unmatched_overlapping_intervals():
    # GH#34292 the level of an overlapping IntervalIndex cannot use get_indexer
    idx = pd.IntervalIndex.from_tuples([(0, 2), (1, 3), (5, 6)], name="iv")
    midx = pd.MultiIndex.from_arrays([idx[:2], ["x", "y"]], names=["iv", "k"])

    result, lidx, ridx = idx.join(midx, how="left", return_indexers=True)

    expected = pd.MultiIndex.from_arrays([idx, ["x", "y", np.nan]], names=["iv", "k"])
    tm.assert_index_equal(result, expected)
    tm.assert_numpy_array_equal(lidx, np.array([0, 1, 2], dtype=np.intp))
    tm.assert_numpy_array_equal(ridx, np.array([0, 1, -1], dtype=np.intp))


def test_join_multi_wrong_order():
    # GH 25760
    # GH 28956

    midx1 = pd.MultiIndex.from_product([[1, 2], [3, 4]], names=["a", "b"])
    midx2 = pd.MultiIndex.from_product([[1, 2], [3, 4]], names=["b", "a"])

    join_idx, lidx, ridx = midx1.join(midx2, return_indexers=True)

    exp_ridx = np.array([-1, -1, -1, -1], dtype=np.intp)

    tm.assert_index_equal(midx1, join_idx)
    assert lidx is None
    tm.assert_numpy_array_equal(ridx, exp_ridx)


def test_join_multi_return_indexers():
    # GH 34074

    midx1 = pd.MultiIndex.from_product([[1, 2], [3, 4], [5, 6]], names=["a", "b", "c"])
    midx2 = pd.MultiIndex.from_product([[1, 2], [3, 4]], names=["a", "b"])

    result = midx1.join(midx2, return_indexers=False)
    tm.assert_index_equal(result, midx1)


def test_join_overlapping_interval_level():
    # GH 44096
    idx_1 = pd.MultiIndex.from_tuples(
        [
            (1, pd.Interval(0.0, 1.0)),
            (1, pd.Interval(1.0, 2.0)),
            (1, pd.Interval(2.0, 5.0)),
            (2, pd.Interval(0.0, 1.0)),
            (2, pd.Interval(1.0, 3.0)),  # interval limit is here at 3.0, not at 2.0
            (2, pd.Interval(3.0, 5.0)),
        ],
        names=["num", "interval"],
    )

    idx_2 = pd.MultiIndex.from_tuples(
        [
            (1, pd.Interval(2.0, 5.0)),
            (1, pd.Interval(0.0, 1.0)),
            (1, pd.Interval(1.0, 2.0)),
            (2, pd.Interval(3.0, 5.0)),
            (2, pd.Interval(0.0, 1.0)),
            (2, pd.Interval(1.0, 3.0)),
        ],
        names=["num", "interval"],
    )

    expected = pd.MultiIndex.from_tuples(
        [
            (1, pd.Interval(0.0, 1.0)),
            (1, pd.Interval(1.0, 2.0)),
            (1, pd.Interval(2.0, 5.0)),
            (2, pd.Interval(0.0, 1.0)),
            (2, pd.Interval(1.0, 3.0)),
            (2, pd.Interval(3.0, 5.0)),
        ],
        names=["num", "interval"],
    )
    result = idx_1.join(idx_2, how="outer")

    tm.assert_index_equal(result, expected)


def test_join_midx_ea():
    # GH#49277
    midx = pd.MultiIndex.from_arrays(
        [pd.Series([1, 1, 3], dtype="Int64"), pd.Series([1, 2, 3], dtype="Int64")],
        names=["a", "b"],
    )
    midx2 = pd.MultiIndex.from_arrays(
        [pd.Series([1], dtype="Int64"), pd.Series([3], dtype="Int64")], names=["a", "c"]
    )
    result = midx.join(midx2, how="inner")
    expected = pd.MultiIndex.from_arrays(
        [
            pd.Series([1, 1], dtype="Int64"),
            pd.Series([1, 2], dtype="Int64"),
            pd.Series([3, 3], dtype="Int64"),
        ],
        names=["a", "b", "c"],
    )
    tm.assert_index_equal(result, expected)


def test_join_midx_string():
    # GH#49277
    midx = pd.MultiIndex.from_arrays(
        [
            pd.Series(["a", "a", "c"], dtype=pd.StringDtype()),
            pd.Series(["a", "b", "c"], dtype=pd.StringDtype()),
        ],
        names=["a", "b"],
    )
    midx2 = pd.MultiIndex.from_arrays(
        [
            pd.Series(["a"], dtype=pd.StringDtype()),
            pd.Series(["c"], dtype=pd.StringDtype()),
        ],
        names=["a", "c"],
    )
    result = midx.join(midx2, how="inner")
    expected = pd.MultiIndex.from_arrays(
        [
            pd.Series(["a", "a"], dtype=pd.StringDtype()),
            pd.Series(["a", "b"], dtype=pd.StringDtype()),
            pd.Series(["c", "c"], dtype=pd.StringDtype()),
        ],
        names=["a", "b", "c"],
    )
    tm.assert_index_equal(result, expected)


def test_join_multi_with_nan():
    # GH29252
    df1 = pd.DataFrame(
        data={"col1": [1.1, 1.2]},
        index=pd.MultiIndex.from_product([["A"], [1.0, 2.0]], names=["id1", "id2"]),
    )
    df2 = pd.DataFrame(
        data={"col2": [2.1, 2.2]},
        index=pd.MultiIndex.from_product([["A"], [np.nan, 2.0]], names=["id1", "id2"]),
    )
    result = df1.join(df2)
    expected = pd.DataFrame(
        data={"col1": [1.1, 1.2], "col2": [np.nan, 2.2]},
        index=pd.MultiIndex.from_product([["A"], [1.0, 2.0]], names=["id1", "id2"]),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("val", [0, 5])
def test_join_dtypes(any_numeric_ea_dtype, val):
    # GH#49830
    midx = pd.MultiIndex.from_arrays(
        [pd.Series([1, 2], dtype=any_numeric_ea_dtype), [3, 4]]
    )
    midx2 = pd.MultiIndex.from_arrays(
        [pd.Series([1, val, val], dtype=any_numeric_ea_dtype), [3, 4, 4]]
    )
    result = midx.join(midx2, how="outer")
    expected = pd.MultiIndex.from_arrays(
        [pd.Series([val, val, 1, 2], dtype=any_numeric_ea_dtype), [4, 4, 3, 4]]
    ).sort_values()
    tm.assert_index_equal(result, expected)


def test_join_dtypes_all_nan(any_numeric_ea_dtype):
    # GH#49830
    midx = pd.MultiIndex.from_arrays(
        [pd.Series([1, 2], dtype=any_numeric_ea_dtype), [np.nan, np.nan]]
    )
    midx2 = pd.MultiIndex.from_arrays(
        [pd.Series([1, 0, 0], dtype=any_numeric_ea_dtype), [np.nan, np.nan, np.nan]]
    )
    result = midx.join(midx2, how="outer")
    expected = pd.MultiIndex.from_arrays(
        [
            pd.Series([0, 0, 1, 2], dtype=any_numeric_ea_dtype),
            [np.nan, np.nan, np.nan, np.nan],
        ]
    )
    tm.assert_index_equal(result, expected)


def test_join_index_levels():
    # GH#53093
    midx = pd.MultiIndex.from_tuples([("a", "2019-02-01"), ("a", "2019-02-01")])
    midx2 = pd.MultiIndex.from_tuples([("a", "2019-01-31")])
    result = midx.join(midx2, how="outer")
    expected = pd.MultiIndex.from_tuples(
        [("a", "2019-01-31"), ("a", "2019-02-01"), ("a", "2019-02-01")]
    )
    tm.assert_index_equal(result.levels[1], expected.levels[1])
    tm.assert_index_equal(result, expected)
