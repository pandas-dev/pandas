"""
Tests for values coercion in setitem-like operations on DataFrame.

For the most part, these should be multi-column DataFrames, otherwise
we would share the tests with Series.
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestDataFrameSetitemCoercion:
    @pytest.mark.parametrize("consolidate", [True, False])
    def test_loc_setitem_multiindex_columns(self, consolidate):
        # GH#18415 Setting values in a single column preserves dtype,
        #  while setting them in multiple columns did unwanted cast.

        # Note that A here has 2 blocks, below we do the same thing
        #  with a consolidated frame.
        A = pd.DataFrame(np.zeros((6, 5), dtype=np.float32))
        A = pd.concat([A, A], axis=1, keys=[1, 2])
        if consolidate:
            A = A._consolidate()

        A.loc[2:3, (1, slice(2, 3))] = np.ones((2, 2), dtype=np.float32)
        assert (A.dtypes == np.float32).all()

        A.loc[0:5, (1, slice(2, 3))] = np.ones((6, 2), dtype=np.float32)

        assert (A.dtypes == np.float32).all()

        A.loc[:, (1, slice(2, 3))] = np.ones((6, 2), dtype=np.float32)
        assert (A.dtypes == np.float32).all()

        # TODO: i think this isn't about MultiIndex and could be done with iloc?


def test_37477():
    # fixed by GH#45121
    orig = pd.DataFrame({"A": [1, 2, 3], "B": [3, 4, 5]})

    df = orig.copy()
    with pytest.raises(TypeError, match="Invalid value"):
        df.at[1, "B"] = 1.2

    with pytest.raises(TypeError, match="Invalid value"):
        df.loc[1, "B"] = 1.2

    with pytest.raises(TypeError, match="Invalid value"):
        df.iat[1, 1] = 1.2

    with pytest.raises(TypeError, match="Invalid value"):
        df.iloc[1, 1] = 1.2


def test_6942(indexer_al):
    # check that the .at __setitem__ after setting "Live" actually sets the data
    start = pd.Timestamp("2014-04-01")
    t1 = pd.Timestamp("2014-04-23 12:42:38.883082")
    t2 = pd.Timestamp("2014-04-24 01:33:30.040039")

    dti = pd.date_range(start, periods=1)
    orig = pd.DataFrame(index=dti, columns=["timenow", "Live"])

    df = orig.copy()
    indexer_al(df)[start, "timenow"] = t1

    df["Live"] = True

    df.at[start, "timenow"] = t2
    assert df.iloc[0, 0] == t2


def test_26395(indexer_al):
    # .at case fixed by GH#45121 (best guess)
    df = pd.DataFrame(index=["A", "B", "C"])
    df["D"] = 0

    indexer_al(df)["C", "D"] = 2
    expected = pd.DataFrame({"D": [0, 0, 2]}, index=["A", "B", "C"], dtype=np.int64)
    tm.assert_frame_equal(df, expected)

    with pytest.raises(TypeError, match="Invalid value"):
        indexer_al(df)["C", "D"] = 44.5

    with pytest.raises(TypeError, match="Invalid value"):
        indexer_al(df)["C", "D"] = "hello"


def test_15231():
    df = pd.DataFrame([[1, 2], [3, 4]], columns=["a", "b"])
    df.loc[2] = pd.Series({"a": 5, "b": 6})
    assert (df.dtypes == np.int64).all()

    df.loc[3] = pd.Series({"a": 7})

    # df["a"] doesn't have any NaNs, should not have been cast
    exp_dtypes = pd.Series([np.int64, np.float64], dtype=object, index=["a", "b"])
    tm.assert_series_equal(df.dtypes, exp_dtypes)


def test_iloc_setitem_unnecesssary_float_upcasting():
    # GH#12255
    df = pd.DataFrame(
        {
            0: np.array([1, 3], dtype=np.float32),
            1: np.array([2, 4], dtype=np.float32),
            2: ["a", "b"],
        }
    )
    orig = df.copy()

    values = df[0].values.reshape(2, 1)
    df.iloc[:, 0:1] = values

    tm.assert_frame_equal(df, orig)


@pytest.mark.parametrize("indexer", [1, (1, slice(None))])
@pytest.mark.parametrize(
    "box",
    [list, tuple, lambda x: np.array(x, dtype=object)],
    ids=["list", "tuple", "array"],
)
def test_12499(indexer, box):
    # GH#12499 setting a new row must not infer one dtype for the whole row,
    #  which turned the NaN into NaT
    ts = pd.Timestamp("2016-03-01 03:13:22.98986", tz="UTC")
    df = pd.DataFrame({"one": [0], "two": [ts]})
    df.loc[indexer] = box([np.nan, pd.NaT])

    expected = pd.DataFrame(
        {"one": [0, np.nan], "two": pd.Series([ts, pd.NaT], dtype="M8[us, UTC]")}
    )
    tm.assert_frame_equal(df, expected)


def test_12499_extension_columns():
    # GH#12499 the row infers as period, which must not turn NaN into NaT
    per = pd.Period("2020", "D")
    df = pd.DataFrame({"a": [per], "b": [pd.Interval(0.0, 1.0)]})
    df.loc[1] = [per, np.nan]

    expected = pd.DataFrame(
        {
            "a": [per, per],
            "b": pd.array([pd.Interval(0.0, 1.0), np.nan], dtype="interval[float64]"),
        }
    )
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize(
    "row, expected",
    [
        ([2**60 + 1, 1.5], [1, 2**60 + 1]),
        ([np.nan, pd.Timestamp("2021")], [1, np.nan]),
    ],
)
def test_12499_object_column(row, expected):
    # GH#12499 an object column would keep a value converted by whole-row
    #  inference
    df = pd.DataFrame({"a": pd.Series([1], dtype=object), "b": [row[1]]})
    df.loc[1] = row

    tm.assert_numpy_array_equal(
        df["a"].to_numpy(), np.array(expected, dtype=object), strict_nan=True
    )


def test_12499_large_int_not_rounded():
    # GH#12499 inferring float64 for the whole row would round the int
    df = pd.DataFrame({"a": [0], "b": [0.5]})
    df.loc[1] = [2**60 + 1, 0.5]

    expected = pd.DataFrame({"a": [0, 2**60 + 1], "b": [0.5, 0.5]})
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize(
    "box", [list, lambda x: dict(enumerate(x))], ids=["list", "dict"]
)
def test_12499_row_expansion_not_fragmented(box):
    # GH#12499 an object row must not leave one block per column
    df = pd.DataFrame({0: [1.5], 1: [2.5], 2: [1], 3: [2]})
    df.loc[1] = box([3.5, 4.5, 3, 4])

    expected = pd.DataFrame({0: [1.5, 3.5], 1: [2.5, 4.5], 2: [1, 3], 3: [2, 4]})
    tm.assert_frame_equal(df, expected)
    assert df._mgr.nblocks == 2


def test_20476():
    mi = pd.MultiIndex.from_product([["A", "B"], ["a", "b", "c"]])
    df = pd.DataFrame(-1, index=range(3), columns=mi)
    filler = pd.DataFrame([[1, 2, 3.0]] * 3, index=range(3), columns=["a", "b", "c"])
    df["A"] = filler

    expected = pd.DataFrame(
        {
            0: [1, 1, 1],
            1: [2, 2, 2],
            2: [3.0, 3.0, 3.0],
            3: [-1, -1, -1],
            4: [-1, -1, -1],
            5: [-1, -1, -1],
        }
    )
    expected.columns = mi
    exp_dtypes = pd.Series(
        [np.dtype(np.int64)] * 2 + [np.dtype(np.float64)] + [np.dtype(np.int64)] * 3,
        index=mi,
    )
    tm.assert_series_equal(df.dtypes, exp_dtypes)


def test_setitem_midx_subframe_keeps_mixed_dtypes():
    # GH#13433
    df = pd.DataFrame({("test", "it"): [23.0, 24.0], ("test", "not"): ["a", "a"]})
    expected = df.copy()
    df["test"] = df["test"]
    tm.assert_frame_equal(df, expected)
