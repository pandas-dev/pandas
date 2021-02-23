import numpy as np
import pytest

from pandas.core.dtypes.base import registry as ea_registry
from pandas.core.dtypes.dtypes import DatetimeTZDtype, IntervalDtype, PeriodDtype

from pandas import (
    Categorical,
    DataFrame,
    Index,
    Interval,
    NaT,
    Period,
    PeriodIndex,
    Series,
    Timestamp,
    date_range,
    notna,
    period_range,
)
import pandas._testing as tm
from pandas.core.arrays import SparseArray


class TestDataFrameSetItem:
    @pytest.mark.parametrize("dtype", ["int32", "int64", "float32", "float64"])
    def test_setitem_dtype(self, dtype, float_frame):
        arr = np.random.randn(len(float_frame))

        float_frame[dtype] = np.array(arr, dtype=dtype)
        assert float_frame[dtype].dtype.name == dtype

    def test_setitem_list_not_dataframe(self, float_frame):
        data = np.random.randn(len(float_frame), 2)
        float_frame[["A", "B"]] = data
        tm.assert_almost_equal(float_frame[["A", "B"]].values, data)

    def test_setitem_error_msmgs(self):

        # GH 7432
        df = DataFrame(
            {"bar": [1, 2, 3], "baz": ["d", "e", "f"]},
            index=Index(["a", "b", "c"], name="foo"),
        )
        ser = Series(
            ["g", "h", "i", "j"],
            index=Index(["a", "b", "c", "a"], name="foo"),
            name="fiz",
        )
        msg = "cannot reindex from a duplicate axis"
        with pytest.raises(ValueError, match=msg):
            df["newcol"] = ser

        # GH 4107, more descriptive error message
        df = DataFrame(np.random.randint(0, 2, (4, 4)), columns=["a", "b", "c", "d"])

        msg = "incompatible index of inserted column with frame index"
        with pytest.raises(TypeError, match=msg):
            df["gr"] = df.groupby(["b", "c"]).count()

    def test_setitem_benchmark(self):
        # from the vb_suite/frame_methods/frame_insert_columns
        N = 10
        K = 5
        df = DataFrame(index=range(N))
        new_col = np.random.randn(N)
        for i in range(K):
            df[i] = new_col
        expected = DataFrame(np.repeat(new_col, K).reshape(N, K), index=range(N))
        tm.assert_frame_equal(df, expected)

    def test_setitem_different_dtype(self):
        df = DataFrame(
            np.random.randn(5, 3), index=np.arange(5), columns=["c", "b", "a"]
        )
        df.insert(0, "foo", df["a"])
        df.insert(2, "bar", df["c"])

        # diff dtype

        # new item
        df["x"] = df["a"].astype("float32")
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")] * 5 + [np.dtype("float32")],
            index=["foo", "c", "bar", "b", "a", "x"],
        )
        tm.assert_series_equal(result, expected)

        # replacing current (in different block)
        df["a"] = df["a"].astype("float32")
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")] * 4 + [np.dtype("float32")] * 2,
            index=["foo", "c", "bar", "b", "a", "x"],
        )
        tm.assert_series_equal(result, expected)

        df["y"] = df["a"].astype("int32")
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")] * 4 + [np.dtype("float32")] * 2 + [np.dtype("int32")],
            index=["foo", "c", "bar", "b", "a", "x", "y"],
        )
        tm.assert_series_equal(result, expected)

    def test_setitem_empty_columns(self):
        # GH 13522
        df = DataFrame(index=["A", "B", "C"])
        df["X"] = df.index
        df["X"] = ["x", "y", "z"]
        exp = DataFrame(data={"X": ["x", "y", "z"]}, index=["A", "B", "C"])
        tm.assert_frame_equal(df, exp)

    def test_setitem_dt64_index_empty_columns(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        df = DataFrame(index=np.arange(len(rng)))

        df["A"] = rng
        assert df["A"].dtype == np.dtype("M8[ns]")

    def test_setitem_timestamp_empty_columns(self):
        # GH#19843
        df = DataFrame(index=range(3))
        df["now"] = Timestamp("20130101", tz="UTC")

        expected = DataFrame(
            [[Timestamp("20130101", tz="UTC")]] * 3, index=[0, 1, 2], columns=["now"]
        )
        tm.assert_frame_equal(df, expected)

    def test_setitem_wrong_length_categorical_dtype_raises(self):
        # GH#29523
        cat = Categorical.from_codes([0, 1, 1, 0, 1, 2], ["a", "b", "c"])
        df = DataFrame(range(10), columns=["bar"])

        msg = (
            rf"Length of values \({len(cat)}\) "
            rf"does not match length of index \({len(df)}\)"
        )
        with pytest.raises(ValueError, match=msg):
            df["foo"] = cat

    def test_setitem_with_sparse_value(self):
        # GH#8131
        df = DataFrame({"c_1": ["a", "b", "c"], "n_1": [1.0, 2.0, 3.0]})
        sp_array = SparseArray([0, 0, 1])
        df["new_column"] = sp_array

        expected = Series(sp_array, name="new_column")
        tm.assert_series_equal(df["new_column"], expected)

    def test_setitem_with_unaligned_sparse_value(self):
        df = DataFrame({"c_1": ["a", "b", "c"], "n_1": [1.0, 2.0, 3.0]})
        sp_series = Series(SparseArray([0, 0, 1]), index=[2, 1, 0])

        df["new_column"] = sp_series
        expected = Series(SparseArray([1, 0, 0]), name="new_column")
        tm.assert_series_equal(df["new_column"], expected)

    def test_setitem_dict_preserves_dtypes(self):
        # https://github.com/pandas-dev/pandas/issues/34573
        expected = DataFrame(
            {
                "a": Series([0, 1, 2], dtype="int64"),
                "b": Series([1, 2, 3], dtype=float),
                "c": Series([1, 2, 3], dtype=float),
            }
        )
        df = DataFrame(
            {
                "a": Series([], dtype="int64"),
                "b": Series([], dtype=float),
                "c": Series([], dtype=float),
            }
        )
        for idx, b in enumerate([1, 2, 3]):
            df.loc[df.shape[0]] = {"a": int(idx), "b": float(b), "c": float(b)}
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "obj,dtype",
        [
            (Period("2020-01"), PeriodDtype("M")),
            (Interval(left=0, right=5), IntervalDtype("int64")),
            (
                Timestamp("2011-01-01", tz="US/Eastern"),
                DatetimeTZDtype(tz="US/Eastern"),
            ),
        ],
    )
    def test_setitem_extension_types(self, obj, dtype):
        # GH: 34832
        expected = DataFrame({"idx": [1, 2, 3], "obj": Series([obj] * 3, dtype=dtype)})

        df = DataFrame({"idx": [1, 2, 3]})
        df["obj"] = obj

        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "ea_name",
        [
            dtype.name
            for dtype in ea_registry.dtypes
            # property would require instantiation
            if not isinstance(dtype.name, property)
        ]
        # mypy doesn't allow adding lists of different types
        # https://github.com/python/mypy/issues/5492
        + ["datetime64[ns, UTC]", "period[D]"],  # type: ignore[list-item]
    )
    def test_setitem_with_ea_name(self, ea_name):
        # GH 38386
        result = DataFrame([0])
        result[ea_name] = [1]
        expected = DataFrame({0: [0], ea_name: [1]})
        tm.assert_frame_equal(result, expected)

    def test_setitem_dt64_ndarray_with_NaT_and_diff_time_units(self):
        # GH#7492
        data_ns = np.array([1, "nat"], dtype="datetime64[ns]")
        result = Series(data_ns).to_frame()
        result["new"] = data_ns
        expected = DataFrame({0: [1, None], "new": [1, None]}, dtype="datetime64[ns]")
        tm.assert_frame_equal(result, expected)

        # OutOfBoundsDatetime error shouldn't occur
        data_s = np.array([1, "nat"], dtype="datetime64[s]")
        result["new"] = data_s
        expected = DataFrame({0: [1, None], "new": [1e9, None]}, dtype="datetime64[ns]")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["h", "m", "s", "ms", "D", "M", "Y"])
    def test_frame_setitem_datetime64_col_other_units(self, unit):
        # Check that non-nano dt64 values get cast to dt64 on setitem
        #  into a not-yet-existing column
        n = 100

        dtype = np.dtype(f"M8[{unit}]")
        vals = np.arange(n, dtype=np.int64).view(dtype)
        ex_vals = vals.astype("datetime64[ns]")

        df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
        df[unit] = vals

        assert df[unit].dtype == np.dtype("M8[ns]")
        assert (df[unit].values == ex_vals).all()

    @pytest.mark.parametrize("unit", ["h", "m", "s", "ms", "D", "M", "Y"])
    def test_frame_setitem_existing_datetime64_col_other_units(self, unit):
        # Check that non-nano dt64 values get cast to dt64 on setitem
        #  into an already-existing dt64 column
        n = 100

        dtype = np.dtype(f"M8[{unit}]")
        vals = np.arange(n, dtype=np.int64).view(dtype)
        ex_vals = vals.astype("datetime64[ns]")

        df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
        df["dates"] = np.arange(n, dtype=np.int64).view("M8[ns]")

        # We overwrite existing dt64 column with new, non-nano dt64 vals
        df["dates"] = vals
        assert (df["dates"].values == ex_vals).all()

    def test_setitem_dt64tz(self, timezone_frame):

        df = timezone_frame
        idx = df["B"].rename("foo")

        # setitem
        df["C"] = idx
        tm.assert_series_equal(df["C"], Series(idx, name="C"))

        df["D"] = "foo"
        df["D"] = idx
        tm.assert_series_equal(df["D"], Series(idx, name="D"))
        del df["D"]

        # assert that A & C are not sharing the same base (e.g. they
        # are copies)
        b1 = df._mgr.blocks[1]
        b2 = df._mgr.blocks[2]
        tm.assert_extension_array_equal(b1.values, b2.values)
        b1base = b1.values._data.base
        b2base = b2.values._data.base
        assert b1base is None or (id(b1base) != id(b2base))

        # with nan
        df2 = df.copy()
        df2.iloc[1, 1] = NaT
        df2.iloc[1, 2] = NaT
        result = df2["B"]
        tm.assert_series_equal(notna(result), Series([True, False, True], name="B"))
        tm.assert_series_equal(df2.dtypes, df.dtypes)

    def test_setitem_periodindex(self):
        rng = period_range("1/1/2000", periods=5, name="index")
        df = DataFrame(np.random.randn(5, 3), index=rng)

        df["Index"] = rng
        rs = Index(df["Index"])
        tm.assert_index_equal(rs, rng, check_names=False)
        assert rs.name == "Index"
        assert rng.name == "index"

        rs = df.reset_index().set_index("index")
        assert isinstance(rs.index, PeriodIndex)
        tm.assert_index_equal(rs.index, rng)

    def test_setitem_complete_column_with_array(self):
        # GH#37954
        df = DataFrame({"a": ["one", "two", "three"], "b": [1, 2, 3]})
        arr = np.array([[1, 1], [3, 1], [5, 1]])
        df[["c", "d"]] = arr
        expected = DataFrame(
            {
                "a": ["one", "two", "three"],
                "b": [1, 2, 3],
                "c": [1, 3, 5],
                "d": [1, 1, 1],
            }
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("dtype", ["f8", "i8", "u8"])
    def test_setitem_bool_with_numeric_index(self, dtype):
        # GH#36319
        cols = Index([1, 2, 3], dtype=dtype)
        df = DataFrame(np.random.randn(3, 3), columns=cols)

        df[False] = ["a", "b", "c"]

        expected_cols = Index([1, 2, 3, False], dtype=object)
        if dtype == "f8":
            expected_cols = Index([1.0, 2.0, 3.0, False], dtype=object)

        tm.assert_index_equal(df.columns, expected_cols)


class TestDataFrameSetItemWithExpansion:
    def test_setitem_listlike_views(self):
        # GH#38148
        df = DataFrame({"a": [1, 2, 3], "b": [4, 4, 6]})

        # get one column as a view of df
        ser = df["a"]

        # add columns with list-like indexer
        df[["c", "d"]] = np.array([[0.1, 0.2], [0.3, 0.4], [0.4, 0.5]])

        # edit in place the first column to check view semantics
        df.iloc[0, 0] = 100

        expected = Series([100, 2, 3], name="a")
        tm.assert_series_equal(ser, expected)

    def test_setitem_string_column_numpy_dtype_raising(self):
        # GH#39010
        df = DataFrame([[1, 2], [3, 4]])
        df["0 - Name"] = [5, 6]
        expected = DataFrame([[1, 2, 5], [3, 4, 6]], columns=[0, 1, "0 - Name"])
        tm.assert_frame_equal(df, expected)


class TestDataFrameSetItemSlicing:
    def test_setitem_slice_position(self):
        # GH#31469
        df = DataFrame(np.zeros((100, 1)))
        df[-4:] = 1
        arr = np.zeros((100, 1))
        arr[-4:] = 1
        expected = DataFrame(arr)
        tm.assert_frame_equal(df, expected)


class TestDataFrameSetItemCallable:
    def test_setitem_callable(self):
        # GH#12533
        df = DataFrame({"A": [1, 2, 3, 4], "B": [5, 6, 7, 8]})
        df[lambda x: "A"] = [11, 12, 13, 14]

        exp = DataFrame({"A": [11, 12, 13, 14], "B": [5, 6, 7, 8]})
        tm.assert_frame_equal(df, exp)


class TestDataFrameSetItemBooleanMask:
    @pytest.mark.parametrize(
        "mask_type",
        [lambda df: df > np.abs(df) / 2, lambda df: (df > np.abs(df) / 2).values],
        ids=["dataframe", "array"],
    )
    def test_setitem_boolean_mask(self, mask_type, float_frame):

        # Test for issue #18582
        df = float_frame.copy()
        mask = mask_type(df)

        # index with boolean mask
        result = df.copy()
        result[mask] = np.nan

        expected = df.copy()
        expected.values[np.array(mask)] = np.nan
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("indexer", [lambda x: x, lambda x: x.loc])
    def test_setitem_boolean_mask_aligning(self, indexer):
        # GH#39931
        df = DataFrame({"a": [1, 4, 2, 3], "b": [5, 6, 7, 8]})
        expected = df.copy()
        mask = df["a"] >= 3
        indexer(df)[mask] = indexer(df)[mask].sort_values("a")
        tm.assert_frame_equal(df, expected)
