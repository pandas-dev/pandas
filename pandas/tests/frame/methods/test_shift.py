import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    CategoricalIndex,
    DataFrame,
    Index,
    NaT,
    Series,
    date_range,
    offsets,
)
import pandas._testing as tm


class TestDataFrameShift:
    def test_shift_int(self, datetime_frame, frame_or_series):
        ts = tm.get_obj(datetime_frame, frame_or_series).astype(int)
        shifted = ts.shift(1)
        expected = ts.astype(float).shift(1)
        tm.assert_equal(shifted, expected)

    def test_shift_32bit_take(self, frame_or_series):
        # 32-bit taking
        # GH#8129
        index = date_range("2000-01-01", periods=5)
        for dtype in ["int32", "int64"]:
            arr = np.arange(5, dtype=dtype)
            s1 = frame_or_series(arr, index=index)
            p = arr[1]
            result = s1.shift(periods=p)
            expected = frame_or_series([np.nan, 0, 1, 2, 3], index=index)
            tm.assert_equal(result, expected)

    @pytest.mark.parametrize("periods", [1, 2, 3, 4])
    def test_shift_preserve_freqstr(self, periods, frame_or_series):
        # GH#21275
        obj = frame_or_series(
            range(periods),
            index=date_range("2016-1-1 00:00:00", periods=periods, freq="H"),
        )

        result = obj.shift(1, "2H")

        expected = frame_or_series(
            range(periods),
            index=date_range("2016-1-1 02:00:00", periods=periods, freq="H"),
        )
        tm.assert_equal(result, expected)

    def test_shift_dst(self, frame_or_series):
        # GH#13926
        dates = date_range("2016-11-06", freq="H", periods=10, tz="US/Eastern")
        obj = frame_or_series(dates)

        res = obj.shift(0)
        tm.assert_equal(res, obj)
        assert tm.get_dtype(res) == "datetime64[ns, US/Eastern]"

        res = obj.shift(1)
        exp_vals = [NaT] + dates.astype(object).values.tolist()[:9]
        exp = frame_or_series(exp_vals)
        tm.assert_equal(res, exp)
        assert tm.get_dtype(res) == "datetime64[ns, US/Eastern]"

        res = obj.shift(-2)
        exp_vals = dates.astype(object).values.tolist()[2:] + [NaT, NaT]
        exp = frame_or_series(exp_vals)
        tm.assert_equal(res, exp)
        assert tm.get_dtype(res) == "datetime64[ns, US/Eastern]"

        for ex in [10, -10, 20, -20]:
            res = obj.shift(ex)
            exp = frame_or_series([NaT] * 10, dtype="datetime64[ns, US/Eastern]")
            tm.assert_equal(res, exp)
            assert tm.get_dtype(res) == "datetime64[ns, US/Eastern]"

    def test_shift(self, datetime_frame, frame_or_series):
        # naive shift
        shiftedFrame = datetime_frame.shift(5)
        tm.assert_index_equal(shiftedFrame.index, datetime_frame.index)

        shiftedSeries = datetime_frame["A"].shift(5)
        tm.assert_series_equal(shiftedFrame["A"], shiftedSeries)

        shiftedFrame = datetime_frame.shift(-5)
        tm.assert_index_equal(shiftedFrame.index, datetime_frame.index)

        shiftedSeries = datetime_frame["A"].shift(-5)
        tm.assert_series_equal(shiftedFrame["A"], shiftedSeries)

        # shift by 0
        unshifted = datetime_frame.shift(0)
        tm.assert_frame_equal(unshifted, datetime_frame)

        # shift by DateOffset
        shiftedFrame = datetime_frame.shift(5, freq=offsets.BDay())
        assert len(shiftedFrame) == len(datetime_frame)

        shiftedFrame2 = datetime_frame.shift(5, freq="B")
        tm.assert_frame_equal(shiftedFrame, shiftedFrame2)

        d = datetime_frame.index[0]
        shifted_d = d + offsets.BDay(5)
        tm.assert_series_equal(
            datetime_frame.xs(d), shiftedFrame.xs(shifted_d), check_names=False
        )

    def test_shift_with_periodindex(self, frame_or_series):
        # Shifting with PeriodIndex
        ps = tm.makePeriodFrame()
        ps = tm.get_obj(ps, frame_or_series)

        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        tm.assert_index_equal(shifted.index, ps.index)
        tm.assert_index_equal(unshifted.index, ps.index)
        if frame_or_series is DataFrame:
            tm.assert_numpy_array_equal(
                unshifted.iloc[:, 0].dropna().values, ps.iloc[:-1, 0].values
            )
        else:
            tm.assert_numpy_array_equal(unshifted.dropna().values, ps.values[:-1])

        shifted2 = ps.shift(1, "B")
        shifted3 = ps.shift(1, offsets.BDay())
        tm.assert_equal(shifted2, shifted3)
        tm.assert_equal(ps, shifted2.shift(-1, "B"))

        msg = "does not match PeriodIndex freq"
        with pytest.raises(ValueError, match=msg):
            ps.shift(freq="D")

        # legacy support
        shifted4 = ps.shift(1, freq="B")
        tm.assert_equal(shifted2, shifted4)

        shifted5 = ps.shift(1, freq=offsets.BDay())
        tm.assert_equal(shifted5, shifted4)

    def test_shift_other_axis(self):
        # shift other axis
        # GH#6371
        df = DataFrame(np.random.rand(10, 5))
        expected = pd.concat(
            [DataFrame(np.nan, index=df.index, columns=[0]), df.iloc[:, 0:-1]],
            ignore_index=True,
            axis=1,
        )
        result = df.shift(1, axis=1)
        tm.assert_frame_equal(result, expected)

    def test_shift_named_axis(self):
        # shift named axis
        df = DataFrame(np.random.rand(10, 5))
        expected = pd.concat(
            [DataFrame(np.nan, index=df.index, columns=[0]), df.iloc[:, 0:-1]],
            ignore_index=True,
            axis=1,
        )
        result = df.shift(1, axis="columns")
        tm.assert_frame_equal(result, expected)

    def test_shift_bool(self):
        df = DataFrame({"high": [True, False], "low": [False, False]})
        rs = df.shift(1)
        xp = DataFrame(
            np.array([[np.nan, np.nan], [True, False]], dtype=object),
            columns=["high", "low"],
        )
        tm.assert_frame_equal(rs, xp)

    def test_shift_categorical(self):
        # GH#9416
        s1 = Series(["a", "b", "c"], dtype="category")
        s2 = Series(["A", "B", "C"], dtype="category")
        df = DataFrame({"one": s1, "two": s2})
        rs = df.shift(1)
        xp = DataFrame({"one": s1.shift(1), "two": s2.shift(1)})
        tm.assert_frame_equal(rs, xp)

    def test_shift_categorical_fill_value(self, frame_or_series):
        ts = frame_or_series(["a", "b", "c", "d"], dtype="category")
        res = ts.shift(1, fill_value="a")
        expected = frame_or_series(
            pd.Categorical(
                ["a", "a", "b", "c"], categories=["a", "b", "c", "d"], ordered=False
            )
        )
        tm.assert_equal(res, expected)

        # check for incorrect fill_value
        msg = r"Cannot setitem on a Categorical with a new category \(f\)"
        with pytest.raises(TypeError, match=msg):
            ts.shift(1, fill_value="f")

    def test_shift_fill_value(self, frame_or_series):
        # GH#24128
        dti = date_range("1/1/2000", periods=5, freq="H")

        ts = frame_or_series([1.0, 2.0, 3.0, 4.0, 5.0], index=dti)
        exp = frame_or_series([0.0, 1.0, 2.0, 3.0, 4.0], index=dti)
        # check that fill value works
        result = ts.shift(1, fill_value=0.0)
        tm.assert_equal(result, exp)

        exp = frame_or_series([0.0, 0.0, 1.0, 2.0, 3.0], index=dti)
        result = ts.shift(2, fill_value=0.0)
        tm.assert_equal(result, exp)

        ts = frame_or_series([1, 2, 3])
        res = ts.shift(2, fill_value=0)
        assert tm.get_dtype(res) == tm.get_dtype(ts)

        # retain integer dtype
        obj = frame_or_series([1, 2, 3, 4, 5], index=dti)
        exp = frame_or_series([0, 1, 2, 3, 4], index=dti)
        result = obj.shift(1, fill_value=0)
        tm.assert_equal(result, exp)

        exp = frame_or_series([0, 0, 1, 2, 3], index=dti)
        result = obj.shift(2, fill_value=0)
        tm.assert_equal(result, exp)

    def test_shift_empty(self):
        # Regression test for GH#8019
        df = DataFrame({"foo": []})
        rs = df.shift(-1)

        tm.assert_frame_equal(df, rs)

    def test_shift_duplicate_columns(self):
        # GH#9092; verify that position-based shifting works
        # in the presence of duplicate columns
        column_lists = [list(range(5)), [1] * 5, [1, 1, 2, 2, 1]]
        data = np.random.randn(20, 5)

        shifted = []
        for columns in column_lists:
            df = DataFrame(data.copy(), columns=columns)
            for s in range(5):
                df.iloc[:, s] = df.iloc[:, s].shift(s + 1)
            df.columns = range(5)
            shifted.append(df)

        # sanity check the base case
        nulls = shifted[0].isna().sum()
        tm.assert_series_equal(nulls, Series(range(1, 6), dtype="int64"))

        # check all answers are the same
        tm.assert_frame_equal(shifted[0], shifted[1])
        tm.assert_frame_equal(shifted[0], shifted[2])

    def test_shift_axis1_multiple_blocks(self, using_array_manager):
        # GH#35488
        df1 = DataFrame(np.random.randint(1000, size=(5, 3)))
        df2 = DataFrame(np.random.randint(1000, size=(5, 2)))
        df3 = pd.concat([df1, df2], axis=1)
        if not using_array_manager:
            assert len(df3._mgr.blocks) == 2

        result = df3.shift(2, axis=1)

        expected = df3.take([-1, -1, 0, 1, 2], axis=1)
        expected.iloc[:, :2] = np.nan
        expected.columns = df3.columns

        tm.assert_frame_equal(result, expected)

        # Case with periods < 0
        # rebuild df3 because `take` call above consolidated
        df3 = pd.concat([df1, df2], axis=1)
        if not using_array_manager:
            assert len(df3._mgr.blocks) == 2
        result = df3.shift(-2, axis=1)

        expected = df3.take([2, 3, 4, -1, -1], axis=1)
        expected.iloc[:, -2:] = np.nan
        expected.columns = df3.columns

        tm.assert_frame_equal(result, expected)

    @td.skip_array_manager_not_yet_implemented  # TODO(ArrayManager) axis=1 support
    def test_shift_axis1_multiple_blocks_with_int_fill(self):
        # GH#42719
        df1 = DataFrame(np.random.randint(1000, size=(5, 3)))
        df2 = DataFrame(np.random.randint(1000, size=(5, 2)))
        df3 = pd.concat([df1.iloc[:4, 1:3], df2.iloc[:4, :]], axis=1)
        result = df3.shift(2, axis=1, fill_value=np.int_(0))
        assert len(df3._mgr.blocks) == 2

        expected = df3.take([-1, -1, 0, 1], axis=1)
        expected.iloc[:, :2] = np.int_(0)
        expected.columns = df3.columns

        tm.assert_frame_equal(result, expected)

        # Case with periods < 0
        df3 = pd.concat([df1.iloc[:4, 1:3], df2.iloc[:4, :]], axis=1)
        result = df3.shift(-2, axis=1, fill_value=np.int_(0))
        assert len(df3._mgr.blocks) == 2

        expected = df3.take([2, 3, -1, -1], axis=1)
        expected.iloc[:, -2:] = np.int_(0)
        expected.columns = df3.columns

        tm.assert_frame_equal(result, expected)

    @pytest.mark.filterwarnings("ignore:tshift is deprecated:FutureWarning")
    def test_tshift(self, datetime_frame, frame_or_series):
        # TODO(2.0): remove this test when tshift deprecation is enforced

        # PeriodIndex
        ps = tm.makePeriodFrame()
        ps = tm.get_obj(ps, frame_or_series)
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        tm.assert_equal(unshifted, ps)

        shifted2 = ps.tshift(freq="B")
        tm.assert_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=offsets.BDay())
        tm.assert_equal(shifted, shifted3)

        msg = "Given freq M does not match PeriodIndex freq B"
        with pytest.raises(ValueError, match=msg):
            ps.tshift(freq="M")

        # DatetimeIndex
        dtobj = tm.get_obj(datetime_frame, frame_or_series)
        shifted = dtobj.tshift(1)
        unshifted = shifted.tshift(-1)

        tm.assert_equal(dtobj, unshifted)

        shifted2 = dtobj.tshift(freq=dtobj.index.freq)
        tm.assert_equal(shifted, shifted2)

        inferred_ts = DataFrame(
            datetime_frame.values,
            Index(np.asarray(datetime_frame.index)),
            columns=datetime_frame.columns,
        )
        inferred_ts = tm.get_obj(inferred_ts, frame_or_series)
        shifted = inferred_ts.tshift(1)

        expected = dtobj.tshift(1)
        expected.index = expected.index._with_freq(None)
        tm.assert_equal(shifted, expected)

        unshifted = shifted.tshift(-1)
        tm.assert_equal(unshifted, inferred_ts)

        no_freq = dtobj.iloc[[0, 5, 7]]
        msg = "Freq was not set in the index hence cannot be inferred"
        with pytest.raises(ValueError, match=msg):
            no_freq.tshift()

    def test_tshift_deprecated(self, datetime_frame, frame_or_series):
        # GH#11631
        dtobj = tm.get_obj(datetime_frame, frame_or_series)
        with tm.assert_produces_warning(FutureWarning):
            dtobj.tshift()

    def test_period_index_frame_shift_with_freq(self, frame_or_series):
        ps = tm.makePeriodFrame()
        ps = tm.get_obj(ps, frame_or_series)

        shifted = ps.shift(1, freq="infer")
        unshifted = shifted.shift(-1, freq="infer")
        tm.assert_equal(unshifted, ps)

        shifted2 = ps.shift(freq="B")
        tm.assert_equal(shifted, shifted2)

        shifted3 = ps.shift(freq=offsets.BDay())
        tm.assert_equal(shifted, shifted3)

    def test_datetime_frame_shift_with_freq(self, datetime_frame, frame_or_series):
        dtobj = tm.get_obj(datetime_frame, frame_or_series)
        shifted = dtobj.shift(1, freq="infer")
        unshifted = shifted.shift(-1, freq="infer")
        tm.assert_equal(dtobj, unshifted)

        shifted2 = dtobj.shift(freq=dtobj.index.freq)
        tm.assert_equal(shifted, shifted2)

        inferred_ts = DataFrame(
            datetime_frame.values,
            Index(np.asarray(datetime_frame.index)),
            columns=datetime_frame.columns,
        )
        inferred_ts = tm.get_obj(inferred_ts, frame_or_series)
        shifted = inferred_ts.shift(1, freq="infer")
        expected = dtobj.shift(1, freq="infer")
        expected.index = expected.index._with_freq(None)
        tm.assert_equal(shifted, expected)

        unshifted = shifted.shift(-1, freq="infer")
        tm.assert_equal(unshifted, inferred_ts)

    def test_period_index_frame_shift_with_freq_error(self, frame_or_series):
        ps = tm.makePeriodFrame()
        ps = tm.get_obj(ps, frame_or_series)
        msg = "Given freq M does not match PeriodIndex freq B"
        with pytest.raises(ValueError, match=msg):
            ps.shift(freq="M")

    def test_datetime_frame_shift_with_freq_error(
        self, datetime_frame, frame_or_series
    ):
        dtobj = tm.get_obj(datetime_frame, frame_or_series)
        no_freq = dtobj.iloc[[0, 5, 7]]
        msg = "Freq was not set in the index hence cannot be inferred"
        with pytest.raises(ValueError, match=msg):
            no_freq.shift(freq="infer")

    @td.skip_array_manager_not_yet_implemented  # TODO(ArrayManager) axis=1 support
    def test_shift_dt64values_int_fill_deprecated(self):
        # GH#31971
        ser = Series([pd.Timestamp("2020-01-01"), pd.Timestamp("2020-01-02")])

        with tm.assert_produces_warning(FutureWarning):
            result = ser.shift(1, fill_value=0)
        expected = Series([pd.Timestamp(0), ser[0]])
        tm.assert_series_equal(result, expected)

        df = ser.to_frame()
        with tm.assert_produces_warning(FutureWarning):
            result = df.shift(1, fill_value=0)
        expected = expected.to_frame()
        tm.assert_frame_equal(result, expected)

        # axis = 1
        df2 = DataFrame({"A": ser, "B": ser})
        df2._consolidate_inplace()

        with tm.assert_produces_warning(FutureWarning):
            result = df2.shift(1, axis=1, fill_value=0)

        expected = DataFrame({"A": [pd.Timestamp(0), pd.Timestamp(0)], "B": df2["A"]})
        tm.assert_frame_equal(result, expected)

        # same thing but not consolidated
        # This isn't great that we get different behavior, but
        #  that will go away when the deprecation is enforced
        df3 = DataFrame({"A": ser})
        df3["B"] = ser
        assert len(df3._mgr.arrays) == 2
        result = df3.shift(1, axis=1, fill_value=0)
        expected = DataFrame({"A": [0, 0], "B": df2["A"]})
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "as_cat",
        [
            pytest.param(
                True,
                marks=pytest.mark.xfail(
                    reason="_can_hold_element incorrectly always returns True"
                ),
            ),
            False,
        ],
    )
    @pytest.mark.parametrize(
        "vals",
        [
            date_range("2020-01-01", periods=2),
            date_range("2020-01-01", periods=2, tz="US/Pacific"),
            pd.period_range("2020-01-01", periods=2, freq="D"),
            pd.timedelta_range("2020 Days", periods=2, freq="D"),
            pd.interval_range(0, 3, periods=2),
            pytest.param(
                pd.array([1, 2], dtype="Int64"),
                marks=pytest.mark.xfail(
                    reason="_can_hold_element incorrectly always returns True"
                ),
            ),
            pytest.param(
                pd.array([1, 2], dtype="Float32"),
                marks=pytest.mark.xfail(
                    reason="_can_hold_element incorrectly always returns True"
                ),
            ),
        ],
        ids=lambda x: str(x.dtype),
    )
    # TODO(2.0): remove filtering
    @pytest.mark.filterwarnings("ignore:Index.ravel.*:FutureWarning")
    def test_shift_dt64values_axis1_invalid_fill(
        self, vals, as_cat, using_array_manager, request
    ):
        # GH#44564
        if using_array_manager:
            mark = pytest.mark.xfail(raises=NotImplementedError)
            request.node.add_marker(mark)

        ser = Series(vals)
        if as_cat:
            ser = ser.astype("category")

        df = DataFrame({"A": ser})
        result = df.shift(-1, axis=1, fill_value="foo")
        expected = DataFrame({"A": ["foo", "foo"]})
        tm.assert_frame_equal(result, expected)

        # same thing but multiple blocks
        df2 = DataFrame({"A": ser, "B": ser})
        df2._consolidate_inplace()

        result = df2.shift(-1, axis=1, fill_value="foo")
        expected = DataFrame({"A": df2["B"], "B": ["foo", "foo"]})
        tm.assert_frame_equal(result, expected)

        # same thing but not consolidated
        df3 = DataFrame({"A": ser})
        df3["B"] = ser
        assert len(df3._mgr.arrays) == 2
        result = df3.shift(-1, axis=1, fill_value="foo")
        tm.assert_frame_equal(result, expected)

    def test_shift_axis1_categorical_columns(self):
        # GH#38434
        ci = CategoricalIndex(["a", "b", "c"])
        df = DataFrame(
            {"a": [1, 3], "b": [2, 4], "c": [5, 6]}, index=ci[:-1], columns=ci
        )
        result = df.shift(axis=1)

        expected = DataFrame(
            {"a": [np.nan, np.nan], "b": [1, 3], "c": [2, 4]}, index=ci[:-1], columns=ci
        )
        tm.assert_frame_equal(result, expected)

        # periods != 1
        result = df.shift(2, axis=1)
        expected = DataFrame(
            {"a": [np.nan, np.nan], "b": [np.nan, np.nan], "c": [1, 3]},
            index=ci[:-1],
            columns=ci,
        )
        tm.assert_frame_equal(result, expected)
