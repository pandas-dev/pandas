import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Period, Series, period_range
import pandas._testing as tm


class TestPeriodIndex:
    def test_slice_with_negative_step(self):
        ts = Series(np.arange(20), period_range("2014-01", periods=20, freq="M"))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            tm.assert_series_equal(ts[l_slc], ts.iloc[i_slc])
            tm.assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])
            tm.assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])

        assert_slices_equivalent(SLC[Period("2014-10") :: -1], SLC[9::-1])
        assert_slices_equivalent(SLC["2014-10"::-1], SLC[9::-1])

        assert_slices_equivalent(SLC[: Period("2014-10") : -1], SLC[:8:-1])
        assert_slices_equivalent(SLC[:"2014-10":-1], SLC[:8:-1])

        assert_slices_equivalent(SLC["2015-02":"2014-10":-1], SLC[13:8:-1])
        assert_slices_equivalent(
            SLC[Period("2015-02") : Period("2014-10") : -1], SLC[13:8:-1]
        )
        assert_slices_equivalent(SLC["2015-02" : Period("2014-10") : -1], SLC[13:8:-1])
        assert_slices_equivalent(SLC[Period("2015-02") : "2014-10" : -1], SLC[13:8:-1])

        assert_slices_equivalent(SLC["2014-10":"2015-02":-1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        ts = Series(np.arange(20), period_range("2014-01", periods=20, freq="M"))
        with pytest.raises(ValueError, match="slice step cannot be zero"):
            ts[::0]
        with pytest.raises(ValueError, match="slice step cannot be zero"):
            ts.loc[::0]
        with pytest.raises(ValueError, match="slice step cannot be zero"):
            ts.loc[::0]

    def test_slice_keep_name(self):
        idx = period_range("20010101", periods=10, freq="D", name="bob")
        assert idx.name == idx[1:].name

    def test_pindex_slice_index(self):
        pi = period_range(start="1/1/10", end="12/31/12", freq="M")
        s = Series(np.random.rand(len(pi)), index=pi)
        res = s["2010"]
        exp = s[0:12]
        tm.assert_series_equal(res, exp)
        res = s["2011"]
        exp = s[12:24]
        tm.assert_series_equal(res, exp)

    def test_range_slice_day(self):
        # GH#6716
        didx = pd.date_range(start="2013/01/01", freq="D", periods=400)
        pidx = period_range(start="2013/01/01", freq="D", periods=400)

        for idx in [didx, pidx]:
            # slices against index should raise IndexError
            values = [
                "2014",
                "2013/02",
                "2013/01/02",
                "2013/02/01 9H",
                "2013/02/01 09:00",
            ]
            for v in values:
                with pytest.raises(TypeError):
                    idx[v:]

            s = Series(np.random.rand(len(idx)), index=idx)

            tm.assert_series_equal(s["2013/01/02":], s[1:])
            tm.assert_series_equal(s["2013/01/02":"2013/01/05"], s[1:5])
            tm.assert_series_equal(s["2013/02":], s[31:])
            tm.assert_series_equal(s["2014":], s[365:])

            invalid = ["2013/02/01 9H", "2013/02/01 09:00"]
            for v in invalid:
                with pytest.raises(TypeError):
                    idx[v:]

    def test_range_slice_seconds(self):
        # GH#6716
        didx = pd.date_range(start="2013/01/01 09:00:00", freq="S", periods=4000)
        pidx = period_range(start="2013/01/01 09:00:00", freq="S", periods=4000)

        for idx in [didx, pidx]:
            # slices against index should raise IndexError
            values = [
                "2014",
                "2013/02",
                "2013/01/02",
                "2013/02/01 9H",
                "2013/02/01 09:00",
            ]
            for v in values:
                with pytest.raises(TypeError):
                    idx[v:]

            s = Series(np.random.rand(len(idx)), index=idx)

            tm.assert_series_equal(s["2013/01/01 09:05":"2013/01/01 09:10"], s[300:660])
            tm.assert_series_equal(
                s["2013/01/01 10:00":"2013/01/01 10:05"], s[3600:3960]
            )
            tm.assert_series_equal(s["2013/01/01 10H":], s[3600:])
            tm.assert_series_equal(s[:"2013/01/01 09:30"], s[:1860])
            for d in ["2013/01/01", "2013/01", "2013"]:
                tm.assert_series_equal(s[d:], s)

    def test_range_slice_outofbounds(self):
        # GH#5407
        didx = pd.date_range(start="2013/10/01", freq="D", periods=10)
        pidx = period_range(start="2013/10/01", freq="D", periods=10)

        for idx in [didx, pidx]:
            df = DataFrame(dict(units=[100 + i for i in range(10)]), index=idx)
            empty = DataFrame(index=type(idx)([], freq="D"), columns=["units"])
            empty["units"] = empty["units"].astype("int64")

            tm.assert_frame_equal(df["2013/09/01":"2013/09/30"], empty)
            tm.assert_frame_equal(df["2013/09/30":"2013/10/02"], df.iloc[:2])
            tm.assert_frame_equal(df["2013/10/01":"2013/10/02"], df.iloc[:2])
            tm.assert_frame_equal(df["2013/10/02":"2013/09/30"], empty)
            tm.assert_frame_equal(df["2013/10/15":"2013/10/17"], empty)
            tm.assert_frame_equal(df["2013-06":"2013-09"], empty)
            tm.assert_frame_equal(df["2013-11":"2013-12"], empty)

    def test_partial_slice_doesnt_require_monotonicity(self):
        # See also: DatetimeIndex test ofm the same name
        dti = pd.date_range("2014-01-01", periods=30, freq="30D")
        pi = dti.to_period("D")

        ser_montonic = pd.Series(np.arange(30), index=pi)

        shuffler = list(range(0, 30, 2)) + list(range(1, 31, 2))
        ser = ser_montonic[shuffler]
        nidx = ser.index

        # Manually identified locations of year==2014
        indexer_2014 = np.array(
            [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20], dtype=np.intp
        )
        assert (nidx[indexer_2014].year == 2014).all()
        assert not (nidx[~indexer_2014].year == 2014).any()

        result = nidx.get_loc("2014")
        tm.assert_numpy_array_equal(result, indexer_2014)

        expected = ser[indexer_2014]

        result = nidx.get_value(ser, "2014")
        tm.assert_series_equal(result, expected)

        result = ser.loc["2014"]
        tm.assert_series_equal(result, expected)

        result = ser["2014"]
        tm.assert_series_equal(result, expected)

        # Manually identified locations where ser.index is within Mat 2015
        indexer_may2015 = np.array([23], dtype=np.intp)
        assert nidx[23].year == 2015 and nidx[23].month == 5

        result = nidx.get_loc("May 2015")
        tm.assert_numpy_array_equal(result, indexer_may2015)

        expected = ser[indexer_may2015]

        result = nidx.get_value(ser, "May 2015")
        tm.assert_series_equal(result, expected)

        result = ser.loc["May 2015"]
        tm.assert_series_equal(result, expected)

        result = ser["May 2015"]
        tm.assert_series_equal(result, expected)
