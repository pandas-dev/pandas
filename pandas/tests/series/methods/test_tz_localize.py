import itertools

import pytest
import pytz

from pandas._libs.tslibs import timezones

from pandas import (
    CategoricalDtype,
    DataFrame,
    DatetimeIndex,
    NaT,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm


class TestTZLocalize:
    def test_series_tz_localize_ambiguous_bool(self):
        # make sure that we are correctly accepting bool values as ambiguous

        # GH#14402
        ts = Timestamp("2015-11-01 01:00:03")
        expected0 = Timestamp("2015-11-01 01:00:03-0500", tz="US/Central")
        expected1 = Timestamp("2015-11-01 01:00:03-0600", tz="US/Central")

        ser = Series([ts])
        expected0 = Series([expected0])
        expected1 = Series([expected1])

        with tm.external_error_raised(pytz.AmbiguousTimeError):
            ser.dt.tz_localize("US/Central")

        result = ser.dt.tz_localize("US/Central", ambiguous=True)
        tm.assert_series_equal(result, expected0)

        result = ser.dt.tz_localize("US/Central", ambiguous=[True])
        tm.assert_series_equal(result, expected0)

        result = ser.dt.tz_localize("US/Central", ambiguous=False)
        tm.assert_series_equal(result, expected1)

        result = ser.dt.tz_localize("US/Central", ambiguous=[False])
        tm.assert_series_equal(result, expected1)

    def test_series_tz_localize_matching_index(self):
        # Matching the index of the result with that of the original series
        ts_list = date_range(
            start="2021-01-01T02:00:00", periods=3, freq="1D"
        ).to_list()
        tz_list = date_range(
            start="2021-01-01T02:00:00", periods=3, freq="1D", tz="Europe/Berlin"
        ).to_list()
        df = DataFrame(
            data=list(itertools.islice(itertools.cycle(ts_list), 9)),
            index=[2, 4, 5, 6, 7, 8, 11, 14, 15],
            columns=["pure_datetime"],
        )
        df["categorical_datetime"] = df["pure_datetime"].astype(CategoricalDtype())
        ser = df["categorical_datetime"]
        expected = Series(
            data=list(itertools.islice(itertools.cycle(tz_list), 9)),
            index=[2, 4, 5, 6, 7, 8, 11, 14, 15],
            name="categorical_datetime",
        )
        result = ser.dt.tz_localize("Europe/Berlin")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("tz", ["Europe/Warsaw", "dateutil/Europe/Warsaw"])
    @pytest.mark.parametrize(
        "method, exp",
        [
            ["shift_forward", "2015-03-29 03:00:00"],
            ["NaT", NaT],
            ["raise", None],
            ["foo", "invalid"],
        ],
    )
    def test_series_tz_localize_nonexistent(self, tz, method, exp):
        # GH 8917
        n = 60
        dti = date_range(start="2015-03-29 02:00:00", periods=n, freq="min")
        s = Series(1, dti)
        if method == "raise":
            with tm.external_error_raised(pytz.NonExistentTimeError):
                s.tz_localize(tz, nonexistent=method)
        elif exp == "invalid":
            with pytest.raises(ValueError, match="argument must be one of"):
                dti.tz_localize(tz, nonexistent=method)
        else:
            result = s.tz_localize(tz, nonexistent=method)
            expected = Series(1, index=DatetimeIndex([exp] * n, tz=tz))
            tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_series_tz_localize_empty(self, tzstr):
        # GH#2248
        ser = Series(dtype=object)

        ser2 = ser.tz_localize("utc")
        assert ser2.index.tz == pytz.utc

        ser2 = ser.tz_localize(tzstr)
        timezones.tz_compare(ser2.index.tz, timezones.maybe_get_tz(tzstr))
