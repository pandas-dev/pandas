""" Test Apple numbers read and write """
import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.apple_numbers import read_apple_numbers, to_apple_numbers  # isort:skip


@pytest.mark.single_cpu
class TestFeather:
    def check_error_on_write(self, df, exc, err_msg):
        # check that we are raising the exception
        # on writing

        with pytest.raises(exc, match=err_msg):
            with tm.ensure_clean() as path:
                to_apple_numbers(df, path)

    def check_external_error_on_write(self, df):
        # check that we are raising the exception
        # on writing

        with tm.external_error_raised(Exception):
            with tm.ensure_clean() as path:
                to_apple_numbers(df, path)

    def check_round_trip(self, df, expected=None, write_kwargs={}, **read_kwargs):
        if expected is None:
            expected = df.copy()

        with tm.ensure_clean() as path:
            to_apple_numbers(df, path, **write_kwargs)

            result = read_apple_numbers(path, **read_kwargs)

            tm.assert_frame_equal(result, expected)

    def test_error(self):
        msg = "feather only support IO with DataFrames"
        for obj in [
            pd.Series([1, 2, 3]),
            1,
            "foo",
            pd.Timestamp("20130101"),
            np.array([1, 2, 3]),
        ]:
            self.check_error_on_write(obj, ValueError, msg)

    def test_basic(self):
        df = pd.DataFrame(
            {
                "string": list("abc"),
                "int": list(range(1, 4)),
                "uint": np.arange(3, 6).astype("u1"),
                "float": np.arange(4.0, 7.0, dtype="float64"),
                "float_with_null": [1.0, np.nan, 3],
                "bool": [True, False, True],
                "bool_with_null": [True, np.nan, False],
                "cat": pd.Categorical(list("abc")),
                "dt": pd.DatetimeIndex(
                    list(pd.date_range("20130101", periods=3)), freq=None
                ),
                "dttz": pd.DatetimeIndex(
                    list(pd.date_range("20130101", periods=3, tz="US/Eastern")),
                    freq=None,
                ),
                "dt_with_null": [
                    pd.Timestamp("20130101"),
                    pd.NaT,
                    pd.Timestamp("20130103"),
                ],
                "dtns": pd.DatetimeIndex(
                    list(pd.date_range("20130101", periods=3, freq="ns")), freq=None
                ),
            }
        )
        df["periods"] = pd.period_range("2013", freq="M", periods=3)
        df["timedeltas"] = pd.timedelta_range("1 day", periods=3)
        df["intervals"] = pd.interval_range(0, 3, 3)

        assert df.dttz.dtype.tz.zone == "US/Eastern"

        expected = df.copy()
        expected.loc[1, "bool_with_null"] = None
        self.check_round_trip(df, expected=expected)
