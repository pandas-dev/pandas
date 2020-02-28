from itertools import product

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestSeriesConvertDtypes:
    # The answerdict has keys that have 4 tuples, corresponding to the arguments
    # infer_objects, convert_string, convert_integer, convert_boolean
    # This allows all 16 possible combinations to be tested.  Since common
    # combinations expect the same answer, this provides an easy way to list
    # all the possibilities
    @pytest.mark.parametrize(
        "data, maindtype, answerdict",
        [
            (
                [1, 2, 3],
                np.dtype("int32"),
                {
                    ((True, False), (True, False), (True,), (True, False)): "Int32",
                    ((True, False), (True, False), (False,), (True, False)): np.dtype(
                        "int32"
                    ),
                },
            ),
            (
                [1, 2, 3],
                np.dtype("int64"),
                {
                    ((True, False), (True, False), (True,), (True, False)): "Int64",
                    ((True, False), (True, False), (False,), (True, False)): np.dtype(
                        "int64"
                    ),
                },
            ),
            (
                ["x", "y", "z"],
                np.dtype("O"),
                {
                    (
                        (True, False),
                        (True,),
                        (True, False),
                        (True, False),
                    ): pd.StringDtype(),
                    ((True, False), (False,), (True, False), (True, False)): np.dtype(
                        "O"
                    ),
                },
            ),
            (
                [True, False, np.nan],
                np.dtype("O"),
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True,),
                    ): pd.BooleanDtype(),
                    ((True, False), (True, False), (True, False), (False,)): np.dtype(
                        "O"
                    ),
                },
            ),
            (
                ["h", "i", np.nan],
                np.dtype("O"),
                {
                    (
                        (True, False),
                        (True,),
                        (True, False),
                        (True, False),
                    ): pd.StringDtype(),
                    ((True, False), (False,), (True, False), (True, False)): np.dtype(
                        "O"
                    ),
                },
            ),
            (  # GH32117
                ["h", "i", 1],
                np.dtype("O"),
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): np.dtype("O"),
                },
            ),
            (
                [10, np.nan, 20],
                np.dtype("float"),
                {
                    ((True, False), (True, False), (True,), (True, False)): "Int64",
                    ((True, False), (True, False), (False,), (True, False)): np.dtype(
                        "float"
                    ),
                },
            ),
            (
                [np.nan, 100.5, 200],
                np.dtype("float"),
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): np.dtype("float"),
                },
            ),
            (
                [3, 4, 5],
                "Int8",
                {((True, False), (True, False), (True, False), (True, False)): "Int8"},
            ),
            (
                [[1, 2], [3, 4], [5]],
                None,
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): np.dtype("O"),
                },
            ),
            (
                [4, 5, 6],
                np.dtype("uint32"),
                {
                    ((True, False), (True, False), (True,), (True, False)): "UInt32",
                    ((True, False), (True, False), (False,), (True, False)): np.dtype(
                        "uint32"
                    ),
                },
            ),
            (
                [-10, 12, 13],
                np.dtype("i1"),
                {
                    ((True, False), (True, False), (True,), (True, False)): "Int8",
                    ((True, False), (True, False), (False,), (True, False)): np.dtype(
                        "i1"
                    ),
                },
            ),
            (
                [1, 2.0],
                object,
                {
                    ((True,), (True, False), (True,), (True, False)): "Int64",
                    ((True,), (True, False), (False,), (True, False)): np.dtype(
                        "float"
                    ),
                    ((False,), (True, False), (True, False), (True, False)): np.dtype(
                        "object"
                    ),
                },
            ),
            (
                [1, 2.5],
                object,
                {
                    ((True,), (True, False), (True, False), (True, False)): np.dtype(
                        "float"
                    ),
                    ((False,), (True, False), (True, False), (True, False)): np.dtype(
                        "object"
                    ),
                },
            ),
            (
                ["a", "b"],
                pd.CategoricalDtype(),
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): pd.CategoricalDtype(),
                },
            ),
            (
                pd.to_datetime(["2020-01-14 10:00", "2020-01-15 11:11"]),
                pd.DatetimeTZDtype(tz="UTC"),
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): pd.DatetimeTZDtype(tz="UTC"),
                },
            ),
            (
                pd.to_datetime(["2020-01-14 10:00", "2020-01-15 11:11"]),
                "datetime64[ns]",
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): np.dtype("datetime64[ns]"),
                },
            ),
            (
                pd.to_datetime(["2020-01-14 10:00", "2020-01-15 11:11"]),
                object,
                {
                    ((True,), (True, False), (True, False), (True, False),): np.dtype(
                        "datetime64[ns]"
                    ),
                    ((False,), (True, False), (True, False), (True, False),): np.dtype(
                        "O"
                    ),
                },
            ),
            (
                pd.period_range("1/1/2011", freq="M", periods=3),
                None,
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): pd.PeriodDtype("M"),
                },
            ),
            (
                pd.arrays.IntervalArray([pd.Interval(0, 1), pd.Interval(1, 5)]),
                None,
                {
                    (
                        (True, False),
                        (True, False),
                        (True, False),
                        (True, False),
                    ): pd.IntervalDtype("int64"),
                },
            ),
        ],
    )
    @pytest.mark.parametrize("params", product(*[(True, False)] * 4))
    def test_convert_dtypes(self, data, maindtype, params, answerdict):
        if maindtype is not None:
            series = pd.Series(data, dtype=maindtype)
        else:
            series = pd.Series(data)
        answers = {k: a for (kk, a) in answerdict.items() for k in product(*kk)}

        ns = series.convert_dtypes(*params)
        expected_dtype = answers[tuple(params)]
        expected = pd.Series(series.values, dtype=expected_dtype)
        tm.assert_series_equal(ns, expected)

        # Test that it is a copy
        copy = series.copy(deep=True)
        ns[ns.notna()] = np.nan

        # Make sure original not changed
        tm.assert_series_equal(series, copy)

    def test_convert_string_dtype(self):
        # https://github.com/pandas-dev/pandas/issues/31731 -> converting columns
        # that are already string dtype
        df = pd.DataFrame(
            {"A": ["a", "b", pd.NA], "B": ["ä", "ö", "ü"]}, dtype="string"
        )
        result = df.convert_dtypes()
        tm.assert_frame_equal(df, result)
