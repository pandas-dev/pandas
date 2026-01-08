import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    date_range,
)
import pandas._testing as tm


class TestDataFrameUpdate:
    def test_update_nan(self):
        # #15593 #15617
        # test 1
        df1 = DataFrame({"A": [1.0, 2, 3], "B": date_range("2000", periods=3)})
        df2 = DataFrame({"A": [None, 2, 3]})
        expected = df1.copy()
        df1.update(df2, overwrite=False)

        tm.assert_frame_equal(df1, expected)

        # test 2
        df1 = DataFrame({"A": [1.0, None, 3], "B": date_range("2000", periods=3)})
        df2 = DataFrame({"A": [None, 2, 3]})
        expected = DataFrame({"A": [1.0, 2, 3], "B": date_range("2000", periods=3)})
        df1.update(df2, overwrite=False)

        tm.assert_frame_equal(df1, expected)

    def test_update(self):
        df = DataFrame(
            [[1.5, np.nan, 3.0], [1.5, np.nan, 3.0], [1.5, np.nan, 3], [1.5, np.nan, 3]]
        )

        other = DataFrame([[3.6, 2.0, np.nan], [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other)

        expected = DataFrame(
            [[1.5, np.nan, 3], [3.6, 2, 3], [1.5, np.nan, 3], [1.5, np.nan, 7.0]]
        )
        tm.assert_frame_equal(df, expected)

    def test_update_dtypes(self):
        # gh 3016
        df = DataFrame(
            [[1.0, 2.0, False, True], [4.0, 5.0, True, False]],
            columns=["A", "B", "bool1", "bool2"],
        )

        other = DataFrame([[45, 45]], index=[0], columns=["A", "B"])
        df.update(other)

        expected = DataFrame(
            [[45.0, 45.0, False, True], [4.0, 5.0, True, False]],
            columns=["A", "B", "bool1", "bool2"],
        )
        tm.assert_frame_equal(df, expected)

    def test_update_nooverwrite(self):
        df = DataFrame(
            [[1.5, np.nan, 3.0], [1.5, np.nan, 3.0], [1.5, np.nan, 3], [1.5, np.nan, 3]]
        )

        other = DataFrame([[3.6, 2.0, np.nan], [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, overwrite=False)

        expected = DataFrame(
            [[1.5, np.nan, 3], [1.5, 2, 3], [1.5, np.nan, 3], [1.5, np.nan, 3.0]]
        )
        tm.assert_frame_equal(df, expected)

    def test_update_filtered(self):
        df = DataFrame(
            [[1.5, np.nan, 3.0], [1.5, np.nan, 3.0], [1.5, np.nan, 3], [1.5, np.nan, 3]]
        )

        other = DataFrame([[3.6, 2.0, np.nan], [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, filter_func=lambda x: x > 2)

        expected = DataFrame(
            [[1.5, np.nan, 3], [1.5, np.nan, 3], [1.5, np.nan, 3], [1.5, np.nan, 7.0]]
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "bad_kwarg, exception, msg",
        [
            # errors must be 'ignore' or 'raise'
            ({"errors": "something"}, ValueError, "The parameter errors must.*"),
            ({"join": "inner"}, NotImplementedError, "Only left join is supported"),
        ],
    )
    def test_update_raise_bad_parameter(self, bad_kwarg, exception, msg):
        df = DataFrame([[1.5, 1, 3.0]])
        with pytest.raises(exception, match=msg):
            df.update(df, **bad_kwarg)

    def test_update_raise_on_overlap(self):
        df = DataFrame(
            [[1.5, 1, 3.0], [1.5, np.nan, 3.0], [1.5, np.nan, 3], [1.5, np.nan, 3]]
        )

        other = DataFrame([[2.0, np.nan], [np.nan, 7]], index=[1, 3], columns=[1, 2])
        with pytest.raises(ValueError, match="Data overlaps"):
            df.update(other, errors="raise")

    def test_update_from_non_df(self):
        d = {"a": Series([1, 2, 3, 4]), "b": Series([5, 6, 7, 8])}
        df = DataFrame(d)

        d["a"] = Series([5, 6, 7, 8])
        df.update(d)

        expected = DataFrame(d)

        tm.assert_frame_equal(df, expected)

        d = {"a": [1, 2, 3, 4], "b": [5, 6, 7, 8]}
        df = DataFrame(d)

        d["a"] = [5, 6, 7, 8]
        df.update(d)

        expected = DataFrame(d)

        tm.assert_frame_equal(df, expected)

    def test_update_datetime_tz(self):
        # GH 25807
        result = DataFrame([pd.Timestamp("2019", tz="UTC")])
        with tm.assert_produces_warning(None):
            result.update(result)
        expected = DataFrame([pd.Timestamp("2019", tz="UTC")])
        tm.assert_frame_equal(result, expected)

    def test_update_datetime_tz_in_place(self):
        # https://github.com/pandas-dev/pandas/issues/56227
        result = DataFrame([pd.Timestamp("2019", tz="UTC")])
        orig = result.copy()
        view = result[:]
        result.update(result + pd.Timedelta(days=1))
        expected = DataFrame([pd.Timestamp("2019-01-02", tz="UTC")])
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(view, orig)

    def test_update_with_different_dtype(self):
        # GH#3217
        df = DataFrame({"a": [1, 3], "b": [np.nan, 2]})
        df["c"] = np.nan
        with pytest.raises(TypeError, match="Invalid value"):
            df.update({"c": Series(["foo"], index=[0])})

    @pytest.mark.parametrize("dtype", ["str", object])
    def test_update_modify_view(self, dtype):
        # GH#47188
        df = DataFrame({"A": ["1", np.nan], "B": ["100", np.nan]}, dtype=dtype)
        df2 = DataFrame({"A": ["a", "x"], "B": ["100", "200"]}, dtype=dtype)
        df2_orig = df2.copy()
        result_view = df2[:]
        df2.update(df)
        expected = DataFrame({"A": ["1", "x"], "B": ["100", "200"]}, dtype=dtype)
        tm.assert_frame_equal(df2, expected)
        tm.assert_frame_equal(result_view, df2_orig)

    def test_update_dt_column_with_NaT_create_column(self):
        # GH#16713
        df = DataFrame({"A": [1, None], "B": [pd.NaT, pd.to_datetime("2016-01-01")]})
        df2 = DataFrame({"A": [2, 3]})
        df.update(df2, overwrite=False)
        expected = DataFrame(
            {"A": [1.0, 3.0], "B": [pd.NaT, pd.to_datetime("2016-01-01")]}
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "value_df, value_other, dtype",
        [
            (True, False, bool),
            (1, 2, int),
            (1.0, 2.0, float),
            (1.0 + 1j, 2.0 + 2j, complex),
            (np.uint64(1), np.uint(2), np.dtype("ubyte")),
            (np.uint64(1), np.uint(2), np.dtype("intc")),
            ("a", "b", pd.StringDtype()),
            (
                pd.to_timedelta("1 ms"),
                pd.to_timedelta("2 ms"),
                np.dtype("timedelta64[ns]"),
            ),
            (
                np.datetime64("2000-01-01T00:00:00"),
                np.datetime64("2000-01-02T00:00:00"),
                np.dtype("datetime64[ns]"),
            ),
            (1, 2, pd.Int64Dtype()),
        ],
    )
    def test_update_preserve_dtype(self, value_df, value_other, dtype):
        # GH#55509
        df = DataFrame({"a": [value_df] * 2}, index=[1, 2], dtype=dtype)
        other = DataFrame({"a": [value_other]}, index=[1], dtype=dtype)
        expected = DataFrame({"a": [value_other, value_df]}, index=[1, 2], dtype=dtype)
        df.update(other)
        tm.assert_frame_equal(df, expected)

    def test_update_raises_on_duplicate_argument_index(self):
        # GH#55509
        df = DataFrame({"a": [1, 1]}, index=[1, 2])
        other = DataFrame({"a": [2, 3]}, index=[1, 1])
        with pytest.raises(ValueError, match="duplicate index"):
            df.update(other)

    def test_update_without_intersection(self):
        # GH#63452
        orig = DataFrame({"a": [1]}, index=[1])
        df = orig.copy()
        other = DataFrame({"a": [2]}, index=[2])
        df.update(other)
        tm.assert_frame_equal(df, orig)

    def test_update_on_duplicate_frame_unique_argument_index(self):
        # GH#55509
        df = DataFrame({"a": [1, 1, 1]}, index=[1, 1, 2], dtype=np.dtype("intc"))
        other = DataFrame({"a": [2, 3]}, index=[1, 2], dtype=np.dtype("intc"))
        expected = DataFrame({"a": [2, 2, 3]}, index=[1, 1, 2], dtype=np.dtype("intc"))
        df.update(other)
        tm.assert_frame_equal(df, expected)

    def test_update_preserve_mixed_dtypes(self):
        # GH#44104
        dtype1 = pd.Int64Dtype()
        dtype2 = pd.StringDtype()
        df = DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]})
        df = df.astype({"a": dtype1, "b": dtype2})

        other = DataFrame({"a": [4, 5], "b": ["a", "b"]})
        other = other.astype({"a": dtype1, "b": dtype2})

        expected = DataFrame({"a": [4, 5, 3], "b": ["a", "b", "z"]})
        expected = expected.astype({"a": dtype1, "b": dtype2})

        df.update(other)
        tm.assert_frame_equal(df, expected)
