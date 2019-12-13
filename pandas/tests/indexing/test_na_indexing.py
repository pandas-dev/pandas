import pytest

import pandas as pd
import pandas.util.testing as tm


@pytest.mark.parametrize(
    "values, expected, dtype",
    [
        ([1, 2, 3], [1, None], "Int64"),
        (["a", "b", "c"], ["a", pd.NA], "string"),
        ([None, "b", "c"], [pd.NA, pd.NA], "string"),
        ([True, False, None], [True, pd.NA], "boolean"),
        ([None, False, None], [pd.NA, pd.NA], "boolean"),
        pytest.param(
            [pd.Timestamp("2000"), pd.Timestamp("2001"), pd.Timestamp("2002")],
            [pd.Timestamp("2000"), pd.NaT],
            "datetime64[ns]",
            marks=pytest.mark.xfail(reson="TODO. Change DatetimeBlock._slice"),
        ),
        (
            [
                pd.Timestamp("2000", tz="CET"),
                pd.Timestamp("2001", tz="CET"),
                pd.Timestamp("2002", tz="CET"),
            ],
            [pd.Timestamp("2000", tz="CET"), pd.NaT],
            "datetime64[ns, cet]",
        ),
        pytest.param(
            [pd.Timedelta("1H"), pd.Timedelta("2H"), pd.Timedelta("3H")],
            [pd.Timedelta("1H"), pd.NaT],
            "timedelta64[ns]",
            marks=pytest.mark.xfail(reson="TODO. Change TimeDeltaBlock._slice"),
        ),
        (
            [pd.Period("2000"), pd.Period("2001"), pd.Period("2002")],
            [pd.Period("2000"), pd.NaT],
            "period[A-DEC]",
        ),
        (["a", "b", "c"], ["a", None], pd.CategoricalDtype(["a", "b", "c"])),
    ],
)
def test_mask_series(values, expected, dtype):
    s = pd.Series(values, dtype=dtype, name="name")
    orig = pd.Series(values, dtype=dtype, name="name")
    mask = pd.array([True, False, None], dtype="boolean")
    result = s[mask]
    # expected = pd.Series([1, pd.NA], index=[0, 2], dtype="Int64")
    expected = pd.Series(expected, dtype=dtype, name="name", index=[0, 2])
    tm.assert_series_equal(result, expected)
    assert result.dtype == s.dtype
    assert pd.isna(result.iloc[-1])

    tm.assert_equal(s.array[mask], expected.array)

    # ensure no mutation
    tm.assert_series_equal(s, orig)
