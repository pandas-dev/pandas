import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "a": [1.0, 2.5, np.nan],
            "bb": [4, 5, 6],
            "c": pd.Series(["x", None, "z"], dtype=object),
            "d": np.array(["2020-01-01", "NaT", "2020-01-03"], dtype="M8[ns]"),
        }
    )


EXPECTED = """\
Rows: 3
Columns: 4
$ a         <float64> 1.0, 2.5, NaN
$ bb          <int64> 4, 5, 6
$ c          <object> 'x', None, 'z'
$ d  <datetime64[ns]> 2020-01-01 00:00:00, NaT, 2020-01-03 00:00:00
"""


def test_glimpse_prints(df, capsys):
    # GH#51668
    result = df.glimpse()
    assert result is None
    assert capsys.readouterr().out == EXPECTED


def test_glimpse_return_string(df, capsys):
    result = df.glimpse(return_type="string")
    assert result == EXPECTED
    assert capsys.readouterr().out == ""


def test_glimpse_return_self(df, capsys):
    result = df.glimpse(return_type="self")
    assert result is df
    assert capsys.readouterr().out == EXPECTED


def test_glimpse_return_frame(df, capsys):
    result = df.glimpse(return_type="frame", max_items_per_column=2)
    expected = pd.DataFrame(
        {
            "column": ["a", "bb", "c", "d"],
            "dtype": ["float64", "int64", "object", "datetime64[ns]"],
            "values": [
                ["1.0", "2.5"],
                ["4", "5"],
                ["'x'", None],
                ["2020-01-01 00:00:00", None],
            ],
        }
    )
    tm.assert_frame_equal(result, expected)
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize(
    "values, expected",
    [
        (pd.array([1, None], dtype="Int64"), "1, <NA>"),
        (np.array([0.1, np.nan], dtype="float32"), "0.1, NaN"),
        (np.array([1 + 2j, np.nan]), "(1+2j), NaN"),
        (np.array(["a ", "1"], dtype=object), "'a ', '1'"),
        (pd.Categorical(["a", None]), "'a', NaN"),
        (pd.array([True, None], dtype="boolean"), "True, <NA>"),
        (pd.to_timedelta(["1D", None]), "1 days 00:00:00, NaT"),
        (pd.period_range("2020-01", periods=2, freq="M"), "2020-01, 2020-02"),
        (np.array([[1, 2], None], dtype=object), "[1, 2], None"),
        (pd.array(["x", None], dtype=pd.StringDtype(na_value=np.nan)), "'x', NaN"),
        (
            np.array([np.zeros((2, 2)), "a\nb"], dtype=object),
            r"[[0. 0.]\n [0. 0.]], 'a\nb'",
        ),
    ],
)
def test_glimpse_values(values, expected):
    result = pd.DataFrame({"a": values}).glimpse(return_type="string")
    assert result.splitlines()[2].split("> ", 1)[1] == expected


def test_glimpse_arrow_nan_not_missing():
    # a non-missing NaN shows as NaN, distinct from the missing <NA>
    pa = pytest.importorskip("pyarrow")
    arr = pd.array(pa.array([np.nan, 1.0, None]), dtype=pd.ArrowDtype(pa.float64()))
    result = pd.DataFrame({"a": arr}).glimpse(return_type="string")
    assert result.splitlines()[2] == "$ a <double[pyarrow]> NaN, 1.0, <NA>"


def test_glimpse_max_items_per_column(df):
    result = df[["a", "bb"]].glimpse(max_items_per_column=1, return_type="string")
    assert result.splitlines()[2:] == ["$ a  <float64> 1.0", "$ bb   <int64> 4"]

    result = df[["a"]].glimpse(max_items_per_column=0, return_type="string")
    assert result == "Rows: 3\nColumns: 1\n$ a <float64>\n"


def test_glimpse_max_colname_length():
    df = pd.DataFrame({"abcdef": [1], "ab": [2]})
    result = df.glimpse(max_colname_length=4, return_type="string")
    assert result.splitlines()[2:] == ["$ abc… <int64> 1", "$ ab   <int64> 2"]

    # only the printed labels are truncated
    result = df.glimpse(max_colname_length=4, return_type="frame")
    assert result["column"].tolist() == ["abcdef", "ab"]


def test_glimpse_east_asian_width():
    df = pd.DataFrame({"あい": [1], "b": [2]})
    with pd.option_context("display.unicode.east_asian_width", True):
        result = df.glimpse(return_type="string")
    assert result.splitlines()[2:] == ["$ あい <int64> 1", "$ b    <int64> 2"]


def test_glimpse_labels():
    # duplicate, multiline and multi-level labels each get their own line
    df = pd.DataFrame([[1, 2]], columns=["a", "a"])
    result = df.glimpse(return_type="string")
    assert result.splitlines()[2:] == ["$ a <int64> 1", "$ a <int64> 2"]

    df = pd.DataFrame([[1]], columns=["a\nb"])
    result = df.glimpse(return_type="string")
    assert result.splitlines()[2:] == [r"$ a\nb <int64> 1"]

    df = pd.DataFrame([[1, 2]], columns=pd.MultiIndex.from_tuples([("a", 1), ("b", 2)]))
    result = df.glimpse(return_type="string")
    assert result.splitlines()[2:] == ["$ (a, 1) <int64> 1", "$ (b, 2) <int64> 2"]

    result = df.glimpse(return_type="frame")
    assert result["column"].tolist() == [("a", 1), ("b", 2)]

    # labels are not re-inferred
    df = pd.DataFrame([[1, 2]], columns=pd.Index([1, None], dtype=object))
    result = df.glimpse(return_type="frame")
    tm.assert_series_equal(
        result["column"], pd.Series([1, None], name="column", dtype=object)
    )


def test_glimpse_empty():
    assert pd.DataFrame().glimpse(return_type="string") == "Rows: 0\nColumns: 0\n"
    df = pd.DataFrame(columns=["a"], dtype="int64")
    assert df.glimpse(return_type="string") == "Rows: 0\nColumns: 1\n$ a <int64>\n"

    result = df.glimpse(return_type="frame")
    assert result["values"].tolist() == [[]]

    result = pd.DataFrame(index=range(3)).glimpse(return_type="frame")
    assert result.columns.tolist() == ["column", "dtype", "values"]
    assert result.empty
    assert result["values"].dtype == object
    assert result["dtype"].dtype == pd.Series([], dtype="str").dtype


@pytest.mark.parametrize(
    "kwargs, msg",
    [
        ({"return_type": "str"}, "return_type must be one of"),
        ({"max_items_per_column": -1}, "max_items_per_column must be"),
        ({"max_items_per_column": 1.5}, "max_items_per_column must be"),
        ({"max_colname_length": 0}, "max_colname_length must be"),
        ({"max_colname_length": 5.0}, "max_colname_length must be"),
    ],
)
def test_glimpse_invalid(df, kwargs, msg):
    with pytest.raises(ValueError, match=msg):
        df.glimpse(**kwargs)
