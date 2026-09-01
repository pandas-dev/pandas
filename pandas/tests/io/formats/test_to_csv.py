import contextlib
import csv
from datetime import (
    UTC,
    datetime,
)
import io
import os
import sys
from zipfile import ZipFile

import numpy as np
import pytest

from pandas.compat import is_platform_windows

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import FloatingArray
from pandas.core.indexes.base import get_values_for_csv


@pytest.fixture(params=["python", "pyarrow", "auto"])
def engine(request):
    if request.param in ("pyarrow", "auto"):
        pytest.importorskip("pyarrow")
    return request.param


def check_raises_if_pyarrow(option, engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception for unsupported options. The "auto" engine falls back to the
    python engine (silently) for the same options, so it does not raise.
    """
    if engine == "pyarrow":
        raises_if_pyarrow = pytest.raises(
            ValueError,
            match=f"The {option} option is not supported with the pyarrow engine.",
        )
    else:
        raises_if_pyarrow = contextlib.nullcontext()
    return raises_if_pyarrow


def uses_pyarrow(engine: str) -> bool:
    """
    Whether the given engine actually uses pyarrow to render *compatible*
    data (no bool/timedelta64/timezone-aware/unsupported-option content):
    true for "pyarrow" (forced) and "auto" (pyarrow is available and there
    is nothing for it to fall back from). Used to predict pyarrow's
    cosmetic-only differences from the python engine, e.g. always quoting
    string values/headers, omitting a redundant trailing ".0" on floats,
    or (when no explicit lineterminator is requested) always using "\\n"
    regardless of platform -- see csv_str_for_engine.
    """
    return engine in ("pyarrow", "auto")


def csv_str_for_engine(rows: list[str], pyarrow_used: bool) -> str:
    """
    Like tm.convert_rows_list_to_csv_str, but accounts for pyarrow always
    writing "\\n" as its line terminator regardless of platform, unlike
    the python engine's platform-native os.linesep default. Only accurate
    for a to_csv call that does not pass an explicit `lineterminator`.
    `pyarrow_used` should usually be `uses_pyarrow(engine)`, except where
    some other condition (in addition to engine) also determines whether
    pyarrow actually ends up being used for this specific call.
    """
    sep = "\n" if pyarrow_used else os.linesep
    return sep.join(rows) + sep


def check_warns_if_pyarrow_renders_differently(engine):
    """
    Returns a context manager that ensures that explicit engine="pyarrow"
    warns (and "python"/"auto" do not) when the data contains a
    bool/timedelta64/timezone-aware-datetime64/whole-number-float column,
    which pyarrow renders differently than the python engine.
    """
    if engine == "pyarrow":
        return tm.assert_produces_warning(
            UserWarning,
            match="renders bool, timedelta64",
            # pyarrow's own Table.from_pandas may itself emit unrelated
            # incidental warnings (e.g. converting a tz-aware column)
            raise_on_extra_warnings=False,
        )
    return tm.assert_produces_warning(None)


def check_raises_if_pyarrow_all_na_row(engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception when a row is entirely composed of missing values, since it
    would otherwise write a blank line that read_csv silently drops.
    """
    if engine == "pyarrow":
        raises_if_pyarrow = pytest.raises(
            ValueError,
            match="cannot represent a row that is entirely missing values",
        )
    else:
        raises_if_pyarrow = contextlib.nullcontext()
    return raises_if_pyarrow


def check_raises_if_pyarrow_lineterminator(lineterminator, engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception for a `lineterminator` other than "\\n". pyarrow always writes
    "\\n" itself, regardless of platform, so pass the *effective*
    lineterminator being used (os.linesep, if the call being wrapped does
    not pass one explicitly) -- not necessarily the platform default.
    """
    if lineterminator == "\n":
        return contextlib.nullcontext()
    return check_raises_if_pyarrow("lineterminator", engine)


def check_raises_if_pyarrow_binary_only(engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception when asked to write to something that is not a binary buffer
    or file path (e.g. an already-open text file object, or sys.stdout).
    """
    if engine == "pyarrow":
        raises_if_pyarrow = pytest.raises(
            ValueError,
            match="The pyarrow engine can only write",
        )
    else:
        raises_if_pyarrow = contextlib.nullcontext()
    return raises_if_pyarrow


def check_raises_if_pyarrow_mi_columns(engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception when the columns are a MultiIndex, which it cannot represent
    with pandas' usual multi-row header.
    """
    if engine == "pyarrow":
        raises_if_pyarrow = pytest.raises(
            ValueError,
            match="does not support a MultiIndex on the columns",
        )
    else:
        raises_if_pyarrow = contextlib.nullcontext()
    return raises_if_pyarrow


def check_raises_if_pyarrow_unsupported_data(engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception when the data cannot be represented by pyarrow at all (e.g.
    Period, Interval, complex, or Sparse dtype columns).
    """
    if engine == "pyarrow":
        raises_if_pyarrow = pytest.raises(
            ValueError,
            match="cannot write one or more columns",
        )
    else:
        raises_if_pyarrow = contextlib.nullcontext()
    return raises_if_pyarrow


class TestToCSV:
    def test_to_csv_with_single_column(self, temp_file, engine):
        # see gh-18676, https://bugs.python.org/issue32255
        #
        # Python's CSV library adds an extraneous '""'
        # before the newline when the NaN-value is in
        # the first row. Otherwise, only the newline
        # character is added. This behavior is inconsistent
        # and was patched in https://bugs.python.org/pull_request4672.
        #
        # A single (unheadered, unindexed) column with a missing value is
        # exactly a row that is entirely missing values, which the pyarrow
        # engine cannot represent (see test_to_csv_pyarrow_all_na_row_raises).
        # The column is also whole-number-float (the only non-null value is
        # 1.0), so explicit engine="pyarrow" also warns before raising.
        raises_if_pyarrow = check_raises_if_pyarrow_all_na_row(engine)

        df1 = pd.DataFrame([None, 1])
        expected1 = """\
""
1.0
"""
        with raises_if_pyarrow, check_warns_if_pyarrow_renders_differently(engine):
            df1.to_csv(temp_file, header=None, index=None, engine=engine)
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected1

        df2 = pd.DataFrame([1, None])
        expected2 = """\
1.0
""
"""
        with raises_if_pyarrow, check_warns_if_pyarrow_renders_differently(engine):
            df2.to_csv(temp_file, header=None, index=None, engine=engine)
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected2

    def test_to_csv_default_encoding(self, temp_file, engine):
        # GH17097
        df = pd.DataFrame({"col": ["AAAAA", "ÄÄÄÄÄ", "ßßßßß", "聞聞聞聞聞"]})

        # the default to_csv encoding is utf-8.
        df.to_csv(temp_file, engine=engine)
        tm.assert_frame_equal(pd.read_csv(temp_file, index_col=0), df)

    @pytest.mark.parametrize("encoding", ["utf-16", "utf-16-le", "utf-16-be"])
    def test_to_csv_utf16_encoding(self, encoding, temp_file, engine):
        # GH#10755
        raises_if_pyarrow = check_raises_if_pyarrow("encoding", engine)
        df = pd.DataFrame({"col": ["abc", "déf", "日本語"]})
        with raises_if_pyarrow:
            df.to_csv(temp_file, encoding=encoding, engine=engine)
            result = pd.read_csv(temp_file, index_col=0, encoding=encoding)
            tm.assert_frame_equal(result, df)

    def test_to_csv_quotechar(self, temp_file, engine):
        df = pd.DataFrame({"col": [1, 2]})
        expected = """\
"","col"
"0","1"
"1","2"
"""

        df.to_csv(temp_file, quoting=1, engine=engine)  # 1=QUOTE_ALL
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected

        expected = """\
$$,$col$
$0$,$1$
$1$,$2$
"""

        if engine == "pyarrow":
            raises_if_pyarrow = pytest.raises(
                ValueError,
                match='The pyarrow engine only supports " as a quotechar.',
            )
        else:
            raises_if_pyarrow = contextlib.nullcontext()
        with raises_if_pyarrow:
            df.to_csv(temp_file, quoting=1, quotechar="$", engine=engine)
        if engine != "pyarrow":
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected

        with pytest.raises(TypeError, match="quotechar"):
            df.to_csv(temp_file, quoting=1, quotechar=None, engine=engine)

    def test_to_csv_doublequote(self, temp_file, engine):
        df = pd.DataFrame({"col": ['a"a', '"bb"']})
        expected = '''\
"","col"
"0","a""a"
"1","""bb"""
'''

        df.to_csv(temp_file, quoting=1, doublequote=True, engine=engine)  # QUOTE_ALL
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected

        # no escapechar set: the python engine (used directly, or via auto's
        # fallback, since pyarrow cannot honor doublequote=False either way)
        # raises trying to escape the embedded quote; pyarrow raises its own
        # ValueError before ever reaching a writer
        if engine == "pyarrow":
            with check_raises_if_pyarrow("doublequote", engine):
                df.to_csv(temp_file, doublequote=False, engine=engine)
        else:
            with pytest.raises(csv.Error, match="escapechar"):
                df.to_csv(temp_file, doublequote=False, engine=engine)

    def test_to_csv_escapechar(self, temp_file, engine):
        raises_if_pyarrow = check_raises_if_pyarrow("escapechar", engine)
        df = pd.DataFrame({"col": ['a"a', '"bb"']})
        expected = """\
"","col"
"0","a\\"a"
"1","\\"bb\\""
"""

        with raises_if_pyarrow:
            df.to_csv(
                temp_file, quoting=1, doublequote=False, escapechar="\\", engine=engine
            )
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected

        df = pd.DataFrame({"col": ["a,a", ",bb,"]})
        expected = """\
,col
0,a\\,a
1,\\,bb\\,
"""

        with raises_if_pyarrow:
            df.to_csv(
                temp_file, quoting=3, escapechar="\\", engine=engine
            )  # QUOTE_NONE
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected

    def test_csv_to_string(self, engine):
        df = pd.DataFrame({"col": [1, 2]})
        if uses_pyarrow(engine):
            # pyarrow always quotes string-typed values, including headers
            expected_rows = ['"","col"', "0,1", "1,2"]
        else:
            expected_rows = [",col", "0,1", "1,2"]
        expected = csv_str_for_engine(expected_rows, uses_pyarrow(engine))
        assert df.to_csv(engine=engine) == expected

    def test_to_csv_decimal(self, engine):
        # see gh-781
        raises_if_pyarrow = check_raises_if_pyarrow("decimal", engine)
        df = pd.DataFrame({"col1": [1], "col2": ["a"], "col3": [10.1]})

        if uses_pyarrow(engine):
            # pyarrow always quotes string-typed values, including headers
            expected_rows = ['"","col1","col2","col3"', '0,1,"a",10.1']
        else:
            expected_rows = [",col1,col2,col3", "0,1,a,10.1"]
        expected_default = csv_str_for_engine(expected_rows, uses_pyarrow(engine))
        assert df.to_csv(engine=engine) == expected_default

        expected_rows = [";col1;col2;col3", "0;1;a;10,1"]
        expected_european_excel = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert (
                df.to_csv(engine=engine, decimal=",", sep=";")
                == expected_european_excel
            )

        expected_rows = [",col1,col2,col3", "0,1,a,10.10"]
        expected_float_format_default = tm.convert_rows_list_to_csv_str(expected_rows)
        with check_raises_if_pyarrow("float_format", engine):
            assert (
                df.to_csv(engine=engine, float_format="%.2f")
                == expected_float_format_default
            )

        expected_rows = [";col1;col2;col3", "0;1;a;10,10"]
        expected_float_format = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert (
                df.to_csv(engine=engine, decimal=",", sep=";", float_format="%.2f")
                == expected_float_format
            )

        # see gh-11553: testing if decimal is taken into account for '0.0'
        df = pd.DataFrame({"a": [0, 1.1], "b": [2.2, 3.3], "c": 1})

        expected_rows = ["a,b,c", "0^0,2^2,1", "1^1,3^3,1"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert df.to_csv(engine=engine, index=False, decimal="^") == expected

        # same but for an index
        with raises_if_pyarrow:
            assert df.set_index("a").to_csv(engine=engine, decimal="^") == expected

        # same for a multi-index
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, decimal="^") == expected
            )

    def test_to_csv_float_format(self, engine):
        # testing if float_format is taken into account for the index
        # GH 11553
        raises_if_pyarrow = check_raises_if_pyarrow("float_format", engine)
        df = pd.DataFrame({"a": [0, 1], "b": [2.2, 3.3], "c": 1})

        expected_rows = ["a,b,c", "0,2.20,1", "1,3.30,1"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert (
                df.set_index("a").to_csv(engine=engine, float_format="%.2f") == expected
            )

        # same for a multi-index
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, float_format="%.2f")
                == expected
            )

    def test_to_csv_na_rep(self, engine):
        # see gh-11553
        #
        # Testing if NaN values are correctly represented in the index.
        raises_if_pyarrow = check_raises_if_pyarrow("na_rep", engine)
        df = pd.DataFrame({"a": [0, np.nan], "b": [0, 1], "c": [2, 3]})
        expected_rows = ["a,b,c", "0.0,0,2", "_,1,3"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            assert df.set_index("a").to_csv(engine=engine, na_rep="_") == expected
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, na_rep="_") == expected
            )

        # now with an index containing only NaNs
        df = pd.DataFrame({"a": np.nan, "b": [0, 1], "c": [2, 3]})
        expected_rows = ["a,b,c", "_,0,2", "_,1,3"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            assert df.set_index("a").to_csv(engine=engine, na_rep="_") == expected
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, na_rep="_") == expected
            )

        # check if na_rep parameter does not break anything when no NaN
        df = pd.DataFrame({"a": 0, "b": [0, 1], "c": [2, 3]})
        expected_rows = ["a,b,c", "0,0,2", "0,1,3"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            assert df.set_index("a").to_csv(engine=engine, na_rep="_") == expected
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, na_rep="_") == expected
            )

        with raises_if_pyarrow:
            csv = pd.Series(["a", pd.NA, "c"]).to_csv(engine=engine, na_rep="ZZZZZ")
            expected = tm.convert_rows_list_to_csv_str([",0", "0,a", "1,ZZZZZ", "2,c"])
            assert expected == csv

    def test_to_csv_na_rep_nullable_string(self, nullable_string_dtype, engine):
        # GH 29975
        # Make sure full na_rep shows up when a dtype is provided
        raises_if_pyarrow = check_raises_if_pyarrow("na_rep", engine)
        expected = tm.convert_rows_list_to_csv_str([",0", "0,a", "1,ZZZZZ", "2,c"])
        with raises_if_pyarrow:
            csv = pd.Series(["a", pd.NA, "c"], dtype=nullable_string_dtype).to_csv(
                engine=engine, na_rep="ZZZZZ"
            )
            assert expected == csv

    def test_to_csv_date_format(self, engine):
        # GH 10209
        raises_if_pyarrow = check_raises_if_pyarrow("date_format", engine)
        df_sec = pd.DataFrame({"A": pd.date_range("20130101", periods=5, freq="s")})
        df_day = pd.DataFrame({"A": pd.date_range("20130101", periods=5, freq="D")})

        if uses_pyarrow(engine):
            # The python engine drops the (all-zero) fractional seconds for
            # this column (GH#62111); pyarrow's writer always shows the
            # column's full declared (here, microsecond) precision.
            expected_rows = [
                '"","A"',
                "0,2013-01-01 00:00:00.000000",
                "1,2013-01-01 00:00:01.000000",
                "2,2013-01-01 00:00:02.000000",
                "3,2013-01-01 00:00:03.000000",
                "4,2013-01-01 00:00:04.000000",
            ]
        else:
            expected_rows = [
                ",A",
                "0,2013-01-01 00:00:00",
                "1,2013-01-01 00:00:01",
                "2,2013-01-01 00:00:02",
                "3,2013-01-01 00:00:03",
                "4,2013-01-01 00:00:04",
            ]
        expected_default_sec = csv_str_for_engine(expected_rows, uses_pyarrow(engine))
        assert df_sec.to_csv(engine=engine) == expected_default_sec

        expected_rows = [
            ",A",
            "0,2013-01-01 00:00:00",
            "1,2013-01-02 00:00:00",
            "2,2013-01-03 00:00:00",
            "3,2013-01-04 00:00:00",
            "4,2013-01-05 00:00:00",
        ]
        expected_ymdhms_day = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert (
                df_day.to_csv(date_format="%Y-%m-%d %H:%M:%S", engine=engine)
                == expected_ymdhms_day
            )

        expected_rows = [
            ",A",
            "0,2013-01-01",
            "1,2013-01-01",
            "2,2013-01-01",
            "3,2013-01-01",
            "4,2013-01-01",
        ]
        expected_ymd_sec = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert (
                df_sec.to_csv(date_format="%Y-%m-%d", engine=engine) == expected_ymd_sec
            )

        expected_short_day = tm.convert_rows_list_to_csv_str(
            [
                ",A",
                "0,2013-01-01",
                "1,2013-01-02",
                "2,2013-01-03",
                "3,2013-01-04",
                "4,2013-01-05",
            ]
        )
        if uses_pyarrow(engine):
            # see the freq="s" case above: the python engine drops the
            # (all-zero) time-of-day here entirely, pyarrow does not.
            expected_default_day = csv_str_for_engine(
                [
                    '"","A"',
                    "0,2013-01-01 00:00:00.000000",
                    "1,2013-01-02 00:00:00.000000",
                    "2,2013-01-03 00:00:00.000000",
                    "3,2013-01-04 00:00:00.000000",
                    "4,2013-01-05 00:00:00.000000",
                ],
                uses_pyarrow(engine),
            )
        else:
            expected_default_day = expected_short_day
        # no date_format passed: reflects whichever engine actually runs
        assert df_day.to_csv(engine=engine) == expected_default_day
        # date_format="%Y-%m-%d" passed explicitly: unsupported by pyarrow
        # (raises), and produces the same short-date text on both python
        # and auto (which falls back to python for this call)
        with raises_if_pyarrow:
            assert (
                df_day.to_csv(date_format="%Y-%m-%d", engine=engine)
                == expected_short_day
            )

        # see gh-7791
        #
        # Testing if date_format parameter is taken into account
        # for multi-indexed DataFrames.
        df_sec["B"] = 0
        df_sec["C"] = 1

        expected_rows = ["A,B,C", "2013-01-01,0,1.0"]
        expected_ymd_sec = tm.convert_rows_list_to_csv_str(expected_rows)

        df_sec_grouped = df_sec.groupby([pd.Grouper(key="A", freq="1h"), "B"])
        with raises_if_pyarrow:
            assert (
                df_sec_grouped.mean().to_csv(date_format="%Y-%m-%d", engine=engine)
                == expected_ymd_sec
            )

    def test_to_csv_different_datetime_formats(self, engine):
        # GH#21734
        df = pd.DataFrame(
            {
                "date": pd.to_datetime("1970-01-01"),
                "datetime": pd.date_range("1970-01-01", periods=2, freq="h"),
            }
        )
        if uses_pyarrow(engine):
            # the "date" column is entirely midnight, so the python engine
            # drops the (all-zero) time-of-day (GH#62111); pyarrow does not,
            # and also always quotes the (string) header
            expected_rows = [
                '"date","datetime"',
                "1970-01-01 00:00:00.000000,1970-01-01 00:00:00.000000",
                "1970-01-01 00:00:00.000000,1970-01-01 01:00:00.000000",
            ]
        else:
            expected_rows = [
                "date,datetime",
                "1970-01-01,1970-01-01 00:00:00",
                "1970-01-01,1970-01-01 01:00:00",
            ]
        expected = csv_str_for_engine(expected_rows, uses_pyarrow(engine))
        assert df.to_csv(index=False, engine=engine) == expected

    def test_to_csv_period_columns(self, engine):
        # GH#55426 - exercise the column path for PeriodArray
        # Period dtype can't be represented by pyarrow at all
        raises_if_pyarrow = check_raises_if_pyarrow_unsupported_data(engine)
        df = pd.DataFrame(
            {
                "a": pd.period_range("2020-01-01", periods=2, freq="D"),
                "b": pd.period_range("2020-03-01", periods=2, freq="D"),
            }
        )
        expected_rows = [
            "a,b",
            "2020-01-01,2020-03-01",
            "2020-01-02,2020-03-02",
        ]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert df.to_csv(index=False, engine=engine) == expected

    def test_to_csv_period_columns_date_format(self, engine):
        # GH#51621 - date_format is honored for a PeriodArray column, not just
        # a PeriodIndex
        # date_format is unsupported by pyarrow regardless of the Period
        # dtype also being unsupported
        raises_if_pyarrow = check_raises_if_pyarrow("date_format", engine)
        df = pd.DataFrame({"a": pd.period_range("2000-01-01", periods=2, freq="h")})
        expected_rows = [
            ",a",
            "0,2000-01-01___00:00:00",
            "1,2000-01-01___01:00:00",
        ]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert (
                df.to_csv(date_format="%Y-%m-%d___%H:%M:%S", engine=engine) == expected
            )

    def test_to_csv_interval_columns(self, engine):
        # GH#55426 - exercise the column path for IntervalArray
        # Interval dtype can't be represented by pyarrow at all
        raises_if_pyarrow = check_raises_if_pyarrow_unsupported_data(engine)
        df = pd.DataFrame(
            {
                "a": pd.arrays.IntervalArray.from_breaks([0, 1, 2]),
                "b": pd.arrays.IntervalArray.from_breaks([3, 4, 5]),
            }
        )
        expected_rows = [
            "a,b",
            '"(0, 1]","(3, 4]"',
            '"(1, 2]","(4, 5]"',
        ]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert df.to_csv(index=False, engine=engine) == expected

    def test_to_csv_date_format_in_categorical(self, engine):
        # GH#40754
        raises_if_pyarrow = check_raises_if_pyarrow("date_format", engine)
        # a single unindexed column with a missing value in it is a row
        # that is entirely missing values, which pyarrow cannot represent
        raises_if_pyarrow_na_row = check_raises_if_pyarrow_all_na_row(engine)
        ser = pd.Series(pd.to_datetime(["2021-03-27", pd.NaT], format="%Y-%m-%d"))
        ser = ser.astype("category")
        expected = tm.convert_rows_list_to_csv_str(["0", "2021-03-27", '""'])
        with raises_if_pyarrow_na_row:
            assert ser.to_csv(index=False, engine=engine) == expected

        ser = pd.Series(
            pd.date_range(
                start="2021-03-27", freq="D", periods=1, tz="Europe/Berlin"
            ).append(pd.DatetimeIndex([pd.NaT]))
        )
        ser = ser.astype("category")
        with raises_if_pyarrow:
            assert (
                ser.to_csv(index=False, engine=engine, date_format="%Y-%m-%d")
                == expected
            )

    def test_to_csv_float_ea_float_format(self, engine):
        # GH#45991
        raises_if_pyarrow = check_raises_if_pyarrow("float_format", engine)
        df = pd.DataFrame({"a": [1.1, 2.02, pd.NA, 6.000006], "b": "c"})
        df["a"] = df["a"].astype("Float64")
        with raises_if_pyarrow:
            result = df.to_csv(index=False, engine=engine, float_format="%.5f")
            expected = tm.convert_rows_list_to_csv_str(
                ["a,b", "1.10000,c", "2.02000,c", ",c", "6.00001,c"]
            )
            assert result == expected

    def test_to_csv_float_ea_no_float_format(self, engine):
        # GH#45991
        df = pd.DataFrame({"a": [1.1, 2.02, pd.NA, 6.000006], "b": "c"})
        df["a"] = df["a"].astype("Float64")
        result = df.to_csv(index=False, engine=engine)
        if uses_pyarrow(engine):
            # pyarrow always quotes string-typed values, including headers
            expected = csv_str_for_engine(
                ['"a","b"', '1.1,"c"', '2.02,"c"', ',"c"', '6.000006,"c"'],
                uses_pyarrow(engine),
            )
        else:
            expected = csv_str_for_engine(
                ["a,b", "1.1,c", "2.02,c", ",c", "6.000006,c"], uses_pyarrow(engine)
            )
        assert result == expected

    def test_to_csv_float_ea_nan_distinguish(self, using_nan_is_na, engine):
        # GH#61617, GH#65227 - to_csv should not crash when FloatingArray
        # contains unmasked NaN (with distinguish_nan_and_na=True)
        df = pd.DataFrame(
            {"a": pd.array([np.nan, pd.NA, 3.0], dtype="Float64"), "b": "c"}
        )
        # pyarrow always quotes string-typed values (incl. headers), and
        # drops the redundant trailing ".0" on whole-number floats; "auto"
        # only actually uses pyarrow here when the unmasked NaN is visible
        # (using_nan_is_na=False), since an unmasked NaN otherwise makes
        # column "a" look like a whole-number-only float column, which
        # "auto" avoids rendering with pyarrow (and explicit engine="pyarrow"
        # warns about instead, only when using_nan_is_na=True)
        pyarrow_used = engine == "pyarrow" or (engine == "auto" and not using_nan_is_na)
        if engine == "pyarrow" and using_nan_is_na:
            warns_if_pyarrow = check_warns_if_pyarrow_renders_differently(engine)
        else:
            warns_if_pyarrow = tm.assert_produces_warning(None)
        with warns_if_pyarrow:
            result = df.to_csv(index=False, engine=engine)
        col_b = '"c"' if pyarrow_used else "c"
        val_a = "3" if pyarrow_used else "3.0"
        header = '"a","b"' if pyarrow_used else "a,b"
        if using_nan_is_na:
            expected = csv_str_for_engine(
                [header, f",{col_b}", f",{col_b}", f"{val_a},{col_b}"], pyarrow_used
            )
        else:
            expected = csv_str_for_engine(
                [header, f"nan,{col_b}", f",{col_b}", f"{val_a},{col_b}"], pyarrow_used
            )
        assert result == expected

    def test_to_csv_float_ea_nan_distinguish_series(self, using_nan_is_na, engine):
        # GH#65227 - Series.to_csv with FloatingArray containing both NaN and NA
        ser = pd.Series((1, pd.NA, 0), index=["a", "b", "c"], dtype="Float64", name="x")
        ser = ser / ser
        # pyarrow always quotes string-typed values (incl. headers/index),
        # and drops the redundant trailing ".0" on whole-number floats;
        # see test_to_csv_float_ea_nan_distinguish for why "auto" only
        # actually uses pyarrow when using_nan_is_na=False (and explicit
        # engine="pyarrow" warns about it instead, only when True)
        pyarrow_used = engine == "pyarrow" or (engine == "auto" and not using_nan_is_na)
        if engine == "pyarrow" and using_nan_is_na:
            warns_if_pyarrow = check_warns_if_pyarrow_renders_differently(engine)
        else:
            warns_if_pyarrow = tm.assert_produces_warning(None)
        with warns_if_pyarrow:
            result = ser.to_csv(engine=engine)
        if pyarrow_used:
            header, idx_a, idx_b, idx_c, val_a = (
                '"","x"',
                '"a"',
                '"b"',
                '"c"',
                "1",
            )
        else:
            header, idx_a, idx_b, idx_c, val_a = ",x", "a", "b", "c", "1.0"
        if using_nan_is_na:
            expected = csv_str_for_engine(
                [header, f"{idx_a},{val_a}", f"{idx_b},", f"{idx_c},"], pyarrow_used
            )
        else:
            expected = csv_str_for_engine(
                [header, f"{idx_a},{val_a}", f"{idx_b},", f"{idx_c},nan"], pyarrow_used
            )
        assert result == expected

    def test_to_csv_2d_float_ea(self):
        # GH#64634 - get_values_for_csv on a 2D float ExtensionArray should not fail.
        data = np.array([1.0, np.nan, 3.0, 4.0])
        mask = np.array([False, True, False, False])
        arr2d = FloatingArray(data, mask).reshape(2, 2)

        result = get_values_for_csv(
            arr2d,
            date_format=None,
            na_rep="NA",
            float_format=None,
            decimal=".",
        )
        expected = np.array([["1.0", "NA"], ["3.0", "4.0"]], dtype=object)
        tm.assert_numpy_array_equal(result, expected)

    def test_to_csv_multi_index(self, engine):
        # see gh-6618
        raises_if_pyarrow = check_raises_if_pyarrow_mi_columns(engine)
        df = pd.DataFrame([1], columns=pd.MultiIndex.from_arrays([[1], [2]]))

        exp_rows = [",1", ",2", "0,1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        with raises_if_pyarrow:
            assert df.to_csv(engine=engine) == exp

        exp_rows = ["1", "2", "1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        with raises_if_pyarrow:
            assert df.to_csv(index=False, engine=engine) == exp

        df = pd.DataFrame(
            [1],
            columns=pd.MultiIndex.from_arrays([[1], [2]]),
            index=pd.MultiIndex.from_arrays([[1], [2]]),
        )

        exp_rows = [",,1", ",,2", "1,2,1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        with raises_if_pyarrow:
            assert df.to_csv(engine=engine) == exp

        exp_rows = ["1", "2", "1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        with raises_if_pyarrow:
            assert df.to_csv(index=False, engine=engine) == exp

        df = pd.DataFrame([1], columns=pd.MultiIndex.from_arrays([["foo"], ["bar"]]))

        exp_rows = [",foo", ",bar", "0,1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        with raises_if_pyarrow:
            assert df.to_csv(engine=engine) == exp

        exp_rows = ["foo", "bar", "1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        with raises_if_pyarrow:
            assert df.to_csv(index=False, engine=engine) == exp

    @pytest.mark.parametrize(
        "ind,expected,expected_pyarrow",
        [
            (
                pd.MultiIndex(levels=[[1.0]], codes=[[0]], names=["x"]),
                "x,data\n1.0,1\n",
                '"x","data"\n1,1\n',
            ),
            (
                pd.MultiIndex(
                    levels=[[1.0], [2.0]], codes=[[0], [0]], names=["x", "y"]
                ),
                "x,y,data\n1.0,2.0,1\n",
                '"x","y","data"\n1,2,1\n',
            ),
        ],
    )
    def test_to_csv_single_level_multi_index(
        self, ind, expected, expected_pyarrow, frame_or_series, engine
    ):
        # see gh-19589
        # pyarrow always quotes string-typed values/headers, and drops the
        # redundant trailing ".0" on whole-number floats (here, the index);
        # "auto" falls back to python for this whole-number-float data, so
        # only explicit engine="pyarrow" actually renders this differently
        if engine == "pyarrow":
            expected = expected_pyarrow
        obj = frame_or_series(pd.Series([1], ind, name="data"))
        with check_warns_if_pyarrow_renders_differently(engine):
            result = obj.to_csv(lineterminator="\n", header=True, engine=engine)
        assert result == expected

    def test_to_csv_string_array_ascii(self, temp_file, engine):
        # GH 10813
        raises_if_pyarrow = check_raises_if_pyarrow("encoding", engine)
        str_array = [{"names": ["foo", "bar"]}, {"names": ["baz", "qux"]}]
        df = pd.DataFrame(str_array)
        expected_ascii = """\
,names
0,"['foo', 'bar']"
1,"['baz', 'qux']"
"""
        with raises_if_pyarrow:
            df.to_csv(temp_file, encoding="ascii", engine=engine)
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected_ascii

    def test_to_csv_string_array_utf8(self, temp_file, engine):
        # GH 10813
        # explicitly requesting "utf-8" is compatible with the pyarrow
        # engine (it always writes utf-8), but the "names" column here
        # holds actual Python list objects, which pyarrow infers as a list
        # type that its CSV writer cannot serialize
        if engine == "pyarrow":
            raises_if_pyarrow = pytest.raises(
                ValueError,
                match="cannot write one or more columns",
            )
        else:
            raises_if_pyarrow = contextlib.nullcontext()
        str_array = [{"names": ["foo", "bar"]}, {"names": ["baz", "qux"]}]
        df = pd.DataFrame(str_array)
        expected_utf8 = """\
,names
0,"['foo', 'bar']"
1,"['baz', 'qux']"
"""
        with raises_if_pyarrow:
            df.to_csv(temp_file, encoding="utf-8", engine=engine)
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected_utf8

    def test_to_csv_roundtrip_with_newline_in_field(self, temp_file, engine):
        # GH#22497 - embedded newlines in field values should survive roundtrip
        df = pd.DataFrame({"A": ["test", "te\nst"]})  # codespell:ignore te
        df.to_csv(temp_file, index=False, engine=engine)
        result = pd.read_csv(temp_file)
        tm.assert_frame_equal(result, df)

    def test_to_csv_string_with_lf(self, temp_file, engine):
        # GH 20353
        data = {"int": [1, 2, 3], "str_lf": ["abc", "d\nef", "g\nh\n\ni"]}
        df = pd.DataFrame(data)

        # case 1: The default line terminator (os.linesep for python, but
        # always "\n" for pyarrow, regardless of platform -- PR 21406)
        sep = b"\n" if uses_pyarrow(engine) else os.linesep.encode("utf-8")
        if uses_pyarrow(engine):
            # pyarrow always quotes string-typed values, including headers
            expected_noarg = (
                b'"int","str_lf"'
                + sep
                + b'1,"abc"'
                + sep
                + b'2,"d\nef"'
                + sep
                + b'3,"g\nh\n\ni"'
                + sep
            )
        else:
            expected_noarg = (
                b"int,str_lf"
                + sep
                + b"1,abc"
                + sep
                + b'2,"d\nef"'
                + sep
                + b'3,"g\nh\n\ni"'
                + sep
            )
        df.to_csv(temp_file, index=False, engine=engine)
        with open(temp_file, "rb") as f:
            assert f.read() == expected_noarg

        # case 2: LF as line terminator -- always compatible with pyarrow
        if uses_pyarrow(engine):
            expected_lf = b'"int","str_lf"\n1,"abc"\n2,"d\nef"\n3,"g\nh\n\ni"\n'
        else:
            expected_lf = b'int,str_lf\n1,abc\n2,"d\nef"\n3,"g\nh\n\ni"\n'
        with check_raises_if_pyarrow_lineterminator("\n", engine):
            df.to_csv(temp_file, lineterminator="\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_lf

        # case 3: CRLF as line terminator -- never compatible with pyarrow
        # (it can only ever write "\n"), so this is always the plain
        # python-style rendering, regardless of engine
        # 'lineterminator' should not change inner element
        expected_crlf = b'int,str_lf\r\n1,abc\r\n2,"d\nef"\r\n3,"g\nh\n\ni"\r\n'
        with check_raises_if_pyarrow_lineterminator("\r\n", engine):
            df.to_csv(temp_file, lineterminator="\r\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_crlf

    def test_to_csv_string_with_crlf(self, temp_file, engine):
        # GH 20353
        data = {"int": [1, 2, 3], "str_crlf": ["abc", "d\r\nef", "g\r\nh\r\n\r\ni"]}
        df = pd.DataFrame(data)
        # case 1: The default line terminator (os.linesep for python, but
        # always "\n" for pyarrow, regardless of platform -- PR 21406)
        sep = b"\n" if uses_pyarrow(engine) else os.linesep.encode("utf-8")
        if uses_pyarrow(engine):
            # pyarrow always quotes string-typed values, including headers
            expected_noarg = (
                b'"int","str_crlf"'
                + sep
                + b'1,"abc"'
                + sep
                + b'2,"d\r\nef"'
                + sep
                + b'3,"g\r\nh\r\n\r\ni"'
                + sep
            )
        else:
            expected_noarg = (
                b"int,str_crlf"
                + sep
                + b"1,abc"
                + sep
                + b'2,"d\r\nef"'
                + sep
                + b'3,"g\r\nh\r\n\r\ni"'
                + sep
            )
        df.to_csv(temp_file, index=False, engine=engine)
        with open(temp_file, "rb") as f:
            assert f.read() == expected_noarg

        # case 2: LF as line terminator -- always compatible with pyarrow
        if uses_pyarrow(engine):
            expected_lf = (
                b'"int","str_crlf"\n1,"abc"\n2,"d\r\nef"\n3,"g\r\nh\r\n\r\ni"\n'
            )
        else:
            expected_lf = b'int,str_crlf\n1,abc\n2,"d\r\nef"\n3,"g\r\nh\r\n\r\ni"\n'
        with check_raises_if_pyarrow_lineterminator("\n", engine):
            df.to_csv(temp_file, lineterminator="\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_lf

        # case 3: CRLF as line terminator -- never compatible with pyarrow
        # (it can only ever write "\n"), so this is always the plain
        # python-style rendering, regardless of engine
        # 'lineterminator' should not change inner element
        expected_crlf = (
            b'int,str_crlf\r\n1,abc\r\n2,"d\r\nef"\r\n3,"g\r\nh\r\n\r\ni"\r\n'
        )
        with check_raises_if_pyarrow_lineterminator("\r\n", engine):
            df.to_csv(temp_file, lineterminator="\r\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_crlf

    def test_to_csv_stdout_file(self, capsys, engine):
        # GH 21561
        # sys.stdout is a text stream, which the pyarrow engine cannot
        # write to (nor can it honor a non-default `encoding`)
        raises_if_pyarrow = check_raises_if_pyarrow_binary_only(engine)
        df = pd.DataFrame(
            [["foo", "bar"], ["baz", "qux"]], columns=["name_1", "name_2"]
        )
        expected_rows = [",name_1,name_2", "0,foo,bar", "1,baz,qux"]
        expected_ascii = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            df.to_csv(sys.stdout, encoding="ascii", engine=engine)
            captured = capsys.readouterr()

            assert captured.out == expected_ascii
            assert not sys.stdout.closed

    def test_to_csv_write_to_open_file(self, temp_file, engine, request):
        # GH 21696
        # an already-open text file object is another destination the
        # pyarrow engine cannot write to; explicit engine="pyarrow" raises
        # cleanly before ever touching the (Windows-problematic) write path
        # below, so the Windows xfail only applies to python/auto
        if engine != "pyarrow":
            mark = pytest.mark.xfail(
                is_platform_windows(),
                reason=(
                    "Especially in Windows, file stream should not be passed"
                    "to csv writer without newline='' option."
                    "(https://docs.python.org/3/library/csv.html#csv.writer)"
                ),
            )
            request.applymarker(mark)
        raise_if_pyarrow = check_raises_if_pyarrow_binary_only(engine)
        df = pd.DataFrame({"a": ["x", "y", "z"]})
        expected = """\
manual header
x
y
z
"""
        with open(temp_file, "w", encoding="utf-8") as f:
            f.write("manual header\n")
            with raise_if_pyarrow:
                df.to_csv(f, header=None, index=None, engine=engine)
        if engine != "pyarrow":
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected

    def test_to_csv_write_to_open_file_with_newline_py3(self, temp_file, engine):
        # see gh-21696
        # see gh-20353
        raise_if_pyarrow = check_raises_if_pyarrow_binary_only(engine)
        df = pd.DataFrame({"a": ["x", "y", "z"]})
        expected_rows = ["x", "y", "z"]
        expected = "manual header\n" + tm.convert_rows_list_to_csv_str(expected_rows)

        # TODO: Open in bytes mode for pyarrow
        with open(temp_file, "w", newline="", encoding="utf-8") as f:
            f.write("manual header\n")
            with raise_if_pyarrow:
                df.to_csv(f, header=None, index=None, engine=engine)

        if engine != "pyarrow":
            with open(temp_file, "rb") as f:
                assert f.read() == bytes(expected, "utf-8")

    @pytest.mark.parametrize("to_infer", [True, False])
    @pytest.mark.parametrize("read_infer", [True, False])
    def test_to_csv_compression(
        self,
        compression_only,
        read_infer,
        to_infer,
        compression_to_extension,
        temp_file,
        engine,
    ):
        # see gh-15008
        compression = compression_only

        df = pd.DataFrame({"A": [1]})

        to_compression = "infer" if to_infer else compression
        read_compression = "infer" if read_infer else compression

        path_ext = str(temp_file) + "." + compression_to_extension[compression]
        df.to_csv(path_ext, compression=to_compression, engine=engine)
        result = pd.read_csv(path_ext, index_col=0, compression=read_compression)
        tm.assert_frame_equal(result, df)

    def test_to_csv_compression_dict(self, compression_only, temp_file, engine):
        # GH 26023
        method = compression_only
        df = pd.DataFrame({"ABC": [1]})
        extension = {
            "gzip": "gz",
            "zstd": "zst",
        }.get(method, method)

        path = str(temp_file) + "." + extension
        df.to_csv(path, compression={"method": method}, engine=engine)
        read_df = pd.read_csv(path, index_col=0)
        tm.assert_frame_equal(read_df, df)

    def test_to_csv_compression_dict_no_method_raises(self, temp_file, engine):
        # GH 26023
        df = pd.DataFrame({"ABC": [1]})
        compression = {"some_option": True}
        msg = "must have key 'method'"

        with pytest.raises(ValueError, match=msg):
            df.to_csv(temp_file, compression=compression, engine=engine)

    @pytest.mark.parametrize("compression", ["zip", "infer"])
    @pytest.mark.parametrize("archive_name", ["test_to_csv.csv", "test_to_csv.zip"])
    def test_to_csv_zip_arguments(self, compression, archive_name, temp_file, engine):
        # GH 26023
        df = pd.DataFrame({"ABC": [1]})

        path = str(temp_file) + ".zip"
        df.to_csv(
            path,
            compression={"method": compression, "archive_name": archive_name},
            engine=engine,
        )
        with ZipFile(path) as zp:
            assert len(zp.filelist) == 1
            archived_file = zp.filelist[0].filename
            assert archived_file == archive_name

    @pytest.mark.parametrize(
        "filename,expected_arcname",
        [
            ("archive.csv", "archive.csv"),
            ("archive.tsv", "archive.tsv"),
            ("archive.csv.zip", "archive.csv"),
            ("archive.tsv.zip", "archive.tsv"),
            ("archive.zip", "archive"),
        ],
    )
    def test_to_csv_zip_infer_name(self, tmp_path, filename, expected_arcname, engine):
        # GH 39465
        df = pd.DataFrame({"ABC": [1]})
        path = tmp_path / filename
        df.to_csv(path, compression="zip", engine=engine)
        with ZipFile(path) as zp:
            assert len(zp.filelist) == 1
            archived_file = zp.filelist[0].filename
            assert archived_file == expected_arcname

    @pytest.mark.parametrize("df_new_type", ["Int64"])
    def test_to_csv_na_rep_long_string(self, df_new_type, engine):
        # see gh-25099
        raises_if_pyarrow = check_raises_if_pyarrow("na_rep", engine)
        df = pd.DataFrame({"c": [pd.NA] * 3})
        df = df.astype(df_new_type)
        expected_rows = ["c", "mynull", "mynull", "mynull"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            result = df.to_csv(
                index=False, na_rep="mynull", encoding="ascii", engine=engine
            )

            assert expected == result

    def test_to_csv_timedelta_precision(self, engine):
        # GH 6783
        # writing to a StringIO (a text buffer) is itself incompatible with
        # the pyarrow engine, regardless of the timedelta64 dtype
        raises_if_pyarrow = check_raises_if_pyarrow_binary_only(engine)
        s = pd.Series([1, 1]).astype("timedelta64[ns]")
        buf = io.StringIO()
        with raises_if_pyarrow:
            s.to_csv(buf, engine=engine)
            result = buf.getvalue()
            expected_rows = [
                ",0",
                "0,0 days 00:00:00.000000001",
                "1,0 days 00:00:00.000000001",
            ]
            expected = tm.convert_rows_list_to_csv_str(expected_rows)
            assert result == expected

    def test_na_rep_truncated(self, engine):
        # https://github.com/pandas-dev/pandas/issues/31447
        raises_if_pyarrow = check_raises_if_pyarrow("na_rep", engine)
        with raises_if_pyarrow:
            result = pd.Series(range(8, 12)).to_csv(na_rep="-", engine=engine)
            expected = tm.convert_rows_list_to_csv_str(
                [",0", "0,8", "1,9", "2,10", "3,11"]
            )
            assert result == expected

        with raises_if_pyarrow:
            result = pd.Series([True, False]).to_csv(na_rep="nan", engine=engine)
            expected = tm.convert_rows_list_to_csv_str([",0", "0,True", "1,False"])
            assert result == expected

        with raises_if_pyarrow:
            result = pd.Series([1.1, 2.2]).to_csv(na_rep=".", engine=engine)
            expected = tm.convert_rows_list_to_csv_str([",0", "0,1.1", "1,2.2"])
            assert result == expected

    @pytest.mark.parametrize("errors", ["surrogatepass", "ignore", "replace"])
    def test_to_csv_errors(self, errors, temp_file, engine):
        # GH 22610
        raises_if_pyarrow = check_raises_if_pyarrow("errors", engine)
        data = ["\ud800foo"]
        with raises_if_pyarrow:
            ser = pd.Series(data, index=pd.Index(data, dtype=object), dtype=object)
            ser.to_csv(temp_file, errors=errors, engine=engine)
        # No use in reading back the data as it is not the same anymore
        # due to the error handling

    @pytest.mark.parametrize("mode", ["wb", "w"])
    def test_to_csv_binary_handle(self, mode, temp_file, engine):
        """
        Binary file objects should work (if 'mode' contains a 'b') or even without
        it in most cases.

        GH 35058 and GH 19827
        """
        df = pd.DataFrame(
            1.1 * np.arange(120).reshape((30, 4)),
            columns=pd.Index(list("ABCD")),
            index=pd.Index([f"i-{i}" for i in range(30)]),
        )

        with open(temp_file, mode="w+b") as handle:
            if mode == "w":
                raises_if_pyarrow = check_raises_if_pyarrow_binary_only(engine)
            else:
                raises_if_pyarrow = contextlib.nullcontext()
            with raises_if_pyarrow:
                df.to_csv(handle, mode=mode, engine=engine)
        if not engine == "pyarrow" and mode == "w":
            tm.assert_frame_equal(df, pd.read_csv(temp_file, index_col=0))

    @pytest.mark.parametrize("mode", ["wb", "w"])
    def test_to_csv_encoding_binary_handle(self, mode, temp_file, engine, request):
        """
        Binary file objects should honor a specified encoding.

        GH 23854 and GH 13068 with binary handles
        """

        if mode == "w" and engine == "pyarrow":
            mark = pytest.mark.xfail(
                reason="pyarrow doesn't support non-binary handles."
            )
            request.applymarker(mark)

        raises_if_pyarrow = check_raises_if_pyarrow("encoding", engine)
        # example from GH 23854
        content = "a, b, 🐟".encode("utf-8-sig")
        buffer = io.BytesIO(content)
        df = pd.read_csv(buffer, encoding="utf-8-sig")

        buffer = io.BytesIO()
        with raises_if_pyarrow:
            df.to_csv(
                buffer, mode=mode, encoding="utf-8-sig", index=False, engine=engine
            )
            buffer.seek(0)  # tests whether file handle wasn't closed
            assert buffer.getvalue().startswith(content)

        # example from GH 13068
        with open(temp_file, "w+b") as handle:
            with raises_if_pyarrow:
                pd.DataFrame().to_csv(
                    handle, mode=mode, encoding="utf-8-sig", engine=engine
                )

                handle.seek(0)
                assert handle.read().startswith(b'\xef\xbb\xbf""')


def test_to_csv_iterative_compression_name(compression, temp_file, engine):
    # GH 38714
    df = pd.DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=pd.Index(list("ABCD")),
        index=pd.Index([f"i-{i}" for i in range(30)]),
    )
    df.to_csv(temp_file, compression=compression, chunksize=1, engine=engine)
    tm.assert_frame_equal(
        pd.read_csv(temp_file, compression=compression, index_col=0), df
    )


def test_to_csv_iterative_compression_buffer(compression, engine):
    # GH 38714
    df = pd.DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=pd.Index(list("ABCD")),
        index=pd.Index([f"i-{i}" for i in range(30)]),
    )
    with io.BytesIO() as buffer:
        df.to_csv(buffer, compression=compression, chunksize=1, engine=engine)
        buffer.seek(0)
        tm.assert_frame_equal(
            pd.read_csv(buffer, compression=compression, index_col=0), df
        )
        assert not buffer.closed


def test_new_style_float_format_basic():
    df = pd.DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A\n0,1234.57\n1,9876.54\n"
    assert result == expected


def test_new_style_float_format_thousands():
    df = pd.DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format="{:,.2f}", lineterminator="\n")
    expected = ',A\n0,"1,234.57"\n1,"9,876.54"\n'
    assert result == expected


def test_new_style_scientific_format():
    df = pd.DataFrame({"A": [0.000123, 0.000456]})
    result = df.to_csv(float_format="{:.2e}", lineterminator="\n")
    expected = ",A\n0,1.23e-04\n1,4.56e-04\n"
    assert result == expected


def test_new_style_with_nan():
    df = pd.DataFrame({"A": [1.23, np.nan, 4.56]})
    result = df.to_csv(float_format="{:.2f}", na_rep="NA", lineterminator="\n")
    expected = ",A\n0,1.23\n1,NA\n2,4.56\n"
    assert result == expected


def test_new_style_with_mixed_types():
    df = pd.DataFrame({"A": [1.23, 4.56], "B": ["x", "y"]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A,B\n0,1.23,x\n1,4.56,y\n"
    assert result == expected


def test_new_style_with_mixed_types_in_column():
    df = pd.DataFrame({"A": [1.23, "text", 4.56]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A\n0,1.23\n1,text\n2,4.56\n"
    assert result == expected


def test_invalid_new_style_format_missing_brace():
    df = pd.DataFrame({"A": [1.23]})
    with pytest.raises(ValueError, match="Invalid new-style format string '{:.2f"):
        df.to_csv(float_format="{:.2f")


def test_invalid_new_style_format_specifier():
    df = pd.DataFrame({"A": [1.23]})
    with pytest.raises(ValueError, match="Invalid new-style format string '{:.2z}'"):
        df.to_csv(float_format="{:.2z}")


def test_old_style_format_compatibility():
    df = pd.DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format="%.2f", lineterminator="\n")
    expected = ",A\n0,1234.57\n1,9876.54\n"
    assert result == expected


def test_callable_float_format_compatibility():
    df = pd.DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format=lambda x: f"{x:,.2f}", lineterminator="\n")
    expected = ',A\n0,"1,234.57"\n1,"9,876.54"\n'
    assert result == expected


def test_no_float_format():
    # engine="python": float_format=None (the default) does not disqualify
    # the pyarrow engine, but this test is about float_format handling, not
    # engine choice
    df = pd.DataFrame({"A": [1.23, 4.56]})
    result = df.to_csv(float_format=None, lineterminator="\n", engine="python")
    expected = ",A\n0,1.23\n1,4.56\n"
    assert result == expected


def test_large_numbers():
    df = pd.DataFrame({"A": [1e308, 2e308]})
    result = df.to_csv(float_format="{:.2e}", lineterminator="\n")
    expected = ",A\n0,1.00e+308\n1,inf\n"
    assert result == expected


def test_zero_and_negative():
    df = pd.DataFrame({"A": [0.0, -1.23456]})
    result = df.to_csv(float_format="{:+.2f}", lineterminator="\n")
    expected = ",A\n0,+0.00\n1,-1.23\n"
    assert result == expected


def test_unicode_format():
    df = pd.DataFrame({"A": [1.23, 4.56]})
    result = df.to_csv(float_format="{:.2f}€", encoding="utf-8", lineterminator="\n")
    expected = ",A\n0,1.23€\n1,4.56€\n"
    assert result == expected


def test_empty_dataframe():
    df = pd.DataFrame({"A": []})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A\n"
    assert result == expected


def test_multi_column_float():
    df = pd.DataFrame({"A": [1.23, 4.56], "B": [7.89, 0.12]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A,B\n0,1.23,7.89\n1,4.56,0.12\n"
    assert result == expected


def test_invalid_float_format_type():
    df = pd.DataFrame({"A": [1.23]})
    with pytest.raises(ValueError, match="float_format must be a string or callable"):
        df.to_csv(float_format=123)


def test_new_style_with_inf():
    df = pd.DataFrame({"A": [1.23, np.inf, -np.inf]})
    result = df.to_csv(float_format="{:.2f}", na_rep="NA", lineterminator="\n")
    expected = ",A\n0,1.23\n1,inf\n2,-inf\n"
    assert result == expected


def test_new_style_with_precision_edge():
    df = pd.DataFrame({"A": [1.23456789]})
    result = df.to_csv(float_format="{:.10f}", lineterminator="\n")
    expected = ",A\n0,1.2345678900\n"
    assert result == expected


def test_new_style_with_template():
    df = pd.DataFrame({"A": [1234.56789]})
    result = df.to_csv(float_format="Value: {:,.2f}", lineterminator="\n")
    expected = ',A\n0,"Value: 1,234.57"\n'
    assert result == expected


def test_to_csv_null_byte_no_escapechar():
    # GH#47871 a null byte does not require escapechar to be set
    # (was a CPython _csv regression on 3.10, fixed in 3.11)
    # engine="python": this specifically exercises the python csv module
    df = pd.DataFrame({"A": ["\x00"]})
    result = df.to_csv(index=False, lineterminator="\n", engine="python")
    assert result == "A\n\x00\n"


def test_to_csv_escapechar_roundtrip_trailing_backslash():
    # GH#33735 a value ending in the escapechar must remain readable: to_csv
    # has to escape the escapechar itself so read_csv can parse it back
    df = pd.DataFrame({0: ['"key":"value"'], 1: ["mno,"], 2: ["abc\\"], 3: ["ijk"]})
    csv = df.to_csv(header=False, index=False, escapechar="\\")
    result = pd.read_csv(io.StringIO(csv), header=None, escapechar="\\")
    assert result.iloc[0].tolist() == df.iloc[0].tolist()


def test_to_csv_categorical_tz_timestamp_with_na_rep():
    # GH#55945 to_csv with na_rep on a categorical tz-aware timestamp column
    # must not raise
    ser = pd.Series(pd.to_datetime(["2023-11-10 12:00:00+00:00"] * 3)).astype(
        "category"
    )
    df = pd.DataFrame({"ct": ser})
    result = df.to_csv(na_rep=r"\N")
    assert "2023-11-10 12:00:00+00:00" in result

    # with an actual NaT the na_rep placeholder must be emitted (the crash
    # fired even with no missing values, so also exercise the substitution path)
    ser_na = pd.Series(pd.to_datetime(["2023-11-10 12:00:00+00:00", None])).astype(
        "category"
    )
    result_na = pd.DataFrame({"ct": ser_na}).to_csv(na_rep=r"\N")
    assert r"\N" in result_na


def test_to_csv_datetime_tz_consistent_format(engine):
    # GH#62111
    # timezone-aware datetime64 is rendered differently under pyarrow: a
    # "Z" suffix instead of "+00:00", and always at microsecond precision
    # rather than dropping trailing zeros
    pyarrow_used = engine == "pyarrow"
    df = pd.DataFrame(
        {
            "timestamp": [
                datetime(2025, 8, 14, 12, 34, 56, 0, tzinfo=UTC),
                datetime(2025, 8, 14, 12, 34, 56, 1, tzinfo=UTC),
            ]
        }
    )
    with (
        check_warns_if_pyarrow_renders_differently(engine)
        if pyarrow_used
        else tm.assert_produces_warning(None)
    ):
        result = df.to_csv(index=False, engine=engine)
    if pyarrow_used:
        expected_rows = [
            '"timestamp"',
            "2025-08-14 12:34:56.000000Z",
            "2025-08-14 12:34:56.000001Z",
        ]
    else:
        expected_rows = [
            "timestamp",
            "2025-08-14 12:34:56.000000+00:00",
            "2025-08-14 12:34:56.000001+00:00",
        ]
    expected = csv_str_for_engine(expected_rows, pyarrow_used)
    assert result == expected

    df = pd.DataFrame(
        {
            "timestamp": [
                datetime(2025, 8, 14, 12, 34, 56, 0, tzinfo=UTC),
                datetime(2025, 8, 14, 12, 34, 57, 0, tzinfo=UTC),
            ]
        }
    )
    with (
        check_warns_if_pyarrow_renders_differently(engine)
        if pyarrow_used
        else tm.assert_produces_warning(None)
    ):
        result = df.to_csv(index=False, engine=engine)
    if pyarrow_used:
        expected_rows = [
            '"timestamp"',
            "2025-08-14 12:34:56.000000Z",
            "2025-08-14 12:34:57.000000Z",
        ]
    else:
        expected_rows = [
            "timestamp",
            "2025-08-14 12:34:56+00:00",
            "2025-08-14 12:34:57+00:00",
        ]
    expected = csv_str_for_engine(expected_rows, pyarrow_used)
    assert result == expected
