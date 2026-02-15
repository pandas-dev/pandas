import contextlib
import io
import os
import sys
from zipfile import ZipFile

from _csv import Error
import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    compat,
)
import pandas._testing as tm


@pytest.fixture(params=["python", "pyarrow"])
def engine(request):
    if request.param == "pyarrow":
        pytest.importorskip("pyarrow")
    return request.param


@pytest.fixture
def pyarrow_xfail(request):
    """
    Fixture that xfails a test if the engine is pyarrow.
    """
    engine = request.getfixturevalue("engine")
    if engine == "pyarrow":
        mark = pytest.mark.xfail(reason="pyarrow doesn't support this.")
        request.applymarker(mark)


def check_raises_if_pyarrow(option, engine):
    """
    Returns a context manager that ensures that the pyarrow engine raises an
    exception for unsupported options.
    """
    if engine == "pyarrow":
        raises_if_pyarrow = pytest.raises(
            ValueError,
            match=f"The {option} option is not supported with the pyarrow engine.",
        )
    else:
        raises_if_pyarrow = contextlib.nullcontext()
    return raises_if_pyarrow


xfail_pyarrow = pytest.mark.usefixtures("pyarrow_xfail")


class TestToCSV:
    @xfail_pyarrow
    def test_to_csv_with_single_column(self, temp_file, engine):
        # see gh-18676, https://bugs.python.org/issue32255
        #
        # Python's CSV library adds an extraneous '""'
        # before the newline when the NaN-value is in
        # the first row. Otherwise, only the newline
        # character is added. This behavior is inconsistent
        # and was patched in https://bugs.python.org/pull_request4672.
        df1 = DataFrame([None, 1])
        expected1 = """\
""
1.0
"""
        df1.to_csv(temp_file, header=None, index=None, engine=engine)
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected1

        df2 = DataFrame([1, None])
        expected2 = """\
1.0
""
"""
        df2.to_csv(temp_file, header=None, index=None, engine=engine)
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected2

    def test_to_csv_default_encoding(self, temp_file, engine):
        # GH17097
        df = DataFrame({"col": ["AAAAA", "ÄÄÄÄÄ", "ßßßßß", "聞聞聞聞聞"]})

        # the default to_csv encoding is utf-8.
        df.to_csv(temp_file, engine=engine)
        tm.assert_frame_equal(pd.read_csv(temp_file, index_col=0), df)

    def test_to_csv_quotechar(self, temp_file, engine):
        df = DataFrame({"col": [1, 2]})
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

    def test_to_csv_doublequote(self, temp_file):
        df = DataFrame({"col": ['a"a', '"bb"']})
        expected = '''\
"","col"
"0","a""a"
"1","""bb"""
'''

        df.to_csv(temp_file, quoting=1, doublequote=True, engine=engine)  # QUOTE_ALL
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected

        with pytest.raises(Error, match="escapechar"):
            df.to_csv(temp_file, doublequote=False, engine=engine)  # no escapechar set

    def test_to_csv_escapechar(self, temp_file, engine):
        raises_if_pyarrow = check_raises_if_pyarrow("escapechar", engine)
        df = DataFrame({"col": ['a"a', '"bb"']})
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

        df = DataFrame({"col": ["a,a", ",bb,"]})
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

    @xfail_pyarrow
    def test_csv_to_string(self, engine):
        df = DataFrame({"col": [1, 2]})
        expected_rows = [",col", "0,1", "1,2"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.to_csv(engine=engine) == expected

    @xfail_pyarrow
    def test_to_csv_decimal(self, engine):
        # see gh-781
        raises_if_pyarrow = check_raises_if_pyarrow("decimal", engine)
        df = DataFrame({"col1": [1], "col2": ["a"], "col3": [10.1]})

        expected_rows = [",col1,col2,col3", "0,1,a,10.1"]
        expected_default = tm.convert_rows_list_to_csv_str(expected_rows)
        # This assert fails for the pyarrow engine since it quotes strings
        # and the Python engine doesn't
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
        with raises_if_pyarrow:
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
        df = DataFrame({"a": [0, 1.1], "b": [2.2, 3.3], "c": 1})

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
        df = DataFrame({"a": [0, 1], "b": [2.2, 3.3], "c": 1})

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
        df = DataFrame({"a": [0, np.nan], "b": [0, 1], "c": [2, 3]})
        expected_rows = ["a,b,c", "0.0,0,2", "_,1,3"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            assert df.set_index("a").to_csv(engine=engine, na_rep="_") == expected
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, na_rep="_") == expected
            )

        # now with an index containing only NaNs
        df = DataFrame({"a": np.nan, "b": [0, 1], "c": [2, 3]})
        expected_rows = ["a,b,c", "_,0,2", "_,1,3"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            assert df.set_index("a").to_csv(engine=engine, na_rep="_") == expected
        with raises_if_pyarrow:
            assert (
                df.set_index(["a", "b"]).to_csv(engine=engine, na_rep="_") == expected
            )

        # check if na_rep parameter does not break anything when no NaN
        df = DataFrame({"a": 0, "b": [0, 1], "c": [2, 3]})
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

    @xfail_pyarrow
    def test_to_csv_date_format(self, engine):
        # GH 10209
        raises_if_pyarrow = check_raises_if_pyarrow("date_format", engine)
        df_sec = DataFrame({"A": pd.date_range("20130101", periods=5, freq="s")})
        df_day = DataFrame({"A": pd.date_range("20130101", periods=5, freq="D")})

        expected_rows = [
            ",A",
            "0,2013-01-01 00:00:00",
            "1,2013-01-01 00:00:01",
            "2,2013-01-01 00:00:02",
            "3,2013-01-01 00:00:03",
            "4,2013-01-01 00:00:04",
        ]
        expected_default_sec = tm.convert_rows_list_to_csv_str(expected_rows)
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

        expected_rows = [
            ",A",
            "0,2013-01-01",
            "1,2013-01-02",
            "2,2013-01-03",
            "3,2013-01-04",
            "4,2013-01-05",
        ]
        expected_default_day = tm.convert_rows_list_to_csv_str(expected_rows)
        with raises_if_pyarrow:
            assert df_day.to_csv(engine=engine) == expected_default_day
        with raises_if_pyarrow:
            assert (
                df_day.to_csv(date_format="%Y-%m-%d", engine=engine)
                == expected_default_day
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

    @xfail_pyarrow
    def test_to_csv_different_datetime_formats(self, engine):
        # GH#21734
        df = DataFrame(
            {
                "date": pd.to_datetime("1970-01-01"),
                "datetime": pd.date_range("1970-01-01", periods=2, freq="h"),
            }
        )
        expected_rows = [
            "date,datetime",
            "1970-01-01,1970-01-01 00:00:00",
            "1970-01-01,1970-01-01 01:00:00",
        ]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)
        assert df.to_csv(index=False, engine=engine) == expected

    @xfail_pyarrow
    def test_to_csv_date_format_in_categorical(self, engine):
        # GH#40754
        raises_if_pyarrow = check_raises_if_pyarrow("date_format", engine)
        ser = pd.Series(pd.to_datetime(["2021-03-27", pd.NaT], format="%Y-%m-%d"))
        ser = ser.astype("category")
        expected = tm.convert_rows_list_to_csv_str(["0", "2021-03-27", '""'])
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
        df = DataFrame({"a": [1.1, 2.02, pd.NA, 6.000006], "b": "c"})
        df["a"] = df["a"].astype("Float64")
        with raises_if_pyarrow:
            result = df.to_csv(index=False, engine=engine, float_format="%.5f")
            expected = tm.convert_rows_list_to_csv_str(
                ["a,b", "1.10000,c", "2.02000,c", ",c", "6.00001,c"]
            )
            assert result == expected

    @xfail_pyarrow
    def test_to_csv_float_ea_no_float_format(self, engine):
        # GH#45991
        df = DataFrame({"a": [1.1, 2.02, pd.NA, 6.000006], "b": "c"})
        df["a"] = df["a"].astype("Float64")
        result = df.to_csv(index=False, engine=engine)
        expected = tm.convert_rows_list_to_csv_str(
            ["a,b", "1.1,c", "2.02,c", ",c", "6.000006,c"]
        )
        assert result == expected

    @xfail_pyarrow
    def test_to_csv_multi_index(self, engine):
        # see gh-6618
        df = DataFrame([1], columns=pd.MultiIndex.from_arrays([[1], [2]]))

        exp_rows = [",1", ",2", "0,1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(engine=engine) == exp

        exp_rows = ["1", "2", "1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(index=False, engine=engine) == exp

        df = DataFrame(
            [1],
            columns=pd.MultiIndex.from_arrays([[1], [2]]),
            index=pd.MultiIndex.from_arrays([[1], [2]]),
        )

        exp_rows = [",,1", ",,2", "1,2,1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(engine=engine) == exp

        exp_rows = ["1", "2", "1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(index=False, engine=engine) == exp

        df = DataFrame([1], columns=pd.MultiIndex.from_arrays([["foo"], ["bar"]]))

        exp_rows = [",foo", ",bar", "0,1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(engine=engine) == exp

        exp_rows = ["foo", "bar", "1"]
        exp = tm.convert_rows_list_to_csv_str(exp_rows)
        assert df.to_csv(index=False, engine=engine) == exp

    @xfail_pyarrow
    @pytest.mark.parametrize(
        "ind,expected",
        [
            (
                pd.MultiIndex(levels=[[1.0]], codes=[[0]], names=["x"]),
                "x,data\n1.0,1\n",
            ),
            (
                pd.MultiIndex(
                    levels=[[1.0], [2.0]], codes=[[0], [0]], names=["x", "y"]
                ),
                "x,y,data\n1.0,2.0,1\n",
            ),
        ],
    )
    def test_to_csv_single_level_multi_index(
        self, ind, expected, frame_or_series, engine
    ):
        # see gh-19589
        obj = frame_or_series(pd.Series([1], ind, name="data"))
        result = obj.to_csv(lineterminator="\n", header=True, engine=engine)
        assert result == expected

    def test_to_csv_string_array_ascii(self, temp_file, engine):
        # GH 10813
        raises_if_pyarrow = check_raises_if_pyarrow("encoding", engine)
        str_array = [{"names": ["foo", "bar"]}, {"names": ["baz", "qux"]}]
        df = DataFrame(str_array)
        expected_ascii = """\
,names
0,"['foo', 'bar']"
1,"['baz', 'qux']"
"""
        with raises_if_pyarrow:
            df.to_csv(temp_file, encoding="ascii", engine=engine)
            with open(temp_file, encoding="utf-8") as f:
                assert f.read() == expected_ascii

    @xfail_pyarrow
    def test_to_csv_string_array_utf8(self, temp_file, engine):
        # GH 10813
        str_array = [{"names": ["foo", "bar"]}, {"names": ["baz", "qux"]}]
        df = DataFrame(str_array)
        expected_utf8 = """\
,names
0,"['foo', 'bar']"
1,"['baz', 'qux']"
"""
        df.to_csv(temp_file, encoding="utf-8", engine=engine)
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected_utf8

    @xfail_pyarrow
    def test_to_csv_string_with_lf(self, temp_file, engine):
        # GH 20353
        raises_if_pyarrow = check_raises_if_pyarrow("lineterminator", engine)
        data = {"int": [1, 2, 3], "str_lf": ["abc", "d\nef", "g\nh\n\ni"]}
        df = DataFrame(data)

        # case 1: The default line terminator(=os.linesep)(PR 21406)
        os_linesep = os.linesep.encode("utf-8")
        expected_noarg = (
            b"int,str_lf"
            + os_linesep
            + b"1,abc"
            + os_linesep
            + b'2,"d\nef"'
            + os_linesep
            + b'3,"g\nh\n\ni"'
            + os_linesep
        )
        df.to_csv(temp_file, index=False, engine=engine)
        with open(temp_file, "rb") as f:
            assert f.read() == expected_noarg

        # case 2: LF as line terminator
        expected_lf = b'int,str_lf\n1,abc\n2,"d\nef"\n3,"g\nh\n\ni"\n'
        with raises_if_pyarrow:
            df.to_csv(temp_file, lineterminator="\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_lf

        # case 3: CRLF as line terminator
        # 'lineterminator' should not change inner element
        expected_crlf = b'int,str_lf\r\n1,abc\r\n2,"d\nef"\r\n3,"g\nh\n\ni"\r\n'
        with raises_if_pyarrow:
            df.to_csv(temp_file, lineterminator="\r\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_crlf

    @xfail_pyarrow
    def test_to_csv_string_with_crlf(self, temp_file, engine):
        # GH 20353
        raises_if_pyarrow = check_raises_if_pyarrow("lineterminator", engine)
        data = {"int": [1, 2, 3], "str_crlf": ["abc", "d\r\nef", "g\r\nh\r\n\r\ni"]}
        df = DataFrame(data)
        # case 1: The default line terminator(=os.linesep)(PR 21406)
        os_linesep = os.linesep.encode("utf-8")
        expected_noarg = (
            b"int,str_crlf"
            + os_linesep
            + b"1,abc"
            + os_linesep
            + b'2,"d\r\nef"'
            + os_linesep
            + b'3,"g\r\nh\r\n\r\ni"'
            + os_linesep
        )
        df.to_csv(temp_file, index=False, engine=engine)
        with open(temp_file, "rb") as f:
            assert f.read() == expected_noarg

        # case 2: LF as line terminator
        expected_lf = b'int,str_crlf\n1,abc\n2,"d\r\nef"\n3,"g\r\nh\r\n\r\ni"\n'
        with raises_if_pyarrow:
            df.to_csv(temp_file, lineterminator="\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_lf

        # case 3: CRLF as line terminator
        # 'lineterminator' should not change inner element
        expected_crlf = (
            b'int,str_crlf\r\n1,abc\r\n2,"d\r\nef"\r\n3,"g\r\nh\r\n\r\ni"\r\n'
        )
        with raises_if_pyarrow:
            df.to_csv(temp_file, lineterminator="\r\n", index=False, engine=engine)
            with open(temp_file, "rb") as f:
                assert f.read() == expected_crlf

    @xfail_pyarrow
    def test_to_csv_stdout_file(self, capsys, engine):
        # GH 21561
        df = DataFrame([["foo", "bar"], ["baz", "qux"]], columns=["name_1", "name_2"])
        expected_rows = [",name_1,name_2", "0,foo,bar", "1,baz,qux"]
        expected_ascii = tm.convert_rows_list_to_csv_str(expected_rows)

        df.to_csv(sys.stdout, encoding="ascii", engine=engine)
        captured = capsys.readouterr()

        assert captured.out == expected_ascii
        assert not sys.stdout.closed

    @xfail_pyarrow
    @pytest.mark.xfail(
        compat.is_platform_windows(),
        reason=(
            "Especially in Windows, file stream should not be passed"
            "to csv writer without newline='' option."
            "(https://docs.python.org/3/library/csv.html#csv.writer)"
        ),
    )
    def test_to_csv_write_to_open_file(self, temp_file, engine):
        # GH 21696
        df = DataFrame({"a": ["x", "y", "z"]})
        expected = """\
manual header
x
y
z
"""
        with open(temp_file, "w", encoding="utf-8") as f:
            f.write("manual header\n")
            if engine == "pyarrow":
                raise_if_pyarrow = pytest.raises(
                    ValueError,
                    match="The pyarrow engine can only open files in abinary mode.",
                )
            else:
                raise_if_pyarrow = contextlib.nullcontext()
            with raise_if_pyarrow:
                df.to_csv(f, header=None, index=None, engine=engine)
        with open(temp_file, encoding="utf-8") as f:
            assert f.read() == expected

    @xfail_pyarrow
    def test_to_csv_write_to_open_file_with_newline_py3(self, temp_file, engine):
        # see gh-21696
        # see gh-20353
        df = DataFrame({"a": ["x", "y", "z"]})
        expected_rows = ["x", "y", "z"]
        expected = "manual header\n" + tm.convert_rows_list_to_csv_str(expected_rows)

        # TODO: Open in bytes mode for pyarrow
        with open(temp_file, "w", newline="", encoding="utf-8") as f:
            f.write("manual header\n")
            if engine == "pyarrow":
                raise_if_pyarrow = pytest.raises(
                    ValueError,
                    match="The pyarrow engine can only open file in abinary mode.",
                )
            else:
                raise_if_pyarrow = contextlib.nullcontext()
            with raise_if_pyarrow:
                df.to_csv(f, header=None, index=None, engine=engine)

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

        df = DataFrame({"A": [1]})

        to_compression = "infer" if to_infer else compression
        read_compression = "infer" if read_infer else compression

        path_ext = str(temp_file) + "." + compression_to_extension[compression]
        df.to_csv(path_ext, compression=to_compression, engine=engine)
        result = pd.read_csv(path_ext, index_col=0, compression=read_compression)
        tm.assert_frame_equal(result, df)

    def test_to_csv_compression_dict(self, compression_only, temp_file, engine):
        # GH 26023
        method = compression_only
        df = DataFrame({"ABC": [1]})
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
        df = DataFrame({"ABC": [1]})
        compression = {"some_option": True}
        msg = "must have key 'method'"

        with pytest.raises(ValueError, match=msg):
            df.to_csv(temp_file, compression=compression, engine=engine)

    @pytest.mark.parametrize("compression", ["zip", "infer"])
    @pytest.mark.parametrize("archive_name", ["test_to_csv.csv", "test_to_csv.zip"])
    def test_to_csv_zip_arguments(self, compression, archive_name, temp_file, engine):
        # GH 26023
        df = DataFrame({"ABC": [1]})

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
        df = DataFrame({"ABC": [1]})
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
        df = DataFrame({"c": [pd.NA] * 3})
        df = df.astype(df_new_type)
        expected_rows = ["c", "mynull", "mynull", "mynull"]
        expected = tm.convert_rows_list_to_csv_str(expected_rows)

        with raises_if_pyarrow:
            result = df.to_csv(
                index=False, na_rep="mynull", encoding="ascii", engine=engine
            )

            assert expected == result

    @xfail_pyarrow
    def test_to_csv_timedelta_precision(self, engine):
        # GH 6783
        s = pd.Series([1, 1]).astype("timedelta64[ns]")
        buf = io.StringIO()
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
            ser = pd.Series(data, index=Index(data, dtype=object), dtype=object)
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
        df = DataFrame(
            1.1 * np.arange(120).reshape((30, 4)),
            columns=Index(list("ABCD")),
            index=Index([f"i-{i}" for i in range(30)]),
        )

        with open(temp_file, mode="w+b") as handle:
            if engine == "pyarrow" and mode == "w":
                raises_if_pyarrow = pytest.raises(
                    ValueError,
                    match="The pyarrow engine can only open files in binary mode.",
                )
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
                DataFrame().to_csv(
                    handle, mode=mode, encoding="utf-8-sig", engine=engine
                )

                handle.seek(0)
                assert handle.read().startswith(b'\xef\xbb\xbf""')


def test_to_csv_iterative_compression_name(compression, temp_file, engine):
    # GH 38714
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD")),
        index=Index([f"i-{i}" for i in range(30)]),
    )
    df.to_csv(temp_file, compression=compression, chunksize=1, engine=engine)
    tm.assert_frame_equal(
        pd.read_csv(temp_file, compression=compression, index_col=0), df
    )


def test_to_csv_iterative_compression_buffer(compression, engine):
    # GH 38714
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD")),
        index=Index([f"i-{i}" for i in range(30)]),
    )
    with io.BytesIO() as buffer:
        df.to_csv(buffer, compression=compression, chunksize=1, engine=engine)
        buffer.seek(0)
        tm.assert_frame_equal(
            pd.read_csv(buffer, compression=compression, index_col=0), df
        )
        assert not buffer.closed


def test_new_style_float_format_basic():
    df = DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A\n0,1234.57\n1,9876.54\n"
    assert result == expected


def test_new_style_float_format_thousands():
    df = DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format="{:,.2f}", lineterminator="\n")
    expected = ',A\n0,"1,234.57"\n1,"9,876.54"\n'
    assert result == expected


def test_new_style_scientific_format():
    df = DataFrame({"A": [0.000123, 0.000456]})
    result = df.to_csv(float_format="{:.2e}", lineterminator="\n")
    expected = ",A\n0,1.23e-04\n1,4.56e-04\n"
    assert result == expected


def test_new_style_with_nan():
    df = DataFrame({"A": [1.23, np.nan, 4.56]})
    result = df.to_csv(float_format="{:.2f}", na_rep="NA", lineterminator="\n")
    expected = ",A\n0,1.23\n1,NA\n2,4.56\n"
    assert result == expected


def test_new_style_with_mixed_types():
    df = DataFrame({"A": [1.23, 4.56], "B": ["x", "y"]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A,B\n0,1.23,x\n1,4.56,y\n"
    assert result == expected


def test_new_style_with_mixed_types_in_column():
    df = DataFrame({"A": [1.23, "text", 4.56]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A\n0,1.23\n1,text\n2,4.56\n"
    assert result == expected


def test_invalid_new_style_format_missing_brace():
    df = DataFrame({"A": [1.23]})
    with pytest.raises(ValueError, match="Invalid new-style format string '{:.2f"):
        df.to_csv(float_format="{:.2f")


def test_invalid_new_style_format_specifier():
    df = DataFrame({"A": [1.23]})
    with pytest.raises(ValueError, match="Invalid new-style format string '{:.2z}'"):
        df.to_csv(float_format="{:.2z}")


def test_old_style_format_compatibility():
    df = DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format="%.2f", lineterminator="\n")
    expected = ",A\n0,1234.57\n1,9876.54\n"
    assert result == expected


def test_callable_float_format_compatibility():
    df = DataFrame({"A": [1234.56789, 9876.54321]})
    result = df.to_csv(float_format=lambda x: f"{x:,.2f}", lineterminator="\n")
    expected = ',A\n0,"1,234.57"\n1,"9,876.54"\n'
    assert result == expected


def test_no_float_format():
    df = DataFrame({"A": [1.23, 4.56]})
    result = df.to_csv(float_format=None, lineterminator="\n")
    expected = ",A\n0,1.23\n1,4.56\n"
    assert result == expected


def test_large_numbers():
    df = DataFrame({"A": [1e308, 2e308]})
    result = df.to_csv(float_format="{:.2e}", lineterminator="\n")
    expected = ",A\n0,1.00e+308\n1,inf\n"
    assert result == expected


def test_zero_and_negative():
    df = DataFrame({"A": [0.0, -1.23456]})
    result = df.to_csv(float_format="{:+.2f}", lineterminator="\n")
    expected = ",A\n0,+0.00\n1,-1.23\n"
    assert result == expected


def test_unicode_format():
    df = DataFrame({"A": [1.23, 4.56]})
    result = df.to_csv(float_format="{:.2f}€", encoding="utf-8", lineterminator="\n")
    expected = ",A\n0,1.23€\n1,4.56€\n"
    assert result == expected


def test_empty_dataframe():
    df = DataFrame({"A": []})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A\n"
    assert result == expected


def test_multi_column_float():
    df = DataFrame({"A": [1.23, 4.56], "B": [7.89, 0.12]})
    result = df.to_csv(float_format="{:.2f}", lineterminator="\n")
    expected = ",A,B\n0,1.23,7.89\n1,4.56,0.12\n"
    assert result == expected


def test_invalid_float_format_type():
    df = DataFrame({"A": [1.23]})
    with pytest.raises(ValueError, match="float_format must be a string or callable"):
        df.to_csv(float_format=123)


def test_new_style_with_inf():
    df = DataFrame({"A": [1.23, np.inf, -np.inf]})
    result = df.to_csv(float_format="{:.2f}", na_rep="NA", lineterminator="\n")
    expected = ",A\n0,1.23\n1,inf\n2,-inf\n"
    assert result == expected


def test_new_style_with_precision_edge():
    df = DataFrame({"A": [1.23456789]})
    result = df.to_csv(float_format="{:.10f}", lineterminator="\n")
    expected = ",A\n0,1.2345678900\n"
    assert result == expected


def test_new_style_with_template():
    df = DataFrame({"A": [1234.56789]})
    result = df.to_csv(float_format="Value: {:,.2f}", lineterminator="\n")
    expected = ',A\n0,"Value: 1,234.57"\n'
    assert result == expected
