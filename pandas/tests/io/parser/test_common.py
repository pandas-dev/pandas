"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""

import codecs
from collections import OrderedDict
import csv
from datetime import datetime
from io import BytesIO, StringIO
import os
import platform
from tempfile import TemporaryFile

import numpy as np
import pytest

from pandas._libs.tslib import Timestamp
from pandas.errors import DtypeWarning, EmptyDataError, ParserError

from pandas import DataFrame, Index, MultiIndex, Series, compat, concat
import pandas.util.testing as tm

from pandas.io.common import URLError
from pandas.io.parsers import CParserWrapper, TextFileReader, TextParser


def test_override_set_noconvert_columns():
    # see gh-17351
    #
    # Usecols needs to be sorted in _set_noconvert_columns based
    # on the test_usecols_with_parse_dates test from test_usecols.py
    class MyTextFileReader(TextFileReader):
        def __init__(self):
            self._currow = 0
            self.squeeze = False

    class MyCParserWrapper(CParserWrapper):
        def _set_noconvert_columns(self):
            if self.usecols_dtype == "integer":
                # self.usecols is a set, which is documented as unordered
                # but in practice, a CPython set of integers is sorted.
                # In other implementations this assumption does not hold.
                # The following code simulates a different order, which
                # before GH 17351 would cause the wrong columns to be
                # converted via the parse_dates parameter
                self.usecols = list(self.usecols)
                self.usecols.reverse()
            return CParserWrapper._set_noconvert_columns(self)

    data = """a,b,c,d,e
0,1,20140101,0900,4
0,1,20140102,1000,4"""

    parse_dates = [[1, 2]]
    cols = {
        "a": [0, 0],
        "c_d": [Timestamp("2014-01-01 09:00:00"), Timestamp("2014-01-02 10:00:00")],
    }
    expected = DataFrame(cols, columns=["c_d", "a"])

    parser = MyTextFileReader()
    parser.options = {
        "usecols": [0, 2, 3],
        "parse_dates": parse_dates,
        "delimiter": ",",
    }
    parser._engine = MyCParserWrapper(StringIO(data), **parser.options)

    result = parser.read()
    tm.assert_frame_equal(result, expected)


def test_bytes_io_input(all_parsers):
    encoding = "cp1255"
    parser = all_parsers

    data = BytesIO("שלום:1234\n562:123".encode(encoding))
    result = parser.read_csv(data, sep=":", encoding=encoding)

    expected = DataFrame([[562, 123]], columns=["שלום", "1234"])
    tm.assert_frame_equal(result, expected)


def test_empty_decimal_marker(all_parsers):
    data = """A|B|C
1|2,334|5
10|13|10.
"""
    # Parsers support only length-1 decimals
    msg = "Only length-1 decimal markers supported"
    parser = all_parsers

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), decimal="")


def test_bad_stream_exception(all_parsers, csv_dir_path):
    # see gh-13652
    #
    # This test validates that both the Python engine and C engine will
    # raise UnicodeDecodeError instead of C engine raising ParserError
    # and swallowing the exception that caused read to fail.
    path = os.path.join(csv_dir_path, "sauron.SHIFT_JIS.csv")
    codec = codecs.lookup("utf-8")
    utf8 = codecs.lookup("utf-8")
    parser = all_parsers
    msg = "'utf-8' codec can't decode byte"

    # Stream must be binary UTF8.
    with open(path, "rb") as handle, codecs.StreamRecoder(
        handle, utf8.encode, utf8.decode, codec.streamreader, codec.streamwriter
    ) as stream:

        with pytest.raises(UnicodeDecodeError, match=msg):
            parser.read_csv(stream)


def test_read_csv_local(all_parsers, csv1):
    prefix = "file:///" if compat.is_platform_windows() else "file://"
    parser = all_parsers

    fname = prefix + str(os.path.abspath(csv1))
    result = parser.read_csv(fname, index_col=0, parse_dates=True)

    expected = DataFrame(
        [
            [0.980269, 3.685731, -0.364216805298, -1.159738],
            [1.047916, -0.041232, -0.16181208307, 0.212549],
            [0.498581, 0.731168, -0.537677223318, 1.346270],
            [1.120202, 1.567621, 0.00364077397681, 0.675253],
            [-0.487094, 0.571455, -1.6116394093, 0.103469],
            [0.836649, 0.246462, 0.588542635376, 1.062782],
            [-0.157161, 1.340307, 1.1957779562, -1.097007],
        ],
        columns=["A", "B", "C", "D"],
        index=Index(
            [
                datetime(2000, 1, 3),
                datetime(2000, 1, 4),
                datetime(2000, 1, 5),
                datetime(2000, 1, 6),
                datetime(2000, 1, 7),
                datetime(2000, 1, 10),
                datetime(2000, 1, 11),
            ],
            name="index",
        ),
    )
    tm.assert_frame_equal(result, expected)


def test_1000_sep(all_parsers):
    parser = all_parsers
    data = """A|B|C
1|2,334|5
10|13|10.
"""
    expected = DataFrame({"A": [1, 10], "B": [2334, 13], "C": [5, 10.0]})

    result = parser.read_csv(StringIO(data), sep="|", thousands=",")
    tm.assert_frame_equal(result, expected)


def test_squeeze(all_parsers):
    data = """\
a,1
b,2
c,3
"""
    parser = all_parsers
    index = Index(["a", "b", "c"], name=0)
    expected = Series([1, 2, 3], name=1, index=index)

    result = parser.read_csv(StringIO(data), index_col=0, header=None, squeeze=True)
    tm.assert_series_equal(result, expected)

    # see gh-8217
    #
    # Series should not be a view.
    assert not result._is_view


def test_malformed(all_parsers):
    # see gh-6607
    parser = all_parsers
    data = """ignore
A,B,C
1,2,3 # comment
1,2,3,4,5
2,3,4
"""
    msg = "Expected 3 fields in line 4, saw 5"
    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), header=1, comment="#")


@pytest.mark.parametrize("nrows", [5, 3, None])
def test_malformed_chunks(all_parsers, nrows):
    data = """ignore
A,B,C
skip
1,2,3
3,5,10 # comment
1,2,3,4,5
2,3,4
"""
    parser = all_parsers
    msg = "Expected 3 fields in line 6, saw 5"
    reader = parser.read_csv(
        StringIO(data), header=1, comment="#", iterator=True, chunksize=1, skiprows=[2]
    )

    with pytest.raises(ParserError, match=msg):
        reader.read(nrows)


def test_unnamed_columns(all_parsers):
    data = """A,B,C,,
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
    parser = all_parsers
    expected = DataFrame(
        [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]],
        dtype=np.int64,
        columns=["A", "B", "C", "Unnamed: 3", "Unnamed: 4"],
    )
    result = parser.read_csv(StringIO(data))
    tm.assert_frame_equal(result, expected)


def test_csv_mixed_type(all_parsers):
    data = """A,B,C
a,1,2
b,3,4
c,4,5
"""
    parser = all_parsers
    expected = DataFrame({"A": ["a", "b", "c"], "B": [1, 3, 4], "C": [2, 4, 5]})
    result = parser.read_csv(StringIO(data))
    tm.assert_frame_equal(result, expected)


def test_read_csv_low_memory_no_rows_with_index(all_parsers):
    # see gh-21141
    parser = all_parsers

    if not parser.low_memory:
        pytest.skip("This is a low-memory specific test")

    data = """A,B,C
1,1,1,2
2,2,3,4
3,3,4,5
"""
    result = parser.read_csv(StringIO(data), low_memory=True, index_col=0, nrows=0)
    expected = DataFrame(columns=["A", "B", "C"])
    tm.assert_frame_equal(result, expected)


def test_read_csv_dataframe(all_parsers, csv1):
    parser = all_parsers
    result = parser.read_csv(csv1, index_col=0, parse_dates=True)

    expected = DataFrame(
        [
            [0.980269, 3.685731, -0.364216805298, -1.159738],
            [1.047916, -0.041232, -0.16181208307, 0.212549],
            [0.498581, 0.731168, -0.537677223318, 1.346270],
            [1.120202, 1.567621, 0.00364077397681, 0.675253],
            [-0.487094, 0.571455, -1.6116394093, 0.103469],
            [0.836649, 0.246462, 0.588542635376, 1.062782],
            [-0.157161, 1.340307, 1.1957779562, -1.097007],
        ],
        columns=["A", "B", "C", "D"],
        index=Index(
            [
                datetime(2000, 1, 3),
                datetime(2000, 1, 4),
                datetime(2000, 1, 5),
                datetime(2000, 1, 6),
                datetime(2000, 1, 7),
                datetime(2000, 1, 10),
                datetime(2000, 1, 11),
            ],
            name="index",
        ),
    )
    tm.assert_frame_equal(result, expected)


def test_read_csv_no_index_name(all_parsers, csv_dir_path):
    parser = all_parsers
    csv2 = os.path.join(csv_dir_path, "test2.csv")
    result = parser.read_csv(csv2, index_col=0, parse_dates=True)

    expected = DataFrame(
        [
            [0.980269, 3.685731, -0.364216805298, -1.159738, "foo"],
            [1.047916, -0.041232, -0.16181208307, 0.212549, "bar"],
            [0.498581, 0.731168, -0.537677223318, 1.346270, "baz"],
            [1.120202, 1.567621, 0.00364077397681, 0.675253, "qux"],
            [-0.487094, 0.571455, -1.6116394093, 0.103469, "foo2"],
        ],
        columns=["A", "B", "C", "D", "E"],
        index=Index(
            [
                datetime(2000, 1, 3),
                datetime(2000, 1, 4),
                datetime(2000, 1, 5),
                datetime(2000, 1, 6),
                datetime(2000, 1, 7),
            ]
        ),
    )
    tm.assert_frame_equal(result, expected)


def test_read_csv_unicode(all_parsers):
    parser = all_parsers
    data = BytesIO("\u0141aski, Jan;1".encode("utf-8"))

    result = parser.read_csv(data, sep=";", encoding="utf-8", header=None)
    expected = DataFrame([["\u0141aski, Jan", 1]])
    tm.assert_frame_equal(result, expected)


def test_read_csv_wrong_num_columns(all_parsers):
    # Too few columns.
    data = """A,B,C,D,E,F
1,2,3,4,5,6
6,7,8,9,10,11,12
11,12,13,14,15,16
"""
    parser = all_parsers
    msg = "Expected 6 fields in line 3, saw 7"

    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data))


def test_read_duplicate_index_explicit(all_parsers):
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo,12,13,14,15
bar,12,13,14,15
"""
    parser = all_parsers
    result = parser.read_csv(StringIO(data), index_col=0)

    expected = DataFrame(
        [
            [2, 3, 4, 5],
            [7, 8, 9, 10],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
        ],
        columns=["A", "B", "C", "D"],
        index=Index(["foo", "bar", "baz", "qux", "foo", "bar"], name="index"),
    )
    tm.assert_frame_equal(result, expected)


def test_read_duplicate_index_implicit(all_parsers):
    data = """A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo,12,13,14,15
bar,12,13,14,15
"""
    parser = all_parsers
    result = parser.read_csv(StringIO(data))

    expected = DataFrame(
        [
            [2, 3, 4, 5],
            [7, 8, 9, 10],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
        ],
        columns=["A", "B", "C", "D"],
        index=Index(["foo", "bar", "baz", "qux", "foo", "bar"]),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data,kwargs,expected",
    [
        (
            "A,B\nTrue,1\nFalse,2\nTrue,3",
            dict(),
            DataFrame([[True, 1], [False, 2], [True, 3]], columns=["A", "B"]),
        ),
        (
            "A,B\nYES,1\nno,2\nyes,3\nNo,3\nYes,3",
            dict(true_values=["yes", "Yes", "YES"], false_values=["no", "NO", "No"]),
            DataFrame(
                [[True, 1], [False, 2], [True, 3], [False, 3], [True, 3]],
                columns=["A", "B"],
            ),
        ),
        (
            "A,B\nTRUE,1\nFALSE,2\nTRUE,3",
            dict(),
            DataFrame([[True, 1], [False, 2], [True, 3]], columns=["A", "B"]),
        ),
        (
            "A,B\nfoo,bar\nbar,foo",
            dict(true_values=["foo"], false_values=["bar"]),
            DataFrame([[True, False], [False, True]], columns=["A", "B"]),
        ),
    ],
)
def test_parse_bool(all_parsers, data, kwargs, expected):
    parser = all_parsers
    result = parser.read_csv(StringIO(data), **kwargs)
    tm.assert_frame_equal(result, expected)


def test_int_conversion(all_parsers):
    data = """A,B
1.0,1
2.0,2
3.0,3
"""
    parser = all_parsers
    result = parser.read_csv(StringIO(data))

    expected = DataFrame([[1.0, 1], [2.0, 2], [3.0, 3]], columns=["A", "B"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("nrows", [3, 3.0])
def test_read_nrows(all_parsers, nrows):
    # see gh-10476
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    expected = DataFrame(
        [["foo", 2, 3, 4, 5], ["bar", 7, 8, 9, 10], ["baz", 12, 13, 14, 15]],
        columns=["index", "A", "B", "C", "D"],
    )
    parser = all_parsers

    result = parser.read_csv(StringIO(data), nrows=nrows)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("nrows", [1.2, "foo", -1])
def test_read_nrows_bad(all_parsers, nrows):
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    msg = r"'nrows' must be an integer >=0"
    parser = all_parsers

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), nrows=nrows)


@pytest.mark.parametrize("index_col", [0, "index"])
def test_read_chunksize_with_index(all_parsers, index_col):
    parser = all_parsers
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

    reader = parser.read_csv(StringIO(data), index_col=0, chunksize=2)
    expected = DataFrame(
        [
            ["foo", 2, 3, 4, 5],
            ["bar", 7, 8, 9, 10],
            ["baz", 12, 13, 14, 15],
            ["qux", 12, 13, 14, 15],
            ["foo2", 12, 13, 14, 15],
            ["bar2", 12, 13, 14, 15],
        ],
        columns=["index", "A", "B", "C", "D"],
    )
    expected = expected.set_index("index")

    chunks = list(reader)
    tm.assert_frame_equal(chunks[0], expected[:2])
    tm.assert_frame_equal(chunks[1], expected[2:4])
    tm.assert_frame_equal(chunks[2], expected[4:])


@pytest.mark.parametrize("chunksize", [1.3, "foo", 0])
def test_read_chunksize_bad(all_parsers, chunksize):
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    msg = r"'chunksize' must be an integer >=1"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), chunksize=chunksize)


@pytest.mark.parametrize("chunksize", [2, 8])
def test_read_chunksize_and_nrows(all_parsers, chunksize):
    # see gh-15755
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    kwargs = dict(index_col=0, nrows=5)

    reader = parser.read_csv(StringIO(data), chunksize=chunksize, **kwargs)
    expected = parser.read_csv(StringIO(data), **kwargs)
    tm.assert_frame_equal(concat(reader), expected)


def test_read_chunksize_and_nrows_changing_size(all_parsers):
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    kwargs = dict(index_col=0, nrows=5)

    reader = parser.read_csv(StringIO(data), chunksize=8, **kwargs)
    expected = parser.read_csv(StringIO(data), **kwargs)

    tm.assert_frame_equal(reader.get_chunk(size=2), expected.iloc[:2])
    tm.assert_frame_equal(reader.get_chunk(size=4), expected.iloc[2:5])

    with pytest.raises(StopIteration, match=""):
        reader.get_chunk(size=3)


def test_get_chunk_passed_chunksize(all_parsers):
    parser = all_parsers
    data = """A,B,C
1,2,3
4,5,6
7,8,9
1,2,3"""

    reader = parser.read_csv(StringIO(data), chunksize=2)
    result = reader.get_chunk()

    expected = DataFrame([[1, 2, 3], [4, 5, 6]], columns=["A", "B", "C"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("kwargs", [dict(), dict(index_col=0)])
def test_read_chunksize_compat(all_parsers, kwargs):
    # see gh-12185
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    reader = parser.read_csv(StringIO(data), chunksize=2, **kwargs)

    result = parser.read_csv(StringIO(data), **kwargs)
    tm.assert_frame_equal(concat(reader), result)


def test_read_chunksize_jagged_names(all_parsers):
    # see gh-23509
    parser = all_parsers
    data = "\n".join(["0"] * 7 + [",".join(["0"] * 10)])

    expected = DataFrame([[0] + [np.nan] * 9] * 7 + [[0] * 10])
    reader = parser.read_csv(StringIO(data), names=range(10), chunksize=4)

    result = concat(reader)
    tm.assert_frame_equal(result, expected)


def test_read_data_list(all_parsers):
    parser = all_parsers
    kwargs = dict(index_col=0)
    data = "A,B,C\nfoo,1,2,3\nbar,4,5,6"

    data_list = [["A", "B", "C"], ["foo", "1", "2", "3"], ["bar", "4", "5", "6"]]
    expected = parser.read_csv(StringIO(data), **kwargs)

    parser = TextParser(data_list, chunksize=2, **kwargs)
    result = parser.read()

    tm.assert_frame_equal(result, expected)


def test_iterator(all_parsers):
    # see gh-6607
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    kwargs = dict(index_col=0)

    expected = parser.read_csv(StringIO(data), **kwargs)
    reader = parser.read_csv(StringIO(data), iterator=True, **kwargs)

    first_chunk = reader.read(3)
    tm.assert_frame_equal(first_chunk, expected[:3])

    last_chunk = reader.read(5)
    tm.assert_frame_equal(last_chunk, expected[3:])


def test_iterator2(all_parsers):
    parser = all_parsers
    data = """A,B,C
foo,1,2,3
bar,4,5,6
baz,7,8,9
"""

    reader = parser.read_csv(StringIO(data), iterator=True)
    result = list(reader)

    expected = DataFrame(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        index=["foo", "bar", "baz"],
        columns=["A", "B", "C"],
    )
    tm.assert_frame_equal(result[0], expected)


def test_reader_list(all_parsers):
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    kwargs = dict(index_col=0)

    lines = list(csv.reader(StringIO(data)))
    reader = TextParser(lines, chunksize=2, **kwargs)

    expected = parser.read_csv(StringIO(data), **kwargs)
    chunks = list(reader)

    tm.assert_frame_equal(chunks[0], expected[:2])
    tm.assert_frame_equal(chunks[1], expected[2:4])
    tm.assert_frame_equal(chunks[2], expected[4:])


def test_reader_list_skiprows(all_parsers):
    data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    parser = all_parsers
    kwargs = dict(index_col=0)

    lines = list(csv.reader(StringIO(data)))
    reader = TextParser(lines, chunksize=2, skiprows=[1], **kwargs)

    expected = parser.read_csv(StringIO(data), **kwargs)
    chunks = list(reader)

    tm.assert_frame_equal(chunks[0], expected[1:3])


def test_iterator_stop_on_chunksize(all_parsers):
    # gh-3967: stopping iteration when chunksize is specified
    parser = all_parsers
    data = """A,B,C
foo,1,2,3
bar,4,5,6
baz,7,8,9
"""

    reader = parser.read_csv(StringIO(data), chunksize=1)
    result = list(reader)

    assert len(result) == 3
    expected = DataFrame(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        index=["foo", "bar", "baz"],
        columns=["A", "B", "C"],
    )
    tm.assert_frame_equal(concat(result), expected)


@pytest.mark.parametrize(
    "kwargs", [dict(iterator=True, chunksize=1), dict(iterator=True), dict(chunksize=1)]
)
def test_iterator_skipfooter_errors(all_parsers, kwargs):
    msg = "'skipfooter' not supported for 'iteration'"
    parser = all_parsers
    data = "a\n1\n2"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), skipfooter=1, **kwargs)


def test_nrows_skipfooter_errors(all_parsers):
    msg = "'skipfooter' not supported with 'nrows'"
    data = "a\n1\n2\n3\n4\n5\n6"
    parser = all_parsers

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), skipfooter=1, nrows=5)


@pytest.mark.parametrize(
    "data,kwargs,expected",
    [
        (
            """foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
""",
            dict(index_col=0, names=["index", "A", "B", "C", "D"]),
            DataFrame(
                [
                    [2, 3, 4, 5],
                    [7, 8, 9, 10],
                    [12, 13, 14, 15],
                    [12, 13, 14, 15],
                    [12, 13, 14, 15],
                    [12, 13, 14, 15],
                ],
                index=Index(["foo", "bar", "baz", "qux", "foo2", "bar2"], name="index"),
                columns=["A", "B", "C", "D"],
            ),
        ),
        (
            """foo,one,2,3,4,5
foo,two,7,8,9,10
foo,three,12,13,14,15
bar,one,12,13,14,15
bar,two,12,13,14,15
""",
            dict(index_col=[0, 1], names=["index1", "index2", "A", "B", "C", "D"]),
            DataFrame(
                [
                    [2, 3, 4, 5],
                    [7, 8, 9, 10],
                    [12, 13, 14, 15],
                    [12, 13, 14, 15],
                    [12, 13, 14, 15],
                ],
                index=MultiIndex.from_tuples(
                    [
                        ("foo", "one"),
                        ("foo", "two"),
                        ("foo", "three"),
                        ("bar", "one"),
                        ("bar", "two"),
                    ],
                    names=["index1", "index2"],
                ),
                columns=["A", "B", "C", "D"],
            ),
        ),
    ],
)
def test_pass_names_with_index(all_parsers, data, kwargs, expected):
    parser = all_parsers
    result = parser.read_csv(StringIO(data), **kwargs)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("index_col", [[0, 1], [1, 0]])
def test_multi_index_no_level_names(all_parsers, index_col):
    data = """index1,index2,A,B,C,D
foo,one,2,3,4,5
foo,two,7,8,9,10
foo,three,12,13,14,15
bar,one,12,13,14,15
bar,two,12,13,14,15
"""
    headless_data = "\n".join(data.split("\n")[1:])

    names = ["A", "B", "C", "D"]
    parser = all_parsers

    result = parser.read_csv(
        StringIO(headless_data), index_col=index_col, header=None, names=names
    )
    expected = parser.read_csv(StringIO(data), index_col=index_col)

    # No index names in headless data.
    expected.index.names = [None] * 2
    tm.assert_frame_equal(result, expected)


def test_multi_index_no_level_names_implicit(all_parsers):
    parser = all_parsers
    data = """A,B,C,D
foo,one,2,3,4,5
foo,two,7,8,9,10
foo,three,12,13,14,15
bar,one,12,13,14,15
bar,two,12,13,14,15
"""

    result = parser.read_csv(StringIO(data))
    expected = DataFrame(
        [
            [2, 3, 4, 5],
            [7, 8, 9, 10],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
            [12, 13, 14, 15],
        ],
        columns=["A", "B", "C", "D"],
        index=MultiIndex.from_tuples(
            [
                ("foo", "one"),
                ("foo", "two"),
                ("foo", "three"),
                ("bar", "one"),
                ("bar", "two"),
            ]
        ),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data,expected,header",
    [
        ("a,b", DataFrame(columns=["a", "b"]), [0]),
        (
            "a,b\nc,d",
            DataFrame(columns=MultiIndex.from_tuples([("a", "c"), ("b", "d")])),
            [0, 1],
        ),
    ],
)
@pytest.mark.parametrize("round_trip", [True, False])
def test_multi_index_blank_df(all_parsers, data, expected, header, round_trip):
    # see gh-14545
    parser = all_parsers
    data = expected.to_csv(index=False) if round_trip else data

    result = parser.read_csv(StringIO(data), header=header)
    tm.assert_frame_equal(result, expected)


def test_no_unnamed_index(all_parsers):
    parser = all_parsers
    data = """ id c0 c1 c2
0 1 0 a b
1 2 0 c d
2 2 2 e f
"""
    result = parser.read_csv(StringIO(data), sep=" ")
    expected = DataFrame(
        [[0, 1, 0, "a", "b"], [1, 2, 0, "c", "d"], [2, 2, 2, "e", "f"]],
        columns=["Unnamed: 0", "id", "c0", "c1", "c2"],
    )
    tm.assert_frame_equal(result, expected)


def test_read_csv_parse_simple_list(all_parsers):
    parser = all_parsers
    data = """foo
bar baz
qux foo
foo
bar"""

    result = parser.read_csv(StringIO(data), header=None)
    expected = DataFrame(["foo", "bar baz", "qux foo", "foo", "bar"])
    tm.assert_frame_equal(result, expected)


@tm.network
def test_url(all_parsers, csv_dir_path):
    # TODO: FTP testing
    parser = all_parsers
    kwargs = dict(sep="\t")

    url = (
        "https://raw.github.com/pandas-dev/pandas/master/"
        "pandas/tests/io/parser/data/salaries.csv"
    )
    url_result = parser.read_csv(url, **kwargs)

    local_path = os.path.join(csv_dir_path, "salaries.csv")
    local_result = parser.read_csv(local_path, **kwargs)
    tm.assert_frame_equal(url_result, local_result)


@pytest.mark.slow
def test_local_file(all_parsers, csv_dir_path):
    parser = all_parsers
    kwargs = dict(sep="\t")

    local_path = os.path.join(csv_dir_path, "salaries.csv")
    local_result = parser.read_csv(local_path, **kwargs)
    url = "file://localhost/" + local_path

    try:
        url_result = parser.read_csv(url, **kwargs)
        tm.assert_frame_equal(url_result, local_result)
    except URLError:
        # Fails on some systems.
        pytest.skip("Failing on: " + " ".join(platform.uname()))


def test_path_path_lib(all_parsers):
    parser = all_parsers
    df = tm.makeDataFrame()
    result = tm.round_trip_pathlib(df.to_csv, lambda p: parser.read_csv(p, index_col=0))
    tm.assert_frame_equal(df, result)


def test_path_local_path(all_parsers):
    parser = all_parsers
    df = tm.makeDataFrame()
    result = tm.round_trip_localpath(
        df.to_csv, lambda p: parser.read_csv(p, index_col=0)
    )
    tm.assert_frame_equal(df, result)


def test_nonexistent_path(all_parsers):
    # gh-2428: pls no segfault
    # gh-14086: raise more helpful FileNotFoundError
    parser = all_parsers
    path = "{}.csv".format(tm.rands(10))

    msg = "does not exist" if parser.engine == "c" else r"\[Errno 2\]"
    with pytest.raises(FileNotFoundError, match=msg) as e:
        parser.read_csv(path)

        filename = e.value.filename
        filename = filename.decode() if isinstance(filename, bytes) else filename

        assert path == filename


def test_missing_trailing_delimiters(all_parsers):
    parser = all_parsers
    data = """A,B,C,D
1,2,3,4
1,3,3,
1,4,5"""

    result = parser.read_csv(StringIO(data))
    expected = DataFrame(
        [[1, 2, 3, 4], [1, 3, 3, np.nan], [1, 4, 5, np.nan]],
        columns=["A", "B", "C", "D"],
    )
    tm.assert_frame_equal(result, expected)


def test_skip_initial_space(all_parsers):
    data = (
        '"09-Apr-2012", "01:10:18.300", 2456026.548822908, 12849, '
        "1.00361,  1.12551, 330.65659, 0355626618.16711,  73.48821, "
        "314.11625,  1917.09447,   179.71425,  80.000, 240.000, -350,  "
        "70.06056, 344.98370, 1,   1, -0.689265, -0.692787,  "
        "0.212036,    14.7674,   41.605,   -9999.0,   -9999.0,   "
        "-9999.0,   -9999.0,   -9999.0,  -9999.0, 000, 012, 128"
    )
    parser = all_parsers

    result = parser.read_csv(
        StringIO(data),
        names=list(range(33)),
        header=None,
        na_values=["-9999.0"],
        skipinitialspace=True,
    )
    expected = DataFrame(
        [
            [
                "09-Apr-2012",
                "01:10:18.300",
                2456026.548822908,
                12849,
                1.00361,
                1.12551,
                330.65659,
                355626618.16711,
                73.48821,
                314.11625,
                1917.09447,
                179.71425,
                80.0,
                240.0,
                -350,
                70.06056,
                344.9837,
                1,
                1,
                -0.689265,
                -0.692787,
                0.212036,
                14.7674,
                41.605,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                0,
                12,
                128,
            ]
        ]
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("sep", [",", "\t"])
@pytest.mark.parametrize("encoding", ["utf-16", "utf-16le", "utf-16be"])
def test_utf16_bom_skiprows(all_parsers, sep, encoding):
    # see gh-2298
    parser = all_parsers
    data = """skip this
skip this too
A,B,C
1,2,3
4,5,6""".replace(
        ",", sep
    )
    path = "__{}__.csv".format(tm.rands(10))
    kwargs = dict(sep=sep, skiprows=2)
    utf8 = "utf-8"

    with tm.ensure_clean(path) as path:
        from io import TextIOWrapper

        bytes_data = data.encode(encoding)

        with open(path, "wb") as f:
            f.write(bytes_data)

        bytes_buffer = BytesIO(data.encode(utf8))
        bytes_buffer = TextIOWrapper(bytes_buffer, encoding=utf8)

        result = parser.read_csv(path, encoding=encoding, **kwargs)
        expected = parser.read_csv(bytes_buffer, encoding=utf8, **kwargs)

        bytes_buffer.close()
        tm.assert_frame_equal(result, expected)


def test_utf16_example(all_parsers, csv_dir_path):
    path = os.path.join(csv_dir_path, "utf16_ex.txt")
    parser = all_parsers
    result = parser.read_csv(path, encoding="utf-16", sep="\t")
    assert len(result) == 50


def test_unicode_encoding(all_parsers, csv_dir_path):
    path = os.path.join(csv_dir_path, "unicode_series.csv")
    parser = all_parsers

    result = parser.read_csv(path, header=None, encoding="latin-1")
    result = result.set_index(0)
    got = result[1][1632]

    expected = "\xc1 k\xf6ldum klaka (Cold Fever) (1994)"
    assert got == expected


def test_trailing_delimiters(all_parsers):
    # see gh-2442
    data = """A,B,C
1,2,3,
4,5,6,
7,8,9,"""
    parser = all_parsers
    result = parser.read_csv(StringIO(data), index_col=False)

    expected = DataFrame({"A": [1, 4, 7], "B": [2, 5, 8], "C": [3, 6, 9]})
    tm.assert_frame_equal(result, expected)


def test_escapechar(all_parsers):
    # http://stackoverflow.com/questions/13824840/feature-request-for-
    # pandas-read-csv
    data = '''SEARCH_TERM,ACTUAL_URL
"bra tv bord","http://www.ikea.com/se/sv/catalog/categories/departments/living_room/10475/?se%7cps%7cnonbranded%7cvardagsrum%7cgoogle%7ctv_bord"
"tv p\xc3\xa5 hjul","http://www.ikea.com/se/sv/catalog/categories/departments/living_room/10475/?se%7cps%7cnonbranded%7cvardagsrum%7cgoogle%7ctv_bord"
"SLAGBORD, \\"Bergslagen\\", IKEA:s 1700-tals serie","http://www.ikea.com/se/sv/catalog/categories/departments/living_room/10475/?se%7cps%7cnonbranded%7cvardagsrum%7cgoogle%7ctv_bord"'''  # noqa

    parser = all_parsers
    result = parser.read_csv(
        StringIO(data), escapechar="\\", quotechar='"', encoding="utf-8"
    )

    assert result["SEARCH_TERM"][2] == (
        'SLAGBORD, "Bergslagen", ' "IKEA:s 1700-tals serie"
    )
    tm.assert_index_equal(result.columns, Index(["SEARCH_TERM", "ACTUAL_URL"]))


def test_int64_min_issues(all_parsers):
    # see gh-2599
    parser = all_parsers
    data = "A,B\n0,0\n0,"
    result = parser.read_csv(StringIO(data))

    expected = DataFrame({"A": [0, 0], "B": [0, np.nan]})
    tm.assert_frame_equal(result, expected)


def test_parse_integers_above_fp_precision(all_parsers):
    data = """Numbers
17007000002000191
17007000002000191
17007000002000191
17007000002000191
17007000002000192
17007000002000192
17007000002000192
17007000002000192
17007000002000192
17007000002000194"""
    parser = all_parsers
    result = parser.read_csv(StringIO(data))
    expected = DataFrame(
        {
            "Numbers": [
                17007000002000191,
                17007000002000191,
                17007000002000191,
                17007000002000191,
                17007000002000192,
                17007000002000192,
                17007000002000192,
                17007000002000192,
                17007000002000192,
                17007000002000194,
            ]
        }
    )
    tm.assert_frame_equal(result, expected)


def test_chunks_have_consistent_numerical_type(all_parsers):
    parser = all_parsers
    integers = [str(i) for i in range(499999)]
    data = "a\n" + "\n".join(integers + ["1.0", "2.0"] + integers)

    # Coercions should work without warnings.
    with tm.assert_produces_warning(None):
        result = parser.read_csv(StringIO(data))

    assert type(result.a[0]) is np.float64
    assert result.a.dtype == np.float


def test_warn_if_chunks_have_mismatched_type(all_parsers):
    warning_type = None
    parser = all_parsers
    integers = [str(i) for i in range(499999)]
    data = "a\n" + "\n".join(integers + ["a", "b"] + integers)

    # see gh-3866: if chunks are different types and can't
    # be coerced using numerical types, then issue warning.
    if parser.engine == "c" and parser.low_memory:
        warning_type = DtypeWarning

    with tm.assert_produces_warning(warning_type):
        df = parser.read_csv(StringIO(data))
    assert df.a.dtype == np.object


@pytest.mark.parametrize("sep", [" ", r"\s+"])
def test_integer_overflow_bug(all_parsers, sep):
    # see gh-2601
    data = "65248E10 11\n55555E55 22\n"
    parser = all_parsers

    result = parser.read_csv(StringIO(data), header=None, sep=sep)
    expected = DataFrame([[6.5248e14, 11], [5.5555e59, 22]])
    tm.assert_frame_equal(result, expected)


def test_catch_too_many_names(all_parsers):
    # see gh-5156
    data = """\
1,2,3
4,,6
7,8,9
10,11,12\n"""
    parser = all_parsers
    msg = (
        "Too many columns specified: expected 4 and found 3"
        if parser.engine == "c"
        else "Number of passed names did not match "
        "number of header fields in the file"
    )

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), header=0, names=["a", "b", "c", "d"])


def test_ignore_leading_whitespace(all_parsers):
    # see gh-3374, gh-6607
    parser = all_parsers
    data = " a b c\n 1 2 3\n 4 5 6\n 7 8 9"
    result = parser.read_csv(StringIO(data), sep=r"\s+")

    expected = DataFrame({"a": [1, 4, 7], "b": [2, 5, 8], "c": [3, 6, 9]})
    tm.assert_frame_equal(result, expected)


def test_chunk_begins_with_newline_whitespace(all_parsers):
    # see gh-10022
    parser = all_parsers
    data = "\n hello\nworld\n"

    result = parser.read_csv(StringIO(data), header=None)
    expected = DataFrame([" hello", "world"])
    tm.assert_frame_equal(result, expected)


def test_empty_with_index(all_parsers):
    # see gh-10184
    data = "x,y"
    parser = all_parsers
    result = parser.read_csv(StringIO(data), index_col=0)

    expected = DataFrame(columns=["y"], index=Index([], name="x"))
    tm.assert_frame_equal(result, expected)


def test_empty_with_multi_index(all_parsers):
    # see gh-10467
    data = "x,y,z"
    parser = all_parsers
    result = parser.read_csv(StringIO(data), index_col=["x", "y"])

    expected = DataFrame(
        columns=["z"], index=MultiIndex.from_arrays([[]] * 2, names=["x", "y"])
    )
    tm.assert_frame_equal(result, expected)


def test_empty_with_reversed_multi_index(all_parsers):
    data = "x,y,z"
    parser = all_parsers
    result = parser.read_csv(StringIO(data), index_col=[1, 0])

    expected = DataFrame(
        columns=["z"], index=MultiIndex.from_arrays([[]] * 2, names=["y", "x"])
    )
    tm.assert_frame_equal(result, expected)


def test_float_parser(all_parsers):
    # see gh-9565
    parser = all_parsers
    data = "45e-1,4.5,45.,inf,-inf"
    result = parser.read_csv(StringIO(data), header=None)

    expected = DataFrame([[float(s) for s in data.split(",")]])
    tm.assert_frame_equal(result, expected)


def test_scientific_no_exponent(all_parsers):
    # see gh-12215
    df = DataFrame.from_dict(
        OrderedDict([("w", ["2e"]), ("x", ["3E"]), ("y", ["42e"]), ("z", ["632E"])])
    )
    data = df.to_csv(index=False)
    parser = all_parsers

    for precision in parser.float_precision_choices:
        df_roundtrip = parser.read_csv(StringIO(data), float_precision=precision)
        tm.assert_frame_equal(df_roundtrip, df)


@pytest.mark.parametrize("conv", [None, np.int64, np.uint64])
def test_int64_overflow(all_parsers, conv):
    data = """ID
00013007854817840016671868
00013007854817840016749251
00013007854817840016754630
00013007854817840016781876
00013007854817840017028824
00013007854817840017963235
00013007854817840018860166"""
    parser = all_parsers

    if conv is None:
        # 13007854817840016671868 > UINT64_MAX, so this
        # will overflow and return object as the dtype.
        result = parser.read_csv(StringIO(data))
        expected = DataFrame(
            [
                "00013007854817840016671868",
                "00013007854817840016749251",
                "00013007854817840016754630",
                "00013007854817840016781876",
                "00013007854817840017028824",
                "00013007854817840017963235",
                "00013007854817840018860166",
            ],
            columns=["ID"],
        )
        tm.assert_frame_equal(result, expected)
    else:
        # 13007854817840016671868 > UINT64_MAX, so attempts
        # to cast to either int64 or uint64 will result in
        # an OverflowError being raised.
        msg = (
            "(Python int too large to convert to C long)|"
            "(long too big to convert)|"
            "(int too big to convert)"
        )

        with pytest.raises(OverflowError, match=msg):
            parser.read_csv(StringIO(data), converters={"ID": conv})


@pytest.mark.parametrize(
    "val", [np.iinfo(np.uint64).max, np.iinfo(np.int64).max, np.iinfo(np.int64).min]
)
def test_int64_uint64_range(all_parsers, val):
    # These numbers fall right inside the int64-uint64
    # range, so they should be parsed as string.
    parser = all_parsers
    result = parser.read_csv(StringIO(str(val)), header=None)

    expected = DataFrame([val])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "val", [np.iinfo(np.uint64).max + 1, np.iinfo(np.int64).min - 1]
)
def test_outside_int64_uint64_range(all_parsers, val):
    # These numbers fall just outside the int64-uint64
    # range, so they should be parsed as string.
    parser = all_parsers
    result = parser.read_csv(StringIO(str(val)), header=None)

    expected = DataFrame([str(val)])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("exp_data", [[str(-1), str(2 ** 63)], [str(2 ** 63), str(-1)]])
def test_numeric_range_too_wide(all_parsers, exp_data):
    # No numerical dtype can hold both negative and uint64
    # values, so they should be cast as string.
    parser = all_parsers
    data = "\n".join(exp_data)
    expected = DataFrame(exp_data)

    result = parser.read_csv(StringIO(data), header=None)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("iterator", [True, False])
def test_empty_with_nrows_chunksize(all_parsers, iterator):
    # see gh-9535
    parser = all_parsers
    expected = DataFrame(columns=["foo", "bar"])

    nrows = 10
    data = StringIO("foo,bar\n")

    if iterator:
        result = next(iter(parser.read_csv(data, chunksize=nrows)))
    else:
        result = parser.read_csv(data, nrows=nrows)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data,kwargs,expected,msg",
    [
        # gh-10728: WHITESPACE_LINE
        (
            "a,b,c\n4,5,6\n ",
            dict(),
            DataFrame([[4, 5, 6]], columns=["a", "b", "c"]),
            None,
        ),
        # gh-10548: EAT_LINE_COMMENT
        (
            "a,b,c\n4,5,6\n#comment",
            dict(comment="#"),
            DataFrame([[4, 5, 6]], columns=["a", "b", "c"]),
            None,
        ),
        # EAT_CRNL_NOP
        (
            "a,b,c\n4,5,6\n\r",
            dict(),
            DataFrame([[4, 5, 6]], columns=["a", "b", "c"]),
            None,
        ),
        # EAT_COMMENT
        (
            "a,b,c\n4,5,6#comment",
            dict(comment="#"),
            DataFrame([[4, 5, 6]], columns=["a", "b", "c"]),
            None,
        ),
        # SKIP_LINE
        (
            "a,b,c\n4,5,6\nskipme",
            dict(skiprows=[2]),
            DataFrame([[4, 5, 6]], columns=["a", "b", "c"]),
            None,
        ),
        # EAT_LINE_COMMENT
        (
            "a,b,c\n4,5,6\n#comment",
            dict(comment="#", skip_blank_lines=False),
            DataFrame([[4, 5, 6]], columns=["a", "b", "c"]),
            None,
        ),
        # IN_FIELD
        (
            "a,b,c\n4,5,6\n ",
            dict(skip_blank_lines=False),
            DataFrame([["4", 5, 6], [" ", None, None]], columns=["a", "b", "c"]),
            None,
        ),
        # EAT_CRNL
        (
            "a,b,c\n4,5,6\n\r",
            dict(skip_blank_lines=False),
            DataFrame([[4, 5, 6], [None, None, None]], columns=["a", "b", "c"]),
            None,
        ),
        # ESCAPED_CHAR
        (
            "a,b,c\n4,5,6\n\\",
            dict(escapechar="\\"),
            None,
            "(EOF following escape character)|(unexpected end of data)",
        ),
        # ESCAPE_IN_QUOTED_FIELD
        (
            'a,b,c\n4,5,6\n"\\',
            dict(escapechar="\\"),
            None,
            "(EOF inside string starting at row 2)|(unexpected end of data)",
        ),
        # IN_QUOTED_FIELD
        (
            'a,b,c\n4,5,6\n"',
            dict(escapechar="\\"),
            None,
            "(EOF inside string starting at row 2)|(unexpected end of data)",
        ),
    ],
    ids=[
        "whitespace-line",
        "eat-line-comment",
        "eat-crnl-nop",
        "eat-comment",
        "skip-line",
        "eat-line-comment",
        "in-field",
        "eat-crnl",
        "escaped-char",
        "escape-in-quoted-field",
        "in-quoted-field",
    ],
)
def test_eof_states(all_parsers, data, kwargs, expected, msg):
    # see gh-10728, gh-10548
    parser = all_parsers

    if expected is None:
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), **kwargs)
    else:
        result = parser.read_csv(StringIO(data), **kwargs)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("usecols", [None, [0, 1], ["a", "b"]])
def test_uneven_lines_with_usecols(all_parsers, usecols):
    # see gh-12203
    parser = all_parsers
    data = r"""a,b,c
0,1,2
3,4,5,6,7
8,9,10"""

    if usecols is None:
        # Make sure that an error is still raised
        # when the "usecols" parameter is not provided.
        msg = r"Expected \d+ fields in line \d+, saw \d+"
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data))
    else:
        expected = DataFrame({"a": [0, 3, 8], "b": [1, 4, 9]})

        result = parser.read_csv(StringIO(data), usecols=usecols)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data,kwargs,expected",
    [
        # First, check to see that the response of parser when faced with no
        # provided columns raises the correct error, with or without usecols.
        ("", dict(), None),
        ("", dict(usecols=["X"]), None),
        (
            ",,",
            dict(names=["Dummy", "X", "Dummy_2"], usecols=["X"]),
            DataFrame(columns=["X"], index=[0], dtype=np.float64),
        ),
        (
            "",
            dict(names=["Dummy", "X", "Dummy_2"], usecols=["X"]),
            DataFrame(columns=["X"]),
        ),
    ],
)
def test_read_empty_with_usecols(all_parsers, data, kwargs, expected):
    # see gh-12493
    parser = all_parsers

    if expected is None:
        msg = "No columns to parse from file"
        with pytest.raises(EmptyDataError, match=msg):
            parser.read_csv(StringIO(data), **kwargs)
    else:
        result = parser.read_csv(StringIO(data), **kwargs)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs,expected",
    [
        # gh-8661, gh-8679: this should ignore six lines, including
        # lines with trailing whitespace and blank lines.
        (
            dict(
                header=None,
                delim_whitespace=True,
                skiprows=[0, 1, 2, 3, 5, 6],
                skip_blank_lines=True,
            ),
            DataFrame([[1.0, 2.0, 4.0], [5.1, np.nan, 10.0]]),
        ),
        # gh-8983: test skipping set of rows after a row with trailing spaces.
        (
            dict(
                delim_whitespace=True, skiprows=[1, 2, 3, 5, 6], skip_blank_lines=True
            ),
            DataFrame({"A": [1.0, 5.1], "B": [2.0, np.nan], "C": [4.0, 10]}),
        ),
    ],
)
def test_trailing_spaces(all_parsers, kwargs, expected):
    data = "A B C  \nrandom line with trailing spaces    \nskip\n1,2,3\n1,2.,4.\nrandom line with trailing tabs\t\t\t\n   \n5.1,NaN,10.0\n"  # noqa
    parser = all_parsers

    result = parser.read_csv(StringIO(data.replace(",", "  ")), **kwargs)
    tm.assert_frame_equal(result, expected)


def test_raise_on_sep_with_delim_whitespace(all_parsers):
    # see gh-6607
    data = "a b c\n1 2 3"
    parser = all_parsers

    with pytest.raises(ValueError, match="you can only specify one"):
        parser.read_csv(StringIO(data), sep=r"\s", delim_whitespace=True)


@pytest.mark.parametrize("delim_whitespace", [True, False])
def test_single_char_leading_whitespace(all_parsers, delim_whitespace):
    # see gh-9710
    parser = all_parsers
    data = """\
MyColumn
a
b
a
b\n"""

    expected = DataFrame({"MyColumn": list("abab")})
    result = parser.read_csv(
        StringIO(data), skipinitialspace=True, delim_whitespace=delim_whitespace
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "sep,skip_blank_lines,exp_data",
    [
        (",", True, [[1.0, 2.0, 4.0], [5.0, np.nan, 10.0], [-70.0, 0.4, 1.0]]),
        (r"\s+", True, [[1.0, 2.0, 4.0], [5.0, np.nan, 10.0], [-70.0, 0.4, 1.0]]),
        (
            ",",
            False,
            [
                [1.0, 2.0, 4.0],
                [np.nan, np.nan, np.nan],
                [np.nan, np.nan, np.nan],
                [5.0, np.nan, 10.0],
                [np.nan, np.nan, np.nan],
                [-70.0, 0.4, 1.0],
            ],
        ),
    ],
)
def test_empty_lines(all_parsers, sep, skip_blank_lines, exp_data):
    parser = all_parsers
    data = """\
A,B,C
1,2.,4.


5.,NaN,10.0

-70,.4,1
"""

    if sep == r"\s+":
        data = data.replace(",", "  ")

    result = parser.read_csv(StringIO(data), sep=sep, skip_blank_lines=skip_blank_lines)
    expected = DataFrame(exp_data, columns=["A", "B", "C"])
    tm.assert_frame_equal(result, expected)


def test_whitespace_lines(all_parsers):
    parser = all_parsers
    data = """

\t  \t\t
\t
A,B,C
\t    1,2.,4.
5.,NaN,10.0
"""
    expected = DataFrame([[1, 2.0, 4.0], [5.0, np.nan, 10.0]], columns=["A", "B", "C"])
    result = parser.read_csv(StringIO(data))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data,expected",
    [
        (
            """   A   B   C   D
a   1   2   3   4
b   1   2   3   4
c   1   2   3   4
""",
            DataFrame(
                [[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                columns=["A", "B", "C", "D"],
                index=["a", "b", "c"],
            ),
        ),
        (
            "    a b c\n1 2 3 \n4 5  6\n 7 8 9",
            DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["a", "b", "c"]),
        ),
    ],
)
def test_whitespace_regex_separator(all_parsers, data, expected):
    # see gh-6607
    parser = all_parsers
    result = parser.read_csv(StringIO(data), sep=r"\s+")
    tm.assert_frame_equal(result, expected)


def test_verbose_read(all_parsers, capsys):
    parser = all_parsers
    data = """a,b,c,d
one,1,2,3
one,1,2,3
,1,2,3
one,1,2,3
,1,2,3
,1,2,3
one,1,2,3
two,1,2,3"""

    # Engines are verbose in different ways.
    parser.read_csv(StringIO(data), verbose=True)
    captured = capsys.readouterr()

    if parser.engine == "c":
        assert "Tokenization took:" in captured.out
        assert "Parser memory cleanup took:" in captured.out
    else:  # Python engine
        assert captured.out == "Filled 3 NA values in column a\n"


def test_verbose_read2(all_parsers, capsys):
    parser = all_parsers
    data = """a,b,c,d
one,1,2,3
two,1,2,3
three,1,2,3
four,1,2,3
five,1,2,3
,1,2,3
seven,1,2,3
eight,1,2,3"""

    parser.read_csv(StringIO(data), verbose=True, index_col=0)
    captured = capsys.readouterr()

    # Engines are verbose in different ways.
    if parser.engine == "c":
        assert "Tokenization took:" in captured.out
        assert "Parser memory cleanup took:" in captured.out
    else:  # Python engine
        assert captured.out == "Filled 1 NA values in column a\n"


def test_iteration_open_handle(all_parsers):
    parser = all_parsers
    kwargs = dict(squeeze=True, header=None)

    with tm.ensure_clean() as path:
        with open(path, "w") as f:
            f.write("AAA\nBBB\nCCC\nDDD\nEEE\nFFF\nGGG")

        with open(path, "r") as f:
            for line in f:
                if "CCC" in line:
                    break

            result = parser.read_csv(f, **kwargs)
            expected = Series(["DDD", "EEE", "FFF", "GGG"], name=0)
            tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "data,thousands,decimal",
    [
        (
            """A|B|C
1|2,334.01|5
10|13|10.
""",
            ",",
            ".",
        ),
        (
            """A|B|C
1|2.334,01|5
10|13|10,
""",
            ".",
            ",",
        ),
    ],
)
def test_1000_sep_with_decimal(all_parsers, data, thousands, decimal):
    parser = all_parsers
    expected = DataFrame({"A": [1, 10], "B": [2334.01, 13], "C": [5, 10.0]})

    result = parser.read_csv(
        StringIO(data), sep="|", thousands=thousands, decimal=decimal
    )
    tm.assert_frame_equal(result, expected)


def test_euro_decimal_format(all_parsers):
    parser = all_parsers
    data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""

    result = parser.read_csv(StringIO(data), sep=";", decimal=",")
    expected = DataFrame(
        [
            [1, 1521.1541, 187101.9543, "ABC", "poi", 4.738797819],
            [2, 121.12, 14897.76, "DEF", "uyt", 0.377320872],
            [3, 878.158, 108013.434, "GHI", "rez", 2.735694704],
        ],
        columns=["Id", "Number1", "Number2", "Text1", "Text2", "Number3"],
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("na_filter", [True, False])
def test_inf_parsing(all_parsers, na_filter):
    parser = all_parsers
    data = """\
,A
a,inf
b,-inf
c,+Inf
d,-Inf
e,INF
f,-INF
g,+INf
h,-INf
i,inF
j,-inF"""
    expected = DataFrame(
        {"A": [float("inf"), float("-inf")] * 5},
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
    )
    result = parser.read_csv(StringIO(data), index_col=0, na_filter=na_filter)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("nrows", [0, 1, 2, 3, 4, 5])
def test_raise_on_no_columns(all_parsers, nrows):
    parser = all_parsers
    data = "\n" * nrows

    msg = "No columns to parse from file"
    with pytest.raises(EmptyDataError, match=msg):
        parser.read_csv(StringIO(data))


def test_memory_map(all_parsers, csv_dir_path):
    mmap_file = os.path.join(csv_dir_path, "test_mmap.csv")
    parser = all_parsers

    expected = DataFrame(
        {"a": [1, 2, 3], "b": ["one", "two", "three"], "c": ["I", "II", "III"]}
    )

    result = parser.read_csv(mmap_file, memory_map=True)
    tm.assert_frame_equal(result, expected)


def test_null_byte_char(all_parsers):
    # see gh-2741
    data = "\x00,foo"
    names = ["a", "b"]
    parser = all_parsers

    if parser.engine == "c":
        expected = DataFrame([[np.nan, "foo"]], columns=names)
        out = parser.read_csv(StringIO(data), names=names)
        tm.assert_frame_equal(out, expected)
    else:
        if compat.PY38:
            msg = "line contains NUL"
        else:
            msg = "NULL byte detected"
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), names=names)


@pytest.mark.parametrize(
    "data,kwargs,expected",
    [
        # Basic test
        ("a\n1", dict(), DataFrame({"a": [1]})),
        # "Regular" quoting
        ('"a"\n1', dict(quotechar='"'), DataFrame({"a": [1]})),
        # Test in a data row instead of header
        ("b\n1", dict(names=["a"]), DataFrame({"a": ["b", "1"]})),
        # Test in empty data row with skipping
        ("\n1", dict(names=["a"], skip_blank_lines=True), DataFrame({"a": [1]})),
        # Test in empty data row without skipping
        (
            "\n1",
            dict(names=["a"], skip_blank_lines=False),
            DataFrame({"a": [np.nan, 1]}),
        ),
    ],
)
def test_utf8_bom(all_parsers, data, kwargs, expected):
    # see gh-4793
    parser = all_parsers
    bom = "\ufeff"
    utf8 = "utf-8"

    def _encode_data_with_bom(_data):
        bom_data = (bom + _data).encode(utf8)
        return BytesIO(bom_data)

    result = parser.read_csv(_encode_data_with_bom(data), encoding=utf8, **kwargs)
    tm.assert_frame_equal(result, expected)


def test_temporary_file(all_parsers):
    # see gh-13398
    parser = all_parsers
    data = "0 0"

    new_file = TemporaryFile("w+")
    new_file.write(data)
    new_file.flush()
    new_file.seek(0)

    result = parser.read_csv(new_file, sep=r"\s+", header=None)
    new_file.close()

    expected = DataFrame([[0, 0]])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("byte", [8, 16])
@pytest.mark.parametrize("fmt", ["utf-{0}", "utf_{0}", "UTF-{0}", "UTF_{0}"])
def test_read_csv_utf_aliases(all_parsers, byte, fmt):
    # see gh-13549
    expected = DataFrame({"mb_num": [4.8], "multibyte": ["test"]})
    parser = all_parsers

    encoding = fmt.format(byte)
    data = "mb_num,multibyte\n4.8,test".encode(encoding)

    result = parser.read_csv(BytesIO(data), encoding=encoding)
    tm.assert_frame_equal(result, expected)


def test_internal_eof_byte(all_parsers):
    # see gh-5500
    parser = all_parsers
    data = "a,b\n1\x1a,2"

    expected = DataFrame([["1\x1a", 2]], columns=["a", "b"])
    result = parser.read_csv(StringIO(data))
    tm.assert_frame_equal(result, expected)


def test_internal_eof_byte_to_file(all_parsers):
    # see gh-16559
    parser = all_parsers
    data = b'c1,c2\r\n"test \x1a    test", test\r\n'
    expected = DataFrame([["test \x1a    test", " test"]], columns=["c1", "c2"])
    path = "__{}__.csv".format(tm.rands(10))

    with tm.ensure_clean(path) as path:
        with open(path, "wb") as f:
            f.write(data)

        result = parser.read_csv(path)
        tm.assert_frame_equal(result, expected)


def test_sub_character(all_parsers, csv_dir_path):
    # see gh-16893
    filename = os.path.join(csv_dir_path, "sub_char.csv")
    expected = DataFrame([[1, 2, 3]], columns=["a", "\x1ab", "c"])

    parser = all_parsers
    result = parser.read_csv(filename)
    tm.assert_frame_equal(result, expected)


def test_file_handle_string_io(all_parsers):
    # gh-14418
    #
    # Don't close user provided file handles.
    parser = all_parsers
    data = "a,b\n1,2"

    fh = StringIO(data)
    parser.read_csv(fh)
    assert not fh.closed


def test_file_handles_with_open(all_parsers, csv1):
    # gh-14418
    #
    # Don't close user provided file handles.
    parser = all_parsers

    with open(csv1, "r") as f:
        parser.read_csv(f)
        assert not f.closed


def test_invalid_file_buffer_class(all_parsers):
    # see gh-15337
    class InvalidBuffer:
        pass

    parser = all_parsers
    msg = "Invalid file path or buffer object type"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(InvalidBuffer())


def test_invalid_file_buffer_mock(all_parsers):
    # see gh-15337
    parser = all_parsers
    msg = "Invalid file path or buffer object type"

    class Foo:
        pass

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(Foo())


def test_valid_file_buffer_seems_invalid(all_parsers):
    # gh-16135: we want to ensure that "tell" and "seek"
    # aren't actually being used when we call `read_csv`
    #
    # Thus, while the object may look "invalid" (these
    # methods are attributes of the `StringIO` class),
    # it is still a valid file-object for our purposes.
    class NoSeekTellBuffer(StringIO):
        def tell(self):
            raise AttributeError("No tell method")

        def seek(self, pos, whence=0):
            raise AttributeError("No seek method")

    data = "a\n1"
    parser = all_parsers
    expected = DataFrame({"a": [1]})

    result = parser.read_csv(NoSeekTellBuffer(data))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs",
    [dict(), dict(error_bad_lines=True)],  # Default is True.  # Explicitly pass in.
)
@pytest.mark.parametrize(
    "warn_kwargs", [dict(), dict(warn_bad_lines=True), dict(warn_bad_lines=False)]
)
def test_error_bad_lines(all_parsers, kwargs, warn_kwargs):
    # see gh-15925
    parser = all_parsers
    kwargs.update(**warn_kwargs)
    data = "a\n1\n1,2,3\n4\n5,6,7"

    msg = "Expected 1 fields in line 3, saw 3"
    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), **kwargs)


def test_warn_bad_lines(all_parsers, capsys):
    # see gh-15925
    parser = all_parsers
    data = "a\n1\n1,2,3\n4\n5,6,7"
    expected = DataFrame({"a": [1, 4]})

    result = parser.read_csv(StringIO(data), error_bad_lines=False, warn_bad_lines=True)
    tm.assert_frame_equal(result, expected)

    captured = capsys.readouterr()
    assert "Skipping line 3" in captured.err
    assert "Skipping line 5" in captured.err


def test_suppress_error_output(all_parsers, capsys):
    # see gh-15925
    parser = all_parsers
    data = "a\n1\n1,2,3\n4\n5,6,7"
    expected = DataFrame({"a": [1, 4]})

    result = parser.read_csv(
        StringIO(data), error_bad_lines=False, warn_bad_lines=False
    )
    tm.assert_frame_equal(result, expected)

    captured = capsys.readouterr()
    assert captured.err == ""


@pytest.mark.skipif(
    compat.is_platform_windows() and not compat.PY36,
    reason="On Python < 3.6 won't pass on Windows",
)
@pytest.mark.parametrize("filename", ["sé-es-vé.csv", "ru-sй.csv"])
def test_filename_with_special_chars(all_parsers, filename):
    # see gh-15086.
    parser = all_parsers
    df = DataFrame({"a": [1, 2, 3]})

    with tm.ensure_clean(filename) as path:
        df.to_csv(path, index=False)

        result = parser.read_csv(path)
        tm.assert_frame_equal(result, df)


def test_read_csv_memory_growth_chunksize(all_parsers):
    # see gh-24805
    #
    # Let's just make sure that we don't crash
    # as we iteratively process all chunks.
    parser = all_parsers

    with tm.ensure_clean() as path:
        with open(path, "w") as f:
            for i in range(1000):
                f.write(str(i) + "\n")

        result = parser.read_csv(path, chunksize=20)

        for _ in result:
            pass


def test_read_table_equivalency_to_read_csv(all_parsers):
    # see gh-21948
    # As of 0.25.0, read_table is undeprecated
    parser = all_parsers
    data = "a\tb\n1\t2\n3\t4"
    expected = parser.read_csv(StringIO(data), sep="\t")
    result = parser.read_table(StringIO(data))
    tm.assert_frame_equal(result, expected)


def test_first_row_bom(all_parsers):
    # see gh-26545
    parser = all_parsers
    data = '''\ufeff"Head1"	"Head2"	"Head3"'''

    result = parser.read_csv(StringIO(data), delimiter="\t")
    expected = DataFrame(columns=["Head1", "Head2", "Head3"])
    tm.assert_frame_equal(result, expected)
