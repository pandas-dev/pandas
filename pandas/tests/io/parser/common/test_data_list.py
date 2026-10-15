"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""

import csv
from io import StringIO

import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.parsers import TextParser

xfail_pyarrow = pytest.mark.usefixtures("pyarrow_xfail")


@xfail_pyarrow
def test_read_data_list(all_parsers):
    parser = all_parsers
    kwargs = {"index_col": 0}
    data = "A,B,C\nfoo,1,2,3\nbar,4,5,6"

    data_list = [["A", "B", "C"], ["foo", "1", "2", "3"], ["bar", "4", "5", "6"]]
    expected = parser.read_csv(StringIO(data), **kwargs)

    with TextParser(data_list, chunksize=2, **kwargs) as parser:
        result = parser.read()

    tm.assert_frame_equal(result, expected)


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
    kwargs = {"index_col": 0}

    lines = list(csv.reader(StringIO(data)))
    with TextParser(lines, chunksize=2, **kwargs) as reader:
        chunks = list(reader)

    expected = parser.read_csv(StringIO(data), **kwargs)

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
    kwargs = {"index_col": 0}

    lines = list(csv.reader(StringIO(data)))
    with TextParser(lines, chunksize=2, skiprows=[1], **kwargs) as reader:
        chunks = list(reader)

    expected = parser.read_csv(StringIO(data), **kwargs)

    tm.assert_frame_equal(chunks[0], expected[1:3])


def test_read_csv_parse_simple_list(all_parsers):
    parser = all_parsers
    data = """foo
bar baz
qux foo
foo
bar"""

    result = parser.read_csv(StringIO(data), header=None)
    expected = pd.DataFrame(["foo", "bar baz", "qux foo", "foo", "bar"])
    tm.assert_frame_equal(result, expected)


def test_skipfooter_counts_skiprows_lines_data_list():
    # GH#36827 a skiprows line inside the footer still counts towards skipfooter,
    # on the list-backed path read_excel uses
    data_list = [
        ["A", "B"],
        ["1", "2"],
        ["3", "4"],
        ["5", "6"],
        ["7", "8"],
        ["footer"],
    ]
    with TextParser(data_list, skiprows=[4], skipfooter=2) as parser:
        result = parser.read()

    expected = pd.DataFrame({"A": [1, 3, 5], "B": [2, 4, 6]})
    tm.assert_frame_equal(result, expected)


def test_skipfooter_with_explicit_row_count_is_positional():
    # GH#36827 with an explicit row count the parser has not reached the end of the
    # file, so there is no line count to trim against and the positional trim stands
    data_list = [["a", "b"], ["1", "2"], ["3", "4"], ["5", "6"], ["foot"]]
    with TextParser(data_list, skipfooter=1) as parser:
        result = parser.read(1)

    tm.assert_frame_equal(result, pd.DataFrame(columns=["a", "b"], dtype=object))
