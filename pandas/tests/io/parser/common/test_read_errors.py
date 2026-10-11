"""
Tests that work on the Python, C and PyArrow engines but do not have a
specific classification into the other test modules.
"""

import codecs
import csv
from io import StringIO
from pathlib import Path
import re

import numpy as np
import pytest

from pandas.errors import (
    EmptyDataError,
    ParserError,
    ParserWarning,
)

import pandas as pd
import pandas._testing as tm

xfail_pyarrow = pytest.mark.usefixtures("pyarrow_xfail")
skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")


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


@pytest.mark.parametrize("escapechar", ["", "xx"])
def test_bad_escapechar(all_parsers, escapechar):
    # GH#66317
    data = "a,b\n1,2"
    parser = all_parsers

    if parser.engine == "python":
        # message and exception type come from the csv module
        err = TypeError
        msg = '"escapechar" must be'
    else:
        err = ValueError
        msg = "Only length-1 escapes supported"

    with pytest.raises(err, match=msg):
        parser.read_csv(StringIO(data), escapechar=escapechar)


def test_empty_sep(pyarrow_parser_only):
    # GH#66317 previously leaked a cryptic error from pyarrow option
    #  conversion; the C engine has a similarly-cryptic error here and the
    #  python engine treats sep="" as a regex, so pyarrow-only.
    data = "a,b\n1,2"
    parser = pyarrow_parser_only

    msg = "Only length-1 delimiters supported"
    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), sep="")


def test_bad_stream_exception(all_parsers, datapath):
    # see gh-13652
    #
    # This test validates that both the Python engine and C engine will
    # raise UnicodeDecodeError instead of C engine raising ParserError
    # and swallowing the exception that caused read to fail.
    path = datapath("io", "parser", "data", "sauron.SHIFT_JIS.csv")
    codec = codecs.lookup("utf-8")
    utf8 = codecs.lookup("utf-8")
    parser = all_parsers
    msg = "'utf-8' codec can't decode byte"

    # Stream must be binary UTF8.
    with (
        open(path, "rb") as handle,
        codecs.StreamRecoder(
            handle, utf8.encode, utf8.decode, codec.streamreader, codec.streamwriter
        ) as stream,
    ):
        with pytest.raises(UnicodeDecodeError, match=msg):
            parser.read_csv(stream)


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
    err = ParserError
    if parser.engine == "pyarrow":
        msg = "The 'comment' option is not supported with the 'pyarrow' engine"
        err = ValueError
    with pytest.raises(err, match=msg):
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

    if parser.engine == "pyarrow":
        msg = "The 'iterator' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data),
                header=1,
                comment="#",
                iterator=True,
                chunksize=1,
                skiprows=[2],
            )
        return

    msg = "Expected 3 fields in line 6, saw 5"
    with parser.read_csv(
        StringIO(data), header=1, comment="#", iterator=True, chunksize=1, skiprows=[2]
    ) as reader:
        with pytest.raises(ParserError, match=msg):
            reader.read(nrows)


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
        else "Number of passed names did not match number of header fields in the file"
    )

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), header=0, names=["a", "b", "c", "d"])


@skip_pyarrow  # CSV parse error: Empty CSV file or block
@pytest.mark.parametrize("nrows", [0, 1, 2, 3, 4, 5])
def test_raise_on_no_columns(all_parsers, nrows):
    parser = all_parsers
    data = "\n" * nrows

    msg = "No columns to parse from file"
    with pytest.raises(EmptyDataError, match=msg):
        parser.read_csv(StringIO(data))


def test_unexpected_keyword_parameter_exception(all_parsers):
    # GH-34976
    parser = all_parsers

    msg = "{}\\(\\) got an unexpected keyword argument 'foo'"
    with pytest.raises(TypeError, match=msg.format("read_csv")):
        parser.read_csv("foo.csv", foo=1)
    with pytest.raises(TypeError, match=msg.format("read_table")):
        parser.read_table("foo.tsv", foo=1)


def test_suppress_error_output(all_parsers):
    # see gh-15925
    parser = all_parsers
    data = "a\n1\n1,2,3\n4\n5,6,7"
    expected = pd.DataFrame({"a": [1, 4]})

    result = parser.read_csv(StringIO(data), on_bad_lines="skip")
    tm.assert_frame_equal(result, expected)


def test_error_bad_lines(all_parsers):
    # see gh-15925
    parser = all_parsers
    data = "a\n1\n1,2,3\n4\n5,6,7"

    msg = "Expected 1 fields in line 3, saw 3"

    if parser.engine == "pyarrow":
        # "CSV parse error: Expected 1 columns, got 3: 1,2,3"
        pytest.skip(reason="https://github.com/apache/arrow/issues/38676")

    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), on_bad_lines="error")


def test_warn_bad_lines(all_parsers):
    # see gh-15925
    parser = all_parsers
    data = "a\n1\n1,2,3\n4\n5,6,7"
    expected = pd.DataFrame({"a": [1, 4]})
    match_msg = "Skipping line"

    expected_warning = ParserWarning
    if parser.engine == "pyarrow":
        match_msg = "Expected 1 columns, but found 3: 1,2,3"

    with tm.assert_produces_warning(
        expected_warning, match=match_msg, check_stacklevel=False
    ):
        result = parser.read_csv(StringIO(data), on_bad_lines="warn")
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_error_bad_lines_index_col_false(all_parsers):
    # GH#49279
    parser = all_parsers
    data = "a,b\n1,2\n3,4,5\n6,7"
    msg = "Expected 2 fields in line 3, saw 3"

    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), index_col=False, on_bad_lines="error")


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_skip_bad_lines_index_col_false(all_parsers):
    # GH#49279
    parser = all_parsers
    data = "a,b\n1,2\n3,4,5\n6,7"
    expected = pd.DataFrame({"a": [1, 6], "b": [2, 7]})

    result = parser.read_csv(StringIO(data), index_col=False, on_bad_lines="skip")
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_index_col_false_wider_rows_are_not_bad_lines(all_parsers):
    # GH#2442, GH#49279 - wider rows are data to truncate, not bad lines
    parser = all_parsers
    data = "A,B,C\n1,2,3,4\n5,6,7,8"
    expected = pd.DataFrame({"A": [1, 5], "B": [2, 6], "C": [3, 7]})

    result = parser.read_csv_check_warnings(
        ParserWarning,
        "Length of header or names does not match length of data",
        StringIO(data),
        index_col=False,
        on_bad_lines="error",
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_index_col_false_bad_line_counted_from_first_row(all_parsers):
    # GH#49279 - the first data row is wider than the header, so it is the width
    # the bad line is measured against and the one the message reports
    parser = all_parsers
    data = "A,B,C\n1,2,3,4\n5,6,7,8,9"

    with pytest.raises(ParserError, match="Expected 4 fields in line 3, saw 5"):
        parser.read_csv(StringIO(data), index_col=False, on_bad_lines="error")


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_index_col_false_bad_line_with_names(all_parsers):
    # GH#49279 - with names the first data row is already buffered by the time the
    # expected width is taken, so it must not be measured on the second row
    parser = all_parsers
    kwargs = {
        "index_col": False,
        "header": 1,
        "names": ["x", "y"],
        "on_bad_lines": "error",
    }
    expected = pd.DataFrame({"x": [1, 4, 6], "y": [2, 5, 7]})

    result = parser.read_csv_check_warnings(
        ParserWarning,
        "Length of header or names does not match length of data",
        StringIO("hdr\na,b\n1,2,3\n4,5\n6,7,8"),
        **kwargs,
    )
    tm.assert_frame_equal(result, expected)

    with pytest.raises(ParserError, match="Expected 3 fields in line 5, saw 4"):
        parser.read_csv(StringIO("hdr\na,b\n1,2,3\n4,5\n6,7,8,9"), **kwargs)


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_index_col_false_single_data_row_with_names(all_parsers):
    # GH#49279 - the only data row is already buffered and there is no second row
    # to fall back on, so the expected width still has to come from the buffer
    parser = all_parsers
    kwargs = {"index_col": False, "header": 0, "names": ["x", "y"]}
    expected = pd.DataFrame({"x": [1], "y": [2]})

    result = parser.read_csv_check_warnings(
        ParserWarning,
        "Length of header or names does not match length of data",
        StringIO("a,b\n1,2,3\n"),
        **kwargs,
    )
    tm.assert_frame_equal(result, expected)

    # GH#2442 - a trailing delimiter is truncated, not treated as a bad line
    result = parser.read_csv(StringIO("a,b\n1,2,\n"), **kwargs)
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # AttributeError: 'bool' object has no attribute 'copy'
def test_index_col_false_short_names_keep_header_width(all_parsers):
    # GH#49279 - names replaces the column count, but the header line is still what
    # a wider row is measured against, so this row is truncated rather than dropped
    parser = all_parsers
    expected = pd.DataFrame({"x": [1, 3], "y": [2, 4]})

    result = parser.read_csv_check_warnings(
        ParserWarning,
        "Length of header or names does not match length of data",
        StringIO("a,b,c\n1,2\n3,4,5\n"),
        names=["x", "y"],
        header=0,
        index_col=False,
        on_bad_lines="error",
    )
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

    if parser.engine == "pyarrow":
        # Expected 6 columns, got 7: 6,7,8,9,10,11,12
        pytest.skip(reason="https://github.com/apache/arrow/issues/38676")

    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data))


def test_null_byte_char(all_parsers):
    # see gh-2741
    data = "\x00,foo"
    names = ["a", "b"]
    parser = all_parsers

    if parser.engine in ["c", "python"]:
        # GH#19886 a lone NUL is a one-character value, not the empty string,
        # so it is not an na_value. The C engine used to compare the
        # NUL-terminated word and report NaN here.
        expected = pd.DataFrame([["\x00", "foo"]], columns=names)
        out = parser.read_csv(StringIO(data), names=names)
        tm.assert_frame_equal(out, expected)
    else:
        if parser.engine == "pyarrow":
            # CSV parse error: Empty CSV file or block: "
            # cannot infer number of columns"
            pytest.skip(reason="https://github.com/apache/arrow/issues/38676")
        else:
            msg = "NULL byte detected"
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), names=names)


def test_open_file(all_parsers, temp_file):
    # GH 39024
    parser = all_parsers

    msg = "Could not determine delimiter"
    err = csv.Error
    if parser.engine == "pyarrow":
        msg = "the 'pyarrow' engine does not support sep=None"
        err = ValueError

    file = Path(temp_file)
    file.write_bytes(b"\xe4\na\n1")

    with tm.assert_produces_warning(None):
        # should not trigger a ResourceWarning
        with pytest.raises(err, match=msg):
            parser.read_csv(file, sep=None, encoding_errors="replace")


def test_invalid_on_bad_line(all_parsers):
    parser = all_parsers
    data = "a\n1\n1,2,3\n4\n5,6,7"
    with pytest.raises(ValueError, match="Argument abc is invalid for on_bad_lines"):
        parser.read_csv(StringIO(data), on_bad_lines="abc")


def test_bad_header_uniform_error(all_parsers):
    parser = all_parsers
    data = "+++123456789...\ncol1,col2,col3,col4\n1,2,3,4\n"
    msg = "Expected 2 fields in line 2, saw 4"
    if parser.engine == "c":
        msg = (
            "Could not construct index. Requested to use 1 "
            "number of columns, but 3 left to parse."
        )
    elif parser.engine == "pyarrow":
        # "CSV parse error: Expected 1 columns, got 4: col1,col2,col3,col4"
        pytest.skip(reason="https://github.com/apache/arrow/issues/38676")

    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), index_col=0, on_bad_lines="error")


def test_on_bad_lines_warn_correct_formatting(all_parsers):
    # see gh-15925
    parser = all_parsers
    data = """1,2
a,b
a,b,c
a,b,d
a,b
"""
    expected = pd.DataFrame({"1": "a", "2": ["b"] * 2})
    match_msg = "Skipping line"

    expected_warning = ParserWarning
    if parser.engine == "pyarrow":
        match_msg = "Expected 2 columns, but found 3: a,b,c"

    with tm.assert_produces_warning(
        expected_warning, match=match_msg, check_stacklevel=False
    ):
        result = parser.read_csv(StringIO(data), on_bad_lines="warn")
    tm.assert_frame_equal(result, expected)


# GH#71034 the explanation pandas adds to pyarrow's column count error
NAMES_MISMATCH_MSG = re.escape(
    "The 'pyarrow' engine read 2 columns from the first row of the file "
    "but 4 names were passed. The number of columns is taken from the "
    "first row, not from 'names', and rows with a different number of "
    "fields are handled by 'on_bad_lines'. Use engine='c' or "
    "engine='python' to read files with ragged rows."
)


def test_names_wider_than_first_row(all_parsers):
    # GH#71034 the c and python engines take the width from "names" and pad the
    #  rows that are narrower; pyarrow takes it from the first row instead, so
    #  the wider rows are bad lines and the read cannot succeed
    parser = all_parsers
    data = "1,2\n3,4,5,6\n"
    names = ["w", "x", "y", "z"]

    if parser.engine == "pyarrow":
        with pytest.raises(ParserError, match=NAMES_MISMATCH_MSG):
            parser.read_csv(StringIO(data), header=None, names=names)
        return

    result = parser.read_csv(StringIO(data), header=None, names=names)
    expected = pd.DataFrame(
        {"w": [1, 3], "x": [2, 4], "y": [np.nan, 5.0], "z": [np.nan, 6.0]}
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("on_bad_lines", ["error", "skip", "warn"])
def test_names_wider_than_first_row_pyarrow(pyarrow_parser_only, on_bad_lines):
    # GH#71034 every way of handling the bad lines explains the limitation
    #  instead of leaving the user with pyarrow's column count or the
    #  "Length mismatch" error from assigning the names afterwards
    parser = pyarrow_parser_only
    data = "a,b\n1,2,3,4\n"
    names = ["w", "x", "y", "z"]

    if on_bad_lines == "error":
        # pyarrow's own message is kept and followed by the explanation
        msg = (
            r"^CSV parse error: Expected 2 columns, got 4: 1,2,3,4\. "
            f"{NAMES_MISMATCH_MSG}$"
        )
    else:
        # the wider rows were skipped, so pandas reports the mismatch itself
        msg = rf"^{NAMES_MISMATCH_MSG}$"

    if on_bad_lines == "warn":
        with tm.assert_produces_warning(
            ParserWarning,
            match="Expected 2 columns, but found 4: 1,2,3,4",
            check_stacklevel=False,
        ):
            with pytest.raises(ParserError, match=msg):
                parser.read_csv(
                    StringIO(data), header=None, names=names, on_bad_lines=on_bad_lines
                )
    else:
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(
                StringIO(data), header=None, names=names, on_bad_lines=on_bad_lines
            )


def test_names_shorter_than_file_pyarrow_still_forms_an_index(pyarrow_parser_only):
    # GH#71034, GH#65862 fewer names than the file must not be reported as a
    #  mismatch: the leading columns the names do not cover become the index
    parser = pyarrow_parser_only
    data = "1,2,3\n4,5,6"
    result = parser.read_csv(StringIO(data), header=None, names=["b", "c"])
    expected = pd.DataFrame({"b": [2, 5], "c": [3, 6]}, index=[1, 4])
    tm.assert_frame_equal(result, expected)


def test_names_length_matches_first_row_error_unchanged(pyarrow_parser_only):
    # GH#71034 the explanation is only appended when the number of columns
    #  pyarrow read disagrees with "names"
    parser = pyarrow_parser_only
    data = "a,b\n1,2,3,4\n"

    msg = "^CSV parse error: Expected 2 columns, got 4: 1,2,3,4$"
    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), header=None, names=["w", "x"])

    # and with no names at all, as there is nothing to explain
    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data))


def test_uniformly_short_rows_with_longer_names_pyarrow(pyarrow_parser_only):
    # GH#71034 no row is bad here, pyarrow simply reads fewer columns than
    #  "names" has, which used to leak a "Length mismatch" ValueError
    parser = pyarrow_parser_only
    data = "1,2\n3,4"

    msg = re.escape(
        "The 'pyarrow' engine read 2 columns from the first row of the file "
        "but 3 names were passed."
    )
    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), header=None, names=["a", "b", "c"])
