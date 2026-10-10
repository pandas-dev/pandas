"""
Tests that apply specifically to the Python parser. Unless specifically
stated as a Python-specific issue, the goal is to eventually move as many of
these tests out of this module as soon as the C parser can accept further
arguments when parsing.
"""

from __future__ import annotations

import csv
from io import (
    BytesIO,
    StringIO,
    TextIOWrapper,
)
from typing import TYPE_CHECKING

import numpy as np
import pytest

from pandas.errors import (
    ParserError,
    ParserWarning,
)

import pandas as pd
import pandas._testing as tm

if TYPE_CHECKING:
    from collections.abc import Iterator


def test_default_separator(python_parser_only):
    # see gh-17333
    #
    # csv.Sniffer in Python treats "o" as separator.
    data = "aob\n1o2\n3o4"
    parser = python_parser_only
    expected = pd.DataFrame({"a": [1, 3], "b": [2, 4]})

    result = parser.read_csv(StringIO(data), sep=None)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("skipfooter", ["foo", 1.5, True])
def test_invalid_skipfooter_non_int(python_parser_only, skipfooter):
    # see gh-15925 (comment)
    data = "a\n1\n2"
    parser = python_parser_only
    msg = "skipfooter must be an integer"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), skipfooter=skipfooter)


def test_invalid_skipfooter_negative(python_parser_only):
    # see gh-15925 (comment)
    data = "a\n1\n2"
    parser = python_parser_only
    msg = "skipfooter cannot be negative"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), skipfooter=-1)


@pytest.mark.parametrize("kwargs", [{"sep": None}, {"delimiter": "|"}])
def test_sniff_delimiter(python_parser_only, kwargs):
    data = """index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
    parser = python_parser_only
    result = parser.read_csv(StringIO(data), index_col=0, **kwargs)
    expected = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        columns=["A", "B", "C"],
        index=pd.Index(["foo", "bar", "baz"], name="index"),
    )
    tm.assert_frame_equal(result, expected)


def test_sniff_delimiter_comment(python_parser_only):
    data = """# comment line
index|A|B|C
# comment line
foo|1|2|3 # ignore | this
bar|4|5|6
baz|7|8|9
"""
    parser = python_parser_only
    result = parser.read_csv(StringIO(data), index_col=0, sep=None, comment="#")
    expected = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        columns=["A", "B", "C"],
        index=pd.Index(["foo", "bar", "baz"], name="index"),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("encoding", [None, "utf-8"])
def test_sniff_delimiter_encoding(python_parser_only, encoding):
    parser = python_parser_only
    data = """ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""

    if encoding is not None:
        data = data.encode(encoding)
        data = BytesIO(data)
        data = TextIOWrapper(data, encoding=encoding)
    else:
        data = StringIO(data)

    result = parser.read_csv(data, index_col=0, sep=None, skiprows=2, encoding=encoding)
    expected = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        columns=["A", "B", "C"],
        index=pd.Index(["foo", "bar", "baz"], name="index"),
    )
    tm.assert_frame_equal(result, expected)


def test_single_line(python_parser_only):
    # see gh-6607: sniff separator
    parser = python_parser_only
    result = parser.read_csv(StringIO("1,2"), names=["a", "b"], header=None, sep=None)

    expected = pd.DataFrame({"a": [1], "b": [2]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("kwargs", [{"skipfooter": 2}, {"nrows": 3}])
def test_skipfooter(python_parser_only, kwargs):
    # see gh-6607
    data = """A,B,C
1,2,3
4,5,6
7,8,9
want to skip this
also also skip this
"""
    parser = python_parser_only
    result = parser.read_csv(StringIO(data), **kwargs)

    expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["A", "B", "C"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "compression,klass", [("gzip", "GzipFile"), ("bz2", "BZ2File")]
)
def test_decompression_regex_sep(
    temp_file, python_parser_only, csv1, compression, klass
):
    # see gh-6607
    parser = python_parser_only

    with open(csv1, "rb") as f:
        data = f.read()

    data = data.replace(b",", b"::")
    expected = parser.read_csv(csv1)

    module = pytest.importorskip(compression)
    klass = getattr(module, klass)

    with klass(temp_file, mode="wb") as tmp:
        tmp.write(data)

    result = parser.read_csv(temp_file, sep="::", compression=compression)
    tm.assert_frame_equal(result, expected)


def test_read_csv_buglet_4x_multi_index(python_parser_only):
    # see gh-6607
    data = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""
    parser = python_parser_only

    expected = pd.DataFrame(
        [
            [-0.5109, -2.3358, -0.4645, 0.05076, 0.3640],
            [0.4473, 1.4152, 0.2834, 1.00661, 0.1744],
            [-0.6662, -0.5243, -0.3580, 0.89145, 2.5838],
        ],
        columns=["A", "B", "C", "D", "E"],
        index=pd.MultiIndex.from_tuples(
            [("a", "b", 10.0032, 5), ("a", "q", 20, 4), ("x", "q", 30, 3)],
            names=["one", "two", "three", "four"],
        ),
    )
    result = parser.read_csv(StringIO(data), sep=r"\s+")
    tm.assert_frame_equal(result, expected)


def test_read_csv_buglet_4x_multi_index2(python_parser_only):
    # see gh-6893
    data = "      A B C\na b c\n1 3 7 0 3 6\n3 1 4 1 5 9"
    parser = python_parser_only

    expected = pd.DataFrame.from_records(
        [(1, 3, 7, 0, 3, 6), (3, 1, 4, 1, 5, 9)],
        columns=list("abcABC"),
        index=list("abc"),
    )
    result = parser.read_csv(StringIO(data), sep=r"\s+")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("add_footer", [True, False])
def test_skipfooter_with_decimal(python_parser_only, add_footer):
    # see gh-6971
    data = "1#2\n3#4"
    parser = python_parser_only
    expected = pd.DataFrame({"a": [1.2, 3.4]})

    if add_footer:
        # The stray footer line should not mess with the
        # casting of the first two lines if we skip it.
        kwargs = {"skipfooter": 1}
        data += "\nFooter"
    else:
        kwargs = {}

    result = parser.read_csv(StringIO(data), names=["a"], decimal="#", **kwargs)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "sep", ["::", "#####", "!!!", "123", "#1!c5", "%!c!d", "@@#4:2", "_!pd#_"]
)
@pytest.mark.parametrize(
    "encoding", ["utf-16", "utf-16-be", "utf-16-le", "utf-32", "cp037"]
)
def test_encoding_non_utf8_multichar_sep(python_parser_only, sep, encoding):
    # see gh-3404
    expected = pd.DataFrame({"a": [1], "b": [2]})
    parser = python_parser_only

    data = "1" + sep + "2"
    encoded_data = data.encode(encoding)

    result = parser.read_csv(
        BytesIO(encoded_data), sep=sep, names=["a", "b"], encoding=encoding
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("quoting", [csv.QUOTE_MINIMAL, csv.QUOTE_NONE])
def test_multi_char_sep_quotes(python_parser_only, quoting):
    # see gh-13374
    kwargs = {"sep": ",,"}
    parser = python_parser_only

    data = 'a,,b\n1,,a\n2,,"2,,b"'

    if quoting == csv.QUOTE_NONE:
        msg = "Expected 2 fields in line 3, saw 3"
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), quoting=quoting, **kwargs)
    else:
        msg = "ignored when a multi-char delimiter is used"
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), quoting=quoting, **kwargs)


def test_none_delimiter(python_parser_only):
    # see gh-13374 and gh-17465
    parser = python_parser_only
    data = "a,b,c\n0,1,2\n3,4,5,6\n7,8,9"
    expected = pd.DataFrame({"a": [0, 7], "b": [1, 8], "c": [2, 9]})

    # We expect the third line in the data to be
    # skipped because it is malformed, but we do
    # not expect any errors to occur.
    with tm.assert_produces_warning(
        ParserWarning, match="Skipping line 3", check_stacklevel=False
    ):
        result = parser.read_csv(
            StringIO(data), header=0, sep=None, on_bad_lines="warn"
        )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("data", ['a\n1\n"b"a', 'a,b,c\ncat,foo,bar\ndog,foo,"baz'])
@pytest.mark.parametrize("skipfooter", [0, 1])
def test_skipfooter_bad_row(python_parser_only, data, skipfooter):
    # see gh-13879 and gh-15910
    parser = python_parser_only
    if skipfooter:
        msg = "parsing errors in the skipped footer rows"
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), skipfooter=skipfooter)
    else:
        msg = "|".join(["unexpected end of data", "expected after"])
        with pytest.raises(ParserError, match=msg):
            parser.read_csv(StringIO(data), skipfooter=skipfooter)


def test_malformed_skipfooter(python_parser_only):
    parser = python_parser_only
    data = """ignore
A,B,C
1,2,3 # comment
1,2,3,4,5
2,3,4
footer
"""
    msg = "Expected 3 fields in line 4, saw 5"
    with pytest.raises(ParserError, match=msg):
        parser.read_csv(StringIO(data), header=1, comment="#", skipfooter=1)


@pytest.mark.parametrize(
    "data, line",
    [
        ('a,b,c\n1,2,"p\nq"\n3,"x"y,5\n', 4),
        # unclosed quote: the reader stops at the end of the file
        ('a,b\n1,2\n3,"x\n4,5\n6,7\n', 3),
        ('a,b\n"p\nq",2\n3,"x\n4,5\n6,7\n', 4),
        # past the rows read ahead for the header and index
        ('a,b,c\n1,2,3\n4,5,6\n7,8,9\n1,2,"p\nq"\n3,"x"y,5\n', 7),
        ('a,b\n1,2\n3,4\n5,6\n"p\nq",8\n7,"x\n8,9\n', 7),
    ],
)
def test_csv_error_line_number_counts_quoted_newlines(python_parser_only, data, line):
    # GH#16286: csv.Error is reported at the line its record starts on
    parser = python_parser_only
    with tm.assert_produces_warning(
        ParserWarning, match=f"Skipping line {line}:", check_stacklevel=False
    ):
        parser.read_csv(StringIO(data), on_bad_lines="warn")


@pytest.mark.parametrize("chunksize", [1, 2, 3])
def test_bad_line_number_counts_quoted_newlines_chunked(python_parser_only, chunksize):
    # GH#16286: the multi-line row and the bad row can fall in different chunks
    parser = python_parser_only
    data = 'a,b,c\n1,2,"x\ny"\n3,4,5\n6,7,8,9\n'
    with pytest.raises(ParserError, match="Expected 3 fields in line 5, saw 4"):
        with parser.read_csv(StringIO(data), chunksize=chunksize) as reader:
            list(reader)


@pytest.mark.parametrize(
    "data, kwargs, line",
    [
        # csv.Error while reading a row
        ('junk\na,b,c\n1,2,"x\ny"\n3,"x"y,5\n', {"skiprows": 1}, 5),
        ('a,b,c\n1,2,"x\ny"\n3,"x"y,5\n', {"sep": None}, 4),
        # too many fields
        ('junk\na,b,c\n1,2,"x\ny"\n4,5,6,7\n', {"skiprows": 1}, 5),
        ('a,b,c\n1,2,"x\ny"\n4,5,6,7\n', {"sep": None}, 4),
    ],
)
def test_bad_line_number_lines_read_before_reader(
    python_parser_only, data, kwargs, line
):
    # GH#16286: leading skiprows and the sep=None sniff are read before the
    # csv reader exists, so its line count does not include them
    parser = python_parser_only
    with tm.assert_produces_warning(
        ParserWarning, match=f"Skipping line {line}:", check_stacklevel=False
    ):
        parser.read_csv(StringIO(data), on_bad_lines="warn", **kwargs)


def test_bad_line_number_line_dropped_between_buffered_rows(python_parser_only):
    # GH#16286: header=[0, 1] buffers rows ahead; a blank line dropped between
    # them must not shift the line reported for the earlier one
    parser = python_parser_only
    data = 'a,b\nc,"d\ne"\n1,2,3\n\n4,5\n'
    with tm.assert_produces_warning(
        ParserWarning, match="Skipping line 4:", check_stacklevel=False
    ):
        parser.read_csv(StringIO(data), header=[0, 1], on_bad_lines="warn")


def test_bad_line_number_comment_cuts_buffered_row(python_parser_only):
    # GH#16286: the comment cuts the line break out of a row read ahead for
    # the header and index, so its text no longer shows the lines it spanned
    parser = python_parser_only
    data = 'a,b\n1,2\n1,2,3,"#x\ny"\n5,6\n'
    with pytest.raises(ParserError, match="Expected 2 fields in line 3, saw 3"):
        parser.read_csv(StringIO(data), comment="#")


@pytest.mark.parametrize(
    "data, kwargs, line",
    [
        # records dropped from the batch still take up lines
        ('a,b,c\n1,2,3\n1,2,"p\nq"\n\n\n4,5,6,7\n', {}, 7),
        ('a,b,c\n1,2,"p\nq"\n\n4,5,6,7\n', {}, 5),
        ('a,b,c\n1,2,"p\nq"\n# c\n4,5,6,7\n', {"comment": "#"}, 5),
        ('a,b,c\n1,2,3\n1,2,"p #x\nq"\n4,5,6,7\n', {"comment": "#"}, 5),
        ('a,b,c\n1,2,"p\nq"\n4,5,6,7\n"f\no"\n', {"skipfooter": 1}, 4),
        (',a,b\nidx,,\n1,2,"p\nq"\n4,5,6,7\n', {"index_col": 0}, 5),
        ('a,b,c\n1,2,"p\nq"\n\n3,4,5\n4,5,6,7\n', {"chunksize": 2}, 6),
    ],
)
def test_bad_line_number_records_dropped_from_batch(
    python_parser_only, data, kwargs, line
):
    # GH#16286
    parser = python_parser_only
    with pytest.raises(ParserError, match=f"Expected 3 fields in line {line}, saw 4"):
        result = parser.read_csv(StringIO(data), **kwargs)
        if "chunksize" in kwargs:
            list(result)


def test_bad_line_number_warn_records_dropped_from_batch(python_parser_only):
    # GH#16286
    parser = python_parser_only
    data = 'a,b,c\n1,2,"p\nq"\n\n4,5,6,7\n# c\n1,"x\ny",3\n\n1,2,3,4\n'
    with tm.assert_produces_warning(
        ParserWarning, match="Skipping line 5:", check_stacklevel=False
    ) as record:
        parser.read_csv(StringIO(data), comment="#", on_bad_lines="warn")
    assert "Skipping line 10:" in str(record[-1].message)


def test_bad_line_number_lone_carriage_return(python_parser_only, tmp_path):
    # GH#16286: a file path is opened with newline="", so csv.reader counts the
    # lone "\r" as a line break; the row number is reported instead
    parser = python_parser_only
    path = tmp_path / "data.csv"
    path.write_bytes(b'a,b,c\n1,2,"x\ry"\n4,5,6,7\n')
    with pytest.raises(ParserError, match="Expected 3 fields in line 3, saw 4"):
        parser.read_csv(path)


def test_bad_line_number_callable_edits_row(python_parser_only):
    # GH#16286: a callable that edits its row in place must not change the
    # line breaks counted for that row
    parser = python_parser_only
    data = 'a,b,c\n1,2,3\n"x\ny",1,9,8\n1,2,3,4,5\n'

    def drop_first(line):
        del line[0]
        return line

    with tm.assert_produces_warning(
        ParserWarning,
        match="in line 5, saw 4 from bad_lines callable",
        check_stacklevel=False,
    ):
        parser.read_csv(StringIO(data), on_bad_lines=drop_first)


def test_python_engine_file_no_next(python_parser_only):
    parser = python_parser_only

    class NoNextBuffer:
        def __init__(self, csv_data) -> None:
            self.data = csv_data

        def __iter__(self) -> Iterator:
            return self.data.__iter__()

        def read(self):
            return self.data

        def readline(self):
            return self.data

    parser.read_csv(NoNextBuffer("a\n1"))


@pytest.mark.parametrize("bad_line_func", [lambda x: ["2", "3"], lambda x: x[:2]])
def test_on_bad_lines_callable(python_parser_only, bad_line_func):
    # GH 5686
    parser = python_parser_only
    data = """a,b
1,2
2,3,4,5,6
3,4
"""
    bad_sio = StringIO(data)
    result = parser.read_csv(bad_sio, on_bad_lines=bad_line_func)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [2, 3, 4]})
    tm.assert_frame_equal(result, expected)


def test_on_bad_lines_callable_write_to_external_list(python_parser_only):
    # GH 5686
    parser = python_parser_only
    data = """a,b
1,2
2,3,4,5,6
3,4
"""
    bad_sio = StringIO(data)
    lst = []

    def bad_line_func(bad_line: list[str]) -> list[str]:
        lst.append(bad_line)
        return ["2", "3"]

    result = parser.read_csv(bad_sio, on_bad_lines=bad_line_func)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [2, 3, 4]})
    tm.assert_frame_equal(result, expected)
    assert lst == [["2", "3", "4", "5", "6"]]


@pytest.mark.parametrize("bad_line_func", [lambda x: ["foo", "bar"], lambda x: x[:2]])
@pytest.mark.parametrize("sep", [",", "111"])
def test_on_bad_lines_callable_iterator_true(python_parser_only, bad_line_func, sep):
    # GH 5686
    # iterator=True has a separate code path than iterator=False
    parser = python_parser_only
    data = f"""
0{sep}1
hi{sep}there
foo{sep}bar{sep}baz
good{sep}bye
"""
    bad_sio = StringIO(data)
    result_iter = parser.read_csv(
        bad_sio, on_bad_lines=bad_line_func, chunksize=1, iterator=True, sep=sep
    )
    expecteds = [
        {"0": "hi", "1": "there"},
        {"0": "foo", "1": "bar"},
        {"0": "good", "1": "bye"},
    ]
    for i, (result, expected) in enumerate(zip(result_iter, expecteds, strict=True)):
        expected = pd.DataFrame(expected, index=range(i, i + 1))
        tm.assert_frame_equal(result, expected)


def test_on_bad_lines_callable_dont_swallow_errors(python_parser_only):
    # GH 5686
    parser = python_parser_only
    data = """a,b
1,2
2,3,4,5,6
3,4
"""
    bad_sio = StringIO(data)
    msg = "This function is buggy."

    def bad_line_func(bad_line):
        raise ValueError(msg)

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(bad_sio, on_bad_lines=bad_line_func)


def test_on_bad_lines_callable_not_expected_length(python_parser_only):
    # GH 5686
    parser = python_parser_only
    data = """a,b
1,2
2,3,4,5,6
3,4
"""
    bad_sio = StringIO(data)

    result = parser.read_csv_check_warnings(
        ParserWarning, "from bad_lines callable", bad_sio, on_bad_lines=lambda x: x
    )
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [2, 3, 4]})
    tm.assert_frame_equal(result, expected)


def test_on_bad_lines_callable_returns_none(python_parser_only):
    # GH 5686
    parser = python_parser_only
    data = """a,b
1,2
2,3,4,5,6
3,4
"""
    bad_sio = StringIO(data)

    result = parser.read_csv(bad_sio, on_bad_lines=lambda x: None)
    expected = pd.DataFrame({"a": [1, 3], "b": [2, 4]})
    tm.assert_frame_equal(result, expected)


def test_on_bad_lines_index_col_inferred(python_parser_only):
    # GH 5686
    parser = python_parser_only
    data = """a,b
1,2,3
4,5,6
"""
    bad_sio = StringIO(data)

    result = parser.read_csv(bad_sio, on_bad_lines=lambda x: ["99", "99"])
    expected = pd.DataFrame({"a": [2, 5], "b": [3, 6]}, index=[1, 4])
    tm.assert_frame_equal(result, expected)


def test_index_col_false_and_header_none(python_parser_only):
    # GH#46955
    parser = python_parser_only
    data = """
0.5,0.03
0.1,0.2,0.3,2
"""
    result = parser.read_csv_check_warnings(
        ParserWarning,
        "Length of header",
        StringIO(data),
        sep=",",
        header=None,
        index_col=False,
    )
    expected = pd.DataFrame({0: [0.5, 0.1], 1: [0.03, 0.2]})
    tm.assert_frame_equal(result, expected)


def test_header_int_do_not_infer_multiindex_names_on_different_line(python_parser_only):
    # GH#46569
    parser = python_parser_only
    data = StringIO("a\na,b\nc,d,e\nf,g,h")
    result = parser.read_csv_check_warnings(
        ParserWarning, "Length of header", data, engine="python", index_col=False
    )
    expected = pd.DataFrame({"a": ["a", "c", "f"]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype", [{"a": object}, {"a": str, "b": np.int64, "c": np.int64}]
)
def test_no_thousand_convert_with_dot_for_non_numeric_cols(python_parser_only, dtype):
    # GH#50270
    parser = python_parser_only
    data = """\
a;b;c
0000.7995;16.000;0
3.03.001.00514;0;4.000
4923.600.041;23.000;131"""
    result = parser.read_csv(
        StringIO(data),
        sep=";",
        dtype=dtype,
        thousands=".",
    )
    expected = pd.DataFrame(
        {
            "a": ["0000.7995", "3.03.001.00514", "4923.600.041"],
            "b": [16000, 0, 23000],
            "c": [0, 4000, 131],
        }
    )
    if dtype["a"] == object:
        expected["a"] = expected["a"].astype(object)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype,expected",
    [
        (
            {"a": str, "b": np.float64, "c": np.int64},
            {
                "b": [16000.1, 0, 23000],
                "c": [0, 4001, 131],
            },
        ),
        (
            str,
            {
                "b": ["16,000.1", "0", "23,000"],
                "c": ["0", "4,001", "131"],
            },
        ),
    ],
)
def test_no_thousand_convert_for_non_numeric_cols(python_parser_only, dtype, expected):
    # GH#50270
    parser = python_parser_only
    data = """a;b;c
0000,7995;16,000.1;0
3,03,001,00514;0;4,001
4923,600,041;23,000;131
"""
    result = parser.read_csv(
        StringIO(data),
        sep=";",
        dtype=dtype,
        thousands=",",
    )
    expected = pd.DataFrame(expected)
    expected.insert(0, "a", ["0000,7995", "3,03,001,00514", "4923,600,041"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("index_col", [None, 0])
def test_on_bad_lines_callable_warns_and_truncates_with_index_col(
    python_parser_only, index_col
):
    # GH#61837
    parser = python_parser_only
    data = "id,field_1,field_2\n101,A,B\n102,C,D,E\n103,F,G\n"

    def fixer(bad_line):
        return [*list(bad_line), "EXTRA1", "EXTRA2"]

    result = parser.read_csv_check_warnings(
        ParserWarning,
        "from bad_lines callable",
        StringIO(data),
        on_bad_lines=fixer,
        index_col=index_col,
    )

    if index_col is None:
        expected = pd.DataFrame(
            {
                "id": [101, 102, 103],
                "field_1": ["A", "C", "F"],
                "field_2": ["B", "D", "G"],
            }
        )
    else:
        expected = pd.DataFrame(
            {"field_1": ["A", "C", "F"], "field_2": ["B", "D", "G"]},
            index=pd.Index([101, 102, 103], name="id"),
        )

    tm.assert_frame_equal(result, expected)


def test_read_csv_leading_quote_skip(python_parser_only):
    # GH 62739
    tbl = """\
    "
a b
1 3
"""
    parser = python_parser_only
    result = parser.read_csv(
        StringIO(tbl),
        delimiter=" ",
        skiprows=1,
    )
    expected = pd.DataFrame({"a": [1], "b": [3]})
    tm.assert_frame_equal(result, expected)


def test_read_csv_unclosed_double_quote_in_data_still_errors(python_parser_only):
    # GH 62739
    tbl = """\
a b
"
1 3
"""
    parser = python_parser_only
    with pytest.raises(ParserError, match="unexpected end of data"):
        parser.read_csv(StringIO(tbl), delimiter=" ", skiprows=1)


def test_read_csv_skiprows_zero(python_parser_only):
    # GH 62739
    tbl = """\
"
a b
1 3
"""
    parser = python_parser_only
    # don't skip anything
    with pytest.raises(ParserError, match="unexpected end of data"):
        parser.read_csv(StringIO(tbl), delimiter=" ", skiprows=0, engine="python")


def test_memory_map_multichar_sep(python_parser_only, tmp_path):
    # GH#34577 - memory_map with multi-character separator should not raise
    # TypeError about mixing string patterns and bytes-like objects
    parser = python_parser_only
    path = tmp_path / "test.csv"
    path.write_text("key1 - value1\nkey2 - value2\nkey3 - value3\n", encoding="utf-8")

    result = parser.read_csv(
        path, header=None, sep=" - ", names=["key", "value"], memory_map=True
    )
    expected = pd.DataFrame(
        {"key": ["key1", "key2", "key3"], "value": ["value1", "value2", "value3"]}
    )
    tm.assert_frame_equal(result, expected)
