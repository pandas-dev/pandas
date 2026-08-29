"""
Tests that apply specifically to the CParser. Unless specifically stated
as a CParser-specific issue, the goal is to eventually move as many of
these tests out of this module as soon as the Python parser can accept
further arguments when parsing.
"""

from decimal import Decimal
from io import (
    BytesIO,
    StringIO,
    TextIOWrapper,
)
import mmap
import tarfile

import numpy as np
import pytest

from pandas._libs import parsers as libparsers
from pandas.compat import WASM
from pandas.errors import (
    Pandas4Warning,
    ParserError,
    ParserWarning,
)
import pandas.util._test_decorators as td

from pandas import (
    ArrowDtype,
    DataFrame,
    StringDtype,
    concat,
    option_context,
    read_csv,
)
import pandas._testing as tm


@pytest.mark.parametrize(
    "malformed,expected_data",
    [
        ("1\r1\r1\r 1\r 1\r", [[1], [1], [1], [1], [1]]),
        ("1\r1\r1\r 1\r 1\r11\r", [[1], [1], [1], [1], [1], [11]]),
        ("1\r1\r1\r 1\r 1\r11\r1\r", [[1], [1], [1], [1], [1], [11], [1]]),
    ],
    ids=["words pointer", "stream pointer", "lines pointer"],
)
def test_buffer_overflow(c_parser_only, malformed, expected_data):
    # see gh-9205, gh-51141: test certain malformed input files that
    # previously caused buffer overflows in tokenizer.c due to
    # \r characters being treated as line terminators, then triggering
    # an infinite re-parsing loop in the WHITESPACE_LINE state.
    parser = c_parser_only
    result = parser.read_csv(StringIO(malformed), header=None)
    expected = DataFrame(expected_data)
    tm.assert_frame_equal(result, expected)


def test_cr_in_field_with_trailing_space(c_parser_only):
    # GH#51141 - embedded \r followed by space in unquoted field
    # should not cause infinite re-parsing or buffer overflow
    parser = c_parser_only
    data = "a,b,c\n1,2,3\n4,5\r X,6\n"
    result = parser.read_csv(StringIO(data))
    assert len(result) == 3
    assert result.shape == (3, 3)


def test_delim_whitespace_custom_terminator(c_parser_only):
    # See gh-12912
    data = "a b c~1 2 3~4 5 6~7 8 9"
    parser = c_parser_only

    df = parser.read_csv(StringIO(data), lineterminator="~", sep=r"\s+")
    expected = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["a", "b", "c"])
    tm.assert_frame_equal(df, expected)


def test_dtype_and_names_error(c_parser_only):
    # see gh-8833: passing both dtype and names
    # resulting in an error reporting issue
    parser = c_parser_only
    data = """
1.0 1
2.0 2
3.0 3
"""
    # base cases
    result = parser.read_csv(StringIO(data), sep=r"\s+", header=None)
    expected = DataFrame([[1.0, 1], [2.0, 2], [3.0, 3]])
    tm.assert_frame_equal(result, expected)

    result = parser.read_csv(StringIO(data), sep=r"\s+", header=None, names=["a", "b"])
    expected = DataFrame([[1.0, 1], [2.0, 2], [3.0, 3]], columns=["a", "b"])
    tm.assert_frame_equal(result, expected)

    # fallback casting
    result = parser.read_csv(
        StringIO(data), sep=r"\s+", header=None, names=["a", "b"], dtype={"a": np.int32}
    )
    expected = DataFrame([[1, 1], [2, 2], [3, 3]], columns=["a", "b"])
    expected["a"] = expected["a"].astype(np.int32)
    tm.assert_frame_equal(result, expected)

    data = """
1.0 1
nan 2
3.0 3
"""
    # fallback casting, but not castable
    if not WASM:  # no fp exception support in wasm
        with pytest.raises(ValueError, match="cannot safely convert"):
            with tm.assert_produces_warning(RuntimeWarning, check_stacklevel=False):
                parser.read_csv(
                    StringIO(data),
                    sep=r"\s+",
                    header=None,
                    names=["a", "b"],
                    dtype={"a": np.int32},
                )


@pytest.mark.parametrize(
    "match,kwargs",
    [
        # For each of these cases, all of the dtypes are valid, just unsupported.
        (
            (
                "the dtype datetime64 is not supported for parsing, "
                "pass this column using parse_dates instead"
            ),
            {"dtype": {"A": "datetime64", "B": "float64"}},
        ),
        (
            (
                "the dtype datetime64 is not supported for parsing, "
                "pass this column using parse_dates instead"
            ),
            {"dtype": {"A": "datetime64", "B": "float64"}, "parse_dates": ["B"]},
        ),
        (
            "the dtype timedelta64 is not supported for parsing",
            {"dtype": {"A": "timedelta64", "B": "float64"}},
        ),
        (
            f"the dtype {tm.ENDIAN}U8 is not supported for parsing",
            {"dtype": {"A": "U8"}},
        ),
    ],
    ids=["dt64-0", "dt64-1", "td64", f"{tm.ENDIAN}U8"],
)
def test_unsupported_dtype(c_parser_only, match, kwargs, temp_file):
    parser = c_parser_only
    df = DataFrame(
        np.random.default_rng(2).random((5, 2)),
        columns=list("AB"),
        index=["1A", "1B", "1C", "1D", "1E"],
    )

    df.to_csv(temp_file)

    with pytest.raises(TypeError, match=match):
        parser.read_csv(temp_file, index_col=0, **kwargs)


@td.skip_if_32bit
@pytest.mark.slow
# test numbers between 1 and 2
@pytest.mark.parametrize("num", np.linspace(1.0, 2.0, num=21))
def test_precise_conversion(c_parser_only, num):
    parser = c_parser_only

    normal_errors = []
    precise_errors = []

    def error(val: float, actual_val: Decimal) -> Decimal:
        return abs(Decimal(f"{val:.100}") - actual_val)

    # 25 decimal digits of precision
    text = f"a\n{num:.25}"

    depr_msg = "float_precision"
    with tm.assert_produces_warning(
        Pandas4Warning, match=depr_msg, check_stacklevel=False
    ):
        normal_val = float(
            parser.read_csv(StringIO(text), float_precision="legacy")["a"][0]
        )
    with tm.assert_produces_warning(
        Pandas4Warning, match=depr_msg, check_stacklevel=False
    ):
        precise_val = float(
            parser.read_csv(StringIO(text), float_precision="high")["a"][0]
        )
    with tm.assert_produces_warning(
        Pandas4Warning, match=depr_msg, check_stacklevel=False
    ):
        roundtrip_val = float(
            parser.read_csv(StringIO(text), float_precision="round_trip")["a"][0]
        )
    actual_val = Decimal(text[2:])

    normal_errors.append(error(normal_val, actual_val))
    precise_errors.append(error(precise_val, actual_val))

    # round-trip should match float()
    assert roundtrip_val == float(text[2:])

    assert sum(precise_errors) <= sum(normal_errors)
    assert max(precise_errors) <= max(normal_errors)


def test_usecols_dtypes(c_parser_only, using_infer_string):
    parser = c_parser_only
    data = """\
1,2,3
4,5,6
7,8,9
10,11,12"""

    result = parser.read_csv(
        StringIO(data),
        usecols=(0, 1, 2),
        names=("a", "b", "c"),
        header=None,
        converters={"a": str},
        dtype={"b": int, "c": float},
    )
    result2 = parser.read_csv(
        StringIO(data),
        usecols=(0, 2),
        names=("a", "b", "c"),
        header=None,
        converters={"a": str},
        dtype={"b": int, "c": float},
    )

    if using_infer_string:
        assert (result.dtypes == ["string", int, float]).all()
        assert (result2.dtypes == ["string", float]).all()
    else:
        assert (result.dtypes == [object, int, float]).all()
        assert (result2.dtypes == [object, float]).all()


def test_disable_bool_parsing(c_parser_only):
    # see gh-2090

    parser = c_parser_only
    data = """A,B,C
Yes,No,Yes
No,Yes,Yes
Yes,,Yes
No,No,No"""

    result = parser.read_csv(StringIO(data), dtype=object)
    assert (result.dtypes == object).all()

    result = parser.read_csv(StringIO(data), dtype=object, na_filter=False)
    assert result["B"][2] == ""


def test_custom_lineterminator(c_parser_only):
    parser = c_parser_only
    data = "a,b,c~1,2,3~4,5,6"

    result = parser.read_csv(StringIO(data), lineterminator="~")
    expected = parser.read_csv(StringIO(data.replace("~", "\n")))

    tm.assert_frame_equal(result, expected)


def test_parse_ragged_csv(c_parser_only):
    parser = c_parser_only
    data = """1,2,3
1,2,3,4
1,2,3,4,5
1,2
1,2,3,4"""

    nice_data = """1,2,3,,
1,2,3,4,
1,2,3,4,5
1,2,,,
1,2,3,4,"""
    result = parser.read_csv(
        StringIO(data), header=None, names=["a", "b", "c", "d", "e"]
    )

    expected = parser.read_csv(
        StringIO(nice_data), header=None, names=["a", "b", "c", "d", "e"]
    )

    tm.assert_frame_equal(result, expected)

    # too many columns, cause segfault if not careful
    data = "1,2\n3,4,5"

    result = parser.read_csv(StringIO(data), header=None, names=range(50))
    expected = parser.read_csv(StringIO(data), header=None, names=range(3)).reindex(
        columns=range(50)
    )

    tm.assert_frame_equal(result, expected)


def test_tokenize_CR_with_quoting(c_parser_only):
    # see gh-3453
    parser = c_parser_only
    data = ' a,b,c\r"a,b","e,d","f,f"'

    result = parser.read_csv(StringIO(data), header=None)
    expected = parser.read_csv(StringIO(data.replace("\r", "\n")), header=None)
    tm.assert_frame_equal(result, expected)

    result = parser.read_csv(StringIO(data))
    expected = parser.read_csv(StringIO(data.replace("\r", "\n")))
    tm.assert_frame_equal(result, expected)


@pytest.mark.slow
@pytest.mark.parametrize("count", [3 * 2**n for n in range(6)])
def test_grow_boundary_at_cap(c_parser_only, count):
    # See gh-12494
    #
    # Cause of error was that the C parser
    # was not increasing the buffer size when
    # the desired space would fill the buffer
    # to capacity, which would later cause a
    # buffer overflow error when checking the
    # EOF terminator of the CSV stream.
    # 3 * 2^n commas was observed to break the parser
    parser = c_parser_only

    with StringIO("," * count) as s:
        expected = DataFrame(columns=[f"Unnamed: {i}" for i in range(count + 1)])
        df = parser.read_csv(s)
    tm.assert_frame_equal(df, expected)


@pytest.mark.slow
@pytest.mark.parametrize("encoding", [None, "utf-8"])
def test_parse_trim_buffers(c_parser_only, encoding):
    # This test is part of a bugfix for gh-13703. It attempts to
    # to stress the system memory allocator, to cause it to move the
    # stream buffer and either let the OS reclaim the region, or let
    # other memory requests of parser otherwise modify the contents
    # of memory space, where it was formally located.
    # This test is designed to cause a `segfault` with unpatched
    # `tokenizer.c`. Sometimes the test fails on `segfault`, other
    # times it fails due to memory corruption, which causes the
    # loaded DataFrame to differ from the expected one.

    # Also force 'utf-8' encoding, so that `_string_convert` would take
    # a different execution branch.

    parser = c_parser_only

    # Generate a large mixed-type CSV file on-the-fly (one record is
    # approx 1.5KiB).
    record_ = (
        """9999-9,99:99,,,,ZZ,ZZ,,,ZZZ-ZZZZ,.Z-ZZZZ,-9.99,,,9.99,Z"""
        """ZZZZ,,-99,9,ZZZ-ZZZZ,ZZ-ZZZZ,,9.99,ZZZ-ZZZZZ,ZZZ-ZZZZZ,"""
        """ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,9"""
        """99,ZZZ-ZZZZ,,ZZ-ZZZZ,,,,,ZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,,,9,9,"""
        """9,9,99,99,999,999,ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZ,9,ZZ-ZZZZ,9."""
        """99,ZZ-ZZZZ,ZZ-ZZZZ,,,,ZZZZ,,,ZZ,ZZ,,,,,,,,,,,,,9,,,999."""
        """99,999.99,,,ZZZZZ,,,Z9,,,,,,,ZZZ,ZZZ,,,,,,,,,,,ZZZZZ,ZZ"""
        """ZZZ,ZZZ-ZZZZZZ,ZZZ-ZZZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,ZZ-ZZZZ,ZZ-ZZ"""
        """ZZ,,,999999,999999,ZZZ,ZZZ,,,ZZZ,ZZZ,999.99,999.99,,,,Z"""
        """ZZ-ZZZ,ZZZ-ZZZ,-9.99,-9.99,9,9,,99,,9.99,9.99,9,9,9.99,"""
        """9.99,,,,9.99,9.99,,99,,99,9.99,9.99,,,ZZZ,ZZZ,,999.99,,"""
        """999.99,ZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,,,ZZZZZ,ZZZZZ,ZZZ,ZZZ,9,9,"""
        """,,,,,ZZZ-ZZZZ,ZZZ999Z,,,999.99,,999.99,ZZZ-ZZZZ,,,9.999"""
        """,9.999,9.999,9.999,-9.999,-9.999,-9.999,-9.999,9.999,9."""
        """999,9.999,9.999,9.999,9.999,9.999,9.999,99999,ZZZ-ZZZZ,"""
        """,9.99,ZZZ,,,,,,,,ZZZ,,,,,9,,,,9,,,,,,,,,,ZZZ-ZZZZ,ZZZ-Z"""
        """ZZZ,,ZZZZZ,ZZZZZ,ZZZZZ,ZZZZZ,,,9.99,,ZZ-ZZZZ,ZZ-ZZZZ,ZZ"""
        """,999,,,,ZZ-ZZZZ,ZZZ,ZZZ,ZZZ-ZZZZ,ZZZ-ZZZZ,,,99.99,99.99"""
        """,,,9.99,9.99,9.99,9.99,ZZZ-ZZZZ,,,ZZZ-ZZZZZ,,,,,-9.99,-"""
        """9.99,-9.99,-9.99,,,,,,,,,ZZZ-ZZZZ,,9,9.99,9.99,99ZZ,,-9"""
        """.99,-9.99,ZZZ-ZZZZ,,,,,,,ZZZ-ZZZZ,9.99,9.99,9999,,,,,,,"""
        """,,,-9.9,Z/Z-ZZZZ,999.99,9.99,,999.99,ZZ-ZZZZ,ZZ-ZZZZ,9."""
        """99,9.99,9.99,9.99,9.99,9.99,,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ-ZZ"""
        """ZZZ,ZZZ-ZZZZZ,ZZZ-ZZZZZ,ZZZ,ZZZ,ZZZ,ZZZ,9.99,,,-9.99,ZZ"""
        """-ZZZZ,-999.99,,-9999,,999.99,,,,999.99,99.99,,,ZZ-ZZZZZ"""
        """ZZZ,ZZ-ZZZZ-ZZZZZZZ,,,,ZZ-ZZ-ZZZZZZZZ,ZZZZZZZZ,ZZZ-ZZZZ"""
        """,9999,999.99,ZZZ-ZZZZ,-9.99,-9.99,ZZZ-ZZZZ,99:99:99,,99"""
        """,99,,9.99,,-99.99,,,,,,9.99,ZZZ-ZZZZ,-9.99,-9.99,9.99,9"""
        """.99,,ZZZ,,,,,,,ZZZ,ZZZ,,,,,"""
    )

    # Set the number of lines so that a call to `parser_trim_buffers`
    # is triggered: after a couple of full chunks are consumed a
    # relatively small 'residual' chunk would cause reallocation
    # within the parser.
    chunksize, n_lines = 128, 2 * 128 + 15
    csv_data = "\n".join([record_] * n_lines) + "\n"

    # We will use StringIO to load the CSV from this text buffer.
    # pd.read_csv() will iterate over the file in chunks and will
    # finally read a residual chunk of really small size.

    # Generate the expected output: manually create the dataframe
    # by splitting by comma and repeating the `n_lines` times.
    row = tuple(val_ if val_ else np.nan for val_ in record_.split(","))
    expected = DataFrame(
        [row for _ in range(n_lines)], dtype=object, columns=None, index=None
    )

    # Iterate over the CSV file in chunks of `chunksize` lines
    with parser.read_csv(
        StringIO(csv_data),
        header=None,
        dtype=object,
        chunksize=chunksize,
        encoding=encoding,
    ) as chunks_:
        result = concat(chunks_, axis=0, ignore_index=True)

    # Check for data corruption if there was no segfault
    tm.assert_frame_equal(result, expected)


def test_internal_null_byte(c_parser_only):
    # see gh-14012
    #
    # The null byte ('\x00') should not be used as a
    # true line terminator, escape character, or comment
    # character, only as a placeholder to indicate that
    # none was specified.
    #
    # This test should be moved to test_common.py ONLY when
    # Python's csv class supports parsing '\x00'.
    parser = c_parser_only

    names = ["a", "b", "c"]
    data = "1,2,3\n4,\x00,6\n7,8,9"
    # GH#19886 the NUL field is a one-character value, not an na_value, so the
    # column stays a string column rather than becoming float-with-NaN.
    expected = DataFrame(
        {"a": [1, 4, 7], "b": ["2", "\x00", "8"], "c": [3, 6, 9]}, columns=names
    )

    result = parser.read_csv(StringIO(data), names=names)
    tm.assert_frame_equal(result, expected)


def test_read_nrows_large(c_parser_only):
    # gh-7626 - Read only nrows of data in for large inputs (>262144b)
    parser = c_parser_only
    header_narrow = "\t".join(["COL_HEADER_" + str(i) for i in range(10)]) + "\n"
    data_narrow = "\t".join(["somedatasomedatasomedata1" for _ in range(10)]) + "\n"
    header_wide = "\t".join(["COL_HEADER_" + str(i) for i in range(15)]) + "\n"
    data_wide = "\t".join(["somedatasomedatasomedata2" for _ in range(15)]) + "\n"
    test_input = header_narrow + data_narrow * 1050 + header_wide + data_wide * 2

    df = parser.read_csv(StringIO(test_input), sep="\t", nrows=1010)

    assert df.size == 1010 * 10


def test_float_precision_round_trip_with_text(c_parser_only):
    # see gh-15140
    parser = c_parser_only
    with tm.assert_produces_warning(
        Pandas4Warning, match="float_precision", check_stacklevel=False
    ):
        df = parser.read_csv(StringIO("a"), header=None, float_precision="round_trip")
    tm.assert_frame_equal(df, DataFrame({0: ["a"]}))


def test_large_difference_in_columns(c_parser_only):
    # see gh-14125
    parser = c_parser_only

    count = 10000
    large_row = ("X," * count)[:-1] + "\n"
    normal_row = "XXXXXX XXXXXX,111111111111111\n"
    test_input = (large_row + normal_row * 6)[:-1]

    result = parser.read_csv(StringIO(test_input), header=None, usecols=[0])
    rows = test_input.split("\n")

    expected = DataFrame([row.split(",")[0] for row in rows])
    tm.assert_frame_equal(result, expected)


def test_data_after_quote(c_parser_only):
    # see gh-15910
    parser = c_parser_only

    data = 'a\n1\n"b"a'
    result = parser.read_csv(StringIO(data))

    expected = DataFrame({"a": ["1", "ba"]})
    tm.assert_frame_equal(result, expected)


def test_comment_whitespace_delimited(c_parser_only):
    parser = c_parser_only
    test_input = """\
1 2
2 2 3
3 2 3 # 3 fields
4 2 3# 3 fields
5 2 # 2 fields
6 2# 2 fields
7 # 1 field, NaN
8# 1 field, NaN
9 2 3 # skipped line
# comment"""
    with tm.assert_produces_warning(
        ParserWarning, match="Skipping line", check_stacklevel=False
    ):
        df = parser.read_csv(
            StringIO(test_input),
            comment="#",
            header=None,
            delimiter="\\s+",
            skiprows=0,
            on_bad_lines="warn",
        )
    expected = DataFrame([[1, 2], [5, 2], [6, 2], [7, np.nan], [8, np.nan]])
    tm.assert_frame_equal(df, expected)


def test_file_like_no_next(c_parser_only):
    # gh-16530: the file-like need not have a "next" or "__next__"
    # attribute despite having an "__iter__" attribute.
    #
    # NOTE: This is only true for the C engine, not Python engine.
    class NoNextBuffer(StringIO):
        def __next__(self):
            raise AttributeError("No next method")

        next = __next__

    parser = c_parser_only
    data = "a\n1"

    expected = DataFrame({"a": [1]})
    result = parser.read_csv(NoNextBuffer(data))

    tm.assert_frame_equal(result, expected)


def test_buffer_rd_bytes_bad_unicode(c_parser_only):
    # see gh-22748
    t = BytesIO(b"\xb0")
    t = TextIOWrapper(t, encoding="UTF-8", errors="surrogateescape")
    msg = "'utf-8' codec can't encode character"
    with pytest.raises(UnicodeError, match=msg):
        c_parser_only.read_csv(t, encoding="UTF-8")


@pytest.mark.parametrize("tar_suffix", [".tar", ".tar.gz"])
def test_read_tarfile(c_parser_only, datapath, tar_suffix):
    # see gh-16530
    #
    # Unfortunately, Python's CSV library can't handle
    # tarfile objects (expects string, not bytes when
    # iterating through a file-like).
    parser = c_parser_only
    tar_path = datapath("io", "parser", "data", "tar_csv" + tar_suffix)

    with tarfile.open(tar_path, "r") as tar:
        data_file = tar.extractfile("tar_data.csv")

        out = parser.read_csv(data_file)
        expected = DataFrame({"a": [1]})
        tm.assert_frame_equal(out, expected)


def test_chunk_whitespace_on_boundary(c_parser_only):
    # see gh-9735: this issue is C parser-specific (bug when
    # parsing whitespace and characters at chunk boundary)
    #
    # This test case has a field too large for the Python parser / CSV library.
    parser = c_parser_only

    chunk1 = "a" * (1024 * 256 - 2) + "\na"
    chunk2 = "\n a"
    result = parser.read_csv(StringIO(chunk1 + chunk2), header=None)

    expected = DataFrame(["a" * (1024 * 256 - 2), "a", " a"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.skipif(WASM, reason="limited file system access on WASM")
def test_file_handles_mmap(c_parser_only, csv1):
    # gh-14418
    #
    # Don't close user provided file handles.
    parser = c_parser_only

    with open(csv1, encoding="utf-8") as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as m:
            parser.read_csv(m)
            assert not m.closed


def test_file_binary_mode(c_parser_only, temp_file):
    # see gh-23779
    parser = c_parser_only
    expected = DataFrame([[1, 2, 3], [4, 5, 6]])

    path = temp_file
    with open(path, "w", encoding="utf-8") as f:
        f.write("1,2,3\n4,5,6")

    with open(path, "rb") as f:
        result = parser.read_csv(f, header=None)
        tm.assert_frame_equal(result, expected)


def test_binary_file_handle_avoids_text_wrapping(c_parser_only):
    # GH#46823 - binary file-like objects should not be wrapped in
    # TextIOWrapper when using the C engine, as the small internal buffer
    # causes many small reads that are very slow for remote filesystems.
    parser = c_parser_only
    data = BytesIO(b"a,b\n1,2\n3,4\n")
    result = parser.read_csv(data)
    expected = DataFrame({"a": [1, 3], "b": [2, 4]})
    tm.assert_frame_equal(result, expected)

    # Verify the handle was not wrapped in TextIOWrapper
    data = BytesIO(b"a,b\n1,2\n3,4\n")
    with parser.read_csv(data, chunksize=2) as reader:
        assert not isinstance(reader.handles.handle, TextIOWrapper)


def test_unix_style_breaks(c_parser_only, temp_file):
    # GH 11020
    parser = c_parser_only
    path = temp_file
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write("blah\n\ncol_1,col_2,col_3\n\n")
    result = parser.read_csv(path, skiprows=2, encoding="utf-8", engine="c")
    expected = DataFrame(columns=["col_1", "col_2", "col_3"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("float_precision", [None, "legacy", "high", "round_trip"])
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
def test_1000_sep_with_decimal(
    c_parser_only, data, thousands, decimal, float_precision
):
    parser = c_parser_only
    expected = DataFrame({"A": [1, 10], "B": [2334.01, 13], "C": [5, 10.0]})

    warn = Pandas4Warning if float_precision is not None else None
    with tm.assert_produces_warning(
        warn, match="float_precision", check_stacklevel=False
    ):
        result = parser.read_csv(
            StringIO(data),
            sep="|",
            thousands=thousands,
            decimal=decimal,
            float_precision=float_precision,
        )
    tm.assert_frame_equal(result, expected)


def test_float_precision_options(c_parser_only):
    # GH 17154, 36228
    parser = c_parser_only
    s = "foo\n243.164\n"
    df = parser.read_csv(StringIO(s))
    depr_msg = "float_precision"
    with tm.assert_produces_warning(
        Pandas4Warning, match=depr_msg, check_stacklevel=False
    ):
        df2 = parser.read_csv(StringIO(s), float_precision="high")

    tm.assert_frame_equal(df, df2)

    # "legacy" now uses the same precise converter as "high"
    with tm.assert_produces_warning(
        Pandas4Warning, match=depr_msg, check_stacklevel=False
    ):
        df3 = parser.read_csv(StringIO(s), float_precision="legacy")
    tm.assert_frame_equal(df, df3)

    msg = "Unrecognized float_precision option: junk"

    with pytest.raises(ValueError, match=msg):
        with tm.assert_produces_warning(
            Pandas4Warning, match=depr_msg, check_stacklevel=False
        ):
            parser.read_csv(StringIO(s), float_precision="junk")


# The C tokenizer copies runs of ordinary characters in bulk, scanning 16
# bytes at a time for special characters (delimiter, line terminator, quote,
# escape, comment). Fields shorter than 16 bytes never enter that stride, so
# the tests below deliberately use fields that span the 16-byte boundary to
# exercise the vectorized scan and its byte-offset arithmetic.


@pytest.mark.parametrize("length", [14, 15, 16, 17, 18, 31, 32, 33, 47, 48, 49, 64])
def test_bulk_scan_unquoted_field_boundaries(c_parser_only, length):
    # GH#64515: a delimiter or line terminator ending a long unquoted field
    # must be detected at the correct offset within a 16-byte scan chunk.
    parser = c_parser_only
    col_a = "a" * length
    col_b = "b" * length
    col_c = "c" * length
    data = f"A,B,C\n{col_a},{col_b},{col_c}\n{col_a},{col_b},{col_c}\n"
    result = parser.read_csv(StringIO(data))
    expected = DataFrame(
        {"A": [col_a, col_a], "B": [col_b, col_b], "C": [col_c, col_c]}
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("length", [8, 15, 16, 17, 32, 33, 48])
def test_bulk_scan_quoted_field_boundaries(c_parser_only, length):
    # GH#64515: inside a quoted field only the quote/escape characters are
    # special, so an embedded delimiter and an embedded line terminator that
    # fall past the first 16-byte chunk must be copied verbatim rather than
    # ending the field or the record.
    parser = c_parser_only
    inner = ("a" * length) + "," + ("b" * length) + "\n" + ("c" * length)
    data = 'col\n"' + inner + '"\n'
    result = parser.read_csv(StringIO(data))
    expected = DataFrame({"col": [inner]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("length", [15, 16, 17, 32, 40])
def test_bulk_scan_comment_char_boundary(c_parser_only, length):
    # GH#64515: a comment character terminating a long unquoted field must be
    # detected at its offset within a 16-byte scan chunk.
    parser = c_parser_only
    field = "a" * length
    data = f"A\n{field}# a fairly long trailing comment to skip\n"
    result = parser.read_csv(StringIO(data), comment="#")
    expected = DataFrame({"A": [field]})
    tm.assert_frame_equal(result, expected)


def test_invalid_utf8_raises(c_parser_only):
    # GH#65283: the pyarrow string fast path memcpys raw bytes into arrow
    # buffers; ensure invalid UTF-8 still raises like the object path
    parser = c_parser_only
    data = BytesIO(b"col\nabc\n\xff\xfe\n")
    with pytest.raises(UnicodeDecodeError, match="invalid start byte"):
        parser.read_csv(data)


def test_string_storage_python_consistent(c_parser_only):
    # GH#65283: the pyarrow string fast path must not produce an
    # ArrowStringArray when mode.string_storage="python"
    pytest.importorskip("pyarrow")
    parser = c_parser_only
    with option_context("future.infer_string", True, "mode.string_storage", "python"):
        result = parser.read_csv(StringIO("col\nabc\nxyz\n"))
        arr = result["col"].array
        assert isinstance(arr.dtype, StringDtype)
        assert arr.dtype.storage == "python"
        assert type(arr) is arr.dtype.construct_array_type()


@pytest.mark.parametrize("lineterm", ["\n", "\r\n"])
def test_block_lane_blank_and_whitespace_lines(c_parser_only, lineterm):
    # The SIMD block fast lane defers lines led by blanks (or entirely
    # blank) to the state machine; a no-progress lane exit must preserve
    # the post-WHITESPACE_LINE state instead of ping-ponging (GH#66274).
    parser = c_parser_only
    rows = []
    for i in range(200):
        if i % 7 == 0:
            rows.append("")  # blank line, skipped by skip_blank_lines
        elif i % 11 == 0:
            rows.append(f"  {i},{i * 2},{i * 3}")  # leading whitespace
        elif i % 13 == 0:
            # trailing whitespace is inert for the lane (only line *starts*
            # can begin a WHITESPACE_LINE); included to prove it
            rows.append(f"{i},{i * 2},{i * 3}  ")
        else:
            rows.append(f"{i},{i * 2},{i * 3}")
    data = "a,b,c" + lineterm + lineterm.join(rows) + lineterm
    result = parser.read_csv(StringIO(data))
    expected = read_csv(StringIO(data), engine="python")
    tm.assert_frame_equal(result, expected)


def test_block_lane_chunked_reads_match(c_parser_only):
    # A line-limited tokenize call can stop mid-lane; the parser state must
    # be reset so the next chunk resumes cleanly (GH#66274).
    parser = c_parser_only
    n_rows = 500
    data = "a,b\n" + "\n".join(f"value{i:04d},{i}" for i in range(n_rows)) + "\n"
    expected = parser.read_csv(StringIO(data))
    with parser.read_csv(StringIO(data), chunksize=7) as reader:
        result = concat(reader, ignore_index=True)
    tm.assert_frame_equal(result, expected)
    result = parser.read_csv(StringIO(data), low_memory=True)
    tm.assert_frame_equal(result, expected)


def test_block_lane_crlf_pairs_at_block_edges(c_parser_only):
    # \r\n pairs are consumed in-lane; pairs split across a 16-byte block
    # boundary or bare \r must defer to the state machine.  Vary field
    # widths so terminators land on every block offset.
    parser = c_parser_only
    rows = [f"{'x' * (i % 23)},{i}" for i in range(300)]
    data = "a,b\r\n" + "\r\n".join(rows) + "\r\n"
    result = parser.read_csv(StringIO(data))
    expected = parser.read_csv(StringIO(data.replace("\r\n", "\n")))
    tm.assert_frame_equal(result, expected)


def test_block_lane_quoted_specials_mid_block(c_parser_only):
    # Any block containing a quote falls back to the state machine; quoted
    # fields with embedded delimiters/newlines must round-trip regardless
    # of their offset inside a 16-byte block.
    parser = c_parser_only
    rows = []
    for i in range(200):
        pad = "y" * (i % 19)
        rows.append(f'{pad},"emb,{i}\nnext",{i}')
    data = "a,b,c\n" + "\n".join(rows) + "\n"
    result = parser.read_csv(StringIO(data))
    expected = read_csv(StringIO(data), engine="python")
    tm.assert_frame_equal(result, expected)


def test_block_lane_nrows_short_row(c_parser_only):
    # GH#66274: the lane must honor nrows when the row hitting the limit is
    # a short row (synthetic trailing fields); the over-wide final row raises
    # only if the reader ignores nrows.
    # This was originally tuned so that the short rows outran the stream
    # reservation and the lane's capacity guard fired mid-parse.  GH#66657
    # made the padding restore that reservation, so the guard can no longer
    # trip here -- the nrows/short-row path is what this still covers.
    parser = c_parser_only
    data = "a,b,c,d\n" + "2222\n" + "1\n" * 55_000 + "1,2,3,4,5,6,7,8,9,10\n"
    for nrows in range(52_415, 52_431):
        result = parser.read_csv(StringIO(data), nrows=nrows)
        assert len(result) == nrows


def test_block_lane_nrows_checked_before_deferring(c_parser_only):
    # GH#66274: block_after_line must test the line limit before any of its
    # bails, because end_line has already counted the line and the state
    # machine's ==-based limit check would never match again -- so a deferral
    # taken first makes the reader overshoot nrows.  Reaches the deferral via
    # the whitespace-line probe (the leading blank on row 5); the sibling bail
    # for a full token stream can no longer coincide with a limit hit now that
    # GH#66657 keeps slack ahead of the unconsumed input.
    # nrows is pinned: 4 stops exactly on the whitespace row, which is the
    # only value that both reaches the probe and stays clear of the over-wide
    # final row (that row legitimately raises for nrows > 5).
    parser = c_parser_only
    over_wide = ",".join(str(k) for k in range(12))
    data = "a,b,c,d\n" + "1,2,3,4\n" * 4 + " 9,9,9\n" + over_wide + "\n"
    result = parser.read_csv(StringIO(data), nrows=4)
    expected = DataFrame(
        [[1, 2, 3, 4]] * 4, columns=["a", "b", "c", "d"], dtype="int64"
    )
    tm.assert_frame_equal(result, expected)


def test_short_rows_before_wide_field(c_parser_only):
    # GH#66657: closing the missing fields of a short row takes token stream
    # space without consuming input, so a run of short rows used to eat into
    # the reservation the unchecked bulk copies rely on.  The wide field then
    # overran the stream buffer and came back truncated.
    parser = c_parser_only
    wide = "X" * 100_000
    data = "a\n" * 40 + wide + "\n"
    names = [f"c{i}" for i in range(5000)]
    result = parser.read_csv(StringIO(data), header=None, names=names)
    assert result.shape == (41, 5000)
    # length first: the truncation is off-by-one, and a bare inequality on a
    # 100k-char field reports nothing useful
    assert len(result.iloc[40, 0]) == len(wide)
    assert result.iloc[40, 0] == wide


@pytest.mark.parametrize("terminator", ["\n", "\r\n", "\r"])
@pytest.mark.parametrize("n_cols", range(5, 12))
@pytest.mark.parametrize("n_short", [5, 8, 20])
def test_short_rows_do_not_exhaust_token_reservation(
    c_parser_only, n_cols, n_short, terminator
):
    # GH#66657: the padding also consumes one word per closed field, so short
    # rows ate the word-vector reservation the same way they ate the stream's,
    # and end_field has no growth path -- it just reports "Buffer overflow
    # caught - possible malformed input file" on ordinary ragged input.
    # Swept rather than pinned to the handful of shapes that reproduced: which
    # ones do depends on where the buffers' power-of-two growth lands, so a
    # pinned pair would quietly stop covering this if that sizing changed.
    # Each terminator reaches end_line through a call site with its own
    # expression for the unconsumed input, so an off-by-one in one of them
    # survives testing the others.
    parser = c_parser_only
    header = ",".join(f"c{i}" for i in range(n_cols))
    row = ",".join("1" for _ in range(n_cols))
    body = header + "{0}" + "1{0}" * n_short + row + "{0}"
    result = parser.read_csv(StringIO(body.format(terminator)))
    # the terminator must not change the values; a lone \r is not readable
    # through StringIO by the python engine, so oracle on the \n spelling
    expected = read_csv(StringIO(body.format("\n")), engine="python")
    tm.assert_frame_equal(result, expected)


def test_short_rows_block_lane_tail_recopy_defer(c_parser_only):
    # Not a GH#66657 reproducer -- this passes without that fix.  It pins the
    # one shape that reaches the block lane's tail-re-copy deferral, which
    # GH#66657 turned from a rare capacity bail into the routine path: the
    # lane hands the rest of the block back to the state machine mid-block,
    # and that handoff carries a resume position, so getting it wrong
    # duplicates or drops rows rather than raising.  The sweep above never
    # reaches the branch.
    parser = c_parser_only
    data = "c0,c1,c2,c3\n" + "1\n" * 9
    result = parser.read_csv(StringIO(data))
    expected = read_csv(StringIO(data), engine="python")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("kwargs", [{}, {"dtype_backend": "pyarrow"}])
def test_pyarrow_string_non_ascii_eager_wrap_mutable(c_parser_only, kwargs):
    # GH#66277: the eager wrap that bypasses the ExtensionArray constructor
    # for non-ASCII columns did not set _cache, so mutating the result raised
    # AttributeError (reachable with low_memory=False, where the array
    # reaches the frame without an intermediate concat)
    pytest.importorskip("pyarrow")
    parser = c_parser_only
    result = parser.read_csv(StringIO("a\ncafé\nnaïve\n"), **kwargs)
    arr = result["a"].array
    arr[0] = "x"
    arr.sort()
    assert list(arr) == ["naïve", "x"]


@pytest.mark.parametrize("tail", ["", "café\nnaïve\n" * 4])
def test_low_memory_string_chunks_combined(c_parser_only, monkeypatch, tail):
    # GH#66277: with low_memory=True, deferred string chunks combine into a
    # single column at the end of the read; all-ASCII chunks take the
    # zero-copy path while non-ASCII chunks are wrapped eagerly, and either
    # mix must match a single-chunk read
    pytest.importorskip("pyarrow")
    parser = c_parser_only
    data = "a\n" + "foo\nbar\nNA\nbaz\n" * 8 + tail
    expected = parser.read_csv(StringIO(data))
    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", 2**3)
        result = parser.read_csv(StringIO(data))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("kwargs", [{}, {"dtype_backend": "pyarrow"}])
def test_pyarrow_string_fast_path_mutable(kwargs):
    # GH#66619: the fast path builds its result without going through the
    # ExtensionArray constructor, so it must set every attribute the
    # constructor does; omitting _cache made mutating the result raise
    # AttributeError.  low_memory=False is required, not incidental: the
    # low-memory path concatenates its chunks, which rebuilds the array and
    # would hide the omission.
    pytest.importorskip("pyarrow")
    # pinned rather than inherited: the default-kwargs case would otherwise get
    # an object-dtype column, and stop exercising the fast path at all, in the
    # PANDAS_FUTURE_INFER_STRING=0 build.
    with option_context("future.infer_string", True):
        result = read_csv(
            StringIO("a\nfoo\nbar\n"), engine="c", low_memory=False, **kwargs
        )
    arr = result["a"].array
    arr[0] = "zzz"
    arr.sort()
    assert list(arr) == ["bar", "zzz"]


@pytest.mark.parametrize("kwargs", [{}, {"dtype_backend": "pyarrow"}])
def test_pyarrow_string_fast_path_attrs_match_constructor(kwargs):
    # GH#66619: the fast path sets the instance attributes itself instead of
    # calling __init__, so it has to track whatever set the constructor
    # establishes.  Adding an attribute to ArrowStringArray.__init__ or
    # ArrowExtensionArray.__init__ without teaching parsers.pyx about it should
    # fail here rather than silently producing a half-built array.
    pytest.importorskip("pyarrow")
    with option_context("future.infer_string", True):
        result = read_csv(
            StringIO("a\nfoo\nbar\n"), engine="c", low_memory=False, **kwargs
        )
    arr = result["a"].array
    expected = type(arr)(arr._pa_array)
    assert vars(arr).keys() == vars(expected).keys()


def test_pyarrow_string_iterator_dtype_stable_across_chunks():
    # GH#66619: a reader resolves its pyarrow target once, when it converts its
    # first string column, so every chunk of one read gets the same dtype even
    # if the options change mid-iteration.  Previously the target was looked up
    # per chunk and the second chunk here came back object-dtype.
    pytest.importorskip("pyarrow")
    with option_context("future.infer_string", True):
        reader = read_csv(
            StringIO("a\nfoo\nbar\n"), engine="c", chunksize=1, iterator=True
        )
        first = next(reader)
    with option_context("future.infer_string", False):
        second = next(reader)
    assert first["a"].dtype == StringDtype(na_value=np.nan)
    assert second["a"].dtype == first["a"].dtype


@pytest.mark.parametrize("kwargs", [{}, {"dtype_backend": "pyarrow"}])
def test_pyarrow_string_fast_path_token_width_tiers(kwargs):
    # GH#66756: the fast path copies a short token at a compile-time-constant
    # 16 or 32 bytes and lets the copy overshoot into buffer slack, so a token
    # one byte either side of a tier boundary is where a mis-sized copy would
    # truncate the value or trail the following token's bytes into it.
    pa = pytest.importorskip("pyarrow")
    widths = [1, 2, 15, 16, 17, 31, 32, 33, 64]
    values = [
        "".join(chr(ord("a") + pos % 26) for pos in range(width)) for width in widths
    ]
    data = "a,b\n" + "".join(f"{value},{value.upper()}\n" for value in values)
    # The fast path requires infer_string *and* pyarrow storage, so pin both;
    # the dtype check below turns any silent fall-back to the object path
    # (which would satisfy the value assertions) into a loud failure.
    with option_context("future.infer_string", True, "mode.string_storage", "pyarrow"):
        result = read_csv(StringIO(data), engine="c", low_memory=False, **kwargs)
    expected_dtype = (
        ArrowDtype(pa.string()) if kwargs else StringDtype("pyarrow", na_value=np.nan)
    )
    assert result["a"].dtype == expected_dtype
    assert result["a"].tolist() == values
    assert result["b"].tolist() == [value.upper() for value in values]


@pytest.mark.parametrize("kwargs", [{}, {"dtype_backend": "pyarrow"}])
def test_pyarrow_string_fast_path_column_outgrows_size_estimate(kwargs):
    # GH#66756: the fast path sizes its data buffer from the column's leading
    # tokens and grows it mid-pass when that estimate falls short, re-copying
    # what it has already written.  A column whose first rows are far narrower
    # than the rest takes that path repeatedly, where a stale buffer pointer or
    # an undersized grow would corrupt everything written before it.
    pa = pytest.importorskip("pyarrow")
    values = ["ab"] * 20 + [f"{num:x}" * 900 for num in range(1, 200)]
    data = "a\n" + "".join(f"{value}\n" for value in values)
    with option_context("future.infer_string", True, "mode.string_storage", "pyarrow"):
        result = read_csv(StringIO(data), engine="c", low_memory=False, **kwargs)
    expected_dtype = (
        ArrowDtype(pa.string()) if kwargs else StringDtype("pyarrow", na_value=np.nan)
    )
    assert result["a"].dtype == expected_dtype
    assert result["a"].tolist() == values


@pytest.mark.parametrize("kwargs", [{}, {"dtype_backend": "pyarrow"}])
@pytest.mark.parametrize("prefix_len", [1, 200])
def test_embedded_nul_byte_roundtrip(c_parser_only, kwargs, prefix_len):
    # GH#66415: the pyarrow string fast path computed token lengths with
    # strlen, so a quoted field with an embedded NUL byte was truncated at the
    # NUL.  Column "a" takes its length from the next token's start and column
    # "b", last in the stream, from the stream end, covering both branches of
    # the length helper; the two prefix_len values keep a length-capped scan
    # from passing by accident.
    pytest.importorskip("pyarrow")
    parser = c_parser_only
    value = b"x" * prefix_len + b"\x00y"
    data = b'a,b\n"' + value + b'","' + value + b'"'
    result = parser.read_csv(BytesIO(data), **kwargs)
    # engine="python" shares none of the C tokenizer's length arithmetic, so it
    # is an independent reference for what the field should decode to
    expected = read_csv(BytesIO(data), engine="python")
    assert result["a"][0] == value.decode()
    assert result["b"][0] == value.decode()
    assert expected["a"][0] == value.decode()
    assert expected["b"][0] == value.decode()


def test_embedded_nul_fixed_width_bytes(c_parser_only):
    # GH#19886: the fixed-width "S" path copied with strncpy, which stops at an
    # embedded NUL even though numpy "S" can hold one.
    parser = c_parser_only
    data = b'a\n"x\x00y"\n"x\x00z"\n'

    result = parser.read_csv(BytesIO(data), dtype="S5")
    assert result["a"].tolist() == [b"x\x00y", b"x\x00z"]
    tm.assert_frame_equal(result, read_csv(BytesIO(data), dtype="S5", engine="python"))

    # a token longer than the width is still truncated to the width
    assert parser.read_csv(BytesIO(b"a\nabcdef\n"), dtype="S3")["a"].tolist() == [
        b"abc"
    ]


def test_embedded_nul_converter(c_parser_only):
    # GH#19886: converters were handed a value truncated at the NUL, so they
    # disagreed with the object path on the same input.
    parser = c_parser_only
    data = b'a\n"x\x00y"\n"x\x00z"\n'

    result = parser.read_csv(BytesIO(data), converters={"a": str})
    assert result["a"].tolist() == ["x\x00y", "x\x00z"]
    tm.assert_frame_equal(
        result, read_csv(BytesIO(data), converters={"a": str}, engine="python")
    )


def test_embedded_nul_column_name(c_parser_only):
    # GH#19886: a column name was truncated at an embedded NUL, which could
    # also collide two distinct names into one.
    parser = c_parser_only
    data = b'"h\x001","h\x002"\n1,2\n'

    result = parser.read_csv(BytesIO(data))
    assert list(result.columns) == ["h\x001", "h\x002"]
    tm.assert_frame_equal(result, read_csv(BytesIO(data), engine="python"))


@pytest.mark.parametrize("dtype", [object, "str", "string", "category", None])
@pytest.mark.parametrize(
    "data, expected",
    [
        (b'a\n"x\x00y"\n"x\x00z"\n"x"\n', ["x\x00y", "x\x00z", "x"]),
        (b'a\n"x\x00y"\n"x\x00z"\n"x\x00w"\n', ["x\x00y", "x\x00z", "x\x00w"]),
    ],
)
def test_embedded_nul_distinct_values(c_parser_only, dtype, data, expected):
    # GH#66525: the object path interned fields in a hash table keyed on the
    # NUL-terminated word, so two fields differing only past an embedded NUL
    # were silently boxed to the same value.  The second fixture is all
    # equal-length, so distinguishing it needs a comparison that reads past
    # the NUL rather than just a key-length check.
    parser = c_parser_only

    result = parser.read_csv(BytesIO(data), dtype=dtype, keep_default_na=False)
    assert list(result["a"]) == expected


@pytest.mark.parametrize(
    "data",
    [
        b'a\n"q\x00\xff"\n',
        b'a\n"q\xff"\n',
        # two *distinct* undecodable fields: the table dedupes on raw bytes, so
        # these stay separate keys but decode to one label
        b'a\n"q\xff"\n"q\xfe"\n',
        b'a\n"q\xff"\n"q\xfe"\n"z"\n"q\xff"\n',
    ],
)
@pytest.mark.parametrize("encoding_errors", ["replace", "ignore"])
def test_categorical_honors_encoding_errors(c_parser_only, data, encoding_errors):
    # GH#66525: _categorical_convert decoded its category labels with strict
    # errors regardless of encoding_errors, so an undecodable byte raised where
    # dtype=object honored the argument.  Two keys that decode to the same label
    # must merge, or the Categorical is built with duplicate categories.
    parser = c_parser_only
    kwargs = {"keep_default_na": False, "encoding_errors": encoding_errors}

    result = parser.read_csv(BytesIO(data), dtype="category", **kwargs)
    expected = parser.read_csv(BytesIO(data), dtype=object, **kwargs)
    assert list(result["a"]) == list(expected["a"])


def test_categorical_encoding_errors_merge_with_na(c_parser_only):
    # GH#66525: a merged column that also contains NA exercises the sentinel
    # guard in the code remap -- a -1 code must not be looked up in the map.
    parser = c_parser_only
    data = b'a,b\n"q\xff",1\n,2\n"q\xfe",3\n"z",4\n'

    result = parser.read_csv(
        BytesIO(data), dtype={"a": "category"}, encoding_errors="replace"
    )["a"]

    assert list(result.cat.categories) == ["q�", "z"]
    assert list(result.cat.codes) == [0, -1, 0, 1]


@pytest.mark.parametrize(
    "field,other",
    [
        ("1\x00xyz", "2"),
        ("-1\x00xyz", "2"),
        ("18446744073709551615\x00xyz", "2"),
        ("1.5\x00xyz", "2.5"),
        ("1e3\x00xyz", "2.5"),
        ("inf\x00xyz", "2.5"),
        ("infinity\x00xyz", "2.5"),
        ("True\x00xyz", "False"),
        # not a default true_values entry, so only to_boolean can accept it
        ("TRue\x00xyz", "False"),
    ],
)
def test_embedded_nul_is_not_a_numeric_or_boolean_literal(c_parser_only, field, other):
    # GH#66524: the numeric and boolean converters finished on a NUL rather
    # than on the end of the token, so a field was silently accepted at its
    # pre-NUL prefix and the trailing bytes were discarded.  `other` keeps the
    # rest of the column parseable, so the column would convert if the bad
    # field were accepted.
    parser = c_parser_only
    data = f'a\n"{field}"\n{other}\n'.encode()

    result = parser.read_csv(BytesIO(data))
    expected = read_csv(BytesIO(data), engine="python")
    tm.assert_frame_equal(result, expected)
    assert result["a"][0] == field


def test_embedded_nul_with_thousands_separator(c_parser_only):
    # GH#66524: the thousands-separator path strips the separator into a scratch
    # buffer, so it needs the same end-of-token check as the plain path.
    parser = c_parser_only
    data = b'a\n"1,234\x00xyz"\n2\n'

    result = parser.read_csv(BytesIO(data), thousands=",")
    expected = read_csv(BytesIO(data), engine="python", thousands=",")
    tm.assert_frame_equal(result, expected)
    assert result["a"][0] == "1,234\x00xyz"


@pytest.mark.parametrize("na_filter", [True, False])
@pytest.mark.parametrize(
    "good,field",
    [
        ("2.5", "1.5\x00xyz"),
        ("2.5", "inf\x00xyz"),
        # the long spelling takes a different arm of the infinity check
        ("2.5", "-Infinity\x00x"),
        # default spelling: covers the true/false hashset lookup, which must
        # compare the full length rather than the pre-NUL prefix
        ("True", "True\x00xyz"),
        ("True", "TRue\x00xyz"),
        # neither spelling is in the default true_values/false_values, so each
        # reaches a different arm of to_boolean
        ("True", "FAlse\x00xyz"),
    ],
)
def test_embedded_nul_in_later_row(c_parser_only, good, field, na_filter):
    # GH#66524: the float and boolean converters probe only the first non-NA
    # token and bail out before their bulk loop when it rejects, so a NUL field
    # in the first row never reaches the per-row conversion.  Put a good value
    # first so the bulk loop is the code under test; na_filter picks between
    # the two separate loops in each converter.
    parser = c_parser_only
    data = f'a\n{good}\n"{field}"\n'.encode()

    result = parser.read_csv(BytesIO(data), na_filter=na_filter)
    expected = read_csv(BytesIO(data), engine="python", na_filter=na_filter)
    tm.assert_frame_equal(result, expected)
    assert result["a"].tolist() == [good, field]


@pytest.mark.parametrize(
    "field",
    [
        "1\x00xyz",
        "18446744073709551614\x00z",
        # the pre-NUL digits exceed uint64, so this leaves str_to_uint64 by the
        # overflow arm rather than by the end-of-token check the others take
        "18446744073709551616\x00z",
    ],
)
def test_embedded_nul_is_not_a_uint64(c_parser_only, field):
    # GH#66524: str_to_uint64 is reached only after str_to_int64 reports a
    # *clean* overflow, which a NUL-bearing token can never produce, so the
    # leading row has to genuinely exceed int64 for the uint64 path to see the
    # NUL-bearing field at all.
    parser = c_parser_only
    data = f'a\n18446744073709551615\n"{field}"\n'.encode()

    result = parser.read_csv(BytesIO(data))
    expected = read_csv(BytesIO(data), engine="python")
    tm.assert_frame_equal(result, expected)
    assert result["a"].tolist() == ["18446744073709551615", field]


def test_embedded_nul_int64_overflow(c_parser_only):
    # GH#66524: the pre-NUL digits are int64max + 1, so the token leaves
    # str_to_int64 by the overflow arm and its truncation has to be caught
    # before the uint64 retry, where the truncated value would fit.
    parser = c_parser_only
    data = b'a\n1\n"9223372036854775808\x00z"\n'

    result = parser.read_csv(BytesIO(data))
    expected = read_csv(BytesIO(data), engine="python")
    tm.assert_frame_equal(result, expected)
    assert result["a"].tolist() == ["1", "9223372036854775808\x00z"]


def test_embedded_nul_is_not_a_python_int(c_parser_only):
    # GH#66524: a column containing a value too large for uint64 falls back to
    # a Python-int path built on PyLong_FromString, which also stops at a NUL.
    parser = c_parser_only
    data = b'a\n99999999999999999999999999\n"1\x00xyz"\n'

    result = parser.read_csv(BytesIO(data))
    expected = read_csv(BytesIO(data), engine="python")
    tm.assert_frame_equal(result, expected)
    assert result["a"].tolist() == ["99999999999999999999999999", "1\x00xyz"]


def test_embedded_nul_raises_for_explicit_int_dtype(c_parser_only):
    # GH#66524: with the dtype pinned there is no string column to fall back
    # to, so the truncated value has to raise rather than parse.
    parser = c_parser_only
    data = b'a\n"1\x00xyz"\n2\n'

    with pytest.raises(ValueError, match="Unable to parse string"):
        parser.read_csv(BytesIO(data), dtype="Int64")


@pytest.mark.parametrize(
    "field", ["1\x00x", "0\x00x", "1.0\x00x", "0.0\x00x", "True\x00x"]
)
def test_embedded_nul_raises_for_boolean_dtype(c_parser_only, field):
    # GH#66524: "1"/"1.0"/"0"/"0.0" are the numeric spellings dtype="boolean"
    # accepts; none of these are one of them.  The "1.0"/"0.0" spellings take a
    # separate arm of the literal check from the one-character ones.
    parser = c_parser_only
    data = f'a\nTrue\n"{field}"\n'.encode()

    with pytest.raises(ValueError, match="cannot be cast to bool"):
        parser.read_csv(BytesIO(data), dtype="boolean")


@pytest.mark.parametrize("value", [b"NA\x00x", b"nan\x00junk", b"null\x00z"])
def test_default_na_value_prefix_is_not_na(c_parser_only, value):
    # GH#19886: a field was compared against na_values only up to its first
    # NUL, so a value merely *starting* with a default na_value -- needing no
    # custom na_values at all -- was read as NaN.
    parser = c_parser_only
    data = b'a\n"' + value + b'"\n'

    result = parser.read_csv(BytesIO(data))
    assert result["a"][0] == value.decode()
    tm.assert_frame_equal(result, read_csv(BytesIO(data), engine="python"))


@pytest.mark.parametrize("value", [b"\x00y", b"\x00\x00\x00", b"\x00", b"\x00 "])
def test_leading_nul_is_not_na(c_parser_only, value):
    # GH#19886: the na_values lookup compared the NUL-terminated word, so a
    # field starting with a NUL byte matched the empty string -- a default
    # na_value -- and was read as NaN.
    parser = c_parser_only
    data = b'a\n"' + value + b'"\n'

    result = parser.read_csv(BytesIO(data))
    assert result["a"][0] == value.decode()
    tm.assert_frame_equal(result, parser.read_csv(BytesIO(data), na_filter=False))
    tm.assert_frame_equal(result, read_csv(BytesIO(data), engine="python"))


def test_na_values_with_embedded_nul(c_parser_only):
    # GH#19886: an na_value containing a NUL was itself truncated when added to
    # the hashset, so it matched any field sharing its pre-NUL prefix.
    parser = c_parser_only
    data = b'a\n"x\x00y"\n"x\x00z"\n"x"\n'

    result = parser.read_csv(BytesIO(data), na_values=["x\x00y"], keep_default_na=False)
    assert result["a"].isna().tolist() == [True, False, False]
    assert result["a"][1] == "x\x00z"
    assert result["a"][2] == "x"


def test_true_false_values_with_embedded_nul(c_parser_only):
    # GH#19886: a true_values/false_values entry containing a NUL was truncated
    # when added to the hashset, so it also matched a field equal to just the
    # prefix before that NUL. The field here must be prefix-only to be
    # load-bearing -- an exactly-matching field behaves the same either way.
    parser = c_parser_only

    result = parser.read_csv(
        BytesIO(b"a\ny\nno\n"), true_values=["y\x00es"], false_values=["no"]
    )
    assert result["a"].tolist() == ["y", "no"]

    matched = parser.read_csv(
        BytesIO(b'a\n"y\x00es"\nno\n'), true_values=["y\x00es"], false_values=["no"]
    )
    assert matched["a"].tolist() == [True, False]


def test_na_values_leading_nul(c_parser_only):
    # GH#19886: exercises the first-byte prefilter for keys under '\x00' --
    # a leading-NUL na_value must not swallow the empty field or a different
    # leading-NUL value.
    parser = c_parser_only
    data = b'a\n"\x00y"\n"\x00z"\n""\n'

    result = parser.read_csv(BytesIO(data), na_values=["\x00y"], keep_default_na=False)
    assert result["a"].isna().tolist() == [True, False, False]
    assert result["a"][1] == "\x00z"
    assert result["a"][2] == ""


def _raise_on_five(value):
    if value == "5":
        raise RuntimeError("boom")
    return value


@pytest.mark.parametrize(
    "data,kwargs,error,match",
    [
        (
            "a,b\n" + "".join(f"{i},{'oops' if i == 5 else i}\n" for i in range(20)),
            {"dtype": {"b": "int64"}},
            ValueError,
            "invalid literal for int",
        ),
        (
            "a,b\n" + "".join(f"{i},{i}\n" for i in range(20)),
            {"converters": {"b": _raise_on_five}},
            RuntimeError,
            "boom",
        ),
        (
            "a,b\n"
            + "".join(
                (f"{i},{i},{i}\n" if i == 15 else f"{i},{i}\n") for i in range(20)
            ),
            {},
            ParserError,
            "Expected 2 fields",
        ),
    ],
)
def test_read_after_chunk_raised(c_parser_only, data, kwargs, error, match):
    # GH#66622: a chunk that raised closed the reader, and reading again then
    # dereferenced the freed tokenizer buffers instead of raising.
    parser = c_parser_only

    with parser.read_csv(StringIO(data), chunksize=10, **kwargs) as reader:
        with pytest.raises(error, match=match):
            for _ in reader:
                pass

        with pytest.raises(ValueError, match="I/O operation on closed file"):
            next(reader)


def test_read_after_close(c_parser_only):
    # GH#66622: reading from a closed reader crashed the interpreter.
    parser = c_parser_only
    data = "a,b\n" + "".join(f"{i},{i}\n" for i in range(20))

    with parser.read_csv(StringIO(data), chunksize=10) as reader:
        assert len(next(reader)) == 10

    with pytest.raises(ValueError, match="I/O operation on closed file"):
        next(reader)


def test_exhausted_reader_keeps_raising_stop_iteration(c_parser_only):
    # GH#66622: exhausting a reader closes it, but it must keep behaving like a
    # spent iterator rather than reporting a closed file.
    parser = c_parser_only
    data = "a,b\n" + "".join(f"{i},{i}\n" for i in range(20))

    reader = parser.read_csv(StringIO(data), chunksize=10)
    assert len(list(reader)) == 2

    for _ in range(2):
        with pytest.raises(StopIteration):
            next(reader)
