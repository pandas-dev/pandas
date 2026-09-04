"""
Tests that work on both the Python and C engines but do not have a
specific classification into the other test modules.
"""

from io import StringIO

import numpy as np
import pytest

from pandas._libs import parsers as libparsers
from pandas.errors import (
    DtypeWarning,
    ParserError,
)

import pandas as pd
import pandas._testing as tm

skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")


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

    expected = pd.DataFrame(
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

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            with parser.read_csv(StringIO(data), index_col=0, chunksize=2) as reader:
                list(reader)
        return

    with parser.read_csv(StringIO(data), index_col=0, chunksize=2) as reader:
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
    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"

    with pytest.raises(ValueError, match=msg):
        with parser.read_csv(StringIO(data), chunksize=chunksize) as _:
            pass


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
    kwargs = {"index_col": 0, "nrows": 5}

    if parser.engine == "pyarrow":
        msg = "The 'nrows' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), **kwargs)
        return

    expected = parser.read_csv(StringIO(data), **kwargs)
    with parser.read_csv(StringIO(data), chunksize=chunksize, **kwargs) as reader:
        tm.assert_frame_equal(pd.concat(reader), expected)


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
    kwargs = {"index_col": 0, "nrows": 5}

    if parser.engine == "pyarrow":
        msg = "The 'nrows' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), **kwargs)
        return

    expected = parser.read_csv(StringIO(data), **kwargs)
    with parser.read_csv(StringIO(data), chunksize=8, **kwargs) as reader:
        tm.assert_frame_equal(reader.get_chunk(size=2), expected.iloc[:2])
        tm.assert_frame_equal(reader.get_chunk(size=4), expected.iloc[2:5])

        with pytest.raises(StopIteration, match="^$"):
            reader.get_chunk(size=3)


def test_get_chunk_passed_chunksize(all_parsers):
    parser = all_parsers
    data = """A,B,C
1,2,3
4,5,6
7,8,9
1,2,3"""

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            with parser.read_csv(StringIO(data), chunksize=2) as reader:
                reader.get_chunk()
        return

    with parser.read_csv(StringIO(data), chunksize=2) as reader:
        result = reader.get_chunk()

    expected = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["A", "B", "C"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("kwargs", [{}, {"index_col": 0}])
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
    result = parser.read_csv(StringIO(data), **kwargs)

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            with parser.read_csv(StringIO(data), chunksize=2, **kwargs) as reader:
                pd.concat(reader)
        return

    with parser.read_csv(StringIO(data), chunksize=2, **kwargs) as reader:
        via_reader = pd.concat(reader)
    tm.assert_frame_equal(via_reader, result)


def test_read_chunksize_jagged_names(all_parsers):
    # see gh-23509
    parser = all_parsers
    data = "\n".join(["0"] * 7 + [",".join(["0"] * 10)])

    expected = pd.DataFrame([[0] + [np.nan] * 9] * 7 + [[0] * 10])

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            with parser.read_csv(
                StringIO(data), names=range(10), chunksize=4
            ) as reader:
                pd.concat(reader)
        return

    with parser.read_csv(StringIO(data), names=range(10), chunksize=4) as reader:
        result = pd.concat(reader)
    tm.assert_frame_equal(result, expected)


def test_chunk_begins_with_newline_whitespace(all_parsers):
    # see gh-10022
    parser = all_parsers
    data = "\n hello\nworld\n"

    result = parser.read_csv(StringIO(data), header=None)
    expected = pd.DataFrame([" hello", "world"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.slow
def test_chunks_have_consistent_numerical_type(all_parsers, monkeypatch):
    # mainly an issue with the C parser
    heuristic = 2**3
    parser = all_parsers
    integers = [str(i) for i in range(heuristic - 1)]
    data = "a\n" + "\n".join([*integers, "1.0", "2.0", *integers])

    # Coercions should work without warnings.
    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
        result = parser.read_csv(StringIO(data))

    assert type(result.a[0]) is np.float64
    assert result.a.dtype == float


def test_warn_if_chunks_have_mismatched_type(
    all_parsers, using_infer_string, monkeypatch
):
    warning_type = None
    parser = all_parsers
    heuristic = 2**3
    size = 10

    # see gh-3866: if chunks are different types and can't
    # be coerced using numerical types, then issue warning.
    if parser.engine == "c" and parser.low_memory:
        warning_type = DtypeWarning
        # Use a size to hit warning path dictated by DEFAULT_BUFFER_HEURISTIC
        # monkeypatched below
        size = heuristic - 1

    integers = [str(i) for i in range(size)]
    data = "a\n" + "\n".join([*integers, "a", "b", *integers])

    buf = StringIO(data)

    if parser.engine == "pyarrow":
        df = parser.read_csv(
            buf,
        )
    else:
        with monkeypatch.context() as m:
            m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
            df = parser.read_csv_check_warnings(
                warning_type,
                r"Columns \(0: a\) have mixed types. "
                "Specify dtype option on import or set low_memory=False.",
                buf,
            )
    if parser.engine == "c" and parser.low_memory:
        assert df.a.dtype == object
    elif using_infer_string:
        assert df.a.dtype == "str"
    else:
        assert df.a.dtype == object


@pytest.mark.parametrize("iterator", [True, False])
def test_empty_with_nrows_chunksize(all_parsers, iterator):
    # see gh-9535
    parser = all_parsers
    expected = pd.DataFrame(columns=["foo", "bar"])

    nrows = 10
    data = StringIO("foo,bar\n")

    if parser.engine == "pyarrow":
        msg = (
            "The '(nrows|chunksize)' option is not supported with the 'pyarrow' engine"
        )
        with pytest.raises(ValueError, match=msg):
            if iterator:
                with parser.read_csv(data, chunksize=nrows) as reader:
                    next(iter(reader))
            else:
                parser.read_csv(data, nrows=nrows)
        return

    if iterator:
        with parser.read_csv(data, chunksize=nrows) as reader:
            result = next(iter(reader))
    else:
        result = parser.read_csv(data, nrows=nrows)

    tm.assert_frame_equal(result, expected)


def test_chunksize_with_usecols_second_block_shorter(all_parsers):
    # GH#21211
    parser = all_parsers
    data = """1,2,3,4
5,6,7,8
9,10,11
"""

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data),
                names=["a", "b"],
                chunksize=2,
                usecols=[0, 1],
                header=None,
            )
        return

    result_chunks = parser.read_csv(
        StringIO(data),
        names=["a", "b"],
        chunksize=2,
        usecols=[0, 1],
        header=None,
    )

    expected_frames = [
        pd.DataFrame({"a": [1, 5], "b": [2, 6]}),
        pd.DataFrame({"a": [9], "b": [10]}, index=[2]),
    ]

    for i, result in enumerate(result_chunks):
        tm.assert_frame_equal(result, expected_frames[i])


def test_chunksize_second_block_shorter(all_parsers):
    # GH#21211
    parser = all_parsers
    data = """a,b,c,d
1,2,3,4
5,6,7,8
9,10,11
"""

    if parser.engine == "pyarrow":
        msg = "The 'chunksize' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), chunksize=2)
        return

    result_chunks = parser.read_csv(StringIO(data), chunksize=2)

    expected_frames = [
        pd.DataFrame({"a": [1, 5], "b": [2, 6], "c": [3, 7], "d": [4, 8]}),
        pd.DataFrame({"a": [9], "b": [10], "c": [11], "d": [np.nan]}, index=[2]),
    ]

    for i, result in enumerate(result_chunks):
        tm.assert_frame_equal(result, expected_frames[i])


@skip_pyarrow  # MultiIndex columns are not supported by the pyarrow engine
@pytest.mark.parametrize("run_start", [3, 4, 5, 6])
def test_short_line_run_at_chunk_boundary(all_parsers, monkeypatch, run_start):
    # GH#40587: with the C engine and low_memory=True, a run of consecutive
    # lines with too few fields crossing an internal chunk boundary raised
    # ParserError instead of NaN-padding the missing fields; parametrized so
    # the run crosses the 4-row boundary at every possible offset
    parser = all_parsers
    heuristic = 2**5  # 4-row chunks at this table width
    rows = ["0,1,2,3,4"] * 24
    rows[run_start] = rows[run_start + 1] = "0,1,2"
    data = "a,a,a,b,b\nx,y,z,w,v\n" + "\n".join(rows) + "\n"

    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
        result = parser.read_csv(StringIO(data), header=[0, 1])

    expected = pd.DataFrame(
        {
            "0": [0] * 24,
            "1": [1] * 24,
            "2": [2] * 24,
            "3": [3.0] * 24,
            "4": [4.0] * 24,
        }
    )
    expected.iloc[[run_start, run_start + 1], 3:] = np.nan
    expected.columns = pd.MultiIndex.from_arrays(
        [["a", "a", "a", "b", "b"], ["x", "y", "z", "w", "v"]]
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # the pyarrow engine does not NaN-pad short lines
@pytest.mark.parametrize("run_start", [3, 4, 5, 6])
def test_short_line_run_at_chunk_boundary_single_header(
    all_parsers, monkeypatch, run_start
):
    # GH#40587: the single-header variant of the case above never hit the
    # bug (only a multi-row header leaves buffer slots exempt from padding),
    # so this guards the fixed header bookkeeping in the common case instead
    parser = all_parsers
    heuristic = 2**5  # 4-row chunks at this table width
    rows = ["0,1,2,3,4"] * 24
    rows[run_start] = rows[run_start + 1] = "0,1,2"
    data = "a,b,c,d,e\n" + "\n".join(rows) + "\n"

    with monkeypatch.context() as m:
        m.setattr(libparsers, "DEFAULT_BUFFER_HEURISTIC", heuristic)
        result = parser.read_csv(StringIO(data))

    expected = pd.DataFrame(
        {
            "a": [0] * 24,
            "b": [1] * 24,
            "c": [2] * 24,
            "d": [3.0] * 24,
            "e": [4.0] * 24,
        }
    )
    expected.iloc[[run_start, run_start + 1], 3:] = np.nan
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # the pyarrow engine does not support chunksize
@pytest.mark.parametrize("on_bad_lines", ["error", "skip"])
@pytest.mark.parametrize("chunksize", [1, 2])
def test_bad_line_first_in_chunk(all_parsers, on_bad_lines, chunksize):
    # GH#40587: reading with chunksize consumes each chunk's rows exactly, so
    # the next line lands in buffer slot 0 with no predecessor left to
    # compare against; a too-long line landing there escaped on_bad_lines
    # handling and was silently truncated to the table width and kept
    parser = all_parsers
    data = "1,2,3\n4,5,6\n7,8,9,10\n11,12,13\n"

    if on_bad_lines == "error":
        with pytest.raises(ParserError, match="Expected 3 fields"):
            with parser.read_csv(
                StringIO(data),
                header=None,
                on_bad_lines="error",
                chunksize=chunksize,
            ) as reader:
                pd.concat(reader)
    else:
        with parser.read_csv(
            StringIO(data),
            header=None,
            on_bad_lines="skip",
            chunksize=chunksize,
        ) as reader:
            result = pd.concat(reader)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [11, 12, 13]])
        tm.assert_frame_equal(result, expected)
