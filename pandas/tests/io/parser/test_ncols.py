"""
Tests for the ``ncols`` argument of read_csv and read_table (GH#72210).
"""

from io import StringIO

import numpy as np
import pytest

from pandas.errors import (
    Pandas4Warning,
    ParserError,
)

import pandas as pd
import pandas._testing as tm

skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")


@skip_pyarrow  # ncols is not supported by the pyarrow engine
def test_ncols_pads_short_rows(all_parsers):
    parser = all_parsers
    result = parser.read_csv(StringIO("1,2\n3,4,5"), header=None, ncols=3)
    expected = pd.DataFrame({0: [1, 3], 1: [2, 4], 2: [np.nan, 5.0]})
    tm.assert_frame_equal(result, expected)

    result = parser.read_table(StringIO("1\t2\n3\t4\t5"), header=None, ncols=3)
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # ncols is not supported by the pyarrow engine
def test_ncols_keys_are_positions(all_parsers):
    parser = all_parsers
    result = parser.read_csv(
        StringIO("01,2\n03,4,5"),
        header=None,
        ncols=3,
        dtype={0: str},
        index_col=0,
    )
    expected = pd.DataFrame(
        {1: [2, 4], 2: [np.nan, 5.0]}, index=pd.Index(["01", "03"], name=0)
    )
    tm.assert_frame_equal(result, expected)


@skip_pyarrow  # ncols is not supported by the pyarrow engine
@pytest.mark.parametrize("data", ["1,2,3\n4,5", "1,2\n4,5,6"])
def test_ncols_long_row_raises(all_parsers, data):
    # unlike names=range(2), a long first row does not become an index
    parser = all_parsers
    with pytest.raises(ParserError, match="Expected 2 fields"):
        parser.read_csv(StringIO(data), header=None, ncols=2)


@skip_pyarrow  # ncols is not supported by the pyarrow engine
def test_ncols_chunksize(all_parsers):
    parser = all_parsers
    with parser.read_csv(
        StringIO("1,2\n3,4,5"), header=None, ncols=3, chunksize=1
    ) as reader:
        result = pd.concat(reader)
    expected = pd.DataFrame({0: [1, 3], 1: [2, 4], 2: [np.nan, 5.0]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs, msg",
    [
        ({}, "'ncols' requires header=None"),
        ({"header": 0}, "'ncols' requires header=None"),
        ({"header": None, "ncols": 0}, "'ncols' must be an integer >=1"),
        ({"header": None, "ncols": "3"}, "'ncols' must be an integer >=1"),
    ],
)
def test_ncols_invalid(kwargs, msg):
    kwargs = {"ncols": 3} | kwargs
    with pytest.raises(ValueError, match=msg):
        pd.read_csv(StringIO("1,2\n3,4,5"), **kwargs)


def test_ncols_with_names_raises():
    msg = "Cannot specify both 'names' and 'ncols'"
    with tm.assert_produces_warning(Pandas4Warning, match="'names' argument"):
        with pytest.raises(ValueError, match=msg):
            pd.read_csv(
                StringIO("1,2\n3,4,5"), header=None, ncols=3, names=["a", "b", "c"]
            )


def test_ncols_pyarrow_unsupported(pyarrow_parser_only):
    parser = pyarrow_parser_only
    msg = "The 'ncols' option is not supported with the 'pyarrow' engine"
    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO("1,2\n3,4,5"), header=None, ncols=3)


def test_read_fwf_rejects_ncols():
    # read_fwf already pads short rows
    with pytest.raises(TypeError, match="unexpected keyword argument 'ncols'"):
        pd.read_fwf(StringIO("1 2\n3 4 5"), header=None, ncols=3)
