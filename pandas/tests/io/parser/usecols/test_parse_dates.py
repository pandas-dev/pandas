"""
Tests the usecols functionality during parsing
for all of the parsers defined in parsers.py
"""

from io import StringIO

import pytest

from pandas import (
    DataFrame,
    Index,
    Timestamp,
)
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)
xfail_pyarrow = pytest.mark.usefixtures("pyarrow_xfail")
skip_pyarrow = pytest.mark.usefixtures("pyarrow_skip")

_msg_pyarrow_requires_names = (
    "The pyarrow engine does not allow 'usecols' to be integer column "
    "positions. Pass a list of string column names instead."
)


@skip_pyarrow  # pyarrow.lib.ArrowKeyError: Column 'fdate' in include_columns
def test_usecols_with_parse_dates2(all_parsers):
    # see gh-13604
    parser = all_parsers
    data = """2008-02-07 09:40,1032.43
2008-02-07 09:50,1042.54
2008-02-07 10:00,1051.65"""

    names = ["date", "values"]
    usecols = names[:]
    parse_dates = [0]

    index = Index(
        [
            Timestamp("2008-02-07 09:40"),
            Timestamp("2008-02-07 09:50"),
            Timestamp("2008-02-07 10:00"),
        ],
        name="date",
    )
    cols = {"values": [1032.43, 1042.54, 1051.65]}
    expected = DataFrame(cols, index=index)

    result = parser.read_csv(
        StringIO(data),
        parse_dates=parse_dates,
        index_col=0,
        usecols=usecols,
        header=None,
        names=names,
    )
    tm.assert_frame_equal(result, expected)


def test_usecols_with_parse_dates3(all_parsers):
    # see gh-14792
    parser = all_parsers
    data = """a,b,c,d,e,f,g,h,i,j
2016/09/21,1,1,2,3,4,5,6,7,8"""

    usecols = list("abcdefghij")
    parse_dates = [0]

    cols = {
        "a": Timestamp("2016-09-21"),
        "b": [1],
        "c": [1],
        "d": [2],
        "e": [3],
        "f": [4],
        "g": [5],
        "h": [6],
        "i": [7],
        "j": [8],
    }
    expected = DataFrame(cols, columns=usecols)

    result = parser.read_csv(StringIO(data), usecols=usecols, parse_dates=parse_dates)
    tm.assert_frame_equal(result, expected)
