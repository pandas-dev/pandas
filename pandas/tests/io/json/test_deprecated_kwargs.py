"""
Tests for the deprecated keyword arguments for `read_json`.
"""

from io import StringIO

import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm

from pandas.io.json import read_json


def test_good_kwargs():
    df = pd.DataFrame({"A": [2, 4, 6], "B": [3, 6, 9]}, index=[0, 1, 2])

    with tm.assert_produces_warning(None):
        data1 = StringIO(df.to_json(orient="split"))
        tm.assert_frame_equal(df, read_json(data1, orient="split"))
        data2 = StringIO(df.to_json(orient="columns"))
        tm.assert_frame_equal(df, read_json(data2, orient="columns"))
        data3 = StringIO(df.to_json(orient="index"))
        tm.assert_frame_equal(df, read_json(data3, orient="index"))


@pytest.mark.parametrize("val", [True, False, ["A"]])
def test_deprecated_convert_dates(val):
    # GH#59161
    df = pd.DataFrame({"A": [1, 2, 3]})
    json = StringIO(df.to_json())
    msg = "The 'convert_dates' keyword in read_json is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        read_json(json, convert_dates=val)


@pytest.mark.parametrize("val", [True, False])
def test_deprecated_keep_default_dates(val):
    # GH#59161
    df = pd.DataFrame({"A": [1, 2, 3]})
    json = StringIO(df.to_json())
    msg = "The 'keep_default_dates' keyword in read_json is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        read_json(json, keep_default_dates=val)
