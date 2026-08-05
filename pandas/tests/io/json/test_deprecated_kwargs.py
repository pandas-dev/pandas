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


@pytest.mark.parametrize("kwarg", ["convert_dates", "keep_default_dates"])
@pytest.mark.parametrize("value", [True, False])
def test_convert_dates_kwargs_deprecated(kwarg, value):
    # GH#59161 convert_dates and keep_default_dates are deprecated in favor of
    #  controlling type conversion through dtype
    df = pd.DataFrame({"A": [1, 2, 3]})
    data = StringIO(df.to_json())
    msg = f"The '{kwarg}' keyword in read_json is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = read_json(data, **{kwarg: value})
    # behavior is unchanged: the frame still round-trips
    tm.assert_frame_equal(result, df)


def test_convert_dates_list_deprecated():
    # GH#59161 passing a list of columns to convert_dates also warns
    df = pd.DataFrame({"A": [1, 2, 3]})
    data = StringIO(df.to_json())
    msg = "The 'convert_dates' keyword in read_json is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        read_json(data, convert_dates=["A"])
