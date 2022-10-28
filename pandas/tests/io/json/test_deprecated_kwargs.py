"""
Tests for the deprecated keyword arguments for `read_json`.
"""

import pandas as pd
import pandas._testing as tm

from pandas.io.json import read_json


def test_good_kwargs():
    df = pd.DataFrame({"A": [2, 4, 6], "B": [3, 6, 9]}, index=[0, 1, 2])
    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(df, read_json(df.to_json(orient="split"), orient="split"))
        tm.assert_frame_equal(
            df, read_json(df.to_json(orient="columns"), orient="columns")
        )
        tm.assert_frame_equal(df, read_json(df.to_json(orient="index"), orient="index"))
