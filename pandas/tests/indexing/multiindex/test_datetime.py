from datetime import datetime

import numpy as np

import pandas as pd
import pandas._testing as tm


def test_multiindex_period_datetime():
    # GH4861, using datetime in period of multiindex raises exception

    idx1 = pd.Index(["a", "a", "a", "b", "b"])
    idx2 = pd.period_range("2012-01", periods=len(idx1), freq="M")
    s = pd.Series(np.random.default_rng(2).standard_normal(len(idx1)), [idx1, idx2])

    # try Period as index
    expected = s.iloc[0]
    result = s.loc["a", pd.Period("2012-01")]
    assert result == expected

    # try datetime as index
    result = s.loc["a", datetime(2012, 1, 1)]
    assert result == expected


def test_multiindex_datetime_columns():
    # GH35015, using datetime as column indices raises exception

    mi = pd.MultiIndex.from_tuples(
        [(pd.to_datetime("02/29/2020"), pd.to_datetime("03/01/2020"))], names=["a", "b"]
    )

    df = pd.DataFrame([], columns=mi)

    expected_df = pd.DataFrame(
        [],
        columns=pd.MultiIndex.from_arrays(
            [[pd.to_datetime("02/29/2020")], [pd.to_datetime("03/01/2020")]],
            names=["a", "b"],
        ),
    )

    tm.assert_frame_equal(df, expected_df)
