from datetime import datetime

import numpy as np

from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    MultiIndex,
    Period,
    Series,
    period_range,
    to_datetime,
)


def test_multiindex_period_datetime():
    # GH4861, using datetime in period of multiindex raises exception

    idx1 = Index(["a", "a", "a", "b", "b"])
    idx2 = period_range("2012-01", periods=len(idx1), freq="M")
    s = Series(np.random.randn(len(idx1)), [idx1, idx2])

    # try Period as index
    expected = s.iloc[0]
    result = s.loc["a", Period("2012-01")]
    assert result == expected

    # try datetime as index
    result = s.loc["a", datetime(2012, 1, 1)]
    assert result == expected


def test_multiindex_datetime_columns():
    # GH35015, using datetime as column indices raises exception

    mi = MultiIndex.from_tuples(
        [(to_datetime("02/29/2020"), to_datetime("03/01/2020"))], names=["a", "b"]
    )

    df = DataFrame([], columns=mi)
    assert list(df.columns.names) == ["a", "b"]

    assert all(isinstance(mi.get_level_values(i), DatetimeIndex) for i in range(2))
