from datetime import datetime

import numpy as np

from pandas import Index, Period, Series, period_range


def test_multiindex_period_datetime():
    # GH4861, using datetime in period of multiindex raises exception

    idx1 = Index(['a', 'a', 'a', 'b', 'b'])
    idx2 = period_range('2012-01', periods=len(idx1), freq='M')
    s = Series(np.random.randn(len(idx1)), [idx1, idx2])

    # try Period as index
    expected = s.iloc[0]
    result = s.loc['a', Period('2012-01')]
    assert result == expected

    # try datetime as index
    result = s.loc['a', datetime(2012, 1, 1)]
    assert result == expected
