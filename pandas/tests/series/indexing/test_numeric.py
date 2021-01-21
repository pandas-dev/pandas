import numpy as np

from pandas import DataFrame, Index, Series
import pandas._testing as tm


def test_slice_float64():
    values = np.arange(10.0, 50.0, 2)
    index = Index(values)

    start, end = values[[5, 15]]

    s = Series(np.random.randn(20), index=index)

    result = s[start:end]
    expected = s.iloc[5:16]
    tm.assert_series_equal(result, expected)

    result = s.loc[start:end]
    tm.assert_series_equal(result, expected)

    df = DataFrame(np.random.randn(20, 3), index=index)

    result = df[start:end]
    expected = df.iloc[5:16]
    tm.assert_frame_equal(result, expected)

    result = df.loc[start:end]
    tm.assert_frame_equal(result, expected)


def test_getitem_setitem_slice_bug():
    s = Series(range(10), index=list(range(10)))
    result = s[-12:]
    tm.assert_series_equal(result, s)

    result = s[-7:]
    tm.assert_series_equal(result, s[3:])

    result = s[:-12]
    tm.assert_series_equal(result, s[:0])

    s = Series(range(10), index=list(range(10)))
    s[-12:] = 0
    assert (s == 0).all()

    s[:-12] = 5
    assert (s == 0).all()


def test_getitem_setitem_slice_integers():
    s = Series(np.random.randn(8), index=[2, 4, 6, 8, 10, 12, 14, 16])

    result = s[:4]
    expected = s.reindex([2, 4, 6, 8])
    tm.assert_series_equal(result, expected)

    s[:4] = 0
    assert (s[:4] == 0).all()
    assert not (s[4:] == 0).any()
