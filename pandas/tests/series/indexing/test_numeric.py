import numpy as np

from pandas import DataFrame, Index, Series
import pandas._testing as tm


def test_slice_float64(frame_or_series):
    values = np.arange(10.0, 50.0, 2)
    index = Index(values)

    start, end = values[[5, 15]]

    data = np.random.randn(20, 3)
    if frame_or_series is not DataFrame:
        data = data[:, 0]

    obj = frame_or_series(data, index=index)

    result = obj[start:end]
    expected = obj.iloc[5:16]
    tm.assert_equal(result, expected)

    result = obj.loc[start:end]
    tm.assert_equal(result, expected)


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
