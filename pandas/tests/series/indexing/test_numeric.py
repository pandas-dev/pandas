import numpy as np
import pytest

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


def test_getitem_negative_out_of_bounds():
    s = Series(tm.rands_array(5, 10), index=tm.rands_array(10, 10))

    msg = "index -11 is out of bounds for axis 0 with size 10"
    with pytest.raises(IndexError, match=msg):
        s[-11]
    with pytest.raises(IndexError, match=msg):
        s[-11] = "foo"


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


def test_slice_float_get_set(datetime_series):
    msg = (
        "cannot do slice indexing on DatetimeIndex with these indexers "
        r"\[{key}\] of type float"
    )
    with pytest.raises(TypeError, match=msg.format(key=r"4\.0")):
        datetime_series[4.0:10.0]

    with pytest.raises(TypeError, match=msg.format(key=r"4\.0")):
        datetime_series[4.0:10.0] = 0

    with pytest.raises(TypeError, match=msg.format(key=r"4\.5")):
        datetime_series[4.5:10.0]
    with pytest.raises(TypeError, match=msg.format(key=r"4\.5")):
        datetime_series[4.5:10.0] = 0
