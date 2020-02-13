import numpy as np
import pytest

from pandas import DataFrame, Index, Series
import pandas._testing as tm


def test_delitem():
    # GH 5542
    # should delete the item inplace
    s = Series(range(5))
    del s[0]

    expected = Series(range(1, 5), index=range(1, 5))
    tm.assert_series_equal(s, expected)

    del s[1]
    expected = Series(range(2, 5), index=range(2, 5))
    tm.assert_series_equal(s, expected)

    # empty
    s = Series(dtype=object)

    with pytest.raises(KeyError, match=r"^0$"):
        del s[0]

    # only 1 left, del, add, del
    s = Series(1)
    del s[0]
    tm.assert_series_equal(s, Series(dtype="int64", index=Index([], dtype="int64")))
    s[0] = 1
    tm.assert_series_equal(s, Series(1))
    del s[0]
    tm.assert_series_equal(s, Series(dtype="int64", index=Index([], dtype="int64")))

    # Index(dtype=object)
    s = Series(1, index=["a"])
    del s["a"]
    tm.assert_series_equal(s, Series(dtype="int64", index=Index([], dtype="object")))
    s["a"] = 1
    tm.assert_series_equal(s, Series(1, index=["a"]))
    del s["a"]
    tm.assert_series_equal(s, Series(dtype="int64", index=Index([], dtype="object")))


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


def test_getitem_regression():
    s = Series(range(5), index=list(range(5)))
    result = s[list(range(5))]
    tm.assert_series_equal(result, s)


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


def test_setitem_float_labels():
    # note labels are floats
    s = Series(["a", "b", "c"], index=[0, 0.5, 1])
    tmp = s.copy()

    s.loc[1] = "zoo"
    tmp.iloc[2] = "zoo"

    tm.assert_series_equal(s, tmp)


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


def test_slice_floats2():
    s = Series(np.random.rand(10), index=np.arange(10, 20, dtype=float))

    assert len(s.loc[12.0:]) == 8
    assert len(s.loc[12.5:]) == 7

    i = np.arange(10, 20, dtype=float)
    i[2] = 12.2
    s.index = i
    assert len(s.loc[12.0:]) == 8
    assert len(s.loc[12.5:]) == 7


def test_int_indexing():
    s = Series(np.random.randn(6), index=[0, 0, 1, 1, 2, 2])

    with pytest.raises(KeyError, match=r"^5$"):
        s[5]

    with pytest.raises(KeyError, match=r"^'c'$"):
        s["c"]

    # not monotonic
    s = Series(np.random.randn(6), index=[2, 2, 0, 0, 1, 1])

    with pytest.raises(KeyError, match=r"^5$"):
        s[5]

    with pytest.raises(KeyError, match=r"^'c'$"):
        s["c"]


def test_getitem_int64(datetime_series):
    idx = np.int64(5)
    assert datetime_series[idx] == datetime_series[5]
