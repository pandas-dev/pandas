from pandas import (
    DataFrame,
    Index,
    Series,
)
import pandas._testing as tm


def test_index_constructor():
    ser = Series([1, 2, 3])
    idx = Index(ser)
    # assert not np.shares_memory(idx.values, get_array(ser))
    ser.iloc[0] = 0
    tm.assert_index_equal(idx, Index([1, 2, 3]))


def test_series_constructor():
    ser = Series([1, 2, 3])
    ser2 = Series(ser, index=ser)
    ser2.iloc[0] = 0
    tm.assert_index_equal(ser2.index, Index([1, 2, 3]))


def test_dataframe_constructor():
    ser = Series([1, 2, 3])
    df = DataFrame({"a": ser}, index=ser)
    ser.iloc[0] = 0
    tm.assert_index_equal(df.index, Index([1, 2, 3]))
    df.iloc[0, 0] = 10
    tm.assert_index_equal(df.index, Index([1, 2, 3]))


def test_dataframe_set_index_inplace():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df.set_index("a", drop=False, inplace=True)
    df.iloc[0, 0] = 0
    tm.assert_index_equal(df.index, Index([1, 2, 3], name="a"))


def test_index_attribute_assignment():
    ser = Series([1, 2, 3])
    ser.index = ser
    ser.iloc[0] = 10
    tm.assert_index_equal(ser.index, Index([1, 2, 3]))
