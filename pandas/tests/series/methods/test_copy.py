import numpy as np
import pytest

from pandas import (
    Series,
    Timestamp,
)
import pandas._testing as tm


class TestCopy:
    @pytest.mark.parametrize("deep", ["default", None, False, True])
    def test_copy(self, deep):
        ser = Series(np.arange(10), dtype="float64")

        # default deep is True
        if deep == "default":
            ser2 = ser.copy()
        else:
            ser2 = ser.copy(deep=deep)

        # INFO(CoW) a shallow copy doesn't yet copy the data
        # but parent will not be modified (CoW)
        if deep is None or deep is False:
            assert np.may_share_memory(ser.values, ser2.values)
        else:
            assert not np.may_share_memory(ser.values, ser2.values)

        ser2[::2] = np.nan

        # Did not modify original Series
        assert np.isnan(ser2[0])
        assert not np.isnan(ser[0])

    @pytest.mark.filterwarnings("ignore:Setting a value on a view:FutureWarning")
    @pytest.mark.parametrize("deep", ["default", None, False, True])
    def test_copy_tzaware(self, deep):
        # GH#11794
        # copy of tz-aware
        expected = Series([Timestamp("2012/01/01", tz="UTC")])
        expected2 = Series([Timestamp("1999/01/01", tz="UTC")])

        ser = Series([Timestamp("2012/01/01", tz="UTC")])

        if deep == "default":
            ser2 = ser.copy()
        else:
            ser2 = ser.copy(deep=deep)

        # INFO(CoW) a shallow copy doesn't yet copy the data
        # but parent will not be modified (CoW)
        if deep is None or deep is False:
            assert np.may_share_memory(ser.values, ser2.values)
        else:
            assert not np.may_share_memory(ser.values, ser2.values)

        ser2[0] = Timestamp("1999/01/01", tz="UTC")

        # Did not modify original Series
        tm.assert_series_equal(ser2, expected2)
        tm.assert_series_equal(ser, expected)

    def test_copy_name(self, datetime_series):
        result = datetime_series.copy()
        assert result.name == datetime_series.name

    def test_copy_index_name_checking(self, datetime_series):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy

        datetime_series.index.name = None
        assert datetime_series.index.name is None
        assert datetime_series is datetime_series

        cp = datetime_series.copy()
        cp.index.name = "foo"
        assert datetime_series.index.name is None
