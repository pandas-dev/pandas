import numpy as np

import pandas as pd
from pandas import Series
import pandas._testing as tm
from pandas.core.internals.blocks import IntBlock


class TestSeriesInternals:
    def test_constructor_no_pandas_array(self):
        ser = Series([1, 2, 3])
        result = Series(ser.array)
        tm.assert_series_equal(ser, result)
        assert isinstance(result._mgr.blocks[0], IntBlock)

    def test_astype_no_pandas_dtype(self):
        # https://github.com/pandas-dev/pandas/pull/24866
        ser = Series([1, 2], dtype="int64")
        # Don't have PandasDtype in the public API, so we use `.array.dtype`,
        # which is a PandasDtype.
        result = ser.astype(ser.array.dtype)
        tm.assert_series_equal(result, ser)

    def test_from_array(self):
        result = Series(pd.array(["1H", "2H"], dtype="timedelta64[ns]"))
        assert result._mgr.blocks[0].is_extension is False

        result = Series(pd.array(["2015"], dtype="datetime64[ns]"))
        assert result._mgr.blocks[0].is_extension is False

    def test_from_list_dtype(self):
        result = Series(["1H", "2H"], dtype="timedelta64[ns]")
        assert result._mgr.blocks[0].is_extension is False

        result = Series(["2015"], dtype="datetime64[ns]")
        assert result._mgr.blocks[0].is_extension is False


def test_hasnans_uncached_for_series():
    # GH#19700
    idx = pd.Index([0, 1])
    assert idx.hasnans is False
    assert "hasnans" in idx._cache
    ser = idx.to_series()
    assert ser.hasnans is False
    assert not hasattr(ser, "_cache")
    ser.iloc[-1] = np.nan
    assert ser.hasnans is True
    assert Series.hasnans.__doc__ == pd.Index.hasnans.__doc__
