import numpy as np

import pandas as pd
import pandas._testing as tm


class TestUnique:
    def test_unique_uint64(self):
        ser = pd.Series([1, 2, 2**63, 2**63], dtype=np.uint64)
        res = ser.unique()
        exp = np.array([1, 2, 2**63], dtype=np.uint64)
        tm.assert_numpy_array_equal(res, exp)

    def test_unique_data_ownership(self):
        # it works! GH#1807
        pd.Series(pd.Series(["a", "c", "b"]).unique()).sort_values()

    def test_unique(self):
        # GH#714 also, dtype=float
        ser = pd.Series([1.2345] * 100)
        ser[::2] = np.nan
        result = ser.unique()
        assert len(result) == 2

        # explicit f4 dtype
        ser = pd.Series([1.2345] * 100, dtype="f4")
        ser[::2] = np.nan
        result = ser.unique()
        assert len(result) == 2

    def test_unique_nan_object_dtype(self):
        # NAs in object arrays GH#714
        ser = pd.Series(["foo"] * 100, dtype="O")
        ser[::2] = np.nan
        result = ser.unique()
        assert len(result) == 2

    def test_unique_none(self):
        # decision about None
        ser = pd.Series([1, 2, 3, None, None, None], dtype=object)
        result = ser.unique()
        expected = np.array([1, 2, 3, None], dtype=object)
        tm.assert_numpy_array_equal(result, expected)

    def test_unique_categorical(self):
        # GH#18051
        cat = pd.Categorical([])
        ser = pd.Series(cat)
        result = ser.unique()
        tm.assert_categorical_equal(result, cat)

        cat = pd.Categorical([np.nan])
        ser = pd.Series(cat)
        result = ser.unique()
        tm.assert_categorical_equal(result, cat)

    def test_tz_unique(self):
        # GH 46128
        dti1 = pd.date_range("2016-01-01", periods=3)
        ii1 = pd.IntervalIndex.from_breaks(dti1)
        ser1 = pd.Series(ii1)
        uni1 = ser1.unique()
        tm.assert_interval_array_equal(ser1.array, uni1)

        dti2 = pd.date_range("2016-01-01", periods=3, tz="US/Eastern")
        ii2 = pd.IntervalIndex.from_breaks(dti2)
        ser2 = pd.Series(ii2)
        uni2 = ser2.unique()
        tm.assert_interval_array_equal(ser2.array, uni2)

        assert uni1.dtype != uni2.dtype
