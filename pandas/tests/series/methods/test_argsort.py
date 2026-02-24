import numpy as np
import pytest

from pandas import (
    Series,
    Timestamp,
    isna,
)
import pandas._testing as tm


class TestSeriesArgsort:
    def test_argsort_axis(self):
        # GH#54257
        ser = Series(range(3))

        msg = "No axis named 2 for object type Series"
        with pytest.raises(ValueError, match=msg):
            ser.argsort(axis=2)

    def test_argsort_numpy(self, datetime_series):
        ser = datetime_series
        res = np.argsort(ser).values
        expected = np.argsort(np.array(ser))
        tm.assert_numpy_array_equal(res, expected)

    def test_argsort_numpy_missing(self):
        data = [0.1, np.nan, 0.2, np.nan, 0.3]
        ser = Series(data)
        result = np.argsort(ser)
        expected = np.argsort(np.array(data))

        tm.assert_numpy_array_equal(result.values, expected)

    def test_argsort(self, datetime_series):
        argsorted = datetime_series.argsort()
        assert issubclass(argsorted.dtype.type, np.integer)

    def test_argsort_dt64(self, unit):
        # GH#2967 (introduced bug in 0.11-dev I think)
        ser = Series(
            [Timestamp(f"201301{i:02d}") for i in range(1, 6)], dtype=f"M8[{unit}]"
        )
        assert ser.dtype == f"datetime64[{unit}]"
        shifted = ser.shift(-1)
        assert shifted.dtype == f"datetime64[{unit}]"
        assert isna(shifted[4])

        result = ser.argsort()
        expected = Series(range(5), dtype=np.intp)
        tm.assert_series_equal(result, expected)

        result = shifted.argsort()
        expected = Series([*list(range(4)), 4], dtype=np.intp)
        tm.assert_series_equal(result, expected)

    def test_argsort_stable(self):
        ser = Series(np.random.default_rng(2).integers(0, 100, size=10000))
        mindexer = ser.argsort(kind="mergesort")
        qindexer = ser.argsort()

        mexpected = np.argsort(ser.values, kind="mergesort")
        qexpected = np.argsort(ser.values, kind="quicksort")

        tm.assert_series_equal(mindexer.astype(np.intp), Series(mexpected))
        tm.assert_series_equal(qindexer.astype(np.intp), Series(qexpected))
        msg = (
            r"ndarray Expected type <class 'numpy\.ndarray'>, "
            r"found <class 'pandas\.Series'> instead"
        )
        with pytest.raises(AssertionError, match=msg):
            tm.assert_numpy_array_equal(qindexer, mindexer)

    def test_argsort_preserve_name(self, datetime_series):
        result = datetime_series.argsort()
        assert result.name == datetime_series.name

    def test_argsort_stable_behavior():
        ser = Series([3, 1, 2])
        arr = np.array([3, 1, 2])

        # 1. stable=True matches NumPy
        result_stable_true = ser.argsort(stable=True).to_numpy()
        expected_stable_true = np.argsort(arr, stable=True)
        tm.assert_numpy_array_equal(result_stable_true, expected_stable_true)

        # 2. stable=False matches NumPy
        result_stable_false = ser.argsort(stable=False).to_numpy()
        expected_stable_false = np.argsort(arr, stable=False)
        tm.assert_numpy_array_equal(result_stable_false, expected_stable_false)

        # 3. Conflicting kind and stable should raise
        with pytest.raises(ValueError, match="kind.*keyword parameters"):
            ser.argsort(kind="heapsort", stable=True)
