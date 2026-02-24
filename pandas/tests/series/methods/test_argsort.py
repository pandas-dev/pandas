import numpy as np
import pytest

from pandas import (
    Series,
    Timestamp,
    isna,
)
import pandas._testing as tm


class TestSeriesArgsort:
    @pytest.fixture
    def argsort_stability_series(self):
        return Series(np.random.default_rng(2).integers(0, 100, size=10000))

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

    def test_argsort_kind_stability(self, argsort_stability_series):
        ser = argsort_stability_series
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

    def test_argsort_stable(self, argsort_stability_series):
        ser = argsort_stability_series
        true_indexer = ser.argsort(stable=True)
        false_indexer = ser.argsort(stable=False)

        true_expected = np.argsort(ser.values, stable=True)
        false_expected = np.argsort(ser.values, stable=False)


        tm.assert_numpy_array_equal(false_indexer.values, false_expected)
        tm.assert_numpy_array_equal(true_indexer.values, true_expected)

    def test_argsort_numpy_stable(self, argsort_stability_series):
        ser = argsort_stability_series
        result = np.argsort(ser)
        expected = np.argsort(ser.values)
        tm.assert_numpy_array_equal(result.values, expected)

    @pytest.mark.parametrize(
        "kind, stable",
        [
            ("quicksort", True),
            ("quicksort", False),
            ("mergesort", True),
            ("mergesort", False),
        ],
    )
    def test_argsort_kind_and_stable(
        self,
        argsort_stability_series,
        kind,
        stable,
    ):
        ser = argsort_stability_series
        match = "`kind` and keyword parameters can't be provided at the same time."
        with tm.assert_produces_warning(
            FutureWarning, match=match, check_stacklevel=False
        ):
            np.argsort(ser, kind=kind, stable=stable)

    @pytest.mark.parametrize(
        "kind, stable",
        [
            ("quicksort", True),
            ("quicksort", False),
            ("mergesort", True),
            ("mergesort", False),
        ],
    )
    def test_argsort_numpy_kind_and_stable(
        self,
        argsort_stability_series,
        kind,
        stable,
    ):
        ser = argsort_stability_series
        match = "`kind` and `stable` can't be provided at the same time."
        with tm.assert_produces_warning(
            FutureWarning, match=match, check_stacklevel=False
        ):
            np.argsort(ser, kind=kind, stable=stable)

    def test_argsort_order(self, datetime_series):
        ser = datetime_series

        with tm.assert_produces_warning(False):
            ser.argsort(order="ts")
            ser.argsort(order=["ts"])

        with tm.assert_produces_warning(
            FutureWarning,
            check_stacklevel=False,
            match="`order` should match Series.name if specified",
        ):
            ser.argsort(order="something else")
            ser.argsort(order=["something else"])

    def test_argsort_numpy_order(self, datetime_series):
        ser = datetime_series

        with tm.assert_produces_warning(False):
            np.argsort(ser, order="ts")
            np.argsort(ser, order=["ts"])

        with tm.assert_produces_warning(
            FutureWarning,
            check_stacklevel=False,
            match="`order` should match Series.name if specified",
        ):
            np.argsort(ser, order="something else")
            np.argsort(ser, order=["something else"])

    def test_argsort_preserve_name(self, datetime_series):
        result = datetime_series.argsort()
        assert result.name == datetime_series.name
