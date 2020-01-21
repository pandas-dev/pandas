import numpy as np
import pytest

from pandas._libs import groupby, lib, reduction as libreduction

from pandas.core.dtypes.common import ensure_int64

from pandas import Index, Series, isna
import pandas._testing as tm


def test_series_grouper():
    obj = Series(np.random.randn(10))
    dummy = obj[:0]

    labels = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1, 1], dtype=np.int64)

    grouper = libreduction.SeriesGrouper(obj, np.mean, labels, 2, dummy)
    result, counts = grouper.get_result()

    expected = np.array([obj[3:6].mean(), obj[6:].mean()])
    tm.assert_almost_equal(result, expected)

    exp_counts = np.array([3, 4], dtype=np.int64)
    tm.assert_almost_equal(counts, exp_counts)


def test_series_grouper_requires_nonempty_raises():
    # GH#29500
    obj = Series(np.random.randn(10))
    dummy = obj[:0]
    labels = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1, 1], dtype=np.int64)

    with pytest.raises(ValueError, match="SeriesGrouper requires non-empty `series`"):
        libreduction.SeriesGrouper(dummy, np.mean, labels, 2, dummy)


def test_series_bin_grouper():
    obj = Series(np.random.randn(10))
    dummy = obj[:0]

    bins = np.array([3, 6])

    grouper = libreduction.SeriesBinGrouper(obj, np.mean, bins, dummy)
    result, counts = grouper.get_result()

    expected = np.array([obj[:3].mean(), obj[3:6].mean(), obj[6:].mean()])
    tm.assert_almost_equal(result, expected)

    exp_counts = np.array([3, 3, 4], dtype=np.int64)
    tm.assert_almost_equal(counts, exp_counts)


@pytest.mark.parametrize(
    "binner,closed,expected",
    [
        (
            np.array([0, 3, 6, 9], dtype=np.int64),
            "left",
            np.array([2, 5, 6], dtype=np.int64),
        ),
        (
            np.array([0, 3, 6, 9], dtype=np.int64),
            "right",
            np.array([3, 6, 6], dtype=np.int64),
        ),
        (np.array([0, 3, 6], dtype=np.int64), "left", np.array([2, 5], dtype=np.int64)),
        (
            np.array([0, 3, 6], dtype=np.int64),
            "right",
            np.array([3, 6], dtype=np.int64),
        ),
    ],
)
def test_generate_bins(binner, closed, expected):
    values = np.array([1, 2, 3, 4, 5, 6], dtype=np.int64)
    result = lib.generate_bins_dt64(values, binner, closed=closed)
    tm.assert_numpy_array_equal(result, expected)


def test_group_ohlc():
    def _check(dtype):
        obj = np.array(np.random.randn(20), dtype=dtype)

        bins = np.array([6, 12, 20])
        out = np.zeros((3, 4), dtype)
        counts = np.zeros(len(out), dtype=np.int64)
        labels = ensure_int64(np.repeat(np.arange(3), np.diff(np.r_[0, bins])))

        func = getattr(groupby, f"group_ohlc_{dtype}")
        func(out, counts, obj[:, None], labels)

        def _ohlc(group):
            if isna(group).all():
                return np.repeat(np.nan, 4)
            return [group[0], group.max(), group.min(), group[-1]]

        expected = np.array([_ohlc(obj[:6]), _ohlc(obj[6:12]), _ohlc(obj[12:])])

        tm.assert_almost_equal(out, expected)
        tm.assert_numpy_array_equal(counts, np.array([6, 6, 8], dtype=np.int64))

        obj[:6] = np.nan
        func(out, counts, obj[:, None], labels)
        expected[0] = np.nan
        tm.assert_almost_equal(out, expected)

    _check("float32")
    _check("float64")


class TestMoments:
    pass


class TestReducer:
    def test_int_index(self):
        arr = np.random.randn(100, 4)

        msg = "Must pass either dummy and labels, or neither"
        # we must pass either both labels and dummy, or neither
        with pytest.raises(ValueError, match=msg):
            libreduction.compute_reduction(arr, np.sum, labels=Index(np.arange(4)))

        with pytest.raises(ValueError, match=msg):
            libreduction.compute_reduction(
                arr, np.sum, axis=1, labels=Index(np.arange(100))
            )

        dummy = Series(0.0, index=np.arange(100))
        result = libreduction.compute_reduction(
            arr, np.sum, dummy=dummy, labels=Index(np.arange(4))
        )
        expected = arr.sum(0)
        tm.assert_almost_equal(result, expected)

        dummy = Series(0.0, index=np.arange(4))
        result = libreduction.compute_reduction(
            arr, np.sum, axis=1, dummy=dummy, labels=Index(np.arange(100))
        )
        expected = arr.sum(1)
        tm.assert_almost_equal(result, expected)

        result = libreduction.compute_reduction(
            arr, np.sum, axis=1, dummy=dummy, labels=Index(np.arange(100))
        )
        tm.assert_almost_equal(result, expected)
