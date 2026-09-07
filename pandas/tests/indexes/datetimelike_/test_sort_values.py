import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def check_freq_ascending(ordered, orig, ascending):
    """
    Check the expected freq on a PeriodIndex/DatetimeIndex/TimedeltaIndex
    when the original index is generated (or generate-able) with
    period_range/date_range/timedelta_range.
    """
    if isinstance(ordered, pd.PeriodIndex):
        assert ordered.freq == orig.freq
    elif isinstance(ordered, (pd.DatetimeIndex, pd.TimedeltaIndex)):
        if ascending:
            assert ordered.freq.n == orig.freq.n
        else:
            assert ordered.freq.n == -1 * orig.freq.n


def check_freq_nonmonotonic(ordered, orig):
    """
    Check the expected freq on a PeriodIndex/DatetimeIndex/TimedeltaIndex
    when the original index is _not_ generated (or generate-able) with
    period_range/date_range//timedelta_range.
    """
    if isinstance(ordered, pd.PeriodIndex):
        assert ordered.freq == orig.freq
    elif isinstance(ordered, (pd.DatetimeIndex, pd.TimedeltaIndex)):
        assert ordered.freq is None


class TestSortValues:
    @pytest.fixture(params=[pd.DatetimeIndex, pd.TimedeltaIndex, pd.PeriodIndex])
    def non_monotonic_idx(self, request):
        if request.param is pd.DatetimeIndex:
            return pd.DatetimeIndex(["2000-01-04", "2000-01-01", "2000-01-02"])
        elif request.param is pd.PeriodIndex:
            dti = pd.DatetimeIndex(["2000-01-04", "2000-01-01", "2000-01-02"])
            return dti.to_period("D")
        else:
            return pd.TimedeltaIndex(
                ["1 day 00:00:05", "1 day 00:00:01", "1 day 00:00:02"]
            )

    def test_argmin_argmax(self, non_monotonic_idx):
        assert non_monotonic_idx.argmin() == 1
        assert non_monotonic_idx.argmax() == 0

    def test_sort_values(self, non_monotonic_idx):
        idx = non_monotonic_idx
        ordered = idx.sort_values()
        assert ordered.is_monotonic_increasing
        ordered = idx.sort_values(ascending=False)
        assert ordered[::-1].is_monotonic_increasing

        ordered, dexer = idx.sort_values(return_indexer=True)
        assert ordered.is_monotonic_increasing
        tm.assert_numpy_array_equal(dexer, np.array([1, 2, 0], dtype=np.intp))

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        assert ordered[::-1].is_monotonic_increasing
        tm.assert_numpy_array_equal(dexer, np.array([0, 2, 1], dtype=np.intp))

    def check_sort_values_with_freq(self, idx):
        ordered = idx.sort_values()
        tm.assert_index_equal(ordered, idx)
        check_freq_ascending(ordered, idx, True)

        ordered = idx.sort_values(ascending=False)
        expected = idx[::-1]
        tm.assert_index_equal(ordered, expected)
        check_freq_ascending(ordered, idx, False)

        ordered, indexer = idx.sort_values(return_indexer=True)
        tm.assert_index_equal(ordered, idx)
        tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2], dtype=np.intp))
        check_freq_ascending(ordered, idx, True)

        ordered, indexer = idx.sort_values(return_indexer=True, ascending=False)
        expected = idx[::-1]
        tm.assert_index_equal(ordered, expected)
        tm.assert_numpy_array_equal(indexer, np.array([2, 1, 0], dtype=np.intp))
        check_freq_ascending(ordered, idx, False)

    @pytest.mark.parametrize("freq", ["D", "h"])
    def test_sort_values_with_freq_timedeltaindex(self, freq):
        # GH#10295
        idx = pd.timedelta_range(start=f"1{freq}", periods=3, freq=freq).rename("idx")

        self.check_sort_values_with_freq(idx)

    @pytest.mark.parametrize(
        "idx",
        [
            pd.DatetimeIndex(
                ["2011-01-01", "2011-01-02", "2011-01-03"], freq="D", name="idx"
            ),
            pd.DatetimeIndex(
                ["2011-01-01 09:00", "2011-01-01 10:00", "2011-01-01 11:00"],
                freq="h",
                name="tzidx",
                tz="Asia/Tokyo",
            ),
        ],
    )
    def test_sort_values_with_freq_datetimeindex(self, idx):
        self.check_sort_values_with_freq(idx)

    @pytest.mark.parametrize("freq", ["D", "2D", "4D"])
    def test_sort_values_with_freq_periodindex(self, freq):
        # here with_freq refers to being period_range-like
        idx = pd.PeriodIndex(
            ["2011-01-01", "2011-01-02", "2011-01-03"], freq=freq, name="idx"
        )
        self.check_sort_values_with_freq(idx)

    @pytest.mark.parametrize(
        "idx",
        [
            pd.PeriodIndex(["2011", "2012", "2013"], name="pidx", freq="Y"),
            pd.Index([2011, 2012, 2013], name="idx"),  # for compatibility check
        ],
    )
    def test_sort_values_with_freq_periodindex2(self, idx):
        # here with_freq indicates this is period_range-like
        self.check_sort_values_with_freq(idx)

    def check_sort_values_without_freq(self, idx, expected):
        ordered = idx.sort_values(na_position="first")
        tm.assert_index_equal(ordered, expected)
        check_freq_nonmonotonic(ordered, idx)

        if not idx.isna().any():
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, expected)
            check_freq_nonmonotonic(ordered, idx)

        ordered = idx.sort_values(ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])
        check_freq_nonmonotonic(ordered, idx)

        ordered, indexer = idx.sort_values(return_indexer=True, na_position="first")
        tm.assert_index_equal(ordered, expected)

        exp = np.array([0, 4, 3, 1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, exp)
        check_freq_nonmonotonic(ordered, idx)

        if not idx.isna().any():
            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, expected)

            exp = np.array([0, 4, 3, 1, 2], dtype=np.intp)
            tm.assert_numpy_array_equal(indexer, exp)
            check_freq_nonmonotonic(ordered, idx)

        ordered, indexer = idx.sort_values(return_indexer=True, ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])

        exp = np.array([2, 1, 3, 0, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, exp)
        check_freq_nonmonotonic(ordered, idx)

    def test_sort_values_without_freq_timedeltaindex(self):
        # GH#10295

        idx = pd.TimedeltaIndex(
            ["1 hour", "3 hour", "5 hour", "2 hour ", "1 hour"], name="idx1"
        )
        expected = pd.TimedeltaIndex(
            ["1 hour", "1 hour", "2 hour", "3 hour", "5 hour"], name="idx1"
        )
        self.check_sort_values_without_freq(idx, expected)

    @pytest.mark.parametrize(
        "index_dates,expected_dates",
        [
            (
                ["2011-01-01", "2011-01-03", "2011-01-05", "2011-01-02", "2011-01-01"],
                ["2011-01-01", "2011-01-01", "2011-01-02", "2011-01-03", "2011-01-05"],
            ),
            (
                [pd.NaT, "2011-01-03", "2011-01-05", "2011-01-02", pd.NaT],
                [pd.NaT, pd.NaT, "2011-01-02", "2011-01-03", "2011-01-05"],
            ),
        ],
    )
    def test_sort_values_without_freq_datetimeindex(
        self, index_dates, expected_dates, tz_naive_fixture
    ):
        tz = tz_naive_fixture

        # without freq
        idx = pd.DatetimeIndex(index_dates, tz=tz, name="idx")
        expected = pd.DatetimeIndex(expected_dates, tz=tz, name="idx")

        self.check_sort_values_without_freq(idx, expected)

    @pytest.mark.parametrize(
        "idx,expected",
        [
            (
                pd.PeriodIndex(
                    [
                        "2011-01-01",
                        "2011-01-03",
                        "2011-01-05",
                        "2011-01-02",
                        "2011-01-01",
                    ],
                    freq="D",
                    name="idx1",
                ),
                pd.PeriodIndex(
                    [
                        "2011-01-01",
                        "2011-01-01",
                        "2011-01-02",
                        "2011-01-03",
                        "2011-01-05",
                    ],
                    freq="D",
                    name="idx1",
                ),
            ),
            (
                pd.PeriodIndex(
                    [
                        "2011-01-01",
                        "2011-01-03",
                        "2011-01-05",
                        "2011-01-02",
                        "2011-01-01",
                    ],
                    freq="D",
                    name="idx2",
                ),
                pd.PeriodIndex(
                    [
                        "2011-01-01",
                        "2011-01-01",
                        "2011-01-02",
                        "2011-01-03",
                        "2011-01-05",
                    ],
                    freq="D",
                    name="idx2",
                ),
            ),
            (
                pd.PeriodIndex(
                    [pd.NaT, "2011-01-03", "2011-01-05", "2011-01-02", pd.NaT],
                    freq="D",
                    name="idx3",
                ),
                pd.PeriodIndex(
                    [pd.NaT, pd.NaT, "2011-01-02", "2011-01-03", "2011-01-05"],
                    freq="D",
                    name="idx3",
                ),
            ),
            (
                pd.PeriodIndex(
                    ["2011", "2013", "2015", "2012", "2011"], name="pidx", freq="Y"
                ),
                pd.PeriodIndex(
                    ["2011", "2011", "2012", "2013", "2015"], name="pidx", freq="Y"
                ),
            ),
            (
                # For compatibility check
                pd.Index([2011, 2013, 2015, 2012, 2011], name="idx"),
                pd.Index([2011, 2011, 2012, 2013, 2015], name="idx"),
            ),
        ],
    )
    def test_sort_values_without_freq_periodindex(self, idx, expected):
        # here without_freq means not generateable by period_range
        self.check_sort_values_without_freq(idx, expected)

    def test_sort_values_without_freq_periodindex_nat(self):
        # doesn't quite fit into check_sort_values_without_freq
        idx = pd.PeriodIndex(["2011", "2013", "NaT", "2011"], name="pidx", freq="D")
        expected = pd.PeriodIndex(
            ["NaT", "2011", "2011", "2013"], name="pidx", freq="D"
        )

        ordered = idx.sort_values(na_position="first")
        tm.assert_index_equal(ordered, expected)
        check_freq_nonmonotonic(ordered, idx)

        ordered = idx.sort_values(ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])
        check_freq_nonmonotonic(ordered, idx)


def test_order_stability_compat():
    # GH#35922. sort_values is stable both for normal and datetime-like Index
    pidx = pd.PeriodIndex(
        ["2011", "2013", "2015", "2012", "2011"], name="pidx", freq="Y"
    )
    iidx = pd.Index([2011, 2013, 2015, 2012, 2011], name="idx")
    ordered1, indexer1 = pidx.sort_values(return_indexer=True, ascending=False)
    ordered2, indexer2 = iidx.sort_values(return_indexer=True, ascending=False)
    tm.assert_numpy_array_equal(indexer1, indexer2)
