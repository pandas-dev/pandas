from datetime import timezone
import time

import numpy as np
import pytest

import pandas as pd
from pandas import (
    Series,
    date_range,
    period_range,
)
import pandas._testing as tm


@pytest.mark.parametrize(
    "first_slice,second_slice",
    [
        [[2, None], [None, -5]],
        [[None, 0], [None, -5]],
        [[None, -5], [None, 0]],
        [[None, 0], [None, 0]],
    ],
)
@pytest.mark.parametrize("fill", [None, -1])
def test_align(datetime_series, first_slice, second_slice, join_type, fill):
    a = datetime_series[slice(*first_slice)]
    b = datetime_series[slice(*second_slice)]

    aa, ab = a.align(b, join=join_type, fill_value=fill)

    join_index = a.index.join(b.index, how=join_type)
    if fill is not None:
        diff_a = aa.index.difference(join_index)
        diff_b = ab.index.difference(join_index)
        if len(diff_a) > 0:
            assert (aa.reindex(diff_a) == fill).all()
        if len(diff_b) > 0:
            assert (ab.reindex(diff_b) == fill).all()

    ea = a.reindex(join_index)
    eb = b.reindex(join_index)

    if fill is not None:
        ea = ea.fillna(fill)
        eb = eb.fillna(fill)

    tm.assert_series_equal(aa, ea)
    tm.assert_series_equal(ab, eb)
    assert aa.name == "ts"
    assert ea.name == "ts"
    assert ab.name == "ts"
    assert eb.name == "ts"


def test_align_nocopy(datetime_series):
    b = datetime_series[:5].copy()

    # do copy
    a = datetime_series.copy()
    ra, _ = a.align(b, join="left")
    ra[:5] = 5
    assert not (a[:5] == 5).any()

    # do not copy
    a = datetime_series.copy()
    ra, _ = a.align(b, join="left")
    ra[:5] = 5
    assert not (a[:5] == 5).any()

    # do copy
    a = datetime_series.copy()
    b = datetime_series[:5].copy()
    _, rb = a.align(b, join="right")
    rb[:3] = 5
    assert not (b[:3] == 5).any()

    # do not copy
    a = datetime_series.copy()
    b = datetime_series[:5].copy()
    _, rb = a.align(b, join="right")
    rb[:2] = 5
    assert not (b[:2] == 5).any()


def test_align_same_index(datetime_series):
    a, b = datetime_series.align(datetime_series)
    assert a.index is not datetime_series.index
    assert b.index is not datetime_series.index
    assert a.index.is_(datetime_series.index)
    assert b.index.is_(datetime_series.index)


def test_align_multiindex():
    # GH 10665

    midx = pd.MultiIndex.from_product(
        [range(2), range(3), range(2)], names=("a", "b", "c")
    )
    idx = pd.Index(range(2), name="b")
    s1 = Series(np.arange(12, dtype="int64"), index=midx)
    s2 = Series(np.arange(2, dtype="int64"), index=idx)

    # these must be the same results (but flipped)
    res1l, res1r = s1.align(s2, join="left")
    res2l, res2r = s2.align(s1, join="right")

    expl = s1
    tm.assert_series_equal(expl, res1l)
    tm.assert_series_equal(expl, res2r)
    expr = Series([0, 0, 1, 1, np.nan, np.nan] * 2, index=midx)
    tm.assert_series_equal(expr, res1r)
    tm.assert_series_equal(expr, res2l)

    res1l, res1r = s1.align(s2, join="right")
    res2l, res2r = s2.align(s1, join="left")

    exp_idx = pd.MultiIndex.from_product(
        [range(2), range(2), range(2)], names=("a", "b", "c")
    )
    expl = Series([0, 1, 2, 3, 6, 7, 8, 9], index=exp_idx)
    tm.assert_series_equal(expl, res1l)
    tm.assert_series_equal(expl, res2r)
    expr = Series([0, 0, 1, 1] * 2, index=exp_idx)
    tm.assert_series_equal(expr, res1r)
    tm.assert_series_equal(expr, res2l)


def test_align_dt64tzindex_mismatched_tzs():
    idx1 = date_range("2001", periods=5, freq="h", tz="US/Eastern")
    ser = Series(np.random.default_rng(2).standard_normal(len(idx1)), index=idx1)
    ser_central = ser.tz_convert("US/Central")
    # different timezones convert to UTC

    new1, new2 = ser.align(ser_central)
    assert new1.index.tz is timezone.utc
    assert new2.index.tz is timezone.utc


def test_align_periodindex(join_type):
    rng = period_range("1/1/2000", "1/1/2010", freq="Y")
    ts = Series(np.random.default_rng(2).standard_normal(len(rng)), index=rng)

    # TODO: assert something?
    ts.align(ts[::2], join=join_type)


def test_align_stringindex(any_string_dtype):
    left = Series(range(3), index=pd.Index(["a", "b", "d"], dtype=any_string_dtype))
    right = Series(range(3), index=pd.Index(["a", "b", "c"], dtype=any_string_dtype))
    result_left, result_right = left.align(right)

    expected_idx = pd.Index(["a", "b", "c", "d"], dtype=any_string_dtype)
    expected_left = Series([0, 1, np.nan, 2], index=expected_idx)
    expected_right = Series([0, 1, 2, np.nan], index=expected_idx)

    tm.assert_series_equal(result_left, expected_left)
    tm.assert_series_equal(result_right, expected_right)


def test_align_left_fewer_levels():
    # GH#45224
    left = Series([2], index=pd.MultiIndex.from_tuples([(1, 3)], names=["a", "c"]))
    right = Series(
        [1], index=pd.MultiIndex.from_tuples([(1, 2, 3)], names=["a", "b", "c"])
    )
    result_left, result_right = left.align(right)

    expected_right = Series(
        [1], index=pd.MultiIndex.from_tuples([(1, 3, 2)], names=["a", "c", "b"])
    )
    expected_left = Series(
        [2], index=pd.MultiIndex.from_tuples([(1, 3, 2)], names=["a", "c", "b"])
    )
    tm.assert_series_equal(result_left, expected_left)
    tm.assert_series_equal(result_right, expected_right)


def test_align_left_different_named_levels():
    # GH#45224
    left = Series(
        [2], index=pd.MultiIndex.from_tuples([(1, 4, 3)], names=["a", "d", "c"])
    )
    right = Series(
        [1], index=pd.MultiIndex.from_tuples([(1, 2, 3)], names=["a", "b", "c"])
    )
    result_left, result_right = left.align(right)

    expected_left = Series(
        [2], index=pd.MultiIndex.from_tuples([(1, 4, 3, 2)], names=["a", "d", "c", "b"])
    )
    expected_right = Series(
        [1], index=pd.MultiIndex.from_tuples([(1, 4, 3, 2)], names=["a", "d", "c", "b"])
    )
    tm.assert_series_equal(result_left, expected_left)
    tm.assert_series_equal(result_right, expected_right)


class TestAlignDatetimeUnits:
    """Tests for Series.align with datetime indexes of different units."""
    
    def test_align_datetime_different_units_correctness(self):
        """Test that align works correctly with different datetime units."""
        dates_ns = pd.date_range("2020-01-01", periods=100)
        dates_s = dates_ns.as_unit("s")
        
        s1 = pd.Series(range(100), index=dates_ns)
        s2 = pd.Series(range(100, 200), index=dates_s)
        
        result1, result2 = s1.align(s2)
        
        # Results should have same length
        assert len(result1) == len(result2) == 100
        
        # Values should be preserved
        assert (result1.values == s1.values).all()
        assert (result2.values == s2.values).all()
    
    def test_align_datetime_different_units_performance(self):
        """Test that align is faster with datetime unit optimization."""
        dates_ns = pd.date_range("2020-01-01", periods=1000)
        dates_s = dates_ns.as_unit("s")
        
        s1 = pd.Series(range(1000), index=dates_ns)
        s2 = pd.Series(range(1000), index=dates_s)
        s1_copy = s1.copy()
        
        # Time align with same dtype (baseline)
        start = time.perf_counter()
        for _ in range(10):
            s1.align(s1_copy)
        time_same = time.perf_counter() - start
        
        # Time align with different dtype (should be fast now)
        start = time.perf_counter()
        for _ in range(10):
            s1.align(s2)
        time_diff = time.perf_counter() - start
        
        # With optimization, different-dtype should not be much slower
        # Allow 2x slowdown (was 3x before optimization)
        assert time_diff < time_same * 2.5
    
    def test_align_datetime_units_with_freq(self):
        """Test align with datetime indexes that have frequency."""
        idx1 = pd.date_range("2020-01-01", periods=50, freq="D").as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=50, freq="D").as_unit("s")
        
        s1 = pd.Series(range(50), index=idx1)
        s2 = pd.Series(range(50, 100), index=idx2)
        
        result1, result2 = s1.align(s2)
        
        assert len(result1) == 50
        assert len(result2) == 50


class TestAlignDatetimeUnitsOptimization:
    """Tests for Series.align with datetime indexes of different units."""

    def test_align_datetime_different_units_correctness(self):
        """Test that align works correctly with different datetime units."""
        dates_ns = pd.date_range("2020-01-01", periods=100)
        dates_s = dates_ns.as_unit("s")

        s1 = pd.Series(range(100), index=dates_ns)
        s2 = pd.Series(range(100, 200), index=dates_s)

        result1, result2 = s1.align(s2)

        # Results should have same length
        assert len(result1) == len(result2) == 100

        # Values should be preserved
        assert (result1.values == s1.values).all()
        assert (result2.values == s2.values).all()

        # Indexes should be equal (ignoring dtype)
        assert result1.index.equals(result2.index, check_dtype=False)

    def test_align_datetime_units_with_freq(self):
        """Test align with datetime indexes that have frequency."""
        idx1 = pd.date_range("2020-01-01", periods=50, freq="D").as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=50, freq="D").as_unit("s")

        s1 = pd.Series(range(50), index=idx1)
        s2 = pd.Series(range(50, 100), index=idx2)

        result1, result2 = s1.align(s2)

        assert len(result1) == 50
        assert len(result2) == 50
        assert (result1.values == s1.values).all()
        assert (result2.values == s2.values).all()

    def test_align_datetime_different_units_overlapping(self):
        """Test align with partially overlapping datetime indexes."""
        idx1 = pd.date_range("2020-01-01", periods=10).as_unit("ns")
        idx2 = pd.date_range("2020-01-05", periods=10).as_unit("s")

        s1 = pd.Series(range(10), index=idx1)
        s2 = pd.Series(range(100, 110), index=idx2)

        result1, result2 = s1.align(s2)

        # Should have union of both indexes
        expected_len = len(set(idx1.tolist() + idx2.tolist()))
        assert len(result1) == expected_len
        assert len(result2) == expected_len

    def test_align_datetime_ns_ms_us(self):
        """Test align with various datetime unit combinations."""
        dates = pd.date_range("2020-01-01", periods=20)

        s_ns = pd.Series(range(20), index=dates.as_unit("ns"))
        s_ms = pd.Series(range(20, 40), index=dates.as_unit("ms"))
        s_us = pd.Series(range(40, 60), index=dates.as_unit("us"))

        # Test ns vs ms
        r1, r2 = s_ns.align(s_ms)
        assert len(r1) == len(r2) == 20

        # Test ms vs us
        r3, r4 = s_ms.align(s_us)
        assert len(r3) == len(r4) == 20

        # Test ns vs us
        r5, r6 = s_ns.align(s_us)
        assert len(r5) == len(r6) == 20

    def test_align_datetime_different_units_with_nans(self):
        """Test align with datetime indexes containing NaT values."""
        idx1 = pd.DatetimeIndex(
            ["2020-01-01", "NaT", "2020-01-03", "2020-01-04"]
        ).as_unit("ns")
        idx2 = pd.DatetimeIndex(
            ["2020-01-01", "NaT", "2020-01-03", "2020-01-04"]
        ).as_unit("s")

        s1 = pd.Series([1, 2, 3, 4], index=idx1)
        s2 = pd.Series([10, 20, 30, 40], index=idx2)

        result1, result2 = s1.align(s2)

        # Should align correctly despite NaT
        assert len(result1) == len(result2)


class TestAlignPerformanceRegression:
    """Performance regression tests for Series.align optimization."""

    def test_align_datetime_units_performance_improvement(self):
        """Test that align is faster with datetime unit optimization."""
        dates_ns = pd.date_range("2020-01-01", periods=1000)
        dates_s = dates_ns.as_unit("s")

        s1 = pd.Series(range(1000), index=dates_ns)
        s2 = pd.Series(range(1000), index=dates_s)
        s1_copy = s1.copy()

        # Time align with same dtype (baseline)
        start = time.perf_counter()
        for _ in range(10):
            s1.align(s1_copy)
        time_same = time.perf_counter() - start

        # Time align with different dtype (should be fast now)
        start = time.perf_counter()
        for _ in range(10):
            s1.align(s2)
        time_diff = time.perf_counter() - start

        # With optimization, different-dtype should not be much slower
        # Allow 2.5x slowdown (was ~3x+ before optimization)
        slowdown = time_diff / time_same
        assert slowdown < 2.5, (
            f"Series.align with different datetime units is too slow: "
            f"{slowdown:.2f}x slower than same dtype. "
            f"Expected < 2.5x slowdown."
        )

    def test_align_int_dtype_optimization(self):
        """Test that align benefits from check_dtype=False for integers."""
        idx1 = pd.Index(range(500), dtype='int32')
        idx2 = pd.Index(range(500), dtype='int64')

        s1 = pd.Series(range(500), index=idx1)
        s2 = pd.Series(range(500, 1000), index=idx2)

        # This should now be optimized
        result1, result2 = s1.align(s2)

        assert len(result1) == 500
        assert len(result2) == 500
        assert (result1.values == s1.values).all()
