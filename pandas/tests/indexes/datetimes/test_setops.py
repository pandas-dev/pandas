from datetime import (
    datetime,
    timedelta,
    timezone,
)

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    DataFrame,
    DateOffset,
    DatetimeIndex,
    Index,
    Series,
    Timestamp,
    bdate_range,
    date_range,
)
import pandas._testing as tm
from pandas.core.indexes.datetimelike import _base_freq_ordinal_diff

from pandas.tseries.offsets import (
    BMonthEnd,
    Minute,
    MonthEnd,
)

START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)


class TestDatetimeIndexSetOps:
    tz = [
        None,
        "UTC",
        "Asia/Tokyo",
        "US/Eastern",
        "dateutil/Asia/Singapore",
        "dateutil/US/Pacific",
    ]

    # TODO: moved from test_datetimelike; dedup with version below
    def test_union2(self, sort):
        everything = date_range("2020-01-01", periods=10)
        first = everything[:5]
        second = everything[5:]
        union = first.union(second, sort=sort)
        tm.assert_index_equal(union, everything)

    @pytest.mark.parametrize("box", [np.array, Series, list])
    def test_union3(self, sort, box):
        everything = date_range("2020-01-01", periods=10)
        first = everything[:5]
        second = everything[5:]

        # GH 10149 support listlike inputs other than Index objects
        expected = first.union(second, sort=sort)
        case = box(second.values)
        result = first.union(case, sort=sort)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("tz", tz)
    def test_union(self, tz, sort):
        rng1 = date_range("1/1/2000", freq="D", periods=5, tz=tz, unit="ns")
        other1 = date_range("1/6/2000", freq="D", periods=5, tz=tz, unit="ns")
        expected1 = date_range("1/1/2000", freq="D", periods=10, tz=tz, unit="ns")
        expected1_notsorted = DatetimeIndex(list(other1) + list(rng1))

        rng2 = date_range("1/1/2000", freq="D", periods=5, tz=tz, unit="ns")
        other2 = date_range("1/4/2000", freq="D", periods=5, tz=tz, unit="ns")
        expected2 = date_range("1/1/2000", freq="D", periods=8, tz=tz, unit="ns")
        expected2_notsorted = DatetimeIndex(list(other2) + list(rng2[:3]))

        rng3 = date_range("1/1/2000", freq="D", periods=5, tz=tz, unit="ns")
        other3 = DatetimeIndex([], tz=tz).as_unit("ns")
        expected3 = date_range("1/1/2000", freq="D", periods=5, tz=tz, unit="ns")
        expected3_notsorted = rng3

        for rng, other, exp, exp_notsorted in [
            (rng1, other1, expected1, expected1_notsorted),
            (rng2, other2, expected2, expected2_notsorted),
            (rng3, other3, expected3, expected3_notsorted),
        ]:
            result_union = rng.union(other, sort=sort)
            tm.assert_index_equal(result_union, exp)

            result_union = other.union(rng, sort=sort)
            if sort is None:
                tm.assert_index_equal(result_union, exp)
            else:
                tm.assert_index_equal(result_union, exp_notsorted)

    def test_union_coverage(self, sort):
        idx = DatetimeIndex(["2000-01-03", "2000-01-01", "2000-01-02"])
        ordered = DatetimeIndex(idx.sort_values(), freq="infer")
        result = ordered.union(idx, sort=sort)
        tm.assert_index_equal(result, ordered)

        result = ordered[:0].union(ordered, sort=sort)
        tm.assert_index_equal(result, ordered)
        assert result.freq == ordered.freq

    def test_union_bug_1730(self, sort):
        rng_a = date_range("1/1/2012", periods=4, freq="3h")
        rng_b = date_range("1/1/2012", periods=4, freq="4h")

        result = rng_a.union(rng_b, sort=sort)
        exp = list(rng_a) + list(rng_b[1:])
        if sort is None:
            exp = DatetimeIndex(sorted(exp))
        else:
            exp = DatetimeIndex(exp)
        tm.assert_index_equal(result, exp)

    def test_union_bug_1745(self, sort):
        left = DatetimeIndex(["2012-05-11 15:19:49.695000"])
        right = DatetimeIndex(
            [
                "2012-05-29 13:04:21.322000",
                "2012-05-11 15:27:24.873000",
                "2012-05-11 15:31:05.350000",
            ]
        )

        result = left.union(right, sort=sort)
        exp = DatetimeIndex(
            [
                "2012-05-11 15:19:49.695000",
                "2012-05-29 13:04:21.322000",
                "2012-05-11 15:27:24.873000",
                "2012-05-11 15:31:05.350000",
            ]
        )
        if sort is None:
            exp = exp.sort_values()
        tm.assert_index_equal(result, exp)

    def test_union_bug_4564(self, sort):
        left = date_range("2013-01-01", "2013-02-01")
        right = left + DateOffset(minutes=15)

        result = left.union(right, sort=sort)
        exp = list(left) + list(right)
        if sort is None:
            exp = DatetimeIndex(sorted(exp))
        else:
            exp = DatetimeIndex(exp)
        tm.assert_index_equal(result, exp)

    def test_union_freq_both_none(self, sort):
        # GH11086
        expected = bdate_range("20150101", periods=10)
        expected._data.freq = None

        result = expected.union(expected, sort=sort)
        tm.assert_index_equal(result, expected)
        assert result.freq is None

    def test_union_freq_infer(self):
        # When taking the union of two DatetimeIndexes, we infer
        #  a freq even if the arguments don't have freq.  This matches
        #  TimedeltaIndex behavior.
        dti = date_range("2016-01-01", periods=5)
        left = dti[[0, 1, 3, 4]]
        right = dti[[2, 3, 1]]

        assert left.freq is None
        assert right.freq is None

        result = left.union(right)
        tm.assert_index_equal(result, dti)
        assert result.freq == "D"

    def test_union_dataframe_index(self):
        rng1 = date_range("1/1/1999", "1/1/2012", freq="MS")
        s1 = Series(np.random.default_rng(2).standard_normal(len(rng1)), rng1)

        rng2 = date_range("1/1/1980", "12/1/2001", freq="MS")
        s2 = Series(np.random.default_rng(2).standard_normal(len(rng2)), rng2)
        df = DataFrame({"s1": s1, "s2": s2})

        exp = date_range("1/1/1980", "1/1/2012", freq="MS")
        tm.assert_index_equal(df.index, exp)

    def test_union_with_DatetimeIndex(self, sort):
        i1 = Index(np.arange(0, 20, 2, dtype=np.int64))
        i2 = date_range(start="2012-01-03 00:00:00", periods=10, freq="D")
        # Works
        i1.union(i2, sort=sort)
        # Fails with "AttributeError: can't set attribute"
        i2.union(i1, sort=sort)

    def test_union_same_timezone_different_units(self):
        # GH 55238
        idx1 = date_range("2000-01-01", periods=3, tz="UTC").as_unit("ms")
        idx2 = date_range("2000-01-01", periods=3, tz="UTC").as_unit("us")
        result = idx1.union(idx2)
        expected = date_range("2000-01-01", periods=3, tz="UTC").as_unit("us")
        tm.assert_index_equal(result, expected)

    def test_union_same_nonzero_timezone_different_units(self):
        # GH 60080 - fix timezone being changed to UTC when units differ
        # but timezone is the same
        tz = "UTC+05:00"
        idx1 = date_range("2000-01-01", periods=3, tz=tz).as_unit("us")
        idx2 = date_range("2000-01-01", periods=3, tz=tz).as_unit("ns")

        # Check pre-conditions
        assert idx1.tz == idx2.tz
        assert idx1.dtype != idx2.dtype  # Different units

        # Test union preserves timezone when units differ
        result = idx1.union(idx2)
        expected = date_range("2000-01-01", periods=3, tz=tz).as_unit("ns")
        tm.assert_index_equal(result, expected)

    def test_union_different_dates_same_timezone_different_units(self):
        # GH 60080 - fix timezone being changed to UTC when units differ
        # but timezone is the same
        tz = "UTC+05:00"
        idx1 = date_range("2000-01-01", periods=3, tz=tz).as_unit("us")
        idx3 = date_range("2000-01-03", periods=3, tz=tz).as_unit("us")

        # Test with different dates to ensure it's not just returning one of the inputs
        result = idx1.union(idx3)
        expected = DatetimeIndex(
            ["2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04", "2000-01-05"],
            tz=tz,
        ).as_unit("us")
        tm.assert_index_equal(result, expected)

    def test_intersection_same_timezone_different_units(self):
        # GH 60080 - fix timezone being changed to UTC when units differ
        # but timezone is the same
        tz = "UTC+05:00"
        idx1 = date_range("2000-01-01", periods=3, tz=tz).as_unit("us")
        idx2 = date_range("2000-01-01", periods=3, tz=tz).as_unit("ns")

        # Check pre-conditions
        assert idx1.tz == idx2.tz
        assert idx1.dtype != idx2.dtype  # Different units

        # Test intersection
        result = idx1.intersection(idx2)
        expected = date_range("2000-01-01", periods=3, tz=tz)
        tm.assert_index_equal(result, expected)

        # Intersection dtype is not commutative
        result2 = idx2.intersection(idx1)
        tm.assert_index_equal(result2, expected.as_unit("ns"))

    def test_symmetric_difference_same_timezone_different_units(self):
        # GH 60080 - fix timezone being changed to UTC when units differ
        # but timezone is the same
        tz = "UTC+05:00"
        idx1 = date_range("2000-01-01", periods=3, tz=tz).as_unit("us")
        idx4 = date_range("2000-01-02", periods=3, tz=tz).as_unit("ns")

        # Check pre-conditions
        assert idx1.tz == idx4.tz
        assert idx1.dtype != idx4.dtype  # Different units

        # Test symmetric_difference
        result = idx1.symmetric_difference(idx4)
        expected = DatetimeIndex(["2000-01-01", "2000-01-04"], tz=tz).as_unit("ns")
        tm.assert_index_equal(result, expected)

    # TODO: moved from test_datetimelike; de-duplicate with version below
    def test_intersection2(self):
        first = date_range("2020-01-01", periods=10)
        second = first[5:]
        intersect = first.intersection(second)
        tm.assert_index_equal(intersect, second)

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.intersection(case)
            tm.assert_index_equal(result, second)

        third = Index(["a", "b", "c"])
        result = first.intersection(third)
        expected = Index([], dtype=object)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "tz", [None, "Asia/Tokyo", "US/Eastern", "dateutil/US/Pacific"]
    )
    def test_intersection(self, tz, sort):
        # GH 4690 (with tz)
        base = date_range("6/1/2000", "6/30/2000", freq="D", name="idx", unit="ns")

        # if target has the same name, it is preserved
        rng2 = date_range("5/15/2000", "6/20/2000", freq="D", name="idx", unit="ns")
        expected2 = date_range("6/1/2000", "6/20/2000", freq="D", name="idx", unit="ns")

        # if target name is different, it will be reset
        rng3 = date_range("5/15/2000", "6/20/2000", freq="D", name="other", unit="ns")
        expected3 = date_range("6/1/2000", "6/20/2000", freq="D", name=None, unit="ns")

        rng4 = date_range("7/1/2000", "7/31/2000", freq="D", name="idx", unit="ns")
        expected4 = DatetimeIndex([], freq="D", name="idx", dtype="M8[ns]")

        for rng, expected in [
            (rng2, expected2),
            (rng3, expected3),
            (rng4, expected4),
        ]:
            result = base.intersection(rng)
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

        # non-monotonic
        base = DatetimeIndex(
            ["2011-01-05", "2011-01-04", "2011-01-02", "2011-01-03"], tz=tz, name="idx"
        ).as_unit("ns")

        rng2 = DatetimeIndex(
            ["2011-01-04", "2011-01-02", "2011-02-02", "2011-02-03"], tz=tz, name="idx"
        ).as_unit("ns")
        expected2 = DatetimeIndex(
            ["2011-01-04", "2011-01-02"], tz=tz, name="idx"
        ).as_unit("ns")

        rng3 = DatetimeIndex(
            ["2011-01-04", "2011-01-02", "2011-02-02", "2011-02-03"],
            tz=tz,
            name="other",
        ).as_unit("ns")
        expected3 = DatetimeIndex(
            ["2011-01-04", "2011-01-02"], tz=tz, name=None
        ).as_unit("ns")

        # GH 7880
        rng4 = date_range(
            "7/1/2000", "7/31/2000", freq="D", tz=tz, unit="ns", name="idx"
        )
        expected4 = DatetimeIndex([], tz=tz, name="idx").as_unit("ns")
        assert expected4.freq is None

        for rng, expected in [
            (rng2, expected2),
            (rng3, expected3),
            (rng4, expected4),
        ]:
            result = base.intersection(rng, sort=sort)
            if sort is None:
                expected = expected.sort_values()
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

    # parametrize over both anchored and non-anchored freqs, as they
    #  have different code paths
    @pytest.mark.parametrize("freq", ["min", "B"])
    def test_intersection_empty(self, tz_aware_fixture, freq):
        # empty same freq GH2129
        tz = tz_aware_fixture
        rng = date_range("6/1/2000", "6/15/2000", freq=freq, tz=tz)
        result = rng[0:0].intersection(rng)
        assert len(result) == 0
        assert result.freq == rng.freq

        result = rng.intersection(rng[0:0])
        assert len(result) == 0
        assert result.freq == rng.freq

        # no overlap GH#33604
        check_freq = freq != "min"  # We don't preserve freq on non-anchored offsets
        result = rng[:3].intersection(rng[-3:])
        tm.assert_index_equal(result, rng[:0])
        if check_freq:
            # We don't preserve freq on non-anchored offsets
            assert result.freq == rng.freq

        # swapped left and right
        result = rng[-3:].intersection(rng[:3])
        tm.assert_index_equal(result, rng[:0])
        if check_freq:
            # We don't preserve freq on non-anchored offsets
            assert result.freq == rng.freq

    def test_intersection_bug_1708(self):
        index_1 = date_range("1/1/2012", periods=4, freq="12h")
        index_2 = index_1 + DateOffset(hours=1)

        result = index_1.intersection(index_2)
        assert len(result) == 0

    @pytest.mark.parametrize("tz", tz)
    def test_difference(self, tz, sort):
        rng_dates = ["1/2/2000", "1/3/2000", "1/1/2000", "1/4/2000", "1/5/2000"]

        rng1 = DatetimeIndex(rng_dates, tz=tz)
        other1 = date_range("1/6/2000", freq="D", periods=5, tz=tz)
        expected1 = DatetimeIndex(rng_dates, tz=tz)

        rng2 = DatetimeIndex(rng_dates, tz=tz)
        other2 = date_range("1/4/2000", freq="D", periods=5, tz=tz)
        expected2 = DatetimeIndex(rng_dates[:3], tz=tz)

        rng3 = DatetimeIndex(rng_dates, tz=tz)
        other3 = DatetimeIndex([], tz=tz)
        expected3 = DatetimeIndex(rng_dates, tz=tz)

        for rng, other, expected in [
            (rng1, other1, expected1),
            (rng2, other2, expected2),
            (rng3, other3, expected3),
        ]:
            result_diff = rng.difference(other, sort)
            if sort is None and len(other):
                # We dont sort (yet?) when empty GH#24959
                expected = expected.sort_values()
            tm.assert_index_equal(result_diff, expected)

    def test_difference_freq(self, sort):
        # GH14323: difference of DatetimeIndex should not preserve frequency

        index = date_range("20160920", "20160925", freq="D", unit="ns")
        other = date_range("20160921", "20160924", freq="D", unit="ns")
        expected = DatetimeIndex(["20160920", "20160925"], dtype="M8[ns]", freq=None)
        idx_diff = index.difference(other, sort)
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal("freq", idx_diff, expected)

        # preserve frequency when the difference is a contiguous
        # subset of the original range
        other = date_range("20160922", "20160925", freq="D", unit="ns")
        idx_diff = index.difference(other, sort)
        expected = DatetimeIndex(["20160920", "20160921"], dtype="M8[ns]", freq="D")
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal("freq", idx_diff, expected)

    def test_datetimeindex_diff(self, sort):
        dti1 = date_range(freq="QE-JAN", start=datetime(1997, 12, 31), periods=100)
        dti2 = date_range(freq="QE-JAN", start=datetime(1997, 12, 31), periods=98)
        assert len(dti1.difference(dti2, sort)) == 2

    @pytest.mark.parametrize("tz", [None, "Asia/Tokyo", "US/Eastern"])
    def test_setops_preserve_freq(self, tz):
        rng = date_range("1/1/2000", "1/1/2002", name="idx", tz=tz)

        result = rng[:50].union(rng[50:100])
        assert result.name == rng.name
        assert result.freq == rng.freq
        assert result.tz == rng.tz

        result = rng[:50].union(rng[30:100])
        assert result.name == rng.name
        assert result.freq == rng.freq
        assert result.tz == rng.tz

        result = rng[:50].union(rng[60:100])
        assert result.name == rng.name
        assert result.freq is None
        assert result.tz == rng.tz

        result = rng[:50].intersection(rng[25:75])
        assert result.name == rng.name
        assert result.freqstr == "D"
        assert result.tz == rng.tz

        nofreq = DatetimeIndex(list(rng[25:75]), name="other")
        result = rng[:50].union(nofreq)
        assert result.name is None
        assert result.freq == rng.freq
        assert result.tz == rng.tz

        result = rng[:50].intersection(nofreq)
        assert result.name is None
        assert result.freq == rng.freq
        assert result.tz == rng.tz

    def test_intersection_non_tick_no_fastpath(self):
        # GH#42104
        dti = DatetimeIndex(
            [
                "2018-12-31",
                "2019-03-31",
                "2019-06-30",
                "2019-09-30",
                "2019-12-31",
                "2020-03-31",
            ],
            freq="QE-DEC",
        )
        result = dti[::2].intersection(dti[1::2])
        expected = dti[:0]
        tm.assert_index_equal(result, expected)

    def test_dti_intersection(self):
        rng = date_range("1/1/2011", periods=100, freq="h", tz="utc")

        left = rng[10:90][::-1]
        right = rng[20:80][::-1]

        assert left.tz == rng.tz
        result = left.intersection(right)
        assert result.tz == left.tz

    # Note: not difference, as there is no symmetry requirement there
    @pytest.mark.parametrize("setop", ["union", "intersection", "symmetric_difference"])
    def test_dti_setop_aware(self, setop):
        # non-overlapping
        # GH#39328 as of 2.0 we cast these to UTC instead of object
        rng = date_range("2012-11-15 00:00:00", periods=6, freq="h", tz="US/Central")

        rng2 = date_range("2012-11-15 12:00:00", periods=6, freq="h", tz="US/Eastern")

        result = getattr(rng, setop)(rng2)

        left = rng.tz_convert("UTC")
        right = rng2.tz_convert("UTC")
        expected = getattr(left, setop)(right)
        tm.assert_index_equal(result, expected)
        assert result.tz == left.tz
        if len(result):
            assert result[0].tz is timezone.utc
            assert result[-1].tz is timezone.utc

    def test_dti_union_mixed(self):
        # GH#21671
        rng = DatetimeIndex([Timestamp("2011-01-01"), pd.NaT])
        rng2 = DatetimeIndex(["2012-01-01", "2012-01-02"], tz="Asia/Tokyo")
        result = rng.union(rng2)
        expected = Index(
            [
                Timestamp("2011-01-01"),
                pd.NaT,
                Timestamp("2012-01-01", tz="Asia/Tokyo"),
                Timestamp("2012-01-02", tz="Asia/Tokyo"),
            ],
            dtype=object,
        )
        tm.assert_index_equal(result, expected)


class TestBusinessDatetimeIndex:
    def test_union(self, sort):
        rng = bdate_range(START, END)
        # overlapping
        left = rng[:10]
        right = rng[5:10]

        the_union = left.union(right, sort=sort)
        assert isinstance(the_union, DatetimeIndex)

        # non-overlapping, gap in middle
        left = rng[:5]
        right = rng[10:]

        the_union = left.union(right, sort=sort)
        assert isinstance(the_union, Index)

        # non-overlapping, no gap
        left = rng[:5]
        right = rng[5:10]

        the_union = left.union(right, sort=sort)
        assert isinstance(the_union, DatetimeIndex)

        # order does not matter
        if sort is None:
            tm.assert_index_equal(right.union(left, sort=sort), the_union)
        else:
            expected = DatetimeIndex(list(right) + list(left))
            tm.assert_index_equal(right.union(left, sort=sort), expected)

        # overlapping, but different offset
        rng = date_range(START, END, freq=BMonthEnd())

        the_union = rng.union(rng, sort=sort)
        assert isinstance(the_union, DatetimeIndex)

    def test_union_not_cacheable(self, sort):
        rng = date_range("1/1/2000", periods=50, freq=Minute())
        rng1 = rng[10:]
        rng2 = rng[:25]
        the_union = rng1.union(rng2, sort=sort)
        if sort is None:
            tm.assert_index_equal(the_union, rng)
        else:
            expected = DatetimeIndex(list(rng[10:]) + list(rng[:10]))
            tm.assert_index_equal(the_union, expected)

        rng1 = rng[10:]
        rng2 = rng[15:35]
        the_union = rng1.union(rng2, sort=sort)
        expected = rng[10:]
        tm.assert_index_equal(the_union, expected)

    def test_intersection(self):
        rng = date_range("1/1/2000", periods=50, freq=Minute(), unit="ns")
        rng1 = rng[10:]
        rng2 = rng[:25]
        the_int = rng1.intersection(rng2)
        expected = rng[10:25]
        tm.assert_index_equal(the_int, expected)
        assert isinstance(the_int, DatetimeIndex)
        assert the_int.freq == rng.freq

        the_int = rng1.intersection(rng2)
        tm.assert_index_equal(the_int, expected)

        # non-overlapping
        the_int = rng[:10].intersection(rng[10:])
        expected = DatetimeIndex([]).as_unit("ns")
        tm.assert_index_equal(the_int, expected)

    def test_intersection_bug(self):
        # GH #771
        a = bdate_range("11/30/2011", "12/31/2011")
        b = bdate_range("12/10/2011", "12/20/2011")
        result = a.intersection(b)
        tm.assert_index_equal(result, b)
        assert result.freq == b.freq

    def test_intersection_list(self):
        # GH#35876
        # values is not an Index -> no name -> retain "a"
        values = [Timestamp("2020-01-01"), Timestamp("2020-02-01")]
        idx = DatetimeIndex(values, name="a")
        res = idx.intersection(values)
        tm.assert_index_equal(res, idx)

    def test_month_range_union_tz_pytz(self, sort):
        pytz = pytest.importorskip("pytz")
        tz = pytz.timezone("US/Eastern")

        early_start = datetime(2011, 1, 1)
        early_end = datetime(2011, 3, 1)

        late_start = datetime(2011, 3, 1)
        late_end = datetime(2011, 5, 1)

        early_dr = date_range(start=early_start, end=early_end, tz=tz, freq=MonthEnd())
        late_dr = date_range(start=late_start, end=late_end, tz=tz, freq=MonthEnd())

        early_dr.union(late_dr, sort=sort)

    @td.skip_if_windows
    def test_month_range_union_tz_dateutil(self, sort):
        from pandas._libs.tslibs.timezones import dateutil_gettz

        tz = dateutil_gettz("US/Eastern")

        early_start = datetime(2011, 1, 1)
        early_end = datetime(2011, 3, 1)

        late_start = datetime(2011, 3, 1)
        late_end = datetime(2011, 5, 1)

        early_dr = date_range(start=early_start, end=early_end, tz=tz, freq=MonthEnd())
        late_dr = date_range(start=late_start, end=late_end, tz=tz, freq=MonthEnd())

        early_dr.union(late_dr, sort=sort)

    @pytest.mark.parametrize("sort", [False, None])
    def test_intersection_duplicates(self, sort):
        # GH#38196
        idx1 = Index(
            [
                Timestamp("2019-12-13"),
                Timestamp("2019-12-12"),
                Timestamp("2019-12-12"),
            ]
        )
        result = idx1.intersection(idx1, sort=sort)
        expected = Index([Timestamp("2019-12-13"), Timestamp("2019-12-12")])
        tm.assert_index_equal(result, expected)


class TestCustomDatetimeIndex:
    def test_union(self, sort):
        # overlapping
        rng = bdate_range(START, END, freq="C")
        left = rng[:10]
        right = rng[5:10]

        the_union = left.union(right, sort=sort)
        assert isinstance(the_union, DatetimeIndex)

        # non-overlapping, gap in middle
        left = rng[:5]
        right = rng[10:]

        the_union = left.union(right, sort)
        assert isinstance(the_union, Index)

        # non-overlapping, no gap
        left = rng[:5]
        right = rng[5:10]

        the_union = left.union(right, sort=sort)
        assert isinstance(the_union, DatetimeIndex)

        # order does not matter
        if sort is None:
            tm.assert_index_equal(right.union(left, sort=sort), the_union)

        # overlapping, but different offset
        rng = date_range(START, END, freq=BMonthEnd())

        the_union = rng.union(rng, sort=sort)
        assert isinstance(the_union, DatetimeIndex)

    def test_intersection_bug(self):
        # GH #771
        a = bdate_range("11/30/2011", "12/31/2011", freq="C")
        b = bdate_range("12/10/2011", "12/20/2011", freq="C")
        result = a.intersection(b)
        tm.assert_index_equal(result, b)
        assert result.freq == b.freq

    @pytest.mark.parametrize(
        "tz", [None, "UTC", "Europe/Berlin", timezone(timedelta(hours=-1))]
    )
    def test_intersection_dst_transition(self, tz):
        # GH 46702: Europe/Berlin has DST transition
        idx1 = date_range("2020-03-27", periods=5, freq="D", tz=tz)
        idx2 = date_range("2020-03-30", periods=5, freq="D", tz=tz)
        result = idx1.intersection(idx2)
        expected = date_range("2020-03-30", periods=2, freq="D", tz=tz)
        tm.assert_index_equal(result, expected)

        # GH#45863 same problem for union
        index1 = date_range("2021-10-28", periods=3, freq="D", tz="Europe/London")
        index2 = date_range("2021-10-30", periods=4, freq="D", tz="Europe/London")
        result = index1.union(index2)
        expected = date_range("2021-10-28", periods=6, freq="D", tz="Europe/London")
        tm.assert_index_equal(result, expected)


def test_union_non_nano_rangelike():
    # GH 59036
    l1 = DatetimeIndex(
        ["2024-05-11", "2024-05-12"], dtype="datetime64[us]", name="Date", freq="D"
    )
    l2 = DatetimeIndex(["2024-05-13"], dtype="datetime64[us]", name="Date", freq="D")
    result = l1.union(l2)
    expected = DatetimeIndex(
        ["2024-05-11", "2024-05-12", "2024-05-13"],
        dtype="datetime64[us]",
        name="Date",
        freq="D",
    )
    tm.assert_index_equal(result, expected)


def test_intersection_non_nano_rangelike():
    # GH 59271
    l1 = date_range("2024-01-01", "2024-01-03", unit="s")
    l2 = date_range("2024-01-02", "2024-01-04", unit="s")
    result = l1.intersection(l2)
    expected = DatetimeIndex(
        ["2024-01-02", "2024-01-03"],
        dtype="datetime64[s]",
        freq="D",
    )
    tm.assert_index_equal(result, expected)


def test_union_across_dst_boundary():
    # GH#62915 union should work when one index ends at DST boundary
    # and the other extends past it
    index1 = date_range("2025-10-25", "2025-10-26", freq="D", tz="Europe/Helsinki")
    index2 = date_range("2025-10-25", "2025-10-28", freq="D", tz="Europe/Helsinki")
    result = index1.union(index2)
    expected = date_range("2025-10-25", "2025-10-28", freq="D", tz="Europe/Helsinki")
    tm.assert_index_equal(result, expected)


class TestCanFastIntersect:
    """Tests for _can_fast_intersect alignment checks (GH#44025).

    _can_fast_intersect must return False when two DatetimeIndexes with
    matching freq would not "line up", which causes _fast_intersect to
    produce wrong results.  Three categories of misalignment are covered:

    1. Variable-stride offsets (DateOffset / Week(weekday=None))
    2. Non-normalizing offsets with different time-of-day
    3. n > 1 phase misalignment
    """

    # --- Variable-stride offsets (always return False) ---

    @pytest.mark.parametrize(
        "freq",
        [
            DateOffset(days=3),
            DateOffset(months=1),
            DateOffset(months=1, days=5),
            DateOffset(weeks=1),
            DateOffset(years=1),
        ],
    )
    def test_dateoffset_never_fast(self, freq):
        # GH#44025 DateOffset has no fixed stride
        dti1 = date_range("2021-01-01", periods=5, freq=freq)
        dti2 = date_range("2021-01-02", periods=5, freq=freq)
        assert not dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert len(result) == 0

    def test_week_no_weekday_misaligned_not_fast(self):
        # GH#44025 Week(weekday=None) with different day-of-week means
        #  the grids don't align so _fast_intersect can't be used.
        freq = pd.offsets.Week()
        assert freq.weekday is None

        dti1 = date_range("2021-10-04", periods=5, freq=freq)  # Monday
        dti2 = date_range("2021-10-05", periods=5, freq=freq)  # Tuesday
        assert not dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert len(result) == 0

    def test_week_no_weekday_aligned_fast(self):
        # GH#44025 Week(weekday=None) with same day-of-week is fast.
        freq = pd.offsets.Week()
        dti1 = date_range("2021-10-04", periods=5, freq=freq)
        dti2 = date_range("2021-10-18", periods=5, freq=freq)
        assert dti1._can_fast_intersect(dti2)

        # Result is still correct via the slow path
        result = dti1.intersection(dti2)
        assert set(result) == set(dti1).intersection(set(dti2))

    # --- Non-normalizing offsets with time-of-day mismatch ---

    @pytest.mark.parametrize(
        "freq",
        [
            pd.offsets.CDay(1),
            pd.offsets.BDay(1),
            pd.offsets.MonthEnd(1),
            pd.offsets.Week(weekday=0),
        ],
    )
    def test_different_times(self, freq):
        # GH#44025
        ts1 = Timestamp("2021-10-13 09:00")
        ts2 = Timestamp("2021-10-13 10:00")
        dti1 = date_range(start=ts1, periods=10, freq=freq)
        dti2 = date_range(start=ts2, periods=10, freq=freq)
        assert not dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert len(result) == 0

    @pytest.mark.parametrize(
        "freq",
        [
            pd.offsets.CDay(1),
            pd.offsets.BDay(1),
            pd.offsets.MonthEnd(1),
            pd.offsets.Week(weekday=0),
        ],
    )
    def test_same_times_fast_path_preserved(self, freq):
        # GH#44025 fast path should still be used when times match
        ts = Timestamp("2021-10-13 09:00")
        dti1 = date_range(start=ts, periods=10, freq=freq)
        dti2 = date_range(start=dti1[2], periods=10, freq=freq)
        assert dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert set(result) == set(dti1).intersection(set(dti2))
        assert len(result) > 0

    def test_normalize_on_offset_does_not_affect_date_range(self):
        # GH#44025 normalize=True on the offset does not normalize
        # date_range output, so different start times still misalign.
        freq = pd.offsets.CDay(1, normalize=True)
        dti1 = date_range("2021-10-13 09:00", periods=10, freq=freq)
        dti2 = date_range("2021-10-13 10:00", periods=10, freq=freq)
        assert not dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert len(result) == 0

    # --- n > 1 phase alignment ---

    @pytest.mark.parametrize(
        "freq_str, base_start, aligned_start, misaligned_start",
        [
            ("2ME", "2021-01-31", "2021-03-31", "2021-02-28"),
            ("2QE", "2021-03-31", "2021-09-30", "2021-06-30"),
            ("2YE", "2020-12-31", "2022-12-31", "2021-12-31"),
            ("2W-MON", "2021-10-04", "2021-10-18", "2021-10-11"),
        ],
    )
    def test_phase_alignment(
        self, freq_str, base_start, aligned_start, misaligned_start
    ):
        # GH#44025
        dti_base = date_range(base_start, periods=6, freq=freq_str)
        dti_aligned = date_range(aligned_start, periods=6, freq=freq_str)
        dti_misaligned = date_range(misaligned_start, periods=6, freq=freq_str)

        assert dti_base._can_fast_intersect(dti_aligned)
        result = dti_base.intersection(dti_aligned)
        assert set(result) == set(dti_base).intersection(set(dti_aligned))
        assert len(result) > 0

        assert not dti_base._can_fast_intersect(dti_misaligned)
        result = dti_base.intersection(dti_misaligned)
        assert len(result) == 0

    @pytest.mark.parametrize(
        "n, start1, start2, expect_overlap",
        [
            (3, "2021-01-31", "2021-04-30", True),
            (3, "2021-01-31", "2021-02-28", False),
            (4, "2021-01-31", "2021-05-31", True),
            (4, "2021-01-31", "2021-03-31", False),
        ],
    )
    def test_monthend_various_n(self, n, start1, start2, expect_overlap):
        # GH#44025
        freq = pd.offsets.MonthEnd(n)
        dti1 = date_range(start1, periods=6, freq=freq)
        dti2 = date_range(start2, periods=6, freq=freq)

        expected_set = set(dti1).intersection(set(dti2))
        assert (len(expected_set) > 0) == expect_overlap

        assert dti1._can_fast_intersect(dti2) == expect_overlap
        result = dti1.intersection(dti2)
        assert set(result) == expected_set

    def test_semimonthend_phase(self):
        # GH#44025
        freq = pd.offsets.SemiMonthEnd(2)
        dti1 = date_range("2021-01-15", periods=6, freq=freq)
        dti2 = date_range("2021-02-15", periods=6, freq=freq)  # aligned
        dti3 = date_range("2021-01-31", periods=6, freq=freq)  # misaligned

        assert dti1._can_fast_intersect(dti2)
        assert not dti1._can_fast_intersect(dti3)

        result = dti1.intersection(dti3)
        assert len(result) == 0

    def test_halfyearend_phase(self):
        # GH#44025
        freq = pd.offsets.HalfYearEnd(2)
        dti1 = date_range("2021-06-30", periods=4, freq=freq)
        dti2 = date_range("2022-06-30", periods=4, freq=freq)  # aligned
        dti3 = date_range("2021-12-31", periods=4, freq=freq)  # misaligned

        assert dti1._can_fast_intersect(dti2)
        assert not dti1._can_fast_intersect(dti3)

        result = dti1.intersection(dti3)
        assert len(result) == 0

    def test_n_gt1_with_time_mismatch(self):
        # GH#44025 time mismatch caught before ordinal check
        freq = pd.offsets.MonthEnd(2)
        dti1 = date_range("2021-01-31 09:00", periods=5, freq=freq)
        dti2 = date_range("2021-03-31 10:00", periods=5, freq=freq)
        assert not dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert len(result) == 0

    def test_n_gt1_unsupported_offset_falls_back(self):
        # GH#44025 BDay doesn't support ordinal computation for n > 1
        dti1 = date_range("2021-10-04", periods=5, freq="2B")
        dti2 = date_range("2021-10-06", periods=5, freq="2B")
        assert not dti1._can_fast_intersect(dti2)

        result = dti1.intersection(dti2)
        assert set(result) == set(dti1).intersection(set(dti2))


class TestBaseFreqOrdinalDiff:
    """Unit tests for _base_freq_ordinal_diff helper."""

    def test_month_offset(self):
        ts1 = Timestamp("2021-01-31")
        ts2 = Timestamp("2021-04-30")
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.MonthEnd(2)) == -3
        assert _base_freq_ordinal_diff(ts2, ts1, pd.offsets.MonthEnd(2)) == 3

    def test_quarter_offset(self):
        ts1 = Timestamp("2021-03-31")
        ts2 = Timestamp("2021-09-30")
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.QuarterEnd(2)) == -2

    def test_year_offset(self):
        ts1 = Timestamp("2020-12-31")
        ts2 = Timestamp("2023-12-31")
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.YearEnd(2)) == -3

    def test_week_offset(self):
        ts1 = Timestamp("2021-10-04")  # Monday
        ts2 = Timestamp("2021-10-18")  # Monday, 2 weeks later
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.Week(weekday=0)) == -2

    def test_semimonth_offset(self):
        ts1 = Timestamp("2021-01-15")
        ts2 = Timestamp("2021-01-31")
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.SemiMonthEnd(2)) == -1

    def test_halfyear_offset(self):
        ts1 = Timestamp("2021-06-30")
        ts2 = Timestamp("2022-06-30")
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.HalfYearEnd(2)) == -2

    def test_unsupported_returns_none(self):
        ts1 = Timestamp("2021-10-04")
        ts2 = Timestamp("2021-10-06")
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.BDay(2)) is None
        assert _base_freq_ordinal_diff(ts1, ts2, pd.offsets.CDay(2)) is None
