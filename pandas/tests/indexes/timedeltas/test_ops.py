from datetime import timedelta

import numpy as np
import pytest

from pandas.core.dtypes.generic import ABCDateOffset

import pandas as pd
from pandas import Series, TimedeltaIndex, timedelta_range
import pandas._testing as tm

from pandas.tseries.offsets import Day, Hour


class TestTimedeltaIndexOps:
    def test_value_counts_unique(self):
        # GH 7735
        idx = timedelta_range("1 days 09:00:00", freq="H", periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = TimedeltaIndex(np.repeat(idx.values, range(1, len(idx) + 1)))

        exp_idx = timedelta_range("1 days 18:00:00", freq="-1H", periods=10)
        expected = Series(range(10, 0, -1), index=exp_idx, dtype="int64")

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        expected = timedelta_range("1 days 09:00:00", freq="H", periods=10)
        tm.assert_index_equal(idx.unique(), expected)

        idx = TimedeltaIndex(
            [
                "1 days 09:00:00",
                "1 days 09:00:00",
                "1 days 09:00:00",
                "1 days 08:00:00",
                "1 days 08:00:00",
                pd.NaT,
            ]
        )

        exp_idx = TimedeltaIndex(["1 days 09:00:00", "1 days 08:00:00"])
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = TimedeltaIndex(["1 days 09:00:00", "1 days 08:00:00", pd.NaT])
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    def test_nonunique_contains(self):
        # GH 9512
        for idx in map(
            TimedeltaIndex,
            (
                [0, 1, 0],
                [0, 0, -1],
                [0, -1, -1],
                ["00:01:00", "00:01:00", "00:02:00"],
                ["00:01:00", "00:01:00", "00:00:01"],
            ),
        ):
            assert idx[0] in idx

    def test_unknown_attribute(self):
        # see gh-9680
        tdi = pd.timedelta_range(start=0, periods=10, freq="1s")
        ts = pd.Series(np.random.normal(size=10), index=tdi)
        assert "foo" not in ts.__dict__.keys()
        msg = "'Series' object has no attribute 'foo'"
        with pytest.raises(AttributeError, match=msg):
            ts.foo

    def test_order(self):
        # GH 10295
        idx1 = TimedeltaIndex(["1 day", "2 day", "3 day"], freq="D", name="idx")
        idx2 = TimedeltaIndex(["1 hour", "2 hour", "3 hour"], freq="H", name="idx")

        for idx in [idx1, idx2]:
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, idx)
            assert ordered.freq == idx.freq

            ordered = idx.sort_values(ascending=False)
            expected = idx[::-1]
            tm.assert_index_equal(ordered, expected)
            assert ordered.freq == expected.freq
            assert ordered.freq.n == -1

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, idx)
            tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2]), check_dtype=False)
            assert ordered.freq == idx.freq

            ordered, indexer = idx.sort_values(return_indexer=True, ascending=False)
            tm.assert_index_equal(ordered, idx[::-1])
            assert ordered.freq == expected.freq
            assert ordered.freq.n == -1

        idx1 = TimedeltaIndex(
            ["1 hour", "3 hour", "5 hour", "2 hour ", "1 hour"], name="idx1"
        )
        exp1 = TimedeltaIndex(
            ["1 hour", "1 hour", "2 hour", "3 hour", "5 hour"], name="idx1"
        )

        idx2 = TimedeltaIndex(
            ["1 day", "3 day", "5 day", "2 day", "1 day"], name="idx2"
        )

        # TODO(wesm): unused?
        # exp2 = TimedeltaIndex(['1 day', '1 day', '2 day',
        #                        '3 day', '5 day'], name='idx2')

        # idx3 = TimedeltaIndex([pd.NaT, '3 minute', '5 minute',
        #                        '2 minute', pd.NaT], name='idx3')
        # exp3 = TimedeltaIndex([pd.NaT, pd.NaT, '2 minute', '3 minute',
        #                        '5 minute'], name='idx3')

        for idx, expected in [(idx1, exp1), (idx1, exp1), (idx1, exp1)]:
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, expected)
            assert ordered.freq is None

            ordered = idx.sort_values(ascending=False)
            tm.assert_index_equal(ordered, expected[::-1])
            assert ordered.freq is None

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, expected)

            exp = np.array([0, 4, 3, 1, 2])
            tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            assert ordered.freq is None

            ordered, indexer = idx.sort_values(return_indexer=True, ascending=False)
            tm.assert_index_equal(ordered, expected[::-1])

            exp = np.array([2, 1, 3, 4, 0])
            tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            assert ordered.freq is None

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.timedelta_range("1 day", "31 day", freq="D", name="idx")
        result = idx.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert idx.freq == result.freq

        idx_dup = idx.append(idx)
        assert idx_dup.freq is None  # freq is reset
        result = idx_dup.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert result.freq is None

    def test_drop_duplicates(self):
        # to check Index/Series compat
        base = pd.timedelta_range("1 day", "31 day", freq="D", name="idx")
        idx = base.append(base[:5])

        res = idx.drop_duplicates()
        tm.assert_index_equal(res, base)
        res = Series(idx).drop_duplicates()
        tm.assert_series_equal(res, Series(base))

        res = idx.drop_duplicates(keep="last")
        exp = base[5:].append(base[:5])
        tm.assert_index_equal(res, exp)
        res = Series(idx).drop_duplicates(keep="last")
        tm.assert_series_equal(res, Series(exp, index=np.arange(5, 36)))

        res = idx.drop_duplicates(keep=False)
        tm.assert_index_equal(res, base[5:])
        res = Series(idx).drop_duplicates(keep=False)
        tm.assert_series_equal(res, Series(base[5:], index=np.arange(5, 31)))

    @pytest.mark.parametrize(
        "freq", ["D", "3D", "-3D", "H", "2H", "-2H", "T", "2T", "S", "-3S"]
    )
    def test_infer_freq(self, freq):
        # GH#11018
        idx = pd.timedelta_range("1", freq=freq, periods=10)
        result = pd.TimedeltaIndex(idx.asi8, freq="infer")
        tm.assert_index_equal(idx, result)
        assert result.freq == freq

    def test_shift(self):
        pass  # handled in test_arithmetic.py

    def test_repeat(self):
        index = pd.timedelta_range("1 days", periods=2, freq="D")
        exp = pd.TimedeltaIndex(["1 days", "1 days", "2 days", "2 days"])
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

        index = TimedeltaIndex(["1 days", "NaT", "3 days"])
        exp = TimedeltaIndex(
            [
                "1 days",
                "1 days",
                "1 days",
                "NaT",
                "NaT",
                "NaT",
                "3 days",
                "3 days",
                "3 days",
            ]
        )
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

    def test_nat(self):
        assert pd.TimedeltaIndex._na_value is pd.NaT
        assert pd.TimedeltaIndex([])._na_value is pd.NaT

        idx = pd.TimedeltaIndex(["1 days", "2 days"])
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        assert idx.hasnans is False
        tm.assert_numpy_array_equal(idx._nan_idxs, np.array([], dtype=np.intp))

        idx = pd.TimedeltaIndex(["1 days", "NaT"])
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        assert idx.hasnans is True
        tm.assert_numpy_array_equal(idx._nan_idxs, np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        idx = pd.TimedeltaIndex(["1 days", "2 days", "NaT"])
        assert idx.equals(idx)
        assert idx.equals(idx.copy())
        assert idx.equals(idx.astype(object))
        assert idx.astype(object).equals(idx)
        assert idx.astype(object).equals(idx.astype(object))
        assert not idx.equals(list(idx))
        assert not idx.equals(pd.Series(idx))

        idx2 = pd.TimedeltaIndex(["2 days", "1 days", "NaT"])
        assert not idx.equals(idx2)
        assert not idx.equals(idx2.copy())
        assert not idx.equals(idx2.astype(object))
        assert not idx.astype(object).equals(idx2)
        assert not idx.astype(object).equals(idx2.astype(object))
        assert not idx.equals(list(idx2))
        assert not idx.equals(pd.Series(idx2))

        # Check that we dont raise OverflowError on comparisons outside the
        #  implementation range
        oob = pd.Index([timedelta(days=10 ** 6)] * 3, dtype=object)
        assert not idx.equals(oob)
        assert not idx2.equals(oob)

        # FIXME: oob.apply(np.timedelta64) incorrectly overflows
        oob2 = pd.Index([np.timedelta64(x) for x in oob], dtype=object)
        assert not idx.equals(oob2)
        assert not idx2.equals(oob2)

    @pytest.mark.parametrize("values", [["0 days", "2 days", "4 days"], []])
    @pytest.mark.parametrize("freq", ["2D", Day(2), "48H", Hour(48)])
    def test_freq_setter(self, values, freq):
        # GH 20678
        idx = TimedeltaIndex(values)

        # can set to an offset, converting from string if necessary
        idx._data.freq = freq
        assert idx.freq == freq
        assert isinstance(idx.freq, ABCDateOffset)

        # can reset to None
        idx._data.freq = None
        assert idx.freq is None

    def test_freq_setter_errors(self):
        # GH 20678
        idx = TimedeltaIndex(["0 days", "2 days", "4 days"])

        # setting with an incompatible freq
        msg = (
            "Inferred frequency 2D from passed values does not conform to "
            "passed frequency 5D"
        )
        with pytest.raises(ValueError, match=msg):
            idx._data.freq = "5D"

        # setting with a non-fixed frequency
        msg = r"<2 \* BusinessDays> is a non-fixed frequency"
        with pytest.raises(ValueError, match=msg):
            idx._data.freq = "2B"

        # setting with non-freq string
        with pytest.raises(ValueError, match="Invalid frequency"):
            idx._data.freq = "foo"
