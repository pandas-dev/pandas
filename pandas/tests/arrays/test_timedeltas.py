import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays import TimedeltaArray
import pandas.util.testing as tm


class TestTimedeltaArrayConstructor:
    def test_only_1dim_accepted(self):
        # GH#25282
        arr = np.array([0, 1, 2, 3], dtype="m8[h]").astype("m8[ns]")

        with pytest.raises(ValueError, match="Only 1-dimensional"):
            # 3-dim, we allow 2D to sneak in for ops purposes GH#29853
            TimedeltaArray(arr.reshape(2, 2, 1))

        with pytest.raises(ValueError, match="Only 1-dimensional"):
            # 0-dim
            TimedeltaArray(arr[[0]].squeeze())

    def test_freq_validation(self):
        # ensure that the public constructor cannot create an invalid instance
        arr = np.array([0, 0, 1], dtype=np.int64) * 3600 * 10 ** 9

        msg = (
            "Inferred frequency None from passed values does not "
            "conform to passed frequency D"
        )
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray(arr.view("timedelta64[ns]"), freq="D")

    def test_non_array_raises(self):
        with pytest.raises(ValueError, match="list"):
            TimedeltaArray([1, 2, 3])

    def test_other_type_raises(self):
        with pytest.raises(ValueError, match="dtype bool cannot be converted"):
            TimedeltaArray(np.array([1, 2, 3], dtype="bool"))

    def test_incorrect_dtype_raises(self):
        # TODO: why TypeError for 'category' but ValueError for i8?
        with pytest.raises(
            ValueError, match=r"category cannot be converted to timedelta64\[ns\]"
        ):
            TimedeltaArray(np.array([1, 2, 3], dtype="i8"), dtype="category")

        with pytest.raises(
            ValueError, match=r"dtype int64 cannot be converted to timedelta64\[ns\]",
        ):
            TimedeltaArray(np.array([1, 2, 3], dtype="i8"), dtype=np.dtype("int64"))

    def test_copy(self):
        data = np.array([1, 2, 3], dtype="m8[ns]")
        arr = TimedeltaArray(data, copy=False)
        assert arr._data is data

        arr = TimedeltaArray(data, copy=True)
        assert arr._data is not data
        assert arr._data.base is not data


class TestTimedeltaArray:
    def test_np_sum(self):
        # GH#25282
        vals = np.arange(5, dtype=np.int64).view("m8[h]").astype("m8[ns]")
        arr = TimedeltaArray(vals)
        result = np.sum(arr)
        assert result == vals.sum()

        result = np.sum(pd.TimedeltaIndex(arr))
        assert result == vals.sum()

    def test_from_sequence_dtype(self):
        msg = "dtype .*object.* cannot be converted to timedelta64"
        with pytest.raises(ValueError, match=msg):
            TimedeltaArray._from_sequence([], dtype=object)

    def test_abs(self):
        vals = np.array([-3600 * 10 ** 9, "NaT", 7200 * 10 ** 9], dtype="m8[ns]")
        arr = TimedeltaArray(vals)

        evals = np.array([3600 * 10 ** 9, "NaT", 7200 * 10 ** 9], dtype="m8[ns]")
        expected = TimedeltaArray(evals)

        result = abs(arr)
        tm.assert_timedelta_array_equal(result, expected)

    def test_neg(self):
        vals = np.array([-3600 * 10 ** 9, "NaT", 7200 * 10 ** 9], dtype="m8[ns]")
        arr = TimedeltaArray(vals)

        evals = np.array([3600 * 10 ** 9, "NaT", -7200 * 10 ** 9], dtype="m8[ns]")
        expected = TimedeltaArray(evals)

        result = -arr
        tm.assert_timedelta_array_equal(result, expected)

    def test_neg_freq(self):
        tdi = pd.timedelta_range("2 Days", periods=4, freq="H")
        arr = TimedeltaArray(tdi, freq=tdi.freq)

        expected = TimedeltaArray(-tdi._data, freq=-tdi.freq)

        result = -arr
        tm.assert_timedelta_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", [int, np.int32, np.int64, "uint32", "uint64"])
    def test_astype_int(self, dtype):
        arr = TimedeltaArray._from_sequence([pd.Timedelta("1H"), pd.Timedelta("2H")])
        result = arr.astype(dtype)

        if np.dtype(dtype).kind == "u":
            expected_dtype = np.dtype("uint64")
        else:
            expected_dtype = np.dtype("int64")
        expected = arr.astype(expected_dtype)

        assert result.dtype == expected_dtype
        tm.assert_numpy_array_equal(result, expected)

    def test_setitem_clears_freq(self):
        a = TimedeltaArray(pd.timedelta_range("1H", periods=2, freq="H"))
        a[0] = pd.Timedelta("1H")
        assert a.freq is None

    @pytest.mark.parametrize(
        "obj",
        [
            pd.Timedelta(seconds=1),
            pd.Timedelta(seconds=1).to_timedelta64(),
            pd.Timedelta(seconds=1).to_pytimedelta(),
        ],
    )
    def test_setitem_objects(self, obj):
        # make sure we accept timedelta64 and timedelta in addition to Timedelta
        tdi = pd.timedelta_range("2 Days", periods=4, freq="H")
        arr = TimedeltaArray(tdi, freq=tdi.freq)

        arr[0] = obj
        assert arr[0] == pd.Timedelta(seconds=1)


class TestReductions:
    @pytest.mark.parametrize("name", ["sum", "std", "min", "max", "median"])
    @pytest.mark.parametrize("skipna", [True, False])
    def test_reductions_empty(self, name, skipna):
        tdi = pd.TimedeltaIndex([])
        arr = tdi.array

        result = getattr(tdi, name)(skipna=skipna)
        assert result is pd.NaT

        result = getattr(arr, name)(skipna=skipna)
        assert result is pd.NaT

    def test_min_max(self):
        arr = TimedeltaArray._from_sequence(["3H", "3H", "NaT", "2H", "5H", "4H"])

        result = arr.min()
        expected = pd.Timedelta("2H")
        assert result == expected

        result = arr.max()
        expected = pd.Timedelta("5H")
        assert result == expected

        result = arr.min(skipna=False)
        assert result is pd.NaT

        result = arr.max(skipna=False)
        assert result is pd.NaT

    def test_sum(self):
        tdi = pd.TimedeltaIndex(["3H", "3H", "NaT", "2H", "5H", "4H"])
        arr = tdi.array

        result = arr.sum(skipna=True)
        expected = pd.Timedelta(hours=17)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = tdi.sum(skipna=True)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = arr.sum(skipna=False)
        assert result is pd.NaT

        result = tdi.sum(skipna=False)
        assert result is pd.NaT

        result = arr.sum(min_count=9)
        assert result is pd.NaT

        result = tdi.sum(min_count=9)
        assert result is pd.NaT

        result = arr.sum(min_count=1)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = tdi.sum(min_count=1)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

    def test_npsum(self):
        # GH#25335 np.sum should return a Timedelta, not timedelta64
        tdi = pd.TimedeltaIndex(["3H", "3H", "2H", "5H", "4H"])
        arr = tdi.array

        result = np.sum(tdi)
        expected = pd.Timedelta(hours=17)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = np.sum(arr)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

    def test_std(self):
        tdi = pd.TimedeltaIndex(["0H", "4H", "NaT", "4H", "0H", "2H"])
        arr = tdi.array

        result = arr.std(skipna=True)
        expected = pd.Timedelta(hours=2)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = tdi.std(skipna=True)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = arr.std(skipna=False)
        assert result is pd.NaT

        result = tdi.std(skipna=False)
        assert result is pd.NaT

    def test_median(self):
        tdi = pd.TimedeltaIndex(["0H", "3H", "NaT", "5H06m", "0H", "2H"])
        arr = tdi.array

        result = arr.median(skipna=True)
        expected = pd.Timedelta(hours=2)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = tdi.median(skipna=True)
        assert isinstance(result, pd.Timedelta)
        assert result == expected

        result = arr.std(skipna=False)
        assert result is pd.NaT

        result = tdi.std(skipna=False)
        assert result is pd.NaT
