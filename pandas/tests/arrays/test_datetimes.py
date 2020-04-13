"""
Tests for DatetimeArray
"""
import operator

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import DatetimeArray
from pandas.core.arrays.datetimes import sequence_to_dt64ns


class TestDatetimeArrayConstructor:
    def test_from_sequence_invalid_type(self):
        mi = pd.MultiIndex.from_product([np.arange(5), np.arange(5)])
        with pytest.raises(TypeError, match="Cannot create a DatetimeArray"):
            DatetimeArray._from_sequence(mi)

    def test_only_1dim_accepted(self):
        arr = np.array([0, 1, 2, 3], dtype="M8[h]").astype("M8[ns]")

        with pytest.raises(ValueError, match="Only 1-dimensional"):
            # 3-dim, we allow 2D to sneak in for ops purposes GH#29853
            DatetimeArray(arr.reshape(2, 2, 1))

        with pytest.raises(ValueError, match="Only 1-dimensional"):
            # 0-dim
            DatetimeArray(arr[[0]].squeeze())

    def test_freq_validation(self):
        # GH#24623 check that invalid instances cannot be created with the
        #  public constructor
        arr = np.arange(5, dtype=np.int64) * 3600 * 10 ** 9

        msg = (
            "Inferred frequency H from passed values does not "
            "conform to passed frequency W-SUN"
        )
        with pytest.raises(ValueError, match=msg):
            DatetimeArray(arr, freq="W")

    @pytest.mark.parametrize(
        "meth",
        [
            DatetimeArray._from_sequence,
            sequence_to_dt64ns,
            pd.to_datetime,
            pd.DatetimeIndex,
        ],
    )
    def test_mixing_naive_tzaware_raises(self, meth):
        # GH#24569
        arr = np.array([pd.Timestamp("2000"), pd.Timestamp("2000", tz="CET")])

        msg = (
            "Cannot mix tz-aware with tz-naive values|"
            "Tz-aware datetime.datetime cannot be converted "
            "to datetime64 unless utc=True"
        )

        for obj in [arr, arr[::-1]]:
            # check that we raise regardless of whether naive is found
            #  before aware or vice-versa
            with pytest.raises(ValueError, match=msg):
                meth(obj)

    def test_from_pandas_array(self):
        arr = pd.array(np.arange(5, dtype=np.int64)) * 3600 * 10 ** 9

        result = DatetimeArray._from_sequence(arr, freq="infer")

        expected = pd.date_range("1970-01-01", periods=5, freq="H")._data
        tm.assert_datetime_array_equal(result, expected)

    def test_mismatched_timezone_raises(self):
        arr = DatetimeArray(
            np.array(["2000-01-01T06:00:00"], dtype="M8[ns]"),
            dtype=DatetimeTZDtype(tz="US/Central"),
        )
        dtype = DatetimeTZDtype(tz="US/Eastern")
        with pytest.raises(TypeError, match="Timezone of the array"):
            DatetimeArray(arr, dtype=dtype)

    def test_non_array_raises(self):
        with pytest.raises(ValueError, match="list"):
            DatetimeArray([1, 2, 3])

    def test_bool_dtype_raises(self):
        arr = np.array([1, 2, 3], dtype="bool")

        with pytest.raises(
            ValueError, match="The dtype of 'values' is incorrect.*bool"
        ):
            DatetimeArray(arr)

        msg = r"dtype bool cannot be converted to datetime64\[ns\]"
        with pytest.raises(TypeError, match=msg):
            DatetimeArray._from_sequence(arr)

        with pytest.raises(TypeError, match=msg):
            sequence_to_dt64ns(arr)

        with pytest.raises(TypeError, match=msg):
            pd.DatetimeIndex(arr)

        with pytest.raises(TypeError, match=msg):
            pd.to_datetime(arr)

    def test_incorrect_dtype_raises(self):
        with pytest.raises(ValueError, match="Unexpected value for 'dtype'."):
            DatetimeArray(np.array([1, 2, 3], dtype="i8"), dtype="category")

    def test_freq_infer_raises(self):
        with pytest.raises(ValueError, match="Frequency inference"):
            DatetimeArray(np.array([1, 2, 3], dtype="i8"), freq="infer")

    def test_copy(self):
        data = np.array([1, 2, 3], dtype="M8[ns]")
        arr = DatetimeArray(data, copy=False)
        assert arr._data is data

        arr = DatetimeArray(data, copy=True)
        assert arr._data is not data


class TestDatetimeArrayComparisons:
    # TODO: merge this into tests/arithmetic/test_datetime64 once it is
    #  sufficiently robust

    def test_cmp_dt64_arraylike_tznaive(self, all_compare_operators):
        # arbitrary tz-naive DatetimeIndex
        opname = all_compare_operators.strip("_")
        op = getattr(operator, opname)

        dti = pd.date_range("2016-01-1", freq="MS", periods=9, tz=None)
        arr = DatetimeArray(dti)
        assert arr.freq == dti.freq
        assert arr.tz == dti.tz

        right = dti

        expected = np.ones(len(arr), dtype=bool)
        if opname in ["ne", "gt", "lt"]:
            # for these the comparisons should be all-False
            expected = ~expected

        result = op(arr, arr)
        tm.assert_numpy_array_equal(result, expected)
        for other in [right, np.array(right)]:
            # TODO: add list and tuple, and object-dtype once those
            #  are fixed in the constructor
            result = op(arr, other)
            tm.assert_numpy_array_equal(result, expected)

            result = op(other, arr)
            tm.assert_numpy_array_equal(result, expected)


class TestDatetimeArray:
    def test_astype_to_same(self):
        arr = DatetimeArray._from_sequence(["2000"], tz="US/Central")
        result = arr.astype(DatetimeTZDtype(tz="US/Central"), copy=False)
        assert result is arr

    @pytest.mark.parametrize("dtype", ["datetime64[ns]", "datetime64[ns, UTC]"])
    @pytest.mark.parametrize(
        "other", ["datetime64[ns]", "datetime64[ns, UTC]", "datetime64[ns, CET]"]
    )
    def test_astype_copies(self, dtype, other):
        # https://github.com/pandas-dev/pandas/pull/32490
        s = pd.Series([1, 2], dtype=dtype)
        orig = s.copy()
        t = s.astype(other)
        t[:] = pd.NaT
        tm.assert_series_equal(s, orig)

    @pytest.mark.parametrize("dtype", [int, np.int32, np.int64, "uint32", "uint64"])
    def test_astype_int(self, dtype):
        arr = DatetimeArray._from_sequence([pd.Timestamp("2000"), pd.Timestamp("2001")])
        result = arr.astype(dtype)

        if np.dtype(dtype).kind == "u":
            expected_dtype = np.dtype("uint64")
        else:
            expected_dtype = np.dtype("int64")
        expected = arr.astype(expected_dtype)

        assert result.dtype == expected_dtype
        tm.assert_numpy_array_equal(result, expected)

    def test_tz_setter_raises(self):
        arr = DatetimeArray._from_sequence(["2000"], tz="US/Central")
        with pytest.raises(AttributeError, match="tz_localize"):
            arr.tz = "UTC"

    def test_setitem_different_tz_raises(self):
        data = np.array([1, 2, 3], dtype="M8[ns]")
        arr = DatetimeArray(data, copy=False, dtype=DatetimeTZDtype(tz="US/Central"))
        with pytest.raises(TypeError, match="Cannot compare tz-naive and tz-aware"):
            arr[0] = pd.Timestamp("2000")

        with pytest.raises(ValueError, match="US/Central"):
            arr[0] = pd.Timestamp("2000", tz="US/Eastern")

    def test_setitem_clears_freq(self):
        a = DatetimeArray(pd.date_range("2000", periods=2, freq="D", tz="US/Central"))
        a[0] = pd.Timestamp("2000", tz="US/Central")
        assert a.freq is None

    @pytest.mark.parametrize(
        "obj",
        [
            pd.Timestamp.now(),
            pd.Timestamp.now().to_datetime64(),
            pd.Timestamp.now().to_pydatetime(),
        ],
    )
    def test_setitem_objects(self, obj):
        # make sure we accept datetime64 and datetime in addition to Timestamp
        dti = pd.date_range("2000", periods=2, freq="D")
        arr = dti._data

        arr[0] = obj
        assert arr[0] == obj

    def test_repeat_preserves_tz(self):
        dti = pd.date_range("2000", periods=2, freq="D", tz="US/Central")
        arr = DatetimeArray(dti)

        repeated = arr.repeat([1, 1])

        # preserves tz and values, but not freq
        expected = DatetimeArray(arr.asi8, freq=None, dtype=arr.dtype)
        tm.assert_equal(repeated, expected)

    def test_value_counts_preserves_tz(self):
        dti = pd.date_range("2000", periods=2, freq="D", tz="US/Central")
        arr = DatetimeArray(dti).repeat([4, 3])

        result = arr.value_counts()

        # Note: not tm.assert_index_equal, since `freq`s do not match
        assert result.index.equals(dti)

        arr[-2] = pd.NaT
        result = arr.value_counts()
        expected = pd.Series([1, 4, 2], index=[pd.NaT, dti[0], dti[1]])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("method", ["pad", "backfill"])
    def test_fillna_preserves_tz(self, method):
        dti = pd.date_range("2000-01-01", periods=5, freq="D", tz="US/Central")
        arr = DatetimeArray(dti, copy=True)
        arr[2] = pd.NaT

        fill_val = dti[1] if method == "pad" else dti[3]
        expected = DatetimeArray._from_sequence(
            [dti[0], dti[1], fill_val, dti[3], dti[4]], freq=None, tz="US/Central"
        )

        result = arr.fillna(method=method)
        tm.assert_extension_array_equal(result, expected)

        # assert that arr and dti were not modified in-place
        assert arr[2] is pd.NaT
        assert dti[2] == pd.Timestamp("2000-01-03", tz="US/Central")

    def test_array_interface_tz(self):
        tz = "US/Central"
        data = DatetimeArray(pd.date_range("2017", periods=2, tz=tz))
        result = np.asarray(data)

        expected = np.array(
            [
                pd.Timestamp("2017-01-01T00:00:00", tz=tz),
                pd.Timestamp("2017-01-02T00:00:00", tz=tz),
            ],
            dtype=object,
        )
        tm.assert_numpy_array_equal(result, expected)

        result = np.asarray(data, dtype=object)
        tm.assert_numpy_array_equal(result, expected)

        result = np.asarray(data, dtype="M8[ns]")

        expected = np.array(
            ["2017-01-01T06:00:00", "2017-01-02T06:00:00"], dtype="M8[ns]"
        )
        tm.assert_numpy_array_equal(result, expected)

    def test_array_interface(self):
        data = DatetimeArray(pd.date_range("2017", periods=2))
        expected = np.array(
            ["2017-01-01T00:00:00", "2017-01-02T00:00:00"], dtype="datetime64[ns]"
        )

        result = np.asarray(data)
        tm.assert_numpy_array_equal(result, expected)

        result = np.asarray(data, dtype=object)
        expected = np.array(
            [pd.Timestamp("2017-01-01T00:00:00"), pd.Timestamp("2017-01-02T00:00:00")],
            dtype=object,
        )
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("index", [True, False])
    def test_searchsorted_different_tz(self, index):
        data = np.arange(10, dtype="i8") * 24 * 3600 * 10 ** 9
        arr = DatetimeArray(data, freq="D").tz_localize("Asia/Tokyo")
        if index:
            arr = pd.Index(arr)

        expected = arr.searchsorted(arr[2])
        result = arr.searchsorted(arr[2].tz_convert("UTC"))
        assert result == expected

        expected = arr.searchsorted(arr[2:6])
        result = arr.searchsorted(arr[2:6].tz_convert("UTC"))
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("index", [True, False])
    def test_searchsorted_tzawareness_compat(self, index):
        data = np.arange(10, dtype="i8") * 24 * 3600 * 10 ** 9
        arr = DatetimeArray(data, freq="D")
        if index:
            arr = pd.Index(arr)

        mismatch = arr.tz_localize("Asia/Tokyo")

        msg = "Cannot compare tz-naive and tz-aware datetime-like objects"
        with pytest.raises(TypeError, match=msg):
            arr.searchsorted(mismatch[0])
        with pytest.raises(TypeError, match=msg):
            arr.searchsorted(mismatch)

        with pytest.raises(TypeError, match=msg):
            mismatch.searchsorted(arr[0])
        with pytest.raises(TypeError, match=msg):
            mismatch.searchsorted(arr)

    @pytest.mark.parametrize(
        "other",
        [
            1,
            np.int64(1),
            1.0,
            np.timedelta64("NaT"),
            pd.Timedelta(days=2),
            "invalid",
            np.arange(10, dtype="i8") * 24 * 3600 * 10 ** 9,
            np.arange(10).view("timedelta64[ns]") * 24 * 3600 * 10 ** 9,
            pd.Timestamp.now().to_period("D"),
        ],
    )
    @pytest.mark.parametrize("index", [True, False])
    def test_searchsorted_invalid_types(self, other, index):
        data = np.arange(10, dtype="i8") * 24 * 3600 * 10 ** 9
        arr = DatetimeArray(data, freq="D")
        if index:
            arr = pd.Index(arr)

        msg = "|".join(
            [
                "searchsorted requires compatible dtype or scalar",
                "Unexpected type for 'value'",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            arr.searchsorted(other)


class TestSequenceToDT64NS:
    def test_tz_dtype_mismatch_raises(self):
        arr = DatetimeArray._from_sequence(["2000"], tz="US/Central")
        with pytest.raises(TypeError, match="data is already tz-aware"):
            sequence_to_dt64ns(arr, dtype=DatetimeTZDtype(tz="UTC"))

    def test_tz_dtype_matches(self):
        arr = DatetimeArray._from_sequence(["2000"], tz="US/Central")
        result, _, _ = sequence_to_dt64ns(arr, dtype=DatetimeTZDtype(tz="US/Central"))
        tm.assert_numpy_array_equal(arr._data, result)


class TestReductions:
    @pytest.mark.parametrize("tz", [None, "US/Central"])
    def test_min_max(self, tz):
        arr = DatetimeArray._from_sequence(
            [
                "2000-01-03",
                "2000-01-03",
                "NaT",
                "2000-01-02",
                "2000-01-05",
                "2000-01-04",
            ],
            tz=tz,
        )

        result = arr.min()
        expected = pd.Timestamp("2000-01-02", tz=tz)
        assert result == expected

        result = arr.max()
        expected = pd.Timestamp("2000-01-05", tz=tz)
        assert result == expected

        result = arr.min(skipna=False)
        assert result is pd.NaT

        result = arr.max(skipna=False)
        assert result is pd.NaT

    @pytest.mark.parametrize("tz", [None, "US/Central"])
    @pytest.mark.parametrize("skipna", [True, False])
    def test_min_max_empty(self, skipna, tz):
        arr = DatetimeArray._from_sequence([], tz=tz)
        result = arr.min(skipna=skipna)
        assert result is pd.NaT

        result = arr.max(skipna=skipna)
        assert result is pd.NaT
