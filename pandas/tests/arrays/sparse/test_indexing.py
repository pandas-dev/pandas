import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.sparse import SparseArray


@pytest.fixture
def arr_data():
    return np.array([np.nan, np.nan, 1, 2, 3, np.nan, 4, 5, np.nan, 6])


@pytest.fixture
def arr(arr_data):
    return SparseArray(arr_data)


class TestGetitem:
    def test_getitem(self, arr):
        dense = arr.to_dense()
        for i, value in enumerate(arr):
            tm.assert_almost_equal(value, dense[i])
            tm.assert_almost_equal(arr[-i], dense[-i])

    def test_getitem_arraylike_mask(self, arr):
        arr = SparseArray([0, 1, 2])
        result = arr[[True, False, True]]
        expected = SparseArray([0, 2])
        tm.assert_sp_array_equal(result, expected)

    @pytest.mark.parametrize(
        "slc",
        [
            np.s_[:],
            np.s_[1:10],
            np.s_[1:100],
            np.s_[10:1],
            np.s_[:-3],
            np.s_[-5:-4],
            np.s_[:-12],
            np.s_[-12:],
            np.s_[2:],
            np.s_[2::3],
            np.s_[::2],
            np.s_[::-1],
            np.s_[::-2],
            np.s_[1:6:2],
            np.s_[:-6:-2],
        ],
    )
    @pytest.mark.parametrize(
        "as_dense", [[np.nan] * 10, [1] * 10, [np.nan] * 5 + [1] * 5, []]
    )
    def test_getslice(self, slc, as_dense):
        as_dense = np.array(as_dense)
        arr = SparseArray(as_dense)

        result = arr[slc]
        expected = SparseArray(as_dense[slc])

        tm.assert_sp_array_equal(result, expected)

    def test_getslice_tuple(self):
        dense = np.array([np.nan, 0, 3, 4, 0, 5, np.nan, np.nan, 0])

        sparse = SparseArray(dense)
        res = sparse[(slice(4, None),)]
        exp = SparseArray(dense[4:])
        tm.assert_sp_array_equal(res, exp)

        sparse = SparseArray(dense, fill_value=0)
        res = sparse[(slice(4, None),)]
        exp = SparseArray(dense[4:], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        msg = "too many indices for array"
        with pytest.raises(IndexError, match=msg):
            sparse[4:, :]

        with pytest.raises(IndexError, match=msg):
            # check numpy compat
            dense[4:, :]

    def test_boolean_slice_empty(self):
        arr = SparseArray([0, 1, 2])
        res = arr[[False, False, False]]
        assert res.dtype == arr.dtype

    def test_getitem_bool_sparse_array(self, arr):
        # GH 23122
        spar_bool = SparseArray([False, True] * 5, dtype=np.bool_, fill_value=True)
        exp = SparseArray([np.nan, 2, np.nan, 5, 6])
        tm.assert_sp_array_equal(arr[spar_bool], exp)

        spar_bool = ~spar_bool
        res = arr[spar_bool]
        exp = SparseArray([np.nan, 1, 3, 4, np.nan])
        tm.assert_sp_array_equal(res, exp)

        spar_bool = SparseArray(
            [False, True, np.nan] * 3, dtype=np.bool_, fill_value=np.nan
        )
        res = arr[spar_bool]
        exp = SparseArray([np.nan, 3, 5])
        tm.assert_sp_array_equal(res, exp)

    def test_getitem_bool_sparse_array_as_comparison(self):
        # GH 45110
        arr = SparseArray([1, 2, 3, 4, np.nan, np.nan], fill_value=np.nan)
        res = arr[arr > 2]
        exp = SparseArray([3.0, 4.0], fill_value=np.nan)
        tm.assert_sp_array_equal(res, exp)

    def test_get_item(self, arr):
        zarr = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)

        assert np.isnan(arr[1])
        assert arr[2] == 1
        assert arr[7] == 5

        assert zarr[0] == 0
        assert zarr[2] == 1
        assert zarr[7] == 5

        errmsg = "index is out of bounds: must be an integer between -10 and 9"

        with pytest.raises(IndexError, match=errmsg):
            arr[11]

        with pytest.raises(IndexError, match=errmsg):
            arr[-11]

        assert arr[-1] == arr[len(arr) - 1]


class TestSetitem:
    def test_set_item(self, arr_data):
        arr = SparseArray(arr_data).copy()

        def setitem():
            arr[5] = 3

        def setslice():
            arr[1:5] = 2

        with pytest.raises(TypeError, match="assignment via setitem"):
            setitem()

        with pytest.raises(TypeError, match="assignment via setitem"):
            setslice()


class TestTake:
    def test_take_scalar_raises(self, arr):
        msg = "'indices' must be an array, not a scalar '2'."
        with pytest.raises(ValueError, match=msg):
            arr.take(2)

    def test_take(self, arr_data, arr):
        exp = SparseArray(np.take(arr_data, [2, 3]))
        tm.assert_sp_array_equal(arr.take([2, 3]), exp)

        exp = SparseArray(np.take(arr_data, [0, 1, 2]))
        tm.assert_sp_array_equal(arr.take([0, 1, 2]), exp)

    def test_take_all_empty(self):
        sparse = pd.array([0, 0], dtype=pd.SparseDtype("int64"))
        result = sparse.take([0, 1], allow_fill=True, fill_value=np.nan)
        tm.assert_sp_array_equal(sparse, result)

    def test_take_different_fill_value(self):
        # Take with a different fill value shouldn't overwrite the original
        sparse = pd.array([0.0], dtype=pd.SparseDtype("float64", fill_value=0.0))
        result = sparse.take([0, -1], allow_fill=True, fill_value=np.nan)
        expected = pd.array([0, np.nan], dtype=sparse.dtype)
        tm.assert_sp_array_equal(expected, result)

    def test_take_fill_bool_upcasts_to_object(self):
        # GH#32119 numpy bool can't hold NA, so taking a fill position from a
        #  boolean SparseArray upcasts to object (matching dense reindex)
        #  rather than raising.
        sparse = SparseArray([False, False, True], fill_value=False)
        result = sparse.take([0, 2, -1], allow_fill=True)
        expected = SparseArray([False, True, np.nan], fill_value=False)
        assert result.dtype == pd.SparseDtype(object, False)
        tm.assert_sp_array_equal(result, expected)

    def test_take_fill_bool_all_fill_upcasts_to_object(self):
        # GH#32119 same as above for an all-fill (sp_index.npoints == 0) array
        sparse = SparseArray([False, False], fill_value=False)
        result = sparse.take([0, -1], allow_fill=True)
        expected = SparseArray([False, np.nan], fill_value=False)
        assert result.dtype == pd.SparseDtype(object, False)
        tm.assert_sp_array_equal(result, expected)

    def test_take_fill_bool_empty_upcasts_to_object(self):
        # GH#68483 same as above for a length-zero array, which takes a
        #  separate branch and kept upcasting to float
        sparse = SparseArray(np.array([], dtype=bool))
        result = sparse.take([-1, -1], allow_fill=True)
        expected = SparseArray([np.nan, np.nan], dtype=pd.SparseDtype(object, False))
        assert result.dtype == pd.SparseDtype(object, False)
        tm.assert_sp_array_equal(result, expected)

    def test_take_fill_bool_old_fill_upcasts_to_object(self):
        # GH#68483 same for the old-fill arm, reached when the array's own
        #  fill value is NA; it kept upcasting True to 1.0
        sparse = SparseArray(
            np.array([True, np.nan], dtype=object),
            dtype=pd.SparseDtype(bool, np.nan),
        )
        result = sparse.take([0, 1], allow_fill=True)
        expected = SparseArray([True, np.nan], dtype=pd.SparseDtype(object, np.nan))
        assert result.dtype == pd.SparseDtype(object, np.nan)
        tm.assert_sp_array_equal(result, expected)

    def test_take_fill_bool_all_sparse_na_fill_upcasts_to_object(self):
        # GH#68483 the all-sparse arm took the subtype straight, so np.full
        #  wrote the NA fill into a bool array and resolved it to True
        sparse = SparseArray(
            np.array([np.nan, np.nan], dtype=object),
            dtype=pd.SparseDtype(bool, np.nan),
        )
        result = sparse.take([1, 0], allow_fill=True)
        expected = SparseArray([np.nan, np.nan], dtype=pd.SparseDtype(object, np.nan))
        assert result.dtype == pd.SparseDtype(object, np.nan)
        tm.assert_sp_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["uint8", "int8", "int16", "float32"])
    def test_take_all_sparse_preserves_narrow_subtype(self, dtype):
        # GH#68483 the all-sparse arm writes a fill its own subtype already
        #  holds, so it must not promote
        sparse = SparseArray(np.zeros(3, dtype=dtype), fill_value=0)
        result = sparse.take([2, 1, 0], allow_fill=True)
        assert result.dtype == sparse.dtype

    def test_reindex_empty_bool_upcasts_to_object(self):
        # GH#68483 the user-visible path onto the branch above
        ser = pd.Series(SparseArray(np.array([], dtype=bool)))
        result = ser.reindex([0, 1])
        expected = pd.Series(
            SparseArray([np.nan, np.nan], dtype=pd.SparseDtype(object, False))
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "subtype", ["int8", "int32", "uint16", "uint64", "float32"]
    )
    def test_take_fill_preserves_subtype(self, subtype):
        # GH#68469 an old fill position holds self.fill_value, which
        #  SparseDtype._check_fill_value guarantees fits the subtype, so it must
        #  not promote; widening uint64 to float64 would round above 2**53
        big = 2**63 + 12345 if subtype == "uint64" else 3
        sparse = SparseArray(np.array([1, 0, big], dtype=subtype), fill_value=0)
        result = sparse.take(np.array([0, 1, 2]), allow_fill=True)
        tm.assert_sp_array_equal(result, sparse)

    def test_take_fill_all_fill_preserves_subtype(self):
        # GH#68469 same for the all-fill (sp_index.npoints == 0) path, which
        #  promoted whenever an index of -1 came along
        sparse = SparseArray(np.array([0, 0, 0], dtype="int8"), fill_value=0)
        result = sparse.take(np.array([0, 1, -1]), allow_fill=True, fill_value=0)
        assert result.dtype == pd.SparseDtype("int8", 0)
        tm.assert_sp_array_equal(result, sparse)

    @pytest.mark.parametrize("unit", ["M8[s]", "M8[ns]", "m8[s]", "m8[ns]"])
    def test_take_fill_datetimelike_subtype(self, unit):
        # GH#68469 a datetimelike subtype was promoted on type(fill_value):
        #  to object for a nanosecond unit, which numpy renders as raw
        #  integers, and to microseconds for a coarser one
        data = np.array([1, 3, 5], dtype=unit)
        fill = pd.Timestamp(data[1]) if unit[0] == "M" else pd.Timedelta(data[1])
        sparse = SparseArray(data, fill_value=fill)
        result = sparse.take(np.array([0, 1, 2]), allow_fill=True)
        tm.assert_sp_array_equal(result, sparse)

    @pytest.mark.parametrize("unit", ["M8[s]", "M8[ns]", "m8[s]", "m8[ns]"])
    def test_take_fill_datetimelike_new_fill_is_nat(self, unit):
        # GH#68469 a new fill position takes NaT in the subtype, as dense
        #  reindex does, rather than promoting the array
        data = np.array([1, 3, 5], dtype=unit)
        fill = pd.Timestamp(data[1]) if unit[0] == "M" else pd.Timedelta(data[1])
        sparse = SparseArray(data, fill_value=fill)
        result = sparse.take(np.array([0, 1, -1]), allow_fill=True)
        expected = SparseArray(np.array([1, 3, "NaT"], dtype=unit), fill_value=fill)
        tm.assert_sp_array_equal(result, expected)

    @pytest.mark.parametrize("subtype", ["int8", "uint64", "float32"])
    def test_take_fill_empty_preserves_subtype(self, subtype):
        # GH#68469 the len(self) == 0 arm promoted on type(fill_value) too, so an
        #  empty array disagreed with a non-empty one on the same input
        sparse = SparseArray(np.array([], dtype=subtype), fill_value=0)
        result = sparse.take(np.array([-1, -1]), allow_fill=True, fill_value=0)
        assert result.dtype == pd.SparseDtype(subtype, 0)

    @pytest.mark.parametrize("unit", ["M8[s]", "M8[ns]", "m8[s]", "m8[ns]"])
    def test_take_fill_empty_datetimelike(self, unit):
        # GH#68469 the same arm raised DTypePromotionError for a datetimelike
        #  subtype, since it promoted against type(np.nan)
        sparse = SparseArray(np.array([], dtype=unit), fill_value=pd.NaT)
        result = sparse.take(np.array([-1, -1]), allow_fill=True)
        assert result.dtype.subtype == np.dtype(unit)

    def test_take_fill_all_fill_all_new(self):
        # an all-fill array taken entirely at -1 has no old fill position at
        #  all; pins that the promotion stays well defined, see GH#68469
        sparse = SparseArray(np.array([0, 0, 0], dtype="int8"), fill_value=0)
        result = sparse.take(np.array([-1, -1]), allow_fill=True)
        expected = SparseArray(np.array([np.nan, np.nan]), fill_value=0)
        tm.assert_sp_array_equal(result, expected)

    def test_take_fill_object_subtype_keeps_boxed_scalar(self):
        # only a datetimelike subtype needs the numpy scalar; pins that an
        #  object subtype keeps the Timestamp boxed, cf GH#68073
        sparse = SparseArray(np.array(["a", "x", "b"], dtype=object), fill_value="x")
        ts = pd.Timestamp("2000-01-01")
        result = sparse.take(np.array([0, -1]), allow_fill=True, fill_value=ts)
        assert result.dtype == pd.SparseDtype(object, "x")
        assert result[1] == ts and isinstance(result[1], pd.Timestamp)

    @pytest.mark.parametrize("unit", ["M8[s]", "M8[ns]", "m8[s]", "m8[ns]"])
    def test_take_fill_datetimelike_no_gaps(self, unit):
        # GH#68469 with no gaps only the new fill promotes, and promoting a
        #  datetimelike subtype against type(np.nan) raised outright
        data = np.array([1, 3, 5], dtype=unit)
        fill = (
            pd.Timestamp(99, unit=unit[3:-1])
            if unit[0] == "M"
            else pd.Timedelta(99, unit[3:-1])
        )
        sparse = SparseArray(data, fill_value=fill)
        assert sparse.sp_index.ngaps == 0
        result = sparse.take(np.array([0, 1, 2, -1]), allow_fill=True)
        assert result.dtype.subtype == np.dtype(unit)
        assert result[3] is pd.NaT

    def test_take_fill_value(self):
        data = np.array([1, np.nan, 0, 3, 0])
        sparse = SparseArray(data, fill_value=0)

        exp = SparseArray(np.take(data, [0]), fill_value=0)
        tm.assert_sp_array_equal(sparse.take([0]), exp)

        exp = SparseArray(np.take(data, [1, 3, 4]), fill_value=0)
        tm.assert_sp_array_equal(sparse.take([1, 3, 4]), exp)

    def test_take_negative(self, arr_data, arr):
        exp = SparseArray(np.take(arr_data, [-1]))
        tm.assert_sp_array_equal(arr.take([-1]), exp)

        exp = SparseArray(np.take(arr_data, [-4, -3, -2]))
        tm.assert_sp_array_equal(arr.take([-4, -3, -2]), exp)

    def test_bad_take(self, arr):
        with pytest.raises(IndexError, match="bounds"):
            arr.take([11])

    def test_take_filling(self):
        # similar tests as GH 12631
        sparse = SparseArray([np.nan, np.nan, 1, np.nan, 4])
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([np.nan, np.nan, 4])
        tm.assert_sp_array_equal(result, expected)

        # TODO: actionable?
        # Note: test change: fill_value=True -> allow_fill=True
        result = sparse.take(np.array([1, 0, -1]), allow_fill=True)
        expected = SparseArray([np.nan, np.nan, np.nan])
        tm.assert_sp_array_equal(result, expected)

        # allow_fill=False
        result = sparse.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = SparseArray([np.nan, np.nan, 4])
        tm.assert_sp_array_equal(result, expected)

        msg = "Invalid value in 'indices'"
        with pytest.raises(ValueError, match=msg):
            sparse.take(np.array([1, 0, -2]), allow_fill=True)

        with pytest.raises(ValueError, match=msg):
            sparse.take(np.array([1, 0, -5]), allow_fill=True)

        msg = "out of bounds value in 'indices'"
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, -6]))
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, 5]))
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, 5]), allow_fill=True)

    def test_take_filling_fill_value(self):
        # same tests as GH#12631
        sparse = SparseArray([np.nan, 0, 1, 0, 4], fill_value=0)
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([0, np.nan, 4], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        # fill_value
        result = sparse.take(np.array([1, 0, -1]), allow_fill=True)
        # TODO: actionable?
        # Note: behavior change.
        # the old way of filling self.fill_value doesn't follow EA rules.
        # It's supposed to be self.dtype.na_value (nan in this case)
        expected = SparseArray([0, np.nan, np.nan], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        # allow_fill=False
        result = sparse.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = SparseArray([0, np.nan, 4], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        msg = "Invalid value in 'indices'."
        with pytest.raises(ValueError, match=msg):
            sparse.take(np.array([1, 0, -2]), allow_fill=True)
        with pytest.raises(ValueError, match=msg):
            sparse.take(np.array([1, 0, -5]), allow_fill=True)

        msg = "out of bounds value in 'indices'"
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, -6]))
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, 5]))
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, 5]), fill_value=True)

    @pytest.mark.parametrize("kind", ["block", "integer"])
    def test_take_filling_all_nan(self, kind):
        sparse = SparseArray([np.nan, np.nan, np.nan, np.nan, np.nan], kind=kind)
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([np.nan, np.nan, np.nan], kind=kind)
        tm.assert_sp_array_equal(result, expected)

        result = sparse.take(np.array([1, 0, -1]), fill_value=True)
        expected = SparseArray([np.nan, np.nan, np.nan], kind=kind)
        tm.assert_sp_array_equal(result, expected)

        msg = "out of bounds value in 'indices'"
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, -6]))
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, 5]))
        with pytest.raises(IndexError, match=msg):
            sparse.take(np.array([1, 5]), fill_value=True)


class TestWhere:
    def test_where_retain_fill_value(self):
        # GH#45691 don't lose fill_value on _where
        arr = SparseArray([np.nan, 1.0], fill_value=0)

        mask = np.array([True, False])

        res = arr._where(~mask, 1)
        exp = SparseArray([1, 1.0], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        ser = pd.Series(arr)
        res = ser.where(~mask, 1)
        tm.assert_series_equal(res, pd.Series(exp))
