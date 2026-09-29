import numpy as np
import pytest

from pandas._libs.sparse import IntIndex
from pandas.errors import IntCastingNaNError

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.sparse import SparseArray


class TestAstype:
    def test_astype(self):
        # float -> float
        arr = SparseArray([None, None, 0, 2])
        result = arr.astype("Sparse[float32]")
        expected = SparseArray([None, None, 0, 2], dtype=np.dtype("float32"))
        tm.assert_sp_array_equal(result, expected)

        # GH#35795 the NaN gaps are values, so asking for a 0 fill_value stores
        #  them rather than turning them into zeros
        dtype = pd.SparseDtype("float64", fill_value=0)
        result = arr.astype(dtype)
        expected = SparseArray._simple_new(
            np.array([np.nan, np.nan, 2.0], dtype=dtype.subtype),
            IntIndex(4, [0, 1, 3]),
            dtype,
        )
        tm.assert_sp_array_equal(result, expected)

        # GH#35795 and so they block an integer subtype, as they do when stored
        with pytest.raises(IntCastingNaNError, match="non-finite values"):
            arr.astype(pd.SparseDtype("int64", 0))

        arr = SparseArray([0, np.nan, 0, 1], fill_value=0)
        with pytest.raises(ValueError, match="NA"):
            arr.astype("Sparse[i8]")

    def test_astype_bool(self):
        a = SparseArray([1, 0, 0, 1], dtype=pd.SparseDtype(int, 0))
        result = a.astype(bool)
        expected = np.array([1, 0, 0, 1], dtype=bool)
        tm.assert_numpy_array_equal(result, expected)

        # update fill value
        result = a.astype(pd.SparseDtype(bool, False))
        expected = SparseArray(
            [True, False, False, True], dtype=pd.SparseDtype(bool, False)
        )
        tm.assert_sp_array_equal(result, expected)

    def test_astype_all(self, any_real_numpy_dtype):
        vals = np.array([1, 2, 3])
        arr = SparseArray(vals, fill_value=1)
        typ = np.dtype(any_real_numpy_dtype)
        res = arr.astype(typ)
        tm.assert_numpy_array_equal(res, vals.astype(any_real_numpy_dtype))

    @pytest.mark.parametrize(
        "arr, dtype, expected",
        [
            (
                SparseArray([0, 1]),
                "float",
                SparseArray([0.0, 1.0], dtype=pd.SparseDtype(float, 0.0)),
            ),
            (SparseArray([0, 1]), bool, SparseArray([False, True])),
            (
                SparseArray([0, 1], fill_value=1),
                bool,
                SparseArray([False, True], dtype=pd.SparseDtype(bool, True)),
            ),
            pytest.param(
                SparseArray([0, 1]),
                "datetime64[ns]",
                SparseArray(
                    np.array([0, 1], dtype="datetime64[ns]"),
                    dtype=pd.SparseDtype("datetime64[ns]", pd.Timestamp("1970")),
                ),
            ),
            (
                SparseArray([0, 1, 10]),
                np.str_,
                SparseArray(["0", "1", "10"], dtype=pd.SparseDtype(np.str_, "0")),
            ),
            (SparseArray(["10", "20"]), float, SparseArray([10.0, 20.0])),
            (
                SparseArray([0, 1, 0]),
                object,
                SparseArray([0, 1, 0], dtype=pd.SparseDtype(object, 0)),
            ),
        ],
    )
    def test_astype_more(self, arr, dtype, expected):
        result = arr.astype(arr.dtype.update_dtype(dtype))
        tm.assert_sp_array_equal(result, expected)

    def test_astype_nan_raises(self):
        arr = SparseArray([1.0, np.nan])
        with pytest.raises(ValueError, match="Cannot convert non-finite"):
            arr.astype(int)

    def test_astype_copy_false(self):
        # GH#34456 bug caused by using .view instead of .astype in astype_nansafe
        arr = SparseArray([1, 2, 3])

        dtype = pd.SparseDtype(float, 0)

        result = arr.astype(dtype, copy=False)
        expected = SparseArray([1.0, 2.0, 3.0], fill_value=0.0)
        tm.assert_sp_array_equal(result, expected)

    def test_astype_dt64_to_int64(self):
        # GH#49631 match non-sparse behavior
        values = np.array(["NaT", "2016-01-02", "2016-01-03"], dtype="M8[ns]")

        arr = SparseArray(values)
        result = arr.astype("int64")
        expected = values.astype("int64")
        tm.assert_numpy_array_equal(result, expected)

        # we should also be able to cast to equivalent Sparse[int64]
        dtype_int64 = pd.SparseDtype("int64", np.iinfo(np.int64).min)
        result2 = arr.astype(dtype_int64)
        tm.assert_numpy_array_equal(result2.to_numpy(), expected)

        # GH#50087 we should match the non-sparse behavior regardless of
        #  if we have a fill_value other than NaT
        dtype = pd.SparseDtype("datetime64[ns]", values[1])
        arr3 = SparseArray(values, dtype=dtype)
        result3 = arr3.astype("int64")
        tm.assert_numpy_array_equal(result3, expected)

    @pytest.mark.parametrize("dtype", ["Sparse[int64]", pd.SparseDtype("int64")])
    def test_astype_dt64_to_sparse_int64_fill_value(self, dtype):
        # GH#49631 converting to Sparse[int64] should convert the NaT
        # fill_value to iNaT, not silently replace it with 0. This must hold
        # both for the string form and for an explicit SparseDtype object
        # (Series.astype resolves the string to an object before dispatching).
        values = np.array(["NaT", "2016-01-02", "2016-01-03"], dtype="M8[ns]")
        arr = SparseArray(values)
        iNaT = np.iinfo(np.int64).min
        expected = SparseArray(
            values.astype("int64"),
            dtype=pd.SparseDtype("int64", fill_value=iNaT),
        )

        result = arr.astype(dtype)
        tm.assert_sp_array_equal(result, expected)

        # GH#49631 the fill_value conversion must survive the Series.astype path
        result = pd.Series(arr).astype(dtype)
        tm.assert_sp_array_equal(result.array, expected)

    def test_astype_dt64_to_sparse_int64_explicit_fill_value(self):
        # GH#49631 an explicitly-requested non-default fill_value must be
        # preserved, not overridden with iNaT
        values = np.array(["NaT", "2016-01-02", "2016-01-03"], dtype="M8[ns]")
        arr = SparseArray(values)
        dtype = pd.SparseDtype("int64", fill_value=5)
        # GH#35795 the NaT gap keeps its value rather than becoming the
        #  requested fill_value, so nothing is left to compress
        expected = SparseArray(values.astype("int64"), dtype=dtype)

        result = arr.astype(dtype)
        tm.assert_sp_array_equal(result, expected)

        result = pd.Series(arr).astype(dtype)
        tm.assert_sp_array_equal(result.array, expected)

    @pytest.mark.parametrize("dtype", ["Sparse[int64]", pd.SparseDtype("int64")])
    @pytest.mark.parametrize(
        "unit_dtype, fill_type",
        [
            ("M8[ns]", np.datetime64),
            ("M8[ns]", pd.Timestamp),
            ("M8[us]", pd.Timestamp),
            ("m8[ns]", np.timedelta64),
            ("m8[ns]", pd.Timedelta),
            ("m8[s]", np.timedelta64),
        ],
    )
    def test_astype_datetimelike_to_sparse_int64_non_na_fill_value(
        self, dtype, unit_dtype, fill_type
    ):
        # GH#49631 a non-NA datetimelike fill_value must be converted too, not
        # silently replaced with 0; the boxed pandas spelling of that fill value
        # must work as well as the numpy one
        values = np.array([1, 1, 2], dtype=unit_dtype)
        arr = SparseArray(values, fill_value=fill_type(values[0]))
        expected = SparseArray(
            values.astype("int64"),
            dtype=pd.SparseDtype("int64", fill_value=values[0].astype("int64")),
        )

        result = arr.astype(dtype)
        tm.assert_sp_array_equal(result, expected)

        result = pd.Series(arr).astype(dtype)
        tm.assert_sp_array_equal(result.array, expected)

    @pytest.mark.parametrize("unit_dtype", ["M8[ns]", "m8[ns]"])
    def test_astype_datetimelike_nat_fill_value_spelling(self, unit_dtype):
        # GH#49631 a NaT fill value spelled as a pandas scalar must convert like
        # the numpy spelling, not raise
        values = np.array(["NaT", 1, 2], dtype=unit_dtype)
        arr = SparseArray(values, dtype=pd.SparseDtype(unit_dtype, pd.NaT))
        expected = SparseArray(
            values.astype("int64"),
            dtype=pd.SparseDtype("int64", fill_value=np.iinfo(np.int64).min),
        )

        result = arr.astype("Sparse[int64]")
        tm.assert_sp_array_equal(result, expected)

    @pytest.mark.parametrize(
        "unit_dtype, fill_type", [("M8[ns]", np.datetime64), ("m8[ns]", np.timedelta64)]
    )
    def test_astype_datetimelike_fully_dense_fill_value(self, unit_dtype, fill_type):
        # GH#49631 the fill_value belongs to the dtype whether or not it happens
        # to occur in the data, so a fully dense array converts it too
        values = np.array([1, 2, 3], dtype=unit_dtype)
        arr = SparseArray(values, fill_value=fill_type(5, "ns"))
        assert arr.sp_index.npoints == len(arr)  # fill_value absent from the data

        result = arr.astype("Sparse[int64]")
        assert result.dtype == pd.SparseDtype("int64", fill_value=5)
        tm.assert_numpy_array_equal(result.to_dense(), values.astype("int64"))

    @pytest.mark.parametrize("unit_dtype", ["M8[ns]", "m8[ns]"])
    def test_astype_datetimelike_to_sparse_bool_fill_value(self, unit_dtype):
        # GH#49631 the fill_value conversion is not int-specific: a bool target
        # must get the fill value's truthiness rather than False
        values = np.array([1, 1, 2], dtype=unit_dtype)
        arr = SparseArray(values, fill_value=values[0])

        result = arr.astype("Sparse[bool]")
        assert result.dtype == pd.SparseDtype(bool, fill_value=True)
        tm.assert_numpy_array_equal(result.to_dense(), values.astype(bool))

    @pytest.mark.parametrize(
        "arr, dtype",
        [
            (SparseArray([1.0, np.nan, 2.0]), pd.SparseDtype("float64", 1.0)),
            (SparseArray([1.0, 0.0, 2.0], fill_value=0.0), pd.SparseDtype("float64")),
        ],
    )
    def test_astype_changes_only_fill_value(self, arr, dtype):
        # GH#68567 the two dtypes compared equal, so astype returned early and
        # the requested fill_value was dropped. What the gaps end up holding is
        # SparseArray.astype's contract, not this one's — see GH#35795.
        assert arr.dtype.subtype == dtype.subtype

        result = arr.astype(dtype)

        assert result.dtype.subtype == dtype.subtype
        # not `result.dtype == dtype`: that routes through the __eq__ under test,
        #  which passes even when astype hands back the source dtype unchanged
        tm.assert_almost_equal(result.dtype.fill_value, dtype.fill_value)

    def test_astype_changes_fill_value_series_and_frame(self):
        # GH#68567 Series/DataFrame.astype short-circuit on dtype equality of
        #  their own, separately from SparseArray.astype
        arr = SparseArray([1.0, np.nan, 2.0])
        dtype = pd.SparseDtype("float64", 1.0)

        assert pd.Series(arr).astype(dtype).dtype.fill_value == 1.0
        assert pd.DataFrame({"a": arr}).astype(dtype).dtypes["a"].fill_value == 1.0

    def test_astype_different_fill_value_keeps_values(self):
        # GH#35795 the positions holding the old fill_value must be stored, not
        #  relabelled with the requested one
        arr = SparseArray([1, 1, 0, 1], fill_value=1)

        result = arr.astype(pd.SparseDtype("int64", 0))

        assert result.dtype == pd.SparseDtype("int64", 0)
        expected = np.array([1, 1, 0, 1], dtype="int64")
        tm.assert_numpy_array_equal(result.to_dense(), expected)

    def test_astype_different_fill_value_does_not_fillna(self):
        # GH#35795 a NaN fill_value is data, so astype must not silently replace
        #  it with the requested fill_value
        arr = SparseArray([np.nan, 0.0, 2.0], fill_value=np.nan)

        result = arr.astype(pd.SparseDtype("float32", 0.0))

        assert result.dtype == pd.SparseDtype("float32", 0.0)
        expected = np.array([np.nan, 0.0, 2.0], dtype="float32")
        tm.assert_numpy_array_equal(result.to_dense(), expected)

    def test_astype_fill_value_absent_from_values(self):
        # GH#35795 asking for a fill_value the data does not contain is legal and
        #  gives a fully dense result rather than changing the values
        arr = SparseArray([0, 0, 1], fill_value=0)

        result = arr.astype(pd.SparseDtype("int64", 9))

        assert result.sp_index.npoints == 3
        expected = np.array([0, 0, 1], dtype="int64")
        tm.assert_numpy_array_equal(result.to_dense(), expected)

    def test_astype_fully_dense_na_fill_to_int_no_raise(self):
        # GH#49631 a fully dense float SparseArray whose (unused) NaN fill_value
        # cannot be represented as an integer must not raise on astype to int
        arr = SparseArray(np.array([1.0, 2.0, 3.0]), fill_value=np.nan)
        assert arr.sp_index.npoints == len(arr)  # fully dense

        result = arr.astype("Sparse[int64]")
        expected = SparseArray(np.array([1, 2, 3], dtype="int64"))
        tm.assert_sp_array_equal(result, expected)

    def test_astype_uint64_dense_exact(self):
        # GH#68573 densifying for astype must not promote on type(fill_value):
        #  np.result_type(uint64, int) is float64, which rounds above 2**53
        values = np.array([1, 0, 2**63 + 12345], dtype="uint64")
        arr = SparseArray(values, fill_value=0)

        result = arr.astype("uint64")
        tm.assert_numpy_array_equal(result, values)

    @pytest.mark.parametrize(
        "unit_dtype, fill_type", [("M8[ns]", pd.Timestamp), ("m8[ns]", pd.Timedelta)]
    )
    @pytest.mark.parametrize("target", ["coarser_unit", "int64", object])
    def test_astype_dense_boxed_datetimelike_fill_value(
        self, unit_dtype, fill_type, target
    ):
        # GH#68573 a boxed Timestamp/Timedelta fill_value made the densified
        #  values object dtype, so the cast saw the stored values as raw integers
        values = np.array([10**9, 2 * 10**9]).astype(unit_dtype)
        arr = SparseArray(values, fill_value=fill_type(values[0]))
        dtype = unit_dtype.replace("ns", "s") if target == "coarser_unit" else target

        result = arr.astype(dtype)
        expected = pd.array(values).astype(dtype)
        tm.assert_equal(result, expected)


def test_astype_frame_different_fill_value_round_trip():
    # GH#35795 the reported case: a frame round-tripped through two sparse
    #  fill values came back all zeros
    df = pd.DataFrame([[1, 1, 1, 1, 0, 1]]).T

    sdf = df.astype(pd.SparseDtype("int64", fill_value=1))
    result = sdf.astype(pd.SparseDtype("int64", fill_value=0))

    tm.assert_frame_equal(result.astype("int64"), df)


def test_astype_preserves_kind():
    # GH#35795 re-sparsifying for a new fill_value must keep the index kind
    arr = SparseArray([0, 0, 1, 2], fill_value=0, kind="block")

    result = arr.astype(pd.SparseDtype("float64", 9.0))

    assert result.kind == "block"
    tm.assert_numpy_array_equal(result.to_dense(), np.array([0.0, 0.0, 1.0, 2.0]))


@pytest.mark.parametrize("subtype", ["uint64", "int64"])
def test_astype_different_fill_value_does_not_widen(subtype):
    # GH#35795 materializing the gaps must keep the subtype's own dense dtype;
    #  going through float64 rounds values that do not fit its mantissa
    values = np.array([0, 2**62 + 1, 0, 2**62 + 3], dtype=subtype)
    arr = SparseArray(values, fill_value=0)

    result = arr.astype(pd.SparseDtype(subtype, 7))

    tm.assert_numpy_array_equal(result.to_dense(), values)


def test_astype_object_target_distinguishes_type():
    # GH#35795 object data tells 0 and False apart, so a numeric fill value that
    #  merely compares equal to the requested one must still re-sparsify
    arr = SparseArray([0, 1, 2], fill_value=0)

    result = arr.astype(pd.SparseDtype(object, False))

    assert result.sp_index.npoints == 3
    tm.assert_numpy_array_equal(np.asarray(result), np.array([0, 1, 2], dtype=object))
