import numpy as np
import pytest

from pandas.compat.numpy import _np_version_under1p16

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.numpy_ import PandasArray, PandasDtype

from . import base


@pytest.fixture(params=["float", "object"])
def dtype(request):
    return PandasDtype(np.dtype(request.param))


@pytest.fixture
def allow_in_pandas(monkeypatch):
    """
    A monkeypatch to tells pandas to let us in.

    By default, passing a PandasArray to an index / series / frame
    constructor will unbox that PandasArray to an ndarray, and treat
    it as a non-EA column. We don't want people using EAs without
    reason.

    The mechanism for this is a check against ABCPandasArray
    in each constructor.

    But, for testing, we need to allow them in pandas. So we patch
    the _typ of PandasArray, so that we evade the ABCPandasArray
    check.
    """
    with monkeypatch.context() as m:
        m.setattr(PandasArray, "_typ", "extension")
        yield


@pytest.fixture
def data(allow_in_pandas, dtype):
    if dtype.numpy_dtype == "object":
        return pd.Series([(i,) for i in range(100)]).array
    return PandasArray(np.arange(1, 101, dtype=dtype._dtype))


@pytest.fixture
def data_missing(allow_in_pandas, dtype):
    # For NumPy <1.16, np.array([np.nan, (1,)]) raises
    # ValueError: setting an array element with a sequence.
    if dtype.numpy_dtype == "object":
        if _np_version_under1p16:
            raise pytest.skip("Skipping for NumPy <1.16")
        return PandasArray(np.array([np.nan, (1,)], dtype=object))
    return PandasArray(np.array([np.nan, 1.0]))


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def na_cmp():
    def cmp(a, b):
        return np.isnan(a) and np.isnan(b)

    return cmp


@pytest.fixture
def data_for_sorting(allow_in_pandas, dtype):
    """Length-3 array with a known sort order.

    This should be three items [B, C, A] with
    A < B < C
    """
    if dtype.numpy_dtype == "object":
        # Use an empty tuple for first element, then remove,
        # to disable np.array's shape inference.
        return PandasArray(np.array([(), (2,), (3,), (1,)], dtype=object)[1:])
    return PandasArray(np.array([1, 2, 0]))


@pytest.fixture
def data_missing_for_sorting(allow_in_pandas, dtype):
    """Length-3 array with a known sort order.

    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    if dtype.numpy_dtype == "object":
        return PandasArray(np.array([(1,), np.nan, (0,)], dtype=object))
    return PandasArray(np.array([1, np.nan, 0]))


@pytest.fixture
def data_for_grouping(allow_in_pandas, dtype):
    """Data for factorization, grouping, and unique tests.

    Expected to be like [B, B, NA, NA, A, A, B, C]

    Where A < B < C and NA is missing
    """
    if dtype.numpy_dtype == "object":
        a, b, c = (1,), (2,), (3,)
    else:
        a, b, c = np.arange(3)
    return PandasArray(
        np.array([b, b, np.nan, np.nan, a, a, b, c], dtype=dtype.numpy_dtype)
    )


@pytest.fixture
def skip_numpy_object(dtype):
    """
    Tests for PandasArray with nested data. Users typically won't create
    these objects via `pd.array`, but they can show up through `.array`
    on a Series with nested data. Many of the base tests fail, as they aren't
    appropriate for nested data.

    This fixture allows these tests to be skipped when used as a usefixtures
    marker to either an individual test or a test class.
    """
    if dtype == "object":
        raise pytest.skip("Skipping for object dtype.")


skip_nested = pytest.mark.usefixtures("skip_numpy_object")


class BaseNumPyTests:
    pass


class TestCasting(BaseNumPyTests, base.BaseCastingTests):
    @skip_nested
    def test_astype_str(self, data):
        # ValueError: setting an array element with a sequence
        super().test_astype_str(data)


class TestConstructors(BaseNumPyTests, base.BaseConstructorsTests):
    @pytest.mark.skip(reason="We don't register our dtype")
    # We don't want to register. This test should probably be split in two.
    def test_from_dtype(self, data):
        pass

    @skip_nested
    def test_array_from_scalars(self, data):
        # ValueError: PandasArray must be 1-dimensional.
        super().test_array_from_scalars(data)


class TestDtype(BaseNumPyTests, base.BaseDtypeTests):
    @pytest.mark.skip(reason="Incorrect expected.")
    # we unsurprisingly clash with a NumPy name.
    def test_check_dtype(self, data):
        pass


class TestGetitem(BaseNumPyTests, base.BaseGetitemTests):
    @skip_nested
    def test_getitem_scalar(self, data):
        # AssertionError
        super().test_getitem_scalar(data)

    @skip_nested
    def test_take_series(self, data):
        # ValueError: PandasArray must be 1-dimensional.
        super().test_take_series(data)

    @pytest.mark.xfail(reason="astype doesn't recognize data.dtype")
    def test_loc_iloc_frame_single_dtype(self, data):
        super().test_loc_iloc_frame_single_dtype(data)


class TestGroupby(BaseNumPyTests, base.BaseGroupbyTests):
    @skip_nested
    def test_groupby_extension_apply(self, data_for_grouping, groupby_apply_op):
        # ValueError: Names should be list-like for a MultiIndex
        super().test_groupby_extension_apply(data_for_grouping, groupby_apply_op)


class TestInterface(BaseNumPyTests, base.BaseInterfaceTests):
    @skip_nested
    def test_array_interface(self, data):
        # NumPy array shape inference
        super().test_array_interface(data)


class TestMethods(BaseNumPyTests, base.BaseMethodsTests):
    @pytest.mark.skip(reason="TODO: remove?")
    def test_value_counts(self, all_data, dropna):
        pass

    @pytest.mark.skip(reason="Incorrect expected")
    # We have a bool dtype, so the result is an ExtensionArray
    # but expected is not
    def test_combine_le(self, data_repeated):
        super().test_combine_le(data_repeated)

    @skip_nested
    def test_combine_add(self, data_repeated):
        # Not numeric
        super().test_combine_add(data_repeated)

    @skip_nested
    def test_shift_fill_value(self, data):
        # np.array shape inference. Shift implementation fails.
        super().test_shift_fill_value(data)

    @skip_nested
    @pytest.mark.parametrize("box", [pd.Series, lambda x: x])
    @pytest.mark.parametrize("method", [lambda x: x.unique(), pd.unique])
    def test_unique(self, data, box, method):
        # Fails creating expected
        super().test_unique(data, box, method)

    @skip_nested
    def test_fillna_copy_frame(self, data_missing):
        # The "scalar" for this array isn't a scalar.
        super().test_fillna_copy_frame(data_missing)

    @skip_nested
    def test_fillna_copy_series(self, data_missing):
        # The "scalar" for this array isn't a scalar.
        super().test_fillna_copy_series(data_missing)

    @skip_nested
    def test_hash_pandas_object_works(self, data, as_frame):
        # ndarray of tuples not hashable
        super().test_hash_pandas_object_works(data, as_frame)

    @skip_nested
    def test_searchsorted(self, data_for_sorting, as_series):
        # Test setup fails.
        super().test_searchsorted(data_for_sorting, as_series)

    @skip_nested
    def test_where_series(self, data, na_value, as_frame):
        # Test setup fails.
        super().test_where_series(data, na_value, as_frame)

    @skip_nested
    @pytest.mark.parametrize("repeats", [0, 1, 2, [1, 2, 3]])
    def test_repeat(self, data, repeats, as_series, use_numpy):
        # Fails creating expected
        super().test_repeat(data, repeats, as_series, use_numpy)

    @pytest.mark.xfail(reason="PandasArray.diff may fail on dtype")
    def test_diff(self, data, periods):
        return super().test_diff(data, periods)


@skip_nested
class TestArithmetics(BaseNumPyTests, base.BaseArithmeticOpsTests):
    divmod_exc = None
    series_scalar_exc = None
    frame_scalar_exc = None
    series_array_exc = None

    def test_divmod_series_array(self, data):
        s = pd.Series(data)
        self._check_divmod_op(s, divmod, data, exc=None)

    @pytest.mark.skip("We implement ops")
    def test_error(self, data, all_arithmetic_operators):
        pass

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        super().test_arith_series_with_scalar(data, all_arithmetic_operators)

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        super().test_arith_series_with_array(data, all_arithmetic_operators)


class TestPrinting(BaseNumPyTests, base.BasePrintingTests):
    pass


@skip_nested
class TestNumericReduce(BaseNumPyTests, base.BaseNumericReduceTests):
    def check_reduce(self, s, op_name, skipna):
        result = getattr(s, op_name)(skipna=skipna)
        # avoid coercing int -> float. Just cast to the actual numpy type.
        expected = getattr(s.astype(s.dtype._dtype), op_name)(skipna=skipna)
        tm.assert_almost_equal(result, expected)


@skip_nested
class TestBooleanReduce(BaseNumPyTests, base.BaseBooleanReduceTests):
    pass


class TestMissing(BaseNumPyTests, base.BaseMissingTests):
    @skip_nested
    def test_fillna_scalar(self, data_missing):
        # Non-scalar "scalar" values.
        super().test_fillna_scalar(data_missing)

    @skip_nested
    def test_fillna_series_method(self, data_missing, fillna_method):
        # Non-scalar "scalar" values.
        super().test_fillna_series_method(data_missing, fillna_method)

    @skip_nested
    def test_fillna_series(self, data_missing):
        # Non-scalar "scalar" values.
        super().test_fillna_series(data_missing)

    @skip_nested
    def test_fillna_frame(self, data_missing):
        # Non-scalar "scalar" values.
        super().test_fillna_frame(data_missing)


class TestReshaping(BaseNumPyTests, base.BaseReshapingTests):
    @pytest.mark.skip("Incorrect parent test")
    # not actually a mixed concat, since we concat int and int.
    def test_concat_mixed_dtypes(self, data):
        super().test_concat_mixed_dtypes(data)

    @skip_nested
    def test_merge(self, data, na_value):
        # Fails creating expected
        super().test_merge(data, na_value)

    @skip_nested
    def test_merge_on_extension_array(self, data):
        # Fails creating expected
        super().test_merge_on_extension_array(data)

    @skip_nested
    def test_merge_on_extension_array_duplicates(self, data):
        # Fails creating expected
        super().test_merge_on_extension_array_duplicates(data)

    @skip_nested
    def test_transpose(self, data):
        super().test_transpose(data)


class TestSetitem(BaseNumPyTests, base.BaseSetitemTests):
    @skip_nested
    def test_setitem_scalar_series(self, data, box_in_series):
        # AssertionError
        super().test_setitem_scalar_series(data, box_in_series)

    @skip_nested
    def test_setitem_sequence(self, data, box_in_series):
        # ValueError: shape mismatch: value array of shape (2,1) could not
        # be broadcast to indexing result of shape (2,)
        super().test_setitem_sequence(data, box_in_series)

    @skip_nested
    def test_setitem_sequence_mismatched_length_raises(self, data, as_array):
        # ValueError: PandasArray must be 1-dimensional.
        super().test_setitem_sequence_mismatched_length_raises(data, as_array)

    @skip_nested
    def test_setitem_sequence_broadcasts(self, data, box_in_series):
        # ValueError: cannot set using a list-like indexer with a different
        # length than the value
        super().test_setitem_sequence_broadcasts(data, box_in_series)

    @skip_nested
    def test_setitem_loc_scalar_mixed(self, data):
        # AssertionError
        super().test_setitem_loc_scalar_mixed(data)

    @skip_nested
    def test_setitem_loc_scalar_multiple_homogoneous(self, data):
        # AssertionError
        super().test_setitem_loc_scalar_multiple_homogoneous(data)

    @skip_nested
    def test_setitem_iloc_scalar_mixed(self, data):
        # AssertionError
        super().test_setitem_iloc_scalar_mixed(data)

    @skip_nested
    def test_setitem_iloc_scalar_multiple_homogoneous(self, data):
        # AssertionError
        super().test_setitem_iloc_scalar_multiple_homogoneous(data)

    @skip_nested
    @pytest.mark.parametrize("setter", ["loc", None])
    def test_setitem_mask_broadcast(self, data, setter):
        # ValueError: cannot set using a list-like indexer with a different
        # length than the value
        super().test_setitem_mask_broadcast(data, setter)

    @skip_nested
    def test_setitem_scalar_key_sequence_raise(self, data):
        # Failed: DID NOT RAISE <class 'ValueError'>
        super().test_setitem_scalar_key_sequence_raise(data)

    # TODO: there is some issue with PandasArray, therefore,
    #   skip the setitem test for now, and fix it later (GH 31446)

    @skip_nested
    @pytest.mark.parametrize(
        "mask",
        [
            np.array([True, True, True, False, False]),
            pd.array([True, True, True, False, False], dtype="boolean"),
        ],
        ids=["numpy-array", "boolean-array"],
    )
    def test_setitem_mask(self, data, mask, box_in_series):
        super().test_setitem_mask(data, mask, box_in_series)

    @skip_nested
    def test_setitem_mask_raises(self, data, box_in_series):
        super().test_setitem_mask_raises(data, box_in_series)

    @skip_nested
    def test_setitem_mask_boolean_array_raises(self, data, box_in_series):
        super().test_setitem_mask_boolean_array_raises(data, box_in_series)

    @skip_nested
    @pytest.mark.parametrize(
        "idx",
        [[0, 1, 2], pd.array([0, 1, 2], dtype="Int64"), np.array([0, 1, 2])],
        ids=["list", "integer-array", "numpy-array"],
    )
    def test_setitem_integer_array(self, data, idx, box_in_series):
        super().test_setitem_integer_array(data, idx, box_in_series)

    @skip_nested
    @pytest.mark.parametrize(
        "idx, box_in_series",
        [
            ([0, 1, 2, pd.NA], False),
            pytest.param([0, 1, 2, pd.NA], True, marks=pytest.mark.xfail),
            (pd.array([0, 1, 2, pd.NA], dtype="Int64"), False),
            (pd.array([0, 1, 2, pd.NA], dtype="Int64"), False),
        ],
        ids=["list-False", "list-True", "integer-array-False", "integer-array-True"],
    )
    def test_setitem_integer_with_missing_raises(self, data, idx, box_in_series):
        super().test_setitem_integer_with_missing_raises(data, idx, box_in_series)

    @skip_nested
    def test_setitem_slice(self, data, box_in_series):
        super().test_setitem_slice(data, box_in_series)

    @skip_nested
    def test_setitem_loc_iloc_slice(self, data):
        super().test_setitem_loc_iloc_slice(data)


@skip_nested
class TestParsing(BaseNumPyTests, base.BaseParsingTests):
    pass
