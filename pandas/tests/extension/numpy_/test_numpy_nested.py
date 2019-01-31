"""
Tests for PandasArray with nested data. Users typically won't create
these objects via `pd.array`, but they can show up through `.array`
on a Series with nested data.

We partition these tests into their own file, as many of the base
tests fail, as they aren't appropriate for nested data. It is easier
to have a seperate file with its own data generating fixtures, than
trying to skip based upon the value of a fixture.
"""
import pytest

import pandas as pd
from pandas.core.arrays.numpy_ import PandasArray, PandasDtype

from .. import base

# For NumPy <1.16, np.array([np.nan, (1,)]) raises
# ValueError: setting an array element with a sequence.
np = pytest.importorskip('numpy', minversion='1.16.0')


@pytest.fixture
def dtype():
    return PandasDtype(np.dtype('object'))


@pytest.fixture
def data(allow_in_pandas, dtype):
    return pd.Series([(i,) for i in range(100)]).array


@pytest.fixture
def data_missing(allow_in_pandas):
    return PandasArray(np.array([np.nan, (1,)]))


@pytest.fixture
def data_for_sorting(allow_in_pandas):
    """Length-3 array with a known sort order.

    This should be three items [B, C, A] with
    A < B < C
    """
    # Use an empty tuple for first element, then remove,
    # to disable np.array's shape inference.
    return PandasArray(
        np.array([(), (2,), (3,), (1,)])[1:]
    )


@pytest.fixture
def data_missing_for_sorting(allow_in_pandas):
    """Length-3 array with a known sort order.

    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    return PandasArray(
        np.array([(1,), np.nan, (0,)])
    )


@pytest.fixture
def data_for_grouping(allow_in_pandas):
    """Data for factorization, grouping, and unique tests.

    Expected to be like [B, B, NA, NA, A, A, B, C]

    Where A < B < C and NA is missing
    """
    a, b, c = (1,), (2,), (3,)
    return PandasArray(np.array(
        [b, b, np.nan, np.nan, a, a, b, c]
    ))


skip_nested = pytest.mark.skip(reason="Skipping for nested PandasArray")


class BaseNumPyTests(object):
    pass


class TestCasting(BaseNumPyTests, base.BaseCastingTests):

    @skip_nested
    def test_astype_str(self, data):
        pass


class TestConstructors(BaseNumPyTests, base.BaseConstructorsTests):
    @pytest.mark.skip(reason="We don't register our dtype")
    # We don't want to register. This test should probably be split in two.
    def test_from_dtype(self, data):
        pass

    @skip_nested
    def test_array_from_scalars(self, data):
        pass


class TestDtype(BaseNumPyTests, base.BaseDtypeTests):

    @pytest.mark.skip(reason="Incorrect expected.")
    # we unsurprisingly clash with a NumPy name.
    def test_check_dtype(self, data):
        pass


class TestGetitem(BaseNumPyTests, base.BaseGetitemTests):

    @skip_nested
    def test_getitem_scalar(self, data):
        pass

    @skip_nested
    def test_take_series(self, data):
        pass


class TestGroupby(BaseNumPyTests, base.BaseGroupbyTests):
    @skip_nested
    def test_groupby_extension_apply(self, data_for_grouping, op):
        pass


class TestInterface(BaseNumPyTests, base.BaseInterfaceTests):
    @skip_nested
    def test_array_interface(self, data):
        # NumPy array shape inference
        pass


class TestMethods(BaseNumPyTests, base.BaseMethodsTests):

    @pytest.mark.skip(reason="TODO: remove?")
    def test_value_counts(self, all_data, dropna):
        pass

    @pytest.mark.skip(reason="Incorrect expected")
    # We have a bool dtype, so the result is an ExtensionArray
    # but expected is not
    def test_combine_le(self, data_repeated):
        super(TestMethods, self).test_combine_le(data_repeated)

    @skip_nested
    def test_combine_add(self, data_repeated):
        # Not numeric
        pass

    @skip_nested
    def test_shift_fill_value(self, data):
        # np.array shape inference. Shift implementation fails.
        super().test_shift_fill_value(data)

    @skip_nested
    def test_unique(self, data, box, method):
        # Fails creating expected
        pass

    @skip_nested
    def test_fillna_copy_frame(self, data_missing):
        # The "scalar" for this array isn't a scalar.
        pass

    @skip_nested
    def test_fillna_copy_series(self, data_missing):
        # The "scalar" for this array isn't a scalar.
        pass

    @skip_nested
    def test_hash_pandas_object_works(self, data, as_frame):
        # ndarray of tuples not hashable
        pass

    @skip_nested
    def test_searchsorted(self, data_for_sorting, as_series):
        # Test setup fails.
        pass

    @skip_nested
    def test_where_series(self, data, na_value, as_frame):
        # Test setup fails.
        pass

    @skip_nested
    def test_repeat(self, data, repeats, as_series, use_numpy):
        # Fails creating expected
        pass


class TestPrinting(BaseNumPyTests, base.BasePrintingTests):
    pass


class TestMissing(BaseNumPyTests, base.BaseMissingTests):

    @skip_nested
    def test_fillna_scalar(self, data_missing):
        # Non-scalar "scalar" values.
        pass

    @skip_nested
    def test_fillna_series_method(self, data_missing, method):
        # Non-scalar "scalar" values.
        pass

    @skip_nested
    def test_fillna_series(self, data_missing):
        # Non-scalar "scalar" values.
        pass

    @skip_nested
    def test_fillna_frame(self, data_missing):
        # Non-scalar "scalar" values.
        pass


class TestReshaping(BaseNumPyTests, base.BaseReshapingTests):

    @pytest.mark.skip("Incorrect parent test")
    # not actually a mixed concat, since we concat int and int.
    def test_concat_mixed_dtypes(self, data):
        super(TestReshaping, self).test_concat_mixed_dtypes(data)

    @skip_nested
    def test_merge(self, data, na_value):
        # Fails creating expected
        pass

    @skip_nested
    def test_merge_on_extension_array(self, data):
        # Fails creating expected
        pass

    @skip_nested
    def test_merge_on_extension_array_duplicates(self, data):
        # Fails creating expected
        pass


class TestSetitem(BaseNumPyTests, base.BaseSetitemTests):

    @skip_nested
    def test_setitem_scalar_series(self, data, box_in_series):
        pass

    @skip_nested
    def test_setitem_sequence(self, data, box_in_series):
        pass

    @skip_nested
    def test_setitem_sequence_mismatched_length_raises(self, data, as_array):
        pass

    @skip_nested
    def test_setitem_sequence_broadcasts(self, data, box_in_series):
        pass

    @skip_nested
    def test_setitem_loc_scalar_mixed(self, data):
        pass

    @skip_nested
    def test_setitem_loc_scalar_multiple_homogoneous(self, data):
        pass

    @skip_nested
    def test_setitem_iloc_scalar_mixed(self, data):
        pass

    @skip_nested
    def test_setitem_iloc_scalar_multiple_homogoneous(self, data):
        pass

    @skip_nested
    def test_setitem_mask_broadcast(self, data, setter):
        pass

    @skip_nested
    def test_setitem_scalar_key_sequence_raise(self, data):
        pass


# Skip Arithmetics, NumericReduce, BooleanReduce, Parsing
