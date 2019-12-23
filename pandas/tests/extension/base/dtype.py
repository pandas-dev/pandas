import warnings

import numpy as np
import pytest

import pandas as pd

from .base import BaseExtensionTests


class BaseDtypeTests(BaseExtensionTests):
    """Base class for ExtensionDtype classes"""

    def test_name(self, dtype):
        assert isinstance(dtype.name, str)

    def test_kind(self, dtype):
        valid = set("biufcmMOSUV")
        if dtype.kind is not None:
            assert dtype.kind in valid

    def test_construct_from_string_own_name(self, dtype):
        result = dtype.construct_from_string(dtype.name)
        assert type(result) is type(dtype)

        # check OK as classmethod
        result = type(dtype).construct_from_string(dtype.name)
        assert type(result) is type(dtype)

    def test_is_dtype_from_name(self, dtype):
        result = type(dtype).is_dtype(dtype.name)
        assert result is True

    def test_is_dtype_unboxes_dtype(self, data, dtype):
        assert dtype.is_dtype(data) is True

    def test_is_dtype_from_self(self, dtype):
        result = type(dtype).is_dtype(dtype)
        assert result is True

    def test_is_not_string_type(self, dtype):
        return not pd.api.types.is_string_dtype(dtype)

    def test_is_not_object_type(self, dtype):
        return not pd.api.types.is_object_dtype(dtype)

    def test_eq_with_str(self, dtype):
        assert dtype == dtype.name
        assert dtype != dtype.name + "-suffix"

    def test_eq_with_numpy_object(self, dtype):
        assert dtype != np.dtype("object")

    def test_eq_with_self(self, dtype):
        assert dtype == dtype
        assert dtype != object()

    def test_array_type(self, data, dtype):
        assert dtype.construct_array_type() is type(data)

    def test_check_dtype(self, data):
        dtype = data.dtype

        # check equivalency for using .dtypes
        df = pd.DataFrame(
            {"A": pd.Series(data, dtype=dtype), "B": data, "C": "foo", "D": 1}
        )

        # np.dtype('int64') == 'Int64' == 'int64'
        # so can't distinguish
        if dtype.name == "Int64":
            expected = pd.Series([True, True, False, True], index=list("ABCD"))
        else:
            expected = pd.Series([True, True, False, False], index=list("ABCD"))

        # XXX: This should probably be *fixed* not ignored.
        # See libops.scalar_compare
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = df.dtypes == str(dtype)

        self.assert_series_equal(result, expected)

        expected = pd.Series([True, True, False, False], index=list("ABCD"))
        result = df.dtypes.apply(str) == str(dtype)
        self.assert_series_equal(result, expected)

    def test_hashable(self, dtype):
        hash(dtype)  # no error

    def test_str(self, dtype):
        assert str(dtype) == dtype.name

    def test_eq(self, dtype):
        assert dtype == dtype.name
        assert dtype != "anonther_type"

    def test_construct_from_string(self, dtype):
        dtype_instance = type(dtype).construct_from_string(dtype.name)
        assert isinstance(dtype_instance, type(dtype))

    def test_construct_from_string_another_type_raises(self, dtype):
        msg = f"Cannot construct a '{type(dtype).__name__}' from 'another_type'"
        with pytest.raises(TypeError, match=msg):
            type(dtype).construct_from_string("another_type")
