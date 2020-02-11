"""
Tests for the Index constructor conducting inference.
"""
import numpy as np
import pytest

from pandas.core.dtypes.common import is_unsigned_integer_dtype

from pandas import Index, Int64Index, MultiIndex, UInt64Index
import pandas._testing as tm


class TestIndexConstructorInference:
    @pytest.mark.parametrize("na_value", [None, np.nan])
    @pytest.mark.parametrize("vtype", [list, tuple, iter])
    def test_construction_list_tuples_nan(self, na_value, vtype):
        # GH#18505 : valid tuples containing NaN
        values = [(1, "two"), (3.0, na_value)]
        result = Index(vtype(values))
        expected = MultiIndex.from_tuples(values)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [int, "int64", "int32", "int16", "int8", "uint64", "uint32", "uint16", "uint8"],
    )
    def test_constructor_int_dtype_float(self, dtype):
        # GH#18400
        if is_unsigned_integer_dtype(dtype):
            index_type = UInt64Index
        else:
            index_type = Int64Index

        expected = index_type([0, 1, 2, 3])
        result = Index([0.0, 1.0, 2.0, 3.0], dtype=dtype)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("cast_index", [True, False])
    @pytest.mark.parametrize(
        "vals", [[True, False, True], np.array([True, False, True], dtype=bool)]
    )
    def test_constructor_dtypes_to_object(self, cast_index, vals):
        if cast_index:
            index = Index(vals, dtype=bool)
        else:
            index = Index(vals)

        assert type(index) is Index
        assert index.dtype == object
