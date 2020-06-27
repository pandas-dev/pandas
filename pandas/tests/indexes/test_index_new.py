"""
Tests for the Index constructor conducting inference.
"""
import numpy as np
import pytest

from pandas.core.dtypes.common import is_unsigned_integer_dtype

from pandas import (
    NA,
    CategoricalIndex,
    DatetimeIndex,
    Index,
    Int64Index,
    MultiIndex,
    NaT,
    PeriodIndex,
    Series,
    TimedeltaIndex,
    Timestamp,
    UInt64Index,
    period_range,
)
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

    def test_constructor_categorical_to_object(self):
        # GH#32167 Categorical data and dtype=object should return object-dtype
        ci = CategoricalIndex(range(5))
        result = Index(ci, dtype=object)
        assert not isinstance(result, CategoricalIndex)

    def test_constructor_infer_periodindex(self):
        xp = period_range("2012-1-1", freq="M", periods=3)
        rs = Index(xp)
        tm.assert_index_equal(rs, xp)
        assert isinstance(rs, PeriodIndex)

    @pytest.mark.parametrize("pos", [0, 1])
    @pytest.mark.parametrize(
        "klass,dtype,ctor",
        [
            (DatetimeIndex, "datetime64[ns]", np.datetime64("nat")),
            (TimedeltaIndex, "timedelta64[ns]", np.timedelta64("nat")),
        ],
    )
    def test_constructor_infer_nat_dt_like(
        self, pos, klass, dtype, ctor, nulls_fixture, request
    ):
        expected = klass([NaT, NaT])
        assert expected.dtype == dtype
        data = [ctor]
        data.insert(pos, nulls_fixture)

        if nulls_fixture is NA:
            expected = Index([NA, NaT])
            mark = pytest.mark.xfail(reason="Broken with np.NaT ctor; see GH 31884")
            request.node.add_marker(mark)

        result = Index(data)
        tm.assert_index_equal(result, expected)

        result = Index(np.array(data, dtype=object))
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("swap_objs", [True, False])
    def test_constructor_mixed_nat_objs_infers_object(self, swap_objs):
        # mixed np.datetime64/timedelta64 nat results in object
        data = [np.datetime64("nat"), np.timedelta64("nat")]
        if swap_objs:
            data = data[::-1]

        expected = Index(data, dtype=object)
        tm.assert_index_equal(Index(data), expected)
        tm.assert_index_equal(Index(np.array(data, dtype=object)), expected)


class TestIndexConstructorUnwrapping:
    # Test passing different arraylike values to pd.Index

    @pytest.mark.parametrize("klass", [Index, DatetimeIndex])
    def test_constructor_from_series_dt64(self, klass):
        stamps = [Timestamp("20110101"), Timestamp("20120101"), Timestamp("20130101")]
        expected = DatetimeIndex(stamps)
        ser = Series(stamps)
        result = klass(ser)
        tm.assert_index_equal(result, expected)
