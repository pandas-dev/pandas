"""
This file contains a minimal set of tests for compliance with the extension
array interface test suite, and should contain no other tests.
The test suite for the full functionality of the array is located in
`pandas/tests/arrays/`.

The tests in this file are inherited from the BaseExtensionTests, and only
minimal tweaks should be applied to get the tests passing (by overwriting a
parent method).

Additional tests should either be added to one of the BaseExtensionTests
classes (if they are relevant for the extension interface for all dtypes), or
be added to the array-specific tests in `pandas/tests/arrays/`.

"""
import numpy as np
import pytest

from pandas.compat import (
    IS64,
    is_platform_windows,
)

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.boolean import BooleanDtype
from pandas.tests.extension import base

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in divide:RuntimeWarning"
    ),
    pytest.mark.filterwarnings("ignore:Mean of empty slice:RuntimeWarning"),
]


def make_data():
    return [True, False] * 4 + [np.nan] + [True, False] * 44 + [np.nan] + [True, False]


@pytest.fixture
def dtype():
    return BooleanDtype()


@pytest.fixture
def data(dtype):
    return pd.array(make_data(), dtype=dtype)


@pytest.fixture
def data_for_twos(dtype):
    return pd.array(np.ones(100), dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    return pd.array([np.nan, True], dtype=dtype)


@pytest.fixture
def data_for_sorting(dtype):
    return pd.array([True, True, False], dtype=dtype)


@pytest.fixture
def data_missing_for_sorting(dtype):
    return pd.array([True, np.nan, False], dtype=dtype)


@pytest.fixture
def na_cmp():
    # we are pd.NA
    return lambda x, y: x is pd.NA and y is pd.NA


@pytest.fixture
def na_value():
    return pd.NA


@pytest.fixture
def data_for_grouping(dtype):
    b = True
    a = False
    c = b
    na = np.nan
    return pd.array([b, b, na, na, a, a, b, c], dtype=dtype)


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    pass


class TestIndex(base.BaseIndexTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


class TestArithmeticOps(base.BaseArithmeticOpsTests):
    implements = {"__sub__", "__rsub__"}

    def check_opname(self, s, op_name, other, exc=None):
        # overwriting to indicate ops don't raise an error
        exc = None
        if op_name.strip("_").lstrip("r") in ["pow", "truediv", "floordiv"]:
            # match behavior with non-masked bool dtype
            exc = NotImplementedError
        super().check_opname(s, op_name, other, exc=exc)

    def _check_op(self, obj, op, other, op_name, exc=NotImplementedError):
        if exc is None:
            if op_name in self.implements:
                msg = r"numpy boolean subtract"
                with pytest.raises(TypeError, match=msg):
                    op(obj, other)
                return

            result = op(obj, other)
            expected = self._combine(obj, other, op)

            if op_name in (
                "__floordiv__",
                "__rfloordiv__",
                "__pow__",
                "__rpow__",
                "__mod__",
                "__rmod__",
            ):
                # combine keeps boolean type
                expected = expected.astype("Int8")
            elif op_name in ("__truediv__", "__rtruediv__"):
                # combine with bools does not generate the correct result
                #  (numpy behaviour for div is to regard the bools as numeric)
                expected = self._combine(obj.astype(float), other, op)
                expected = expected.astype("Float64")
            if op_name == "__rpow__":
                # for rpow, combine does not propagate NaN
                expected[result.isna()] = np.nan
            self.assert_equal(result, expected)
        else:
            with pytest.raises(exc):
                op(obj, other)

    @pytest.mark.xfail(
        reason="Inconsistency between floordiv and divmod; we raise for floordiv "
        "but not for divmod. This matches what we do for non-masked bool dtype."
    )
    def test_divmod_series_array(self, data, data_for_twos):
        super().test_divmod_series_array(data, data_for_twos)

    @pytest.mark.xfail(
        reason="Inconsistency between floordiv and divmod; we raise for floordiv "
        "but not for divmod. This matches what we do for non-masked bool dtype."
    )
    def test_divmod(self, data):
        super().test_divmod(data)


class TestComparisonOps(base.BaseComparisonOpsTests):
    def check_opname(self, s, op_name, other, exc=None):
        # overwriting to indicate ops don't raise an error
        super().check_opname(s, op_name, other, exc=None)


class TestReshaping(base.BaseReshapingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    _combine_le_expected_dtype = "boolean"


class TestCasting(base.BaseCastingTests):
    pass


class TestGroupby(base.BaseGroupbyTests):
    """
    Groupby-specific tests are overridden because boolean only has 2
    unique values, base tests uses 3 groups.
    """

    @pytest.mark.parametrize("min_count", [0, 10])
    def test_groupby_sum_mincount(self, data_for_grouping, min_count):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": data_for_grouping[:-1]})
        result = df.groupby("A").sum(min_count=min_count)
        if min_count == 0:
            expected = pd.DataFrame(
                {"B": pd.array([3, 0, 0], dtype="Int64")},
                index=pd.Index([1, 2, 3], name="A"),
            )
            tm.assert_frame_equal(result, expected)
        else:
            expected = pd.DataFrame(
                {"B": pd.array([pd.NA] * 3, dtype="Int64")},
                index=pd.Index([1, 2, 3], name="A"),
            )
            tm.assert_frame_equal(result, expected)


class TestNumericReduce(base.BaseNumericReduceTests):
    def check_reduce(self, s, op_name, skipna):
        if op_name == "count":
            result = getattr(s, op_name)()
            expected = getattr(s.astype("float64"), op_name)()
        else:
            result = getattr(s, op_name)(skipna=skipna)
            expected = getattr(s.astype("float64"), op_name)(skipna=skipna)
        # override parent function to cast to bool for min/max
        if np.isnan(expected):
            expected = pd.NA
        elif op_name in ("min", "max"):
            expected = bool(expected)
        tm.assert_almost_equal(result, expected)

    def check_reduce_frame(self, ser: pd.Series, op_name: str, skipna: bool):
        arr = ser.array

        if op_name in ["count", "kurt", "sem"]:
            assert not hasattr(arr, op_name)
            return

        if op_name in ["mean", "median", "var", "std", "skew"]:
            cmp_dtype = "Float64"
        elif op_name in ["min", "max"]:
            cmp_dtype = "boolean"
        elif op_name in ["sum", "prod"]:
            is_windows_or_32bit = is_platform_windows() or not IS64
            cmp_dtype = "Int32" if is_windows_or_32bit else "Int64"
        else:
            raise TypeError("not supposed to reach this")

        result = arr._reduce(op_name, skipna=skipna, keepdims=True)
        if not skipna and ser.isna().any():
            expected = pd.array([pd.NA], dtype=cmp_dtype)
        else:
            exp_value = getattr(ser.dropna().astype(cmp_dtype), op_name)()
            expected = pd.array([exp_value], dtype=cmp_dtype)
        tm.assert_extension_array_equal(result, expected)


class TestBooleanReduce(base.BaseBooleanReduceTests):
    pass


class TestPrinting(base.BasePrintingTests):
    pass


class TestUnaryOps(base.BaseUnaryOpsTests):
    pass


class TestAccumulation(base.BaseAccumulateTests):
    def check_accumulate(self, s, op_name, skipna):
        length = 64
        if not IS64 or is_platform_windows():
            if not s.dtype.itemsize == 8:
                length = 32

        result = getattr(s, op_name)(skipna=skipna)
        expected = getattr(pd.Series(s.astype("float64")), op_name)(skipna=skipna)
        if op_name not in ("cummin", "cummax"):
            expected = expected.astype(f"Int{length}")
        else:
            expected = expected.astype("boolean")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series_raises(self, data, all_numeric_accumulations, skipna):
        pass


class TestParsing(base.BaseParsingTests):
    pass


class Test2DCompat(base.Dim2CompatTests):
    pass
