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

import pandas as pd
from pandas.core.arrays.floating import (
    Float32Dtype,
    Float64Dtype,
)
from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
)
from pandas.tests.extension import (
    base,
    masked_shared,
)

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in divide:RuntimeWarning"
    ),
    pytest.mark.filterwarnings("ignore:Mean of empty slice:RuntimeWarning"),
    # overflow only relevant for Floating dtype cases cases
    pytest.mark.filterwarnings("ignore:overflow encountered in reduce:RuntimeWarning"),
]


def make_data():
    return list(range(1, 9)) + [pd.NA] + list(range(10, 98)) + [pd.NA] + [99, 100]


def make_float_data():
    return (
        list(np.arange(0.1, 0.9, 0.1))
        + [pd.NA]
        + list(np.arange(1, 9.8, 0.1))
        + [pd.NA]
        + [9.9, 10.0]
    )


@pytest.fixture(
    params=[
        Int8Dtype,
        Int16Dtype,
        Int32Dtype,
        Int64Dtype,
        UInt8Dtype,
        UInt16Dtype,
        UInt32Dtype,
        UInt64Dtype,
        Float32Dtype,
        Float64Dtype,
    ]
)
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    if dtype.kind == "f":
        data = make_float_data()
    else:
        data = make_data()
    return pd.array(data, dtype=dtype)


@pytest.fixture
def data_for_twos(dtype):
    return pd.array(np.ones(100) * 2, dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    if dtype.kind == "f":
        return pd.array([pd.NA, 0.1], dtype=dtype)
    return pd.array([pd.NA, 1], dtype=dtype)


@pytest.fixture
def data_for_sorting(dtype):
    if dtype.kind == "f":
        return pd.array([0.1, 0.2, 0.0], dtype=dtype)
    return pd.array([1, 2, 0], dtype=dtype)


@pytest.fixture
def data_missing_for_sorting(dtype):
    if dtype.kind == "f":
        return pd.array([0.1, pd.NA, 0.0], dtype=dtype)
    return pd.array([1, pd.NA, 0], dtype=dtype)


@pytest.fixture
def na_cmp():
    # we are pd.NA
    return lambda x, y: x is pd.NA and y is pd.NA


@pytest.fixture
def na_value():
    return pd.NA


@pytest.fixture
def data_for_grouping(dtype):
    if dtype.kind == "f":
        b = 0.1
        a = 0.0
        c = 0.2
    else:
        b = 1
        a = 0
        c = 2
    na = pd.NA
    return pd.array([b, b, na, na, a, a, b, c], dtype=dtype)


class TestDtype(base.BaseDtypeTests):
    pass


class TestArithmeticOps(masked_shared.Arithmetic):
    pass


class TestComparisonOps(masked_shared.Comparison):
    pass


class TestInterface(base.BaseInterfaceTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass

    # for test_concat_mixed_dtypes test
    # concat of an Integer and Int coerces to object dtype
    # TODO(jreback) once integrated this would


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    pass


class TestIndex(base.BaseIndexTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    _combine_le_expected_dtype = object  # TODO: can we make this boolean?


class TestCasting(base.BaseCastingTests):
    pass


class TestGroupby(base.BaseGroupbyTests):
    pass


class TestNumericReduce(masked_shared.NumericReduce):
    pass


@pytest.mark.skip(reason="Tested in tests/reductions/test_reductions.py")
class TestBooleanReduce(base.BaseBooleanReduceTests):
    pass


class TestAccumulation(masked_shared.Accumulation):
    pass


class TestPrinting(base.BasePrintingTests):
    pass


class TestParsing(base.BaseParsingTests):
    pass


class Test2DCompat(base.Dim2CompatTests):
    pass
