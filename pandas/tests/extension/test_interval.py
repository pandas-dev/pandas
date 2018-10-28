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

import pandas.util.testing as tm
from pandas import Interval
from pandas.core.arrays import IntervalArray
from pandas.core.dtypes.dtypes import IntervalDtype
from pandas.tests.extension import base


def make_data():
    N = 100
    left = np.random.uniform(size=N).cumsum()
    right = left + np.random.uniform(size=N)
    return [Interval(l, r) for l, r in zip(left, right)]


@pytest.fixture
def dtype():
    return IntervalDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    return IntervalArray(make_data())


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return IntervalArray.from_tuples([None, (0, 1)])


@pytest.fixture
def data_for_sorting():
    return IntervalArray.from_tuples([(1, 2), (2, 3), (0, 1)])


@pytest.fixture
def data_missing_for_sorting():
    return IntervalArray.from_tuples([(1, 2), None, (0, 1)])


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping():
    a = (0, 1)
    b = (1, 2)
    c = (2, 3)
    return IntervalArray.from_tuples([b, b, None, None, a, a, b, c])


class BaseInterval(object):
    pass


class TestDtype(BaseInterval, base.BaseDtypeTests):
    pass


class TestCasting(BaseInterval, base.BaseCastingTests):
    pass


class TestConstructors(BaseInterval, base.BaseConstructorsTests):
    pass


class TestGetitem(BaseInterval, base.BaseGetitemTests):
    pass


class TestGrouping(BaseInterval, base.BaseGroupbyTests):
    pass


class TestInterface(BaseInterval, base.BaseInterfaceTests):
    pass


class TestReduce(base.BaseNoReduceTests):
    pass


class TestMethods(BaseInterval, base.BaseMethodsTests):

    @pytest.mark.skip(reason='addition is not defined for intervals')
    def test_combine_add(self, data_repeated):
        pass


class TestMissing(BaseInterval, base.BaseMissingTests):
    # Index.fillna only accepts scalar `value`, so we have to skip all
    # non-scalar fill tests.
    unsupported_fill = pytest.mark.skip("Unsupported fillna option.")

    @unsupported_fill
    def test_fillna_limit_pad(self):
        pass

    @unsupported_fill
    def test_fillna_series_method(self):
        pass

    @unsupported_fill
    def test_fillna_limit_backfill(self):
        pass

    @unsupported_fill
    def test_fillna_series(self):
        pass

    def test_non_scalar_raises(self, data_missing):
        msg = "Got a 'list' instead."
        with tm.assert_raises_regex(TypeError, msg):
            data_missing.fillna([1, 1])


class TestReshaping(BaseInterval, base.BaseReshapingTests):
    pass


class TestSetitem(BaseInterval, base.BaseSetitemTests):
    pass
