import pytest
import numpy as np

import pandas as pd
from pandas._libs.tslib import iNaT
from pandas.tests.extension import base
from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.arrays import PeriodArray


@pytest.fixture
def dtype():
    return PeriodDtype(freq='D')


@pytest.fixture
def data(dtype):
    return PeriodArray(np.arange(1970, 2070), freq=dtype.freq)


@pytest.fixture
def data_for_sorting(dtype):
    return PeriodArray([2018, 2019, 2017], freq=dtype.freq)


@pytest.fixture
def data_missing(dtype):
    return PeriodArray([iNaT, 2017], freq=dtype.freq)


@pytest.fixture
def data_missing_for_sorting(dtype):
    return PeriodArray([2018, iNaT, 2017], freq=dtype.freq)


@pytest.fixture
def data_for_grouping(dtype):
    B = 2018
    NA = iNaT
    A = 2017
    C = 2019
    return PeriodArray([B, B, NA, NA, A, A, B, C], freq=dtype.freq)


@pytest.fixture
def na_value():
    return pd.NaT


class BasePeriodTests(object):
    pass


class TestPeriodDtype(BasePeriodTests, base.BaseDtypeTests):
    pass


class TestConstructors(BasePeriodTests, base.BaseConstructorsTests):
    pass


class TestGetitem(BasePeriodTests, base.BaseGetitemTests):
    pass


class TestMethods(BasePeriodTests, base.BaseMethodsTests):

    def test_combine_add(self, data_repeated):
        # Period + Period is not defined.
        pass


class TestInterface(BasePeriodTests, base.BaseInterfaceTests):

    def test_no_values_attribute(self, data):
        # We have a values attribute.
        pass


class TestArithmeticOps(BasePeriodTests, base.BaseArithmeticOpsTests):

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # we implement substitution...
        op_name = all_arithmetic_operators
        if op_name in ('__sub__', '__rsub__'):
            s = pd.Series(data)
            self.check_opname(s, op_name, s.iloc[0], exc=None)
        else:
            # ... but not the rest.
            super().test_arith_series_with_scalar(data,
                                                  all_arithmetic_operators)

    def test_error(self):
        pass


class TestCasting(BasePeriodTests, base.BaseCastingTests):
    pass


class TestComparisonOps(BasePeriodTests, base.BaseComparisonOpsTests):

    def _compare_other(self):
        # the base test is not appropriate for us. We raise on comparison
        # with (some) integers, depending on the value.
        pass


class TestMissing(BasePeriodTests, base.BaseMissingTests):
    pass



class TestReshaping(BasePeriodTests, base.BaseReshapingTests):
    pass


class TestSetitem(BasePeriodTests, base.BaseSetitemTests):
    pass


class TestGroupby(BasePeriodTests, base.BaseGroupbyTests):
    pass

