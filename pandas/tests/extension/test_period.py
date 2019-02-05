import numpy as np
import pytest

from pandas._libs.tslib import iNaT

from pandas.core.dtypes.dtypes import PeriodDtype

import pandas as pd
from pandas.core.arrays import PeriodArray
from pandas.tests.extension import base


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

    pass


class TestArithmeticOps(BasePeriodTests, base.BaseArithmeticOpsTests):
    implements = {'__sub__', '__rsub__'}

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # we implement substitution...
        if all_arithmetic_operators in self.implements:
            s = pd.Series(data)
            self.check_opname(s, all_arithmetic_operators, s.iloc[0],
                              exc=None)
        else:
            # ... but not the rest.
            super(TestArithmeticOps, self).test_arith_series_with_scalar(
                data, all_arithmetic_operators
            )

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        if all_arithmetic_operators in self.implements:
            s = pd.Series(data)
            self.check_opname(s, all_arithmetic_operators, s.iloc[0],
                              exc=None)
        else:
            # ... but not the rest.
            super(TestArithmeticOps, self).test_arith_series_with_scalar(
                data, all_arithmetic_operators
            )

    def _check_divmod_op(self, s, op, other, exc=NotImplementedError):
        super(TestArithmeticOps, self)._check_divmod_op(
            s, op, other, exc=TypeError
        )

    def test_add_series_with_extension_array(self, data):
        # we don't implement + for Period
        s = pd.Series(data)
        msg = (r"unsupported operand type\(s\) for \+: "
               r"\'PeriodArray\' and \'PeriodArray\'")
        with pytest.raises(TypeError, match=msg):
            s + data

    def test_error(self):
        pass

    def test_direct_arith_with_series_returns_not_implemented(self, data):
        # Override to use __sub__ instead of __add__
        other = pd.Series(data)
        result = data.__sub__(other)
        assert result is NotImplemented


class TestCasting(BasePeriodTests, base.BaseCastingTests):
    pass


class TestComparisonOps(BasePeriodTests, base.BaseComparisonOpsTests):

    def _compare_other(self, s, data, op_name, other):
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


class TestPrinting(BasePeriodTests, base.BasePrintingTests):
    pass


class TestParsing(BasePeriodTests, base.BaseParsingTests):
    @pytest.mark.parametrize('engine', ['c', 'python'])
    def test_EA_types(self, engine, data):
        expected_msg = r'.*must implement _from_sequence_of_strings.*'
        with pytest.raises(NotImplementedError, match=expected_msg):
            super(TestParsing, self).test_EA_types(engine, data)
