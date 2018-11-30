import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray
from pandas.tests.extension import base


# TODO: figure out a way to test non-TZ
@pytest.fixture(params=["US/Central"])
def dtype(request):
    return DatetimeTZDtype(unit="ns", tz=request.param)


@pytest.fixture
def data(dtype):
    data = DatetimeArray(pd.date_range("2000", periods=100, tz=dtype.tz),
                         tz=dtype.tz)
    return data


@pytest.fixture
def data_missing(dtype):
    return DatetimeArray(
        np.array(['NaT', '2000-01-01'], dtype='datetime64[ns]'),
        tz=dtype.tz
    )


@pytest.fixture
def data_for_sorting(dtype):
    a = pd.Timestamp('2000-01-01')
    b = pd.Timestamp('2000-01-02')
    c = pd.Timestamp('2000-01-03')
    return DatetimeArray(np.array([b, c, a], dtype='datetime64[ns]'),
                         tz=dtype.tz)


@pytest.fixture
def data_missing_for_sorting(dtype):
    a = pd.Timestamp('2000-01-01')
    b = pd.Timestamp('2000-01-02')
    return DatetimeArray(np.array([b, 'NaT', a], dtype='datetime64[ns]'),
                         tz=dtype.tz)


@pytest.fixture
def data_for_grouping(dtype):
    """
        Expected to be like [B, B, NA, NA, A, A, B, C]

        Where A < B < C and NA is missing
    """
    a = pd.Timestamp('2000-01-01')
    b = pd.Timestamp('2000-01-02')
    c = pd.Timestamp('2000-01-03')
    na = 'NaT'
    return DatetimeArray(np.array([b, b, na, na, a, a, b, c],
                                  dtype='datetime64[ns]'),
                         tz=dtype.tz)


@pytest.fixture
def na_cmp():
    def cmp(a, b):
        return a is pd.NaT and a is b
    return cmp


@pytest.fixture
def na_value():
    return pd.NaT


# ----------------------------------------------------------------------------
class BaseDatetimeTests(object):
    pass


# ----------------------------------------------------------------------------
# Tests
class TestDatetimeDtype(BaseDatetimeTests, base.BaseDtypeTests):
    pass


class TestConstructors(BaseDatetimeTests, base.BaseConstructorsTests):
    pass


class TestGetitem(BaseDatetimeTests, base.BaseGetitemTests):
    pass


class TestMethods(BaseDatetimeTests, base.BaseMethodsTests):
    @pytest.mark.xfail(reason='GH-22843', strict=True)
    def test_value_counts(self, all_data, dropna):
        # fails without .value_counts
        return super().test_value_counts(all_data, dropna)

    def test_apply_simple_series(self, data):
        if data.tz:
            # fails without .map
            raise pytest.xfail('GH-23179')
        super().test_apply_simple_series(data)

    def test_combine_add(self, data_repeated):
        # Timestamp.__add__(Timestamp) not defined
        pass


class TestInterface(BaseDatetimeTests, base.BaseInterfaceTests):

    @pytest.mark.xfail(reason="Figure out np.array(tz_aware)", strict=False)
    def test_array_interface(self, data):
        # override, because np.array(data)[0] != data[0]
        # since numpy datetime64ns scalars don't compare equal
        # to timestmap objects.
        result = np.array(data)
        # even this fails, since arary(data) is *not* tz aware, and
        # we don't compare tz-aware and tz-naive.
        # this could work if array(data) was object-dtype with timestamps.
        assert data[0] == result[0]


class TestArithmeticOps(BaseDatetimeTests, base.BaseArithmeticOpsTests):
    implements = {'__sub__', '__rsub__'}

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # TODO: move this to the base class?
        # It's duplicated between Period and Datetime now
        if all_arithmetic_operators in self.implements:
            s = pd.Series(data)
            self.check_opname(s, all_arithmetic_operators, s.iloc[0],
                              exc=None)
        else:
            # ... but not the rest.
            super(TestArithmeticOps, self).test_arith_series_with_scalar(
                data, all_arithmetic_operators
            )

    def test_add_series_with_extension_array(self, data):
        # Datetime + Datetime not implemented
        s = pd.Series(data)
        msg = 'cannot add DatetimeArray(Mixin)? and DatetimeArray(Mixin)?'
        with pytest.raises(TypeError, match=msg):
            s + data

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

    def test_error(self, data, all_arithmetic_operators):
        pass

    @pytest.mark.xfail(reason="Not Implemented", strict=False)
    def test_direct_arith_with_series_returns_not_implemented(self, data):
        # Right now, we have trouble with this. Returning NotImplemented
        # fails other tests like
        # tests/arithmetic/test_datetime64::TestTimestampSeriesArithmetic::
        # test_dt64_seris_add_intlike
        return super(
            TestArithmeticOps,
            self
        ).test_direct_arith_with_series_returns_not_implemented(data)


class TestCasting(BaseDatetimeTests, base.BaseCastingTests):
    pass


class TestComparisonOps(BaseDatetimeTests, base.BaseComparisonOpsTests):

    def _compare_other(self, s, data, op_name, other):
        # the base test is not appropriate for us. We raise on comparison
        # with (some) integers, depending on the value.
        pass

    @pytest.mark.xfail(reason="Not Implemented", strict=False)
    def test_direct_arith_with_series_returns_not_implemented(self, data):
        return super(
            TestComparisonOps,
            self
        ).test_direct_arith_with_series_returns_not_implemented(data)


class TestMissing(BaseDatetimeTests, base.BaseMissingTests):
    pass


class TestReshaping(BaseDatetimeTests, base.BaseReshapingTests):

    @pytest.mark.skip(reason="We have DatetimeTZBlock")
    def test_concat(self, data, in_frame):
        pass

    @pytest.mark.xfail(reason="GH-23816", strict=True)
    def test_concat_mixed_dtypes(self, data):
        # concat(Series[datetimetz], Series[category]) uses a
        # plain np.array(values) on the DatetimeArray, which
        # drops the tz.
        super(TestReshaping, self).test_concat_mixed_dtypes(data)

    @pytest.mark.xfail(reason="GH-13287", strict=True)
    def test_unstack(self, data, index, obj):
        # This fails creating the expected.
        # Ahh this is going to always xfail, since we don't have the
        # fixtures...
        return super(TestReshaping, self).test_unstack(data, index, obj)


class TestSetitem(BaseDatetimeTests, base.BaseSetitemTests):
    pass


class TestGroupby(BaseDatetimeTests, base.BaseGroupbyTests):
    pass
