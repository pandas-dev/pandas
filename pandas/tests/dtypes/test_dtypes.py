# -*- coding: utf-8 -*-
import pytest

from itertools import product

import numpy as np
import pandas as pd
from pandas import Series, Categorical, IntervalIndex, date_range

from pandas.core.dtypes.dtypes import (
    DatetimeTZDtype, PeriodDtype,
    IntervalDtype, CategoricalDtype)
from pandas.core.dtypes.common import (
    is_categorical_dtype, is_categorical,
    is_datetime64tz_dtype, is_datetimetz,
    is_period_dtype, is_period,
    is_dtype_equal, is_datetime64_ns_dtype,
    is_datetime64_dtype, is_interval_dtype,
    is_datetime64_any_dtype, is_string_dtype,
    _coerce_to_dtype)
import pandas.util.testing as tm


class Base(object):

    def test_hash(self):
        hash(self.dtype)

    def test_equality_invalid(self):
        assert not self.dtype == 'foo'
        assert not is_dtype_equal(self.dtype, np.int64)

    def test_numpy_informed(self):
        pytest.raises(TypeError, np.dtype, self.dtype)

        assert not self.dtype == np.str_
        assert not np.str_ == self.dtype

    def test_pickle(self):
        result = tm.round_trip_pickle(self.dtype)
        assert result == self.dtype


class TestCategoricalDtype(Base, tm.TestCase):

    def setUp(self):
        self.dtype = CategoricalDtype()

    def test_hash_vs_equality(self):
        # make sure that we satisfy is semantics
        dtype = self.dtype
        dtype2 = CategoricalDtype()
        self.assertTrue(dtype == dtype2)
        self.assertTrue(dtype2 == dtype)
        self.assertTrue(dtype is dtype2)
        self.assertTrue(dtype2 is dtype)
        self.assertTrue(hash(dtype) == hash(dtype2))

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype, 'category'))
        self.assertTrue(is_dtype_equal(self.dtype, CategoricalDtype()))
        assert not is_dtype_equal(self.dtype, 'foo')

    def test_construction_from_string(self):
        result = CategoricalDtype.construct_from_string('category')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        pytest.raises(
            TypeError, lambda: CategoricalDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        self.assertTrue(CategoricalDtype.is_dtype(self.dtype))
        self.assertTrue(CategoricalDtype.is_dtype('category'))
        self.assertTrue(CategoricalDtype.is_dtype(CategoricalDtype()))
        assert not CategoricalDtype.is_dtype('foo')
        assert not CategoricalDtype.is_dtype(np.float64)

    def test_basic(self):

        self.assertTrue(is_categorical_dtype(self.dtype))

        factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])

        s = Series(factor, name='A')

        # dtypes
        self.assertTrue(is_categorical_dtype(s.dtype))
        self.assertTrue(is_categorical_dtype(s))
        assert not is_categorical_dtype(np.dtype('float64'))

        self.assertTrue(is_categorical(s.dtype))
        self.assertTrue(is_categorical(s))
        assert not is_categorical(np.dtype('float64'))
        assert not is_categorical(1.0)


class TestDatetimeTZDtype(Base, tm.TestCase):

    def setUp(self):
        self.dtype = DatetimeTZDtype('ns', 'US/Eastern')

    def test_hash_vs_equality(self):
        # make sure that we satisfy is semantics
        dtype = self.dtype
        dtype2 = DatetimeTZDtype('ns', 'US/Eastern')
        dtype3 = DatetimeTZDtype(dtype2)
        self.assertTrue(dtype == dtype2)
        self.assertTrue(dtype2 == dtype)
        self.assertTrue(dtype3 == dtype)
        self.assertTrue(dtype is dtype2)
        self.assertTrue(dtype2 is dtype)
        self.assertTrue(dtype3 is dtype)
        self.assertTrue(hash(dtype) == hash(dtype2))
        self.assertTrue(hash(dtype) == hash(dtype3))

    def test_construction(self):
        pytest.raises(ValueError,
                      lambda: DatetimeTZDtype('ms', 'US/Eastern'))

    def test_subclass(self):
        a = DatetimeTZDtype('datetime64[ns, US/Eastern]')
        b = DatetimeTZDtype('datetime64[ns, CET]')

        self.assertTrue(issubclass(type(a), type(a)))
        self.assertTrue(issubclass(type(a), type(b)))

    def test_coerce_to_dtype(self):
        self.assertEqual(_coerce_to_dtype('datetime64[ns, US/Eastern]'),
                         DatetimeTZDtype('ns', 'US/Eastern'))
        self.assertEqual(_coerce_to_dtype('datetime64[ns, Asia/Tokyo]'),
                         DatetimeTZDtype('ns', 'Asia/Tokyo'))

    def test_compat(self):
        self.assertTrue(is_datetime64tz_dtype(self.dtype))
        self.assertTrue(is_datetime64tz_dtype('datetime64[ns, US/Eastern]'))
        self.assertTrue(is_datetime64_any_dtype(self.dtype))
        self.assertTrue(is_datetime64_any_dtype('datetime64[ns, US/Eastern]'))
        self.assertTrue(is_datetime64_ns_dtype(self.dtype))
        self.assertTrue(is_datetime64_ns_dtype('datetime64[ns, US/Eastern]'))
        assert not is_datetime64_dtype(self.dtype)
        assert not is_datetime64_dtype('datetime64[ns, US/Eastern]')

    def test_construction_from_string(self):
        result = DatetimeTZDtype('datetime64[ns, US/Eastern]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        result = DatetimeTZDtype.construct_from_string(
            'datetime64[ns, US/Eastern]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        pytest.raises(TypeError,
                      lambda: DatetimeTZDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        assert not DatetimeTZDtype.is_dtype(None)
        self.assertTrue(DatetimeTZDtype.is_dtype(self.dtype))
        self.assertTrue(DatetimeTZDtype.is_dtype('datetime64[ns, US/Eastern]'))
        assert not DatetimeTZDtype.is_dtype('foo')
        self.assertTrue(DatetimeTZDtype.is_dtype(DatetimeTZDtype(
            'ns', 'US/Pacific')))
        assert not DatetimeTZDtype.is_dtype(np.float64)

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype,
                                       'datetime64[ns, US/Eastern]'))
        self.assertTrue(is_dtype_equal(self.dtype, DatetimeTZDtype(
            'ns', 'US/Eastern')))
        assert not is_dtype_equal(self.dtype, 'foo')
        assert not is_dtype_equal(self.dtype, DatetimeTZDtype('ns', 'CET'))
        assert not is_dtype_equal(DatetimeTZDtype('ns', 'US/Eastern'),
                                  DatetimeTZDtype('ns', 'US/Pacific'))

        # numpy compat
        self.assertTrue(is_dtype_equal(np.dtype("M8[ns]"), "datetime64[ns]"))

    def test_basic(self):

        self.assertTrue(is_datetime64tz_dtype(self.dtype))

        dr = date_range('20130101', periods=3, tz='US/Eastern')
        s = Series(dr, name='A')

        # dtypes
        self.assertTrue(is_datetime64tz_dtype(s.dtype))
        self.assertTrue(is_datetime64tz_dtype(s))
        assert not is_datetime64tz_dtype(np.dtype('float64'))
        assert not is_datetime64tz_dtype(1.0)

        self.assertTrue(is_datetimetz(s))
        self.assertTrue(is_datetimetz(s.dtype))
        assert not is_datetimetz(np.dtype('float64'))
        assert not is_datetimetz(1.0)

    def test_dst(self):

        dr1 = date_range('2013-01-01', periods=3, tz='US/Eastern')
        s1 = Series(dr1, name='A')
        self.assertTrue(is_datetimetz(s1))

        dr2 = date_range('2013-08-01', periods=3, tz='US/Eastern')
        s2 = Series(dr2, name='A')
        self.assertTrue(is_datetimetz(s2))
        self.assertEqual(s1.dtype, s2.dtype)

    def test_parser(self):
        # pr #11245
        for tz, constructor in product(('UTC', 'US/Eastern'),
                                       ('M8', 'datetime64')):
            self.assertEqual(
                DatetimeTZDtype('%s[ns, %s]' % (constructor, tz)),
                DatetimeTZDtype('ns', tz),
            )

    def test_empty(self):
        dt = DatetimeTZDtype()
        with pytest.raises(AttributeError):
            str(dt)


class TestPeriodDtype(Base, tm.TestCase):

    def setUp(self):
        self.dtype = PeriodDtype('D')

    def test_construction(self):
        with pytest.raises(ValueError):
            PeriodDtype('xx')

        for s in ['period[D]', 'Period[D]', 'D']:
            dt = PeriodDtype(s)
            self.assertEqual(dt.freq, pd.tseries.offsets.Day())
            self.assertTrue(is_period_dtype(dt))

        for s in ['period[3D]', 'Period[3D]', '3D']:
            dt = PeriodDtype(s)
            self.assertEqual(dt.freq, pd.tseries.offsets.Day(3))
            self.assertTrue(is_period_dtype(dt))

        for s in ['period[26H]', 'Period[26H]', '26H',
                  'period[1D2H]', 'Period[1D2H]', '1D2H']:
            dt = PeriodDtype(s)
            self.assertEqual(dt.freq, pd.tseries.offsets.Hour(26))
            self.assertTrue(is_period_dtype(dt))

    def test_subclass(self):
        a = PeriodDtype('period[D]')
        b = PeriodDtype('period[3D]')

        self.assertTrue(issubclass(type(a), type(a)))
        self.assertTrue(issubclass(type(a), type(b)))

    def test_identity(self):
        assert PeriodDtype('period[D]') == PeriodDtype('period[D]')
        assert PeriodDtype('period[D]') is PeriodDtype('period[D]')

        assert PeriodDtype('period[3D]') == PeriodDtype('period[3D]')
        assert PeriodDtype('period[3D]') is PeriodDtype('period[3D]')

        assert PeriodDtype('period[1S1U]') == PeriodDtype('period[1000001U]')
        assert PeriodDtype('period[1S1U]') is PeriodDtype('period[1000001U]')

    def test_coerce_to_dtype(self):
        self.assertEqual(_coerce_to_dtype('period[D]'),
                         PeriodDtype('period[D]'))
        self.assertEqual(_coerce_to_dtype('period[3M]'),
                         PeriodDtype('period[3M]'))

    def test_compat(self):
        assert not is_datetime64_ns_dtype(self.dtype)
        assert not is_datetime64_ns_dtype('period[D]')
        assert not is_datetime64_dtype(self.dtype)
        assert not is_datetime64_dtype('period[D]')

    def test_construction_from_string(self):
        result = PeriodDtype('period[D]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        result = PeriodDtype.construct_from_string('period[D]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        with pytest.raises(TypeError):
            PeriodDtype.construct_from_string('foo')
        with pytest.raises(TypeError):
            PeriodDtype.construct_from_string('period[foo]')
        with pytest.raises(TypeError):
            PeriodDtype.construct_from_string('foo[D]')

        with pytest.raises(TypeError):
            PeriodDtype.construct_from_string('datetime64[ns]')
        with pytest.raises(TypeError):
            PeriodDtype.construct_from_string('datetime64[ns, US/Eastern]')

    def test_is_dtype(self):
        self.assertTrue(PeriodDtype.is_dtype(self.dtype))
        self.assertTrue(PeriodDtype.is_dtype('period[D]'))
        self.assertTrue(PeriodDtype.is_dtype('period[3D]'))
        self.assertTrue(PeriodDtype.is_dtype(PeriodDtype('3D')))
        self.assertTrue(PeriodDtype.is_dtype('period[U]'))
        self.assertTrue(PeriodDtype.is_dtype('period[S]'))
        self.assertTrue(PeriodDtype.is_dtype(PeriodDtype('U')))
        self.assertTrue(PeriodDtype.is_dtype(PeriodDtype('S')))

        assert not PeriodDtype.is_dtype('D')
        assert not PeriodDtype.is_dtype('3D')
        assert not PeriodDtype.is_dtype('U')
        assert not PeriodDtype.is_dtype('S')
        assert not PeriodDtype.is_dtype('foo')
        assert not PeriodDtype.is_dtype(np.object_)
        assert not PeriodDtype.is_dtype(np.int64)
        assert not PeriodDtype.is_dtype(np.float64)

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype, 'period[D]'))
        self.assertTrue(is_dtype_equal(self.dtype, PeriodDtype('D')))
        self.assertTrue(is_dtype_equal(self.dtype, PeriodDtype('D')))
        self.assertTrue(is_dtype_equal(PeriodDtype('D'), PeriodDtype('D')))

        assert not is_dtype_equal(self.dtype, 'D')
        assert not is_dtype_equal(PeriodDtype('D'), PeriodDtype('2D'))

    def test_basic(self):
        self.assertTrue(is_period_dtype(self.dtype))

        pidx = pd.period_range('2013-01-01 09:00', periods=5, freq='H')

        self.assertTrue(is_period_dtype(pidx.dtype))
        self.assertTrue(is_period_dtype(pidx))
        self.assertTrue(is_period(pidx))

        s = Series(pidx, name='A')
        # dtypes
        # series results in object dtype currently,
        # is_period checks period_arraylike
        assert not is_period_dtype(s.dtype)
        assert not is_period_dtype(s)
        self.assertTrue(is_period(s))

        assert not is_period_dtype(np.dtype('float64'))
        assert not is_period_dtype(1.0)
        assert not is_period(np.dtype('float64'))
        assert not is_period(1.0)

    def test_empty(self):
        dt = PeriodDtype()
        with pytest.raises(AttributeError):
            str(dt)

    def test_not_string(self):
        # though PeriodDtype has object kind, it cannot be string
        assert not is_string_dtype(PeriodDtype('D'))


class TestIntervalDtype(Base, tm.TestCase):

    # TODO: placeholder
    def setUp(self):
        self.dtype = IntervalDtype('int64')

    def test_construction(self):
        with pytest.raises(ValueError):
            IntervalDtype('xx')

        for s in ['interval[int64]', 'Interval[int64]', 'int64']:
            i = IntervalDtype(s)
            self.assertEqual(i.subtype, np.dtype('int64'))
            self.assertTrue(is_interval_dtype(i))

    def test_construction_generic(self):
        # generic
        i = IntervalDtype('interval')
        assert i.subtype is None
        self.assertTrue(is_interval_dtype(i))
        self.assertTrue(str(i) == 'interval')

        i = IntervalDtype()
        assert i.subtype is None
        self.assertTrue(is_interval_dtype(i))
        self.assertTrue(str(i) == 'interval')

    def test_subclass(self):
        a = IntervalDtype('interval[int64]')
        b = IntervalDtype('interval[int64]')

        self.assertTrue(issubclass(type(a), type(a)))
        self.assertTrue(issubclass(type(a), type(b)))

    def test_is_dtype(self):
        self.assertTrue(IntervalDtype.is_dtype(self.dtype))
        self.assertTrue(IntervalDtype.is_dtype('interval'))
        self.assertTrue(IntervalDtype.is_dtype(IntervalDtype('float64')))
        self.assertTrue(IntervalDtype.is_dtype(IntervalDtype('int64')))
        self.assertTrue(IntervalDtype.is_dtype(IntervalDtype(np.int64)))

        assert not IntervalDtype.is_dtype('D')
        assert not IntervalDtype.is_dtype('3D')
        assert not IntervalDtype.is_dtype('U')
        assert not IntervalDtype.is_dtype('S')
        assert not IntervalDtype.is_dtype('foo')
        assert not IntervalDtype.is_dtype(np.object_)
        assert not IntervalDtype.is_dtype(np.int64)
        assert not IntervalDtype.is_dtype(np.float64)

    def test_identity(self):
        self.assertEqual(IntervalDtype('interval[int64]'),
                         IntervalDtype('interval[int64]'))

    def test_coerce_to_dtype(self):
        self.assertEqual(_coerce_to_dtype('interval[int64]'),
                         IntervalDtype('interval[int64]'))

    def test_construction_from_string(self):
        result = IntervalDtype('interval[int64]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        result = IntervalDtype.construct_from_string('interval[int64]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        with pytest.raises(TypeError):
            IntervalDtype.construct_from_string('foo')
        with pytest.raises(TypeError):
            IntervalDtype.construct_from_string('interval[foo]')
        with pytest.raises(TypeError):
            IntervalDtype.construct_from_string('foo[int64]')

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype, 'interval[int64]'))
        self.assertTrue(is_dtype_equal(self.dtype, IntervalDtype('int64')))
        self.assertTrue(is_dtype_equal(self.dtype, IntervalDtype('int64')))
        self.assertTrue(is_dtype_equal(IntervalDtype('int64'),
                                       IntervalDtype('int64')))

        assert not is_dtype_equal(self.dtype, 'int64')
        assert not is_dtype_equal(IntervalDtype('int64'),
                                  IntervalDtype('float64'))

    def test_basic(self):
        self.assertTrue(is_interval_dtype(self.dtype))

        ii = IntervalIndex.from_breaks(range(3))

        self.assertTrue(is_interval_dtype(ii.dtype))
        self.assertTrue(is_interval_dtype(ii))

        s = Series(ii, name='A')

        # dtypes
        # series results in object dtype currently,
        assert not is_interval_dtype(s.dtype)
        assert not is_interval_dtype(s)

    def test_basic_dtype(self):
        self.assertTrue(is_interval_dtype('interval[int64]'))
        self.assertTrue(is_interval_dtype(IntervalIndex.from_tuples([(0, 1)])))
        self.assertTrue(is_interval_dtype
                        (IntervalIndex.from_breaks(np.arange(4))))
        self.assertTrue(is_interval_dtype(
            IntervalIndex.from_breaks(date_range('20130101', periods=3))))
        assert not is_interval_dtype('U')
        assert not is_interval_dtype('S')
        assert not is_interval_dtype('foo')
        assert not is_interval_dtype(np.object_)
        assert not is_interval_dtype(np.int64)
        assert not is_interval_dtype(np.float64)
