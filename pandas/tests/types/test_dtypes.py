# -*- coding: utf-8 -*-
from itertools import product

import numpy as np
import pandas as pd
from pandas import Series, Categorical, date_range

from pandas.types.dtypes import DatetimeTZDtype, PeriodDtype, CategoricalDtype
from pandas.types.common import (is_categorical_dtype, is_categorical,
                                 is_datetime64tz_dtype, is_datetimetz,
                                 is_period_dtype, is_period,
                                 is_dtype_equal, is_datetime64_ns_dtype,
                                 is_datetime64_dtype,
                                 is_datetime64_any_dtype, is_string_dtype,
                                 _coerce_to_dtype)
import pandas.util.testing as tm


class Base(object):

    def test_hash(self):
        hash(self.dtype)

    def test_equality_invalid(self):
        self.assertRaises(self.dtype == 'foo')
        self.assertFalse(is_dtype_equal(self.dtype, np.int64))

    def test_numpy_informed(self):

        # np.dtype doesn't know about our new dtype
        def f():
            np.dtype(self.dtype)

        self.assertRaises(TypeError, f)

        self.assertNotEqual(self.dtype, np.str_)
        self.assertNotEqual(np.str_, self.dtype)

    def test_pickle(self):
        result = self.round_trip_pickle(self.dtype)
        self.assertEqual(result, self.dtype)


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
        self.assertFalse(is_dtype_equal(self.dtype, 'foo'))

    def test_construction_from_string(self):
        result = CategoricalDtype.construct_from_string('category')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        self.assertRaises(
            TypeError, lambda: CategoricalDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        self.assertTrue(CategoricalDtype.is_dtype(self.dtype))
        self.assertTrue(CategoricalDtype.is_dtype('category'))
        self.assertTrue(CategoricalDtype.is_dtype(CategoricalDtype()))
        self.assertFalse(CategoricalDtype.is_dtype('foo'))
        self.assertFalse(CategoricalDtype.is_dtype(np.float64))

    def test_basic(self):

        self.assertTrue(is_categorical_dtype(self.dtype))

        factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])

        s = Series(factor, name='A')

        # dtypes
        self.assertTrue(is_categorical_dtype(s.dtype))
        self.assertTrue(is_categorical_dtype(s))
        self.assertFalse(is_categorical_dtype(np.dtype('float64')))

        self.assertTrue(is_categorical(s.dtype))
        self.assertTrue(is_categorical(s))
        self.assertFalse(is_categorical(np.dtype('float64')))
        self.assertFalse(is_categorical(1.0))


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
        self.assertRaises(ValueError,
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
        self.assertFalse(is_datetime64_dtype(self.dtype))
        self.assertFalse(is_datetime64_dtype('datetime64[ns, US/Eastern]'))

    def test_construction_from_string(self):
        result = DatetimeTZDtype('datetime64[ns, US/Eastern]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        result = DatetimeTZDtype.construct_from_string(
            'datetime64[ns, US/Eastern]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        self.assertRaises(TypeError,
                          lambda: DatetimeTZDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        self.assertTrue(DatetimeTZDtype.is_dtype(self.dtype))
        self.assertTrue(DatetimeTZDtype.is_dtype('datetime64[ns, US/Eastern]'))
        self.assertFalse(DatetimeTZDtype.is_dtype('foo'))
        self.assertTrue(DatetimeTZDtype.is_dtype(DatetimeTZDtype(
            'ns', 'US/Pacific')))
        self.assertFalse(DatetimeTZDtype.is_dtype(np.float64))

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype,
                                       'datetime64[ns, US/Eastern]'))
        self.assertTrue(is_dtype_equal(self.dtype, DatetimeTZDtype(
            'ns', 'US/Eastern')))
        self.assertFalse(is_dtype_equal(self.dtype, 'foo'))
        self.assertFalse(is_dtype_equal(self.dtype, DatetimeTZDtype('ns',
                                                                    'CET')))
        self.assertFalse(is_dtype_equal(
            DatetimeTZDtype('ns', 'US/Eastern'), DatetimeTZDtype(
                'ns', 'US/Pacific')))

        # numpy compat
        self.assertTrue(is_dtype_equal(np.dtype("M8[ns]"), "datetime64[ns]"))

    def test_basic(self):

        self.assertTrue(is_datetime64tz_dtype(self.dtype))

        dr = date_range('20130101', periods=3, tz='US/Eastern')
        s = Series(dr, name='A')

        # dtypes
        self.assertTrue(is_datetime64tz_dtype(s.dtype))
        self.assertTrue(is_datetime64tz_dtype(s))
        self.assertFalse(is_datetime64tz_dtype(np.dtype('float64')))
        self.assertFalse(is_datetime64tz_dtype(1.0))

        self.assertTrue(is_datetimetz(s))
        self.assertTrue(is_datetimetz(s.dtype))
        self.assertFalse(is_datetimetz(np.dtype('float64')))
        self.assertFalse(is_datetimetz(1.0))

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
        with tm.assertRaises(AttributeError):
            str(dt)


class TestPeriodDtype(Base, tm.TestCase):

    def setUp(self):
        self.dtype = PeriodDtype('D')

    def test_construction(self):
        with tm.assertRaises(ValueError):
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
        self.assertEqual(PeriodDtype('period[D]'),
                         PeriodDtype('period[D]'))
        self.assertIs(PeriodDtype('period[D]'),
                      PeriodDtype('period[D]'))

        self.assertEqual(PeriodDtype('period[3D]'),
                         PeriodDtype('period[3D]'))
        self.assertIs(PeriodDtype('period[3D]'),
                      PeriodDtype('period[3D]'))

        self.assertEqual(PeriodDtype('period[1S1U]'),
                         PeriodDtype('period[1000001U]'))
        self.assertIs(PeriodDtype('period[1S1U]'),
                      PeriodDtype('period[1000001U]'))

    def test_coerce_to_dtype(self):
        self.assertEqual(_coerce_to_dtype('period[D]'),
                         PeriodDtype('period[D]'))
        self.assertEqual(_coerce_to_dtype('period[3M]'),
                         PeriodDtype('period[3M]'))

    def test_compat(self):
        self.assertFalse(is_datetime64_ns_dtype(self.dtype))
        self.assertFalse(is_datetime64_ns_dtype('period[D]'))
        self.assertFalse(is_datetime64_dtype(self.dtype))
        self.assertFalse(is_datetime64_dtype('period[D]'))

    def test_construction_from_string(self):
        result = PeriodDtype('period[D]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        result = PeriodDtype.construct_from_string('period[D]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        with tm.assertRaises(TypeError):
            PeriodDtype.construct_from_string('foo')
        with tm.assertRaises(TypeError):
            PeriodDtype.construct_from_string('period[foo]')
        with tm.assertRaises(TypeError):
            PeriodDtype.construct_from_string('foo[D]')

        with tm.assertRaises(TypeError):
            PeriodDtype.construct_from_string('datetime64[ns]')
        with tm.assertRaises(TypeError):
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

        self.assertFalse(PeriodDtype.is_dtype('D'))
        self.assertFalse(PeriodDtype.is_dtype('3D'))
        self.assertFalse(PeriodDtype.is_dtype('U'))
        self.assertFalse(PeriodDtype.is_dtype('S'))
        self.assertFalse(PeriodDtype.is_dtype('foo'))
        self.assertFalse(PeriodDtype.is_dtype(np.object_))
        self.assertFalse(PeriodDtype.is_dtype(np.int64))
        self.assertFalse(PeriodDtype.is_dtype(np.float64))

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype, 'period[D]'))
        self.assertTrue(is_dtype_equal(self.dtype, PeriodDtype('D')))
        self.assertTrue(is_dtype_equal(self.dtype, PeriodDtype('D')))
        self.assertTrue(is_dtype_equal(PeriodDtype('D'), PeriodDtype('D')))

        self.assertFalse(is_dtype_equal(self.dtype, 'D'))
        self.assertFalse(is_dtype_equal(PeriodDtype('D'), PeriodDtype('2D')))

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
        self.assertFalse(is_period_dtype(s.dtype))
        self.assertFalse(is_period_dtype(s))
        self.assertTrue(is_period(s))

        self.assertFalse(is_period_dtype(np.dtype('float64')))
        self.assertFalse(is_period_dtype(1.0))
        self.assertFalse(is_period(np.dtype('float64')))
        self.assertFalse(is_period(1.0))

    def test_empty(self):
        dt = PeriodDtype()
        with tm.assertRaises(AttributeError):
            str(dt)

    def test_not_string(self):
        # though PeriodDtype has object kind, it cannot be string
        self.assertFalse(is_string_dtype(PeriodDtype('D')))
