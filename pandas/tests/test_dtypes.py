# -*- coding: utf-8 -*-
from itertools import product

import nose
import numpy as np
from pandas import Series, Categorical, date_range
import pandas.core.common as com
from pandas.core.common import (CategoricalDtype, is_categorical_dtype, is_categorical,
                                DatetimeTZDtype, is_datetime64tz_dtype, is_datetimetz,
                                is_dtype_equal, is_datetime64_ns_dtype, is_datetime64_dtype)
import pandas.util.testing as tm

_multiprocess_can_split_ = True

class Base(object):

    def test_hash(self):
        hash(self.dtype)

    def test_equality_invalid(self):
        self.assertRaises(self.dtype == 'foo')

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

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype, 'category'))
        self.assertTrue(is_dtype_equal(self.dtype, CategoricalDtype()))
        self.assertFalse(is_dtype_equal(self.dtype, 'foo'))

    def test_construction_from_string(self):
        result = CategoricalDtype.construct_from_string('category')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        self.assertRaises(TypeError, lambda : CategoricalDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        self.assertTrue(CategoricalDtype.is_dtype(self.dtype))
        self.assertTrue(CategoricalDtype.is_dtype('category'))
        self.assertTrue(CategoricalDtype.is_dtype(CategoricalDtype()))
        self.assertFalse(CategoricalDtype.is_dtype('foo'))
        self.assertFalse(CategoricalDtype.is_dtype(np.float64))

    def test_basic(self):

        self.assertTrue(is_categorical_dtype(self.dtype))

        factor = Categorical.from_array(['a', 'b', 'b', 'a',
                                         'a', 'c', 'c', 'c'])

        s = Series(factor,name='A')

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
        self.dtype = DatetimeTZDtype('ns','US/Eastern')

    def test_construction(self):
        self.assertRaises(ValueError, lambda : DatetimeTZDtype('ms','US/Eastern'))

    def test_subclass(self):
        a = DatetimeTZDtype('datetime64[ns, US/Eastern]')
        b = DatetimeTZDtype('datetime64[ns, CET]')

        self.assertTrue(issubclass(type(a), type(a)))
        self.assertTrue(issubclass(type(a), type(b)))

    def test_compat(self):
        self.assertFalse(is_datetime64_ns_dtype(self.dtype))
        self.assertFalse(is_datetime64_ns_dtype('datetime64[ns, US/Eastern]'))
        self.assertFalse(is_datetime64_dtype(self.dtype))
        self.assertFalse(is_datetime64_dtype('datetime64[ns, US/Eastern]'))

    def test_construction_from_string(self):
        result = DatetimeTZDtype('datetime64[ns, US/Eastern]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        result = DatetimeTZDtype.construct_from_string('datetime64[ns, US/Eastern]')
        self.assertTrue(is_dtype_equal(self.dtype, result))
        self.assertRaises(TypeError, lambda : DatetimeTZDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        self.assertTrue(DatetimeTZDtype.is_dtype(self.dtype))
        self.assertTrue(DatetimeTZDtype.is_dtype('datetime64[ns, US/Eastern]'))
        self.assertFalse(DatetimeTZDtype.is_dtype('foo'))
        self.assertTrue(DatetimeTZDtype.is_dtype(DatetimeTZDtype('ns','US/Pacific')))
        self.assertFalse(DatetimeTZDtype.is_dtype(np.float64))

    def test_equality(self):
        self.assertTrue(is_dtype_equal(self.dtype, 'datetime64[ns, US/Eastern]'))
        self.assertTrue(is_dtype_equal(self.dtype, DatetimeTZDtype('ns','US/Eastern')))
        self.assertFalse(is_dtype_equal(self.dtype, 'foo'))
        self.assertFalse(is_dtype_equal(self.dtype, DatetimeTZDtype('ns','CET')))
        self.assertFalse(is_dtype_equal(DatetimeTZDtype('ns','US/Eastern'), DatetimeTZDtype('ns','US/Pacific')))

        # numpy compat
        self.assertTrue(is_dtype_equal(np.dtype("M8[ns]"),"datetime64[ns]"))

    def test_basic(self):

        self.assertTrue(is_datetime64tz_dtype(self.dtype))

        dr = date_range('20130101',periods=3,tz='US/Eastern')
        s = Series(dr,name='A')

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




if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
