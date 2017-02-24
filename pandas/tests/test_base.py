# -*- coding: utf-8 -*-
from __future__ import print_function

import re
import sys
from datetime import datetime, timedelta
import pytest
import numpy as np

import pandas as pd
import pandas.compat as compat
from pandas.types.common import (is_object_dtype, is_datetimetz,
                                 needs_i8_conversion)
import pandas.util.testing as tm
from pandas import (Series, Index, DatetimeIndex, TimedeltaIndex, PeriodIndex,
                    Timedelta)
from pandas.compat import StringIO
from pandas.compat.numpy import np_array_datetime64_compat
from pandas.core.base import PandasDelegate, NoNewAttributesMixin
from pandas.tseries.base import DatetimeIndexOpsMixin


class CheckStringMixin(object):

    def test_string_methods_dont_fail(self):
        repr(self.container)
        str(self.container)
        bytes(self.container)
        if not compat.PY3:
            unicode(self.container)  # noqa

    def test_tricky_container(self):
        if not hasattr(self, 'unicode_container'):
            pytest.skip('Need unicode_container to test with this')
        repr(self.unicode_container)
        str(self.unicode_container)
        bytes(self.unicode_container)
        if not compat.PY3:
            unicode(self.unicode_container)  # noqa


class CheckImmutable(object):
    mutable_regex = re.compile('does not support mutable operations')

    def check_mutable_error(self, *args, **kwargs):
        # pass whatever functions you normally would to assertRaises (after the
        # Exception kind)
        tm.assertRaisesRegexp(TypeError, self.mutable_regex, *args, **kwargs)

    def test_no_mutable_funcs(self):
        def setitem():
            self.container[0] = 5

        self.check_mutable_error(setitem)

        def setslice():
            self.container[1:2] = 3

        self.check_mutable_error(setslice)

        def delitem():
            del self.container[0]

        self.check_mutable_error(delitem)

        def delslice():
            del self.container[0:3]

        self.check_mutable_error(delslice)
        mutable_methods = getattr(self, "mutable_methods", [])
        for meth in mutable_methods:
            self.check_mutable_error(getattr(self.container, meth))

    def test_slicing_maintains_type(self):
        result = self.container[1:2]
        expected = self.lst[1:2]
        self.check_result(result, expected)

    def check_result(self, result, expected, klass=None):
        klass = klass or self.klass
        self.assertIsInstance(result, klass)
        self.assertEqual(result, expected)


class TestPandasDelegate(tm.TestCase):

    class Delegator(object):
        _properties = ['foo']
        _methods = ['bar']

        def _set_foo(self, value):
            self.foo = value

        def _get_foo(self):
            return self.foo

        foo = property(_get_foo, _set_foo, doc="foo property")

        def bar(self, *args, **kwargs):
            """ a test bar method """
            pass

    class Delegate(PandasDelegate):

        def __init__(self, obj):
            self.obj = obj

    def setUp(self):
        pass

    def test_invalida_delgation(self):
        # these show that in order for the delegation to work
        # the _delegate_* methods need to be overriden to not raise a TypeError

        self.Delegate._add_delegate_accessors(
            delegate=self.Delegator,
            accessors=self.Delegator._properties,
            typ='property'
        )
        self.Delegate._add_delegate_accessors(
            delegate=self.Delegator,
            accessors=self.Delegator._methods,
            typ='method'
        )

        delegate = self.Delegate(self.Delegator())

        def f():
            delegate.foo

        self.assertRaises(TypeError, f)

        def f():
            delegate.foo = 5

        self.assertRaises(TypeError, f)

        def f():
            delegate.foo()

        self.assertRaises(TypeError, f)

    def test_memory_usage(self):
        # Delegate does not implement memory_usage.
        # Check that we fall back to in-built `__sizeof__`
        # GH 12924
        delegate = self.Delegate(self.Delegator())
        sys.getsizeof(delegate)


class Ops(tm.TestCase):

    def _allow_na_ops(self, obj):
        """Whether to skip test cases including NaN"""
        if (isinstance(obj, Index) and
                (obj.is_boolean() or not obj._can_hold_na)):
            # don't test boolean / int64 index
            return False
        return True

    def setUp(self):
        self.bool_index = tm.makeBoolIndex(10, name='a')
        self.int_index = tm.makeIntIndex(10, name='a')
        self.float_index = tm.makeFloatIndex(10, name='a')
        self.dt_index = tm.makeDateIndex(10, name='a')
        self.dt_tz_index = tm.makeDateIndex(10, name='a').tz_localize(
            tz='US/Eastern')
        self.period_index = tm.makePeriodIndex(10, name='a')
        self.string_index = tm.makeStringIndex(10, name='a')
        self.unicode_index = tm.makeUnicodeIndex(10, name='a')

        arr = np.random.randn(10)
        self.int_series = Series(arr, index=self.int_index, name='a')
        self.float_series = Series(arr, index=self.float_index, name='a')
        self.dt_series = Series(arr, index=self.dt_index, name='a')
        self.dt_tz_series = self.dt_tz_index.to_series(keep_tz=True)
        self.period_series = Series(arr, index=self.period_index, name='a')
        self.string_series = Series(arr, index=self.string_index, name='a')

        types = ['bool', 'int', 'float', 'dt', 'dt_tz', 'period', 'string',
                 'unicode']
        fmts = ["{0}_{1}".format(t, f)
                for t in types for f in ['index', 'series']]
        self.objs = [getattr(self, f)
                     for f in fmts if getattr(self, f, None) is not None]

    def check_ops_properties(self, props, filter=None, ignore_failures=False):
        for op in props:
            for o in self.is_valid_objs:

                # if a filter, skip if it doesn't match
                if filter is not None:
                    filt = o.index if isinstance(o, Series) else o
                    if not filter(filt):
                        continue

                try:
                    if isinstance(o, Series):
                        expected = Series(
                            getattr(o.index, op), index=o.index, name='a')
                    else:
                        expected = getattr(o, op)
                except (AttributeError):
                    if ignore_failures:
                        continue

                result = getattr(o, op)

                # these couuld be series, arrays or scalars
                if isinstance(result, Series) and isinstance(expected, Series):
                    tm.assert_series_equal(result, expected)
                elif isinstance(result, Index) and isinstance(expected, Index):
                    tm.assert_index_equal(result, expected)
                elif isinstance(result, np.ndarray) and isinstance(expected,
                                                                   np.ndarray):
                    self.assert_numpy_array_equal(result, expected)
                else:
                    self.assertEqual(result, expected)

            # freq raises AttributeError on an Int64Index because its not
            # defined we mostly care about Series hwere anyhow
            if not ignore_failures:
                for o in self.not_valid_objs:

                    # an object that is datetimelike will raise a TypeError,
                    # otherwise an AttributeError
                    if issubclass(type(o), DatetimeIndexOpsMixin):
                        self.assertRaises(TypeError, lambda: getattr(o, op))
                    else:
                        self.assertRaises(AttributeError,
                                          lambda: getattr(o, op))

    def test_binary_ops_docs(self):
        from pandas import DataFrame, Panel
        op_map = {'add': '+',
                  'sub': '-',
                  'mul': '*',
                  'mod': '%',
                  'pow': '**',
                  'truediv': '/',
                  'floordiv': '//'}
        for op_name in ['add', 'sub', 'mul', 'mod', 'pow', 'truediv',
                        'floordiv']:
            for klass in [Series, DataFrame, Panel]:
                operand1 = klass.__name__.lower()
                operand2 = 'other'
                op = op_map[op_name]
                expected_str = ' '.join([operand1, op, operand2])
                self.assertTrue(expected_str in getattr(klass,
                                                        op_name).__doc__)

                # reverse version of the binary ops
                expected_str = ' '.join([operand2, op, operand1])
                self.assertTrue(expected_str in getattr(klass, 'r' +
                                                        op_name).__doc__)


class TestIndexOps(Ops):

    def setUp(self):
        super(TestIndexOps, self).setUp()
        self.is_valid_objs = [o for o in self.objs if o._allow_index_ops]
        self.not_valid_objs = [o for o in self.objs if not o._allow_index_ops]

    def test_none_comparison(self):

        # bug brought up by #1079
        # changed from TypeError in 0.17.0
        for o in self.is_valid_objs:
            if isinstance(o, Series):

                o[0] = np.nan

                # noinspection PyComparisonWithNone
                result = o == None  # noqa
                self.assertFalse(result.iat[0])
                self.assertFalse(result.iat[1])

                # noinspection PyComparisonWithNone
                result = o != None  # noqa
                self.assertTrue(result.iat[0])
                self.assertTrue(result.iat[1])

                result = None == o  # noqa
                self.assertFalse(result.iat[0])
                self.assertFalse(result.iat[1])

                # this fails for numpy < 1.9
                # and oddly for *some* platforms
                # result = None != o  # noqa
                # self.assertTrue(result.iat[0])
                # self.assertTrue(result.iat[1])

                result = None > o
                self.assertFalse(result.iat[0])
                self.assertFalse(result.iat[1])

                result = o < None
                self.assertFalse(result.iat[0])
                self.assertFalse(result.iat[1])

    def test_ndarray_compat_properties(self):

        for o in self.objs:

            # check that we work
            for p in ['shape', 'dtype', 'flags', 'T', 'strides', 'itemsize',
                      'nbytes']:
                self.assertIsNotNone(getattr(o, p, None))
            self.assertTrue(hasattr(o, 'base'))

            # if we have a datetimelike dtype then needs a view to work
            # but the user is responsible for that
            try:
                self.assertIsNotNone(o.data)
            except ValueError:
                pass

            self.assertRaises(ValueError, o.item)  # len > 1
            self.assertEqual(o.ndim, 1)
            self.assertEqual(o.size, len(o))

        self.assertEqual(Index([1]).item(), 1)
        self.assertEqual(Series([1]).item(), 1)

    def test_ops(self):
        for op in ['max', 'min']:
            for o in self.objs:
                result = getattr(o, op)()
                if not isinstance(o, PeriodIndex):
                    expected = getattr(o.values, op)()
                else:
                    expected = pd.Period(ordinal=getattr(o._values, op)(),
                                         freq=o.freq)
                try:
                    self.assertEqual(result, expected)
                except TypeError:
                    # comparing tz-aware series with np.array results in
                    # TypeError
                    expected = expected.astype('M8[ns]').astype('int64')
                    self.assertEqual(result.value, expected)

    def test_nanops(self):
        # GH 7261
        for op in ['max', 'min']:
            for klass in [Index, Series]:

                obj = klass([np.nan, 2.0])
                self.assertEqual(getattr(obj, op)(), 2.0)

                obj = klass([np.nan])
                self.assertTrue(pd.isnull(getattr(obj, op)()))

                obj = klass([])
                self.assertTrue(pd.isnull(getattr(obj, op)()))

                obj = klass([pd.NaT, datetime(2011, 11, 1)])
                # check DatetimeIndex monotonic path
                self.assertEqual(getattr(obj, op)(), datetime(2011, 11, 1))

                obj = klass([pd.NaT, datetime(2011, 11, 1), pd.NaT])
                # check DatetimeIndex non-monotonic path
                self.assertEqual(getattr(obj, op)(), datetime(2011, 11, 1))

        # argmin/max
        obj = Index(np.arange(5, dtype='int64'))
        self.assertEqual(obj.argmin(), 0)
        self.assertEqual(obj.argmax(), 4)

        obj = Index([np.nan, 1, np.nan, 2])
        self.assertEqual(obj.argmin(), 1)
        self.assertEqual(obj.argmax(), 3)

        obj = Index([np.nan])
        self.assertEqual(obj.argmin(), -1)
        self.assertEqual(obj.argmax(), -1)

        obj = Index([pd.NaT, datetime(2011, 11, 1), datetime(2011, 11, 2),
                     pd.NaT])
        self.assertEqual(obj.argmin(), 1)
        self.assertEqual(obj.argmax(), 2)

        obj = Index([pd.NaT])
        self.assertEqual(obj.argmin(), -1)
        self.assertEqual(obj.argmax(), -1)

    def test_value_counts_unique_nunique(self):
        for orig in self.objs:
            o = orig.copy()
            klass = type(o)
            values = o._values

            if isinstance(values, Index):
                # reset name not to affect latter process
                values.name = None

            # create repeated values, 'n'th element is repeated by n+1 times
            # skip boolean, because it only has 2 values at most
            if isinstance(o, Index) and o.is_boolean():
                continue
            elif isinstance(o, Index):
                expected_index = pd.Index(o[::-1])
                expected_index.name = None
                o = o.repeat(range(1, len(o) + 1))
                o.name = 'a'
            else:
                expected_index = pd.Index(values[::-1])
                idx = o.index.repeat(range(1, len(o) + 1))
                rep = np.repeat(values, range(1, len(o) + 1))
                o = klass(rep, index=idx, name='a')

            # check values has the same dtype as the original
            self.assertEqual(o.dtype, orig.dtype)

            expected_s = Series(range(10, 0, -1), index=expected_index,
                                dtype='int64', name='a')

            result = o.value_counts()
            tm.assert_series_equal(result, expected_s)
            self.assertTrue(result.index.name is None)
            self.assertEqual(result.name, 'a')

            result = o.unique()
            if isinstance(o, Index):
                self.assertTrue(isinstance(result, o.__class__))
                self.assert_index_equal(result, orig)
            elif is_datetimetz(o):
                # datetimetz Series returns array of Timestamp
                self.assertEqual(result[0], orig[0])
                for r in result:
                    self.assertIsInstance(r, pd.Timestamp)
                tm.assert_numpy_array_equal(result,
                                            orig._values.asobject.values)
            else:
                tm.assert_numpy_array_equal(result, orig.values)

            self.assertEqual(o.nunique(), len(np.unique(o.values)))

    def test_value_counts_unique_nunique_null(self):

        for null_obj in [np.nan, None]:
            for orig in self.objs:
                o = orig.copy()
                klass = type(o)
                values = o._values

                if not self._allow_na_ops(o):
                    continue

                # special assign to the numpy array
                if is_datetimetz(o):
                    if isinstance(o, DatetimeIndex):
                        v = o.asi8
                        v[0:2] = pd.tslib.iNaT
                        values = o._shallow_copy(v)
                    else:
                        o = o.copy()
                        o[0:2] = pd.tslib.iNaT
                        values = o._values

                elif needs_i8_conversion(o):
                    values[0:2] = pd.tslib.iNaT
                    values = o._shallow_copy(values)
                else:
                    values[0:2] = null_obj
                # check values has the same dtype as the original

                self.assertEqual(values.dtype, o.dtype)

                # create repeated values, 'n'th element is repeated by n+1
                # times
                if isinstance(o, (DatetimeIndex, PeriodIndex)):
                    expected_index = o.copy()
                    expected_index.name = None

                    # attach name to klass
                    o = klass(values.repeat(range(1, len(o) + 1)))
                    o.name = 'a'
                else:
                    if is_datetimetz(o):
                        expected_index = orig._values._shallow_copy(values)
                    else:
                        expected_index = pd.Index(values)
                    expected_index.name = None
                    o = o.repeat(range(1, len(o) + 1))
                    o.name = 'a'

                # check values has the same dtype as the original
                self.assertEqual(o.dtype, orig.dtype)
                # check values correctly have NaN
                nanloc = np.zeros(len(o), dtype=np.bool)
                nanloc[:3] = True
                if isinstance(o, Index):
                    self.assert_numpy_array_equal(pd.isnull(o), nanloc)
                else:
                    exp = pd.Series(nanloc, o.index, name='a')
                    self.assert_series_equal(pd.isnull(o), exp)

                expected_s_na = Series(list(range(10, 2, -1)) + [3],
                                       index=expected_index[9:0:-1],
                                       dtype='int64', name='a')
                expected_s = Series(list(range(10, 2, -1)),
                                    index=expected_index[9:1:-1],
                                    dtype='int64', name='a')

                result_s_na = o.value_counts(dropna=False)
                tm.assert_series_equal(result_s_na, expected_s_na)
                self.assertTrue(result_s_na.index.name is None)
                self.assertEqual(result_s_na.name, 'a')
                result_s = o.value_counts()
                tm.assert_series_equal(o.value_counts(), expected_s)
                self.assertTrue(result_s.index.name is None)
                self.assertEqual(result_s.name, 'a')

                result = o.unique()
                if isinstance(o, Index):
                    tm.assert_index_equal(result,
                                          Index(values[1:], name='a'))
                elif is_datetimetz(o):
                    # unable to compare NaT / nan
                    tm.assert_numpy_array_equal(result[1:],
                                                values[2:].asobject.values)
                    self.assertIs(result[0], pd.NaT)
                else:
                    tm.assert_numpy_array_equal(result[1:], values[2:])

                    self.assertTrue(pd.isnull(result[0]))
                    self.assertEqual(result.dtype, orig.dtype)

                self.assertEqual(o.nunique(), 8)
                self.assertEqual(o.nunique(dropna=False), 9)

    def test_value_counts_inferred(self):
        klasses = [Index, Series]
        for klass in klasses:
            s_values = ['a', 'b', 'b', 'b', 'b', 'c', 'd', 'd', 'a', 'a']
            s = klass(s_values)
            expected = Series([4, 3, 2, 1], index=['b', 'a', 'd', 'c'])
            tm.assert_series_equal(s.value_counts(), expected)

            if isinstance(s, Index):
                exp = Index(np.unique(np.array(s_values, dtype=np.object_)))
                tm.assert_index_equal(s.unique(), exp)
            else:
                exp = np.unique(np.array(s_values, dtype=np.object_))
                tm.assert_numpy_array_equal(s.unique(), exp)

            self.assertEqual(s.nunique(), 4)
            # don't sort, have to sort after the fact as not sorting is
            # platform-dep
            hist = s.value_counts(sort=False).sort_values()
            expected = Series([3, 1, 4, 2], index=list('acbd')).sort_values()
            tm.assert_series_equal(hist, expected)

            # sort ascending
            hist = s.value_counts(ascending=True)
            expected = Series([1, 2, 3, 4], index=list('cdab'))
            tm.assert_series_equal(hist, expected)

            # relative histogram.
            hist = s.value_counts(normalize=True)
            expected = Series([.4, .3, .2, .1], index=['b', 'a', 'd', 'c'])
            tm.assert_series_equal(hist, expected)

    def test_value_counts_bins(self):
        klasses = [Index, Series]
        for klass in klasses:
            s_values = ['a', 'b', 'b', 'b', 'b', 'c', 'd', 'd', 'a', 'a']
            s = klass(s_values)

            # bins
            self.assertRaises(TypeError,
                              lambda bins: s.value_counts(bins=bins), 1)

            s1 = Series([1, 1, 2, 3])
            res1 = s1.value_counts(bins=1)
            exp1 = Series({0.998: 4})
            tm.assert_series_equal(res1, exp1)
            res1n = s1.value_counts(bins=1, normalize=True)
            exp1n = Series({0.998: 1.0})
            tm.assert_series_equal(res1n, exp1n)

            if isinstance(s1, Index):
                tm.assert_index_equal(s1.unique(), Index([1, 2, 3]))
            else:
                exp = np.array([1, 2, 3], dtype=np.int64)
                tm.assert_numpy_array_equal(s1.unique(), exp)

            self.assertEqual(s1.nunique(), 3)

            res4 = s1.value_counts(bins=4)
            exp4 = Series({0.998: 2,
                           1.5: 1,
                           2.0: 0,
                           2.5: 1}, index=[0.998, 2.5, 1.5, 2.0])
            tm.assert_series_equal(res4, exp4)
            res4n = s1.value_counts(bins=4, normalize=True)
            exp4n = Series(
                {0.998: 0.5,
                 1.5: 0.25,
                 2.0: 0.0,
                 2.5: 0.25}, index=[0.998, 2.5, 1.5, 2.0])
            tm.assert_series_equal(res4n, exp4n)

            # handle NA's properly
            s_values = ['a', 'b', 'b', 'b', np.nan, np.nan,
                        'd', 'd', 'a', 'a', 'b']
            s = klass(s_values)
            expected = Series([4, 3, 2], index=['b', 'a', 'd'])
            tm.assert_series_equal(s.value_counts(), expected)

            if isinstance(s, Index):
                exp = Index(['a', 'b', np.nan, 'd'])
                tm.assert_index_equal(s.unique(), exp)
            else:
                exp = np.array(['a', 'b', np.nan, 'd'], dtype=object)
                tm.assert_numpy_array_equal(s.unique(), exp)
            self.assertEqual(s.nunique(), 3)

            s = klass({})
            expected = Series([], dtype=np.int64)
            tm.assert_series_equal(s.value_counts(), expected,
                                   check_index_type=False)
            # returned dtype differs depending on original
            if isinstance(s, Index):
                self.assert_index_equal(s.unique(), Index([]),
                                        exact=False)
            else:
                self.assert_numpy_array_equal(s.unique(), np.array([]),
                                              check_dtype=False)

            self.assertEqual(s.nunique(), 0)

    def test_value_counts_datetime64(self):
        klasses = [Index, Series]
        for klass in klasses:
            # GH 3002, datetime64[ns]
            # don't test names though
            txt = "\n".join(['xxyyzz20100101PIE', 'xxyyzz20100101GUM',
                             'xxyyzz20100101EGG', 'xxyyww20090101EGG',
                             'foofoo20080909PIE', 'foofoo20080909GUM'])
            f = StringIO(txt)
            df = pd.read_fwf(f, widths=[6, 8, 3],
                             names=["person_id", "dt", "food"],
                             parse_dates=["dt"])

            s = klass(df['dt'].copy())
            s.name = None

            idx = pd.to_datetime(['2010-01-01 00:00:00Z',
                                  '2008-09-09 00:00:00Z',
                                  '2009-01-01 00:00:00X'])
            expected_s = Series([3, 2, 1], index=idx)
            tm.assert_series_equal(s.value_counts(), expected_s)

            expected = np_array_datetime64_compat(['2010-01-01 00:00:00Z',
                                                   '2009-01-01 00:00:00Z',
                                                   '2008-09-09 00:00:00Z'],
                                                  dtype='datetime64[ns]')
            if isinstance(s, Index):
                tm.assert_index_equal(s.unique(), DatetimeIndex(expected))
            else:
                tm.assert_numpy_array_equal(s.unique(), expected)

            self.assertEqual(s.nunique(), 3)

            # with NaT
            s = df['dt'].copy()
            s = klass([v for v in s.values] + [pd.NaT])

            result = s.value_counts()
            self.assertEqual(result.index.dtype, 'datetime64[ns]')
            tm.assert_series_equal(result, expected_s)

            result = s.value_counts(dropna=False)
            expected_s[pd.NaT] = 1
            tm.assert_series_equal(result, expected_s)

            unique = s.unique()
            self.assertEqual(unique.dtype, 'datetime64[ns]')

            # numpy_array_equal cannot compare pd.NaT
            if isinstance(s, Index):
                exp_idx = DatetimeIndex(expected.tolist() + [pd.NaT])
                tm.assert_index_equal(unique, exp_idx)
            else:
                tm.assert_numpy_array_equal(unique[:3], expected)
                self.assertTrue(pd.isnull(unique[3]))

            self.assertEqual(s.nunique(), 3)
            self.assertEqual(s.nunique(dropna=False), 4)

            # timedelta64[ns]
            td = df.dt - df.dt + timedelta(1)
            td = klass(td, name='dt')

            result = td.value_counts()
            expected_s = Series([6], index=[Timedelta('1day')], name='dt')
            tm.assert_series_equal(result, expected_s)

            expected = TimedeltaIndex(['1 days'], name='dt')
            if isinstance(td, Index):
                tm.assert_index_equal(td.unique(), expected)
            else:
                tm.assert_numpy_array_equal(td.unique(), expected.values)

            td2 = timedelta(1) + (df.dt - df.dt)
            td2 = klass(td2, name='dt')
            result2 = td2.value_counts()
            tm.assert_series_equal(result2, expected_s)

    def test_factorize(self):
        for orig in self.objs:
            o = orig.copy()

            if isinstance(o, Index) and o.is_boolean():
                exp_arr = np.array([0, 1] + [0] * 8, dtype=np.intp)
                exp_uniques = o
                exp_uniques = Index([False, True])
            else:
                exp_arr = np.array(range(len(o)), dtype=np.intp)
                exp_uniques = o
            labels, uniques = o.factorize()

            self.assert_numpy_array_equal(labels, exp_arr)
            if isinstance(o, Series):
                self.assert_index_equal(uniques, Index(orig),
                                        check_names=False)
            else:
                # factorize explicitly resets name
                self.assert_index_equal(uniques, exp_uniques,
                                        check_names=False)

    def test_factorize_repeated(self):
        for orig in self.objs:
            o = orig.copy()

            # don't test boolean
            if isinstance(o, Index) and o.is_boolean():
                continue

            # sort by value, and create duplicates
            if isinstance(o, Series):
                o = o.sort_values()
                n = o.iloc[5:].append(o)
            else:
                indexer = o.argsort()
                o = o.take(indexer)
                n = o[5:].append(o)

            exp_arr = np.array([5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                               dtype=np.intp)
            labels, uniques = n.factorize(sort=True)

            self.assert_numpy_array_equal(labels, exp_arr)
            if isinstance(o, Series):
                self.assert_index_equal(uniques, Index(orig).sort_values(),
                                        check_names=False)
            else:
                self.assert_index_equal(uniques, o, check_names=False)

            exp_arr = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4],
                               np.intp)
            labels, uniques = n.factorize(sort=False)
            self.assert_numpy_array_equal(labels, exp_arr)

            if isinstance(o, Series):
                expected = Index(o.iloc[5:10].append(o.iloc[:5]))
                self.assert_index_equal(uniques, expected, check_names=False)
            else:
                expected = o[5:10].append(o[:5])
                self.assert_index_equal(uniques, expected, check_names=False)

    def test_duplicated_drop_duplicates_index(self):
        # GH 4060
        for original in self.objs:
            if isinstance(original, Index):

                # special case
                if original.is_boolean():
                    result = original.drop_duplicates()
                    expected = Index([False, True], name='a')
                    tm.assert_index_equal(result, expected)
                    continue

                # original doesn't have duplicates
                expected = np.array([False] * len(original), dtype=bool)
                duplicated = original.duplicated()
                tm.assert_numpy_array_equal(duplicated, expected)
                self.assertTrue(duplicated.dtype == bool)
                result = original.drop_duplicates()
                tm.assert_index_equal(result, original)
                self.assertFalse(result is original)

                # has_duplicates
                self.assertFalse(original.has_duplicates)

                # create repeated values, 3rd and 5th values are duplicated
                idx = original[list(range(len(original))) + [5, 3]]
                expected = np.array([False] * len(original) + [True, True],
                                    dtype=bool)
                duplicated = idx.duplicated()
                tm.assert_numpy_array_equal(duplicated, expected)
                self.assertTrue(duplicated.dtype == bool)
                tm.assert_index_equal(idx.drop_duplicates(), original)

                base = [False] * len(idx)
                base[3] = True
                base[5] = True
                expected = np.array(base)

                duplicated = idx.duplicated(keep='last')
                tm.assert_numpy_array_equal(duplicated, expected)
                self.assertTrue(duplicated.dtype == bool)
                result = idx.drop_duplicates(keep='last')
                tm.assert_index_equal(result, idx[~expected])

                # deprecate take_last
                with tm.assert_produces_warning(FutureWarning):
                    duplicated = idx.duplicated(take_last=True)
                tm.assert_numpy_array_equal(duplicated, expected)
                self.assertTrue(duplicated.dtype == bool)
                with tm.assert_produces_warning(FutureWarning):
                    result = idx.drop_duplicates(take_last=True)
                tm.assert_index_equal(result, idx[~expected])

                base = [False] * len(original) + [True, True]
                base[3] = True
                base[5] = True
                expected = np.array(base)

                duplicated = idx.duplicated(keep=False)
                tm.assert_numpy_array_equal(duplicated, expected)
                self.assertTrue(duplicated.dtype == bool)
                result = idx.drop_duplicates(keep=False)
                tm.assert_index_equal(result, idx[~expected])

                with tm.assertRaisesRegexp(
                        TypeError, r"drop_duplicates\(\) got an unexpected "
                        "keyword argument"):
                    idx.drop_duplicates(inplace=True)

            else:
                expected = Series([False] * len(original),
                                  index=original.index, name='a')
                tm.assert_series_equal(original.duplicated(), expected)
                result = original.drop_duplicates()
                tm.assert_series_equal(result, original)
                self.assertFalse(result is original)

                idx = original.index[list(range(len(original))) + [5, 3]]
                values = original._values[list(range(len(original))) + [5, 3]]
                s = Series(values, index=idx, name='a')

                expected = Series([False] * len(original) + [True, True],
                                  index=idx, name='a')
                tm.assert_series_equal(s.duplicated(), expected)
                tm.assert_series_equal(s.drop_duplicates(), original)

                base = [False] * len(idx)
                base[3] = True
                base[5] = True
                expected = Series(base, index=idx, name='a')

                tm.assert_series_equal(s.duplicated(keep='last'), expected)
                tm.assert_series_equal(s.drop_duplicates(keep='last'),
                                       s[~np.array(base)])

                # deprecate take_last
                with tm.assert_produces_warning(FutureWarning):
                    tm.assert_series_equal(
                        s.duplicated(take_last=True), expected)
                with tm.assert_produces_warning(FutureWarning):
                    tm.assert_series_equal(s.drop_duplicates(take_last=True),
                                           s[~np.array(base)])
                base = [False] * len(original) + [True, True]
                base[3] = True
                base[5] = True
                expected = Series(base, index=idx, name='a')

                tm.assert_series_equal(s.duplicated(keep=False), expected)
                tm.assert_series_equal(s.drop_duplicates(keep=False),
                                       s[~np.array(base)])

                s.drop_duplicates(inplace=True)
                tm.assert_series_equal(s, original)

    def test_drop_duplicates_series_vs_dataframe(self):
        # GH 14192
        df = pd.DataFrame({'a': [1, 1, 1, 'one', 'one'],
                           'b': [2, 2, np.nan, np.nan, np.nan],
                           'c': [3, 3, np.nan, np.nan, 'three'],
                           'd': [1, 2, 3, 4, 4],
                           'e': [datetime(2015, 1, 1), datetime(2015, 1, 1),
                                 datetime(2015, 2, 1), pd.NaT, pd.NaT]
                           })
        for column in df.columns:
            for keep in ['first', 'last', False]:
                dropped_frame = df[[column]].drop_duplicates(keep=keep)
                dropped_series = df[column].drop_duplicates(keep=keep)
                tm.assert_frame_equal(dropped_frame, dropped_series.to_frame())

    def test_fillna(self):
        # # GH 11343
        # though Index.fillna and Series.fillna has separate impl,
        # test here to confirm these works as the same

        for orig in self.objs:

            o = orig.copy()
            values = o.values

            # values will not be changed
            result = o.fillna(o.astype(object).values[0])
            if isinstance(o, Index):
                self.assert_index_equal(o, result)
            else:
                self.assert_series_equal(o, result)
            # check shallow_copied
            self.assertFalse(o is result)

        for null_obj in [np.nan, None]:
            for orig in self.objs:
                o = orig.copy()
                klass = type(o)

                if not self._allow_na_ops(o):
                    continue

                if needs_i8_conversion(o):

                    values = o.astype(object).values
                    fill_value = values[0]
                    values[0:2] = pd.NaT
                else:
                    values = o.values.copy()
                    fill_value = o.values[0]
                    values[0:2] = null_obj

                expected = [fill_value] * 2 + list(values[2:])

                expected = klass(expected)
                o = klass(values)

                # check values has the same dtype as the original
                self.assertEqual(o.dtype, orig.dtype)

                result = o.fillna(fill_value)
                if isinstance(o, Index):
                    self.assert_index_equal(result, expected)
                else:
                    self.assert_series_equal(result, expected)
                # check shallow_copied
                self.assertFalse(o is result)

    def test_memory_usage(self):
        for o in self.objs:
            res = o.memory_usage()
            res_deep = o.memory_usage(deep=True)

            if (is_object_dtype(o) or (isinstance(o, Series) and
                                       is_object_dtype(o.index))):
                # if there are objects, only deep will pick them up
                self.assertTrue(res_deep > res)
            else:
                self.assertEqual(res, res_deep)

            if isinstance(o, Series):
                self.assertEqual(
                    (o.memory_usage(index=False) +
                        o.index.memory_usage()),
                    o.memory_usage(index=True)
                )

            # sys.getsizeof will call the .memory_usage with
            # deep=True, and add on some GC overhead
            diff = res_deep - sys.getsizeof(o)
            self.assertTrue(abs(diff) < 100)

    def test_searchsorted(self):
        # See gh-12238
        for o in self.objs:
            index = np.searchsorted(o, max(o))
            self.assertTrue(0 <= index <= len(o))

            index = np.searchsorted(o, max(o), sorter=range(len(o)))
            self.assertTrue(0 <= index <= len(o))

    def test_validate_bool_args(self):
        invalid_values = [1, "True", [1, 2, 3], 5.0]

        for value in invalid_values:
            with self.assertRaises(ValueError):
                self.int_series.drop_duplicates(inplace=value)


class TestTranspose(Ops):
    errmsg = "the 'axes' parameter is not supported"

    def test_transpose(self):
        for obj in self.objs:
            if isinstance(obj, Index):
                tm.assert_index_equal(obj.transpose(), obj)
            else:
                tm.assert_series_equal(obj.transpose(), obj)

    def test_transpose_non_default_axes(self):
        for obj in self.objs:
            tm.assertRaisesRegexp(ValueError, self.errmsg,
                                  obj.transpose, 1)
            tm.assertRaisesRegexp(ValueError, self.errmsg,
                                  obj.transpose, axes=1)

    def test_numpy_transpose(self):
        for obj in self.objs:
            if isinstance(obj, Index):
                tm.assert_index_equal(np.transpose(obj), obj)
            else:
                tm.assert_series_equal(np.transpose(obj), obj)

            tm.assertRaisesRegexp(ValueError, self.errmsg,
                                  np.transpose, obj, axes=1)


class TestNoNewAttributesMixin(tm.TestCase):

    def test_mixin(self):
        class T(NoNewAttributesMixin):
            pass

        t = T()
        self.assertFalse(hasattr(t, "__frozen"))
        t.a = "test"
        self.assertEqual(t.a, "test")
        t._freeze()
        # self.assertTrue("__frozen" not in dir(t))
        self.assertIs(getattr(t, "__frozen"), True)

        def f():
            t.b = "test"

        self.assertRaises(AttributeError, f)
        self.assertFalse(hasattr(t, "b"))
