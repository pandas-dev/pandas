from __future__ import print_function
import re
from datetime import datetime, timedelta
import numpy as np
import pandas.compat as compat
import pandas as pd
from pandas.compat import u, StringIO
from pandas.core.base import FrozenList, FrozenNDArray, PandasDelegate, DatetimeIndexOpsMixin
from pandas.util.testing import assertRaisesRegexp, assert_isinstance
from pandas.tseries.common import is_datetimelike
from pandas import Series, Index, Int64Index, DatetimeIndex, PeriodIndex
import pandas.tslib as tslib
import nose

import pandas.util.testing as tm

class CheckStringMixin(object):
    def test_string_methods_dont_fail(self):
        repr(self.container)
        str(self.container)
        bytes(self.container)
        if not compat.PY3:
            unicode(self.container)

    def test_tricky_container(self):
        if not hasattr(self, 'unicode_container'):
            raise nose.SkipTest('Need unicode_container to test with this')
        repr(self.unicode_container)
        str(self.unicode_container)
        bytes(self.unicode_container)
        if not compat.PY3:
            unicode(self.unicode_container)


class CheckImmutable(object):
    mutable_regex = re.compile('does not support mutable operations')

    def check_mutable_error(self, *args, **kwargs):
        # pass whatever functions you normally would to assertRaises (after the Exception kind)
        assertRaisesRegexp(TypeError, self.mutable_regex, *args, **kwargs)

    def test_no_mutable_funcs(self):
        def setitem(): self.container[0] = 5

        self.check_mutable_error(setitem)

        def setslice(): self.container[1:2] = 3

        self.check_mutable_error(setslice)

        def delitem(): del self.container[0]

        self.check_mutable_error(delitem)

        def delslice(): del self.container[0:3]

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
        assert_isinstance(result, klass)
        self.assertEqual(result, expected)


class TestFrozenList(CheckImmutable, CheckStringMixin, tm.TestCase):
    mutable_methods = ('extend', 'pop', 'remove', 'insert')
    unicode_container = FrozenList([u("\u05d0"), u("\u05d1"), "c"])

    def setUp(self):
        self.lst = [1, 2, 3, 4, 5]
        self.container = FrozenList(self.lst)
        self.klass = FrozenList

    def test_add(self):
        result = self.container + (1, 2, 3)
        expected = FrozenList(self.lst + [1, 2, 3])
        self.check_result(result, expected)

        result = (1, 2, 3) + self.container
        expected = FrozenList([1, 2, 3] + self.lst)
        self.check_result(result, expected)

    def test_inplace(self):
        q = r = self.container
        q += [5]
        self.check_result(q, self.lst + [5])
        # other shouldn't be mutated
        self.check_result(r, self.lst)


class TestFrozenNDArray(CheckImmutable, CheckStringMixin, tm.TestCase):
    mutable_methods = ('put', 'itemset', 'fill')
    unicode_container = FrozenNDArray([u("\u05d0"), u("\u05d1"), "c"])

    def setUp(self):
        self.lst = [3, 5, 7, -2]
        self.container = FrozenNDArray(self.lst)
        self.klass = FrozenNDArray

    def test_shallow_copying(self):
        original = self.container.copy()
        assert_isinstance(self.container.view(), FrozenNDArray)
        self.assertFalse(isinstance(self.container.view(np.ndarray), FrozenNDArray))
        self.assertIsNot(self.container.view(), self.container)
        self.assert_numpy_array_equal(self.container, original)
        # shallow copy should be the same too
        assert_isinstance(self.container._shallow_copy(), FrozenNDArray)
        # setting should not be allowed
        def testit(container): container[0] = 16

        self.check_mutable_error(testit, self.container)

    def test_values(self):
        original = self.container.view(np.ndarray).copy()
        n = original[0] + 15
        vals = self.container.values()
        self.assert_numpy_array_equal(original, vals)
        self.assertIsNot(original, vals)
        vals[0] = n
        self.assert_numpy_array_equal(self.container, original)
        self.assertEqual(vals[0], n)


class TestPandasDelegate(tm.TestCase):

    def setUp(self):
        pass

    def test_invalida_delgation(self):
        # these show that in order for the delegation to work
        # the _delegate_* methods need to be overriden to not raise a TypeError

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
        Delegate._add_delegate_accessors(delegate=Delegator,
                                         accessors=Delegator._properties,
                                         typ='property')
        Delegate._add_delegate_accessors(delegate=Delegator,
                                         accessors=Delegator._methods,
                                         typ='method')

        delegate = Delegate(Delegator())

        def f():
            delegate.foo
        self.assertRaises(TypeError, f)
        def f():
            delegate.foo = 5
        self.assertRaises(TypeError, f)
        def f():
            delegate.foo()
        self.assertRaises(TypeError, f)


class Ops(tm.TestCase):
    def setUp(self):
        self.int_index     = tm.makeIntIndex(10)
        self.float_index   = tm.makeFloatIndex(10)
        self.dt_index      = tm.makeDateIndex(10)
        self.dt_tz_index   = tm.makeDateIndex(10).tz_localize(tz='US/Eastern')
        self.period_index  = tm.makePeriodIndex(10)
        self.string_index  = tm.makeStringIndex(10)

        arr = np.random.randn(10)
        self.int_series    = Series(arr, index=self.int_index)
        self.float_series  = Series(arr, index=self.int_index)
        self.dt_series     = Series(arr, index=self.dt_index)
        self.dt_tz_series  = self.dt_tz_index.to_series(keep_tz=True)
        self.period_series = Series(arr, index=self.period_index)
        self.string_series = Series(arr, index=self.string_index)

        types = ['int','float','dt', 'dt_tz', 'period','string']
        self.objs = [ getattr(self,"{0}_{1}".format(t,f)) for t in types for f in ['index','series'] ]

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
                        expected = Series(getattr(o.index,op),index=o.index)
                    else:
                        expected = getattr(o,op)
                except (AttributeError):
                    if ignore_failures:
                        continue

                result = getattr(o,op)

                # these couuld be series, arrays or scalars
                if isinstance(result,Series) and isinstance(expected,Series):
                    tm.assert_series_equal(result,expected)
                elif isinstance(result,Index) and isinstance(expected,Index):
                    tm.assert_index_equal(result,expected)
                elif isinstance(result,np.ndarray) and isinstance(expected,np.ndarray):
                    self.assert_numpy_array_equal(result,expected)
                else:
                    self.assertEqual(result, expected)

            # freq raises AttributeError on an Int64Index because its not defined
            # we mostly care about Series hwere anyhow
            if not ignore_failures:
                for o in self.not_valid_objs:

                    # an object that is datetimelike will raise a TypeError, otherwise
                    # an AttributeError
                    if issubclass(type(o), DatetimeIndexOpsMixin):
                        self.assertRaises(TypeError, lambda : getattr(o,op))
                    else:
                        self.assertRaises(AttributeError, lambda : getattr(o,op))


class TestIndexOps(Ops):

    def setUp(self):
        super(TestIndexOps, self).setUp()
        self.is_valid_objs  = [ o for o in self.objs if o._allow_index_ops ]
        self.not_valid_objs = [ o for o in self.objs if not o._allow_index_ops ]

    def test_ndarray_compat_properties(self):

        for o in self.objs:

            # check that we work
            for p in ['shape', 'dtype', 'base', 'flags', 'T',
                      'strides', 'itemsize', 'nbytes']:
                self.assertIsNotNone(getattr(o, p, None))

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
        for op in ['max','min']:
            for o in self.objs:
                result = getattr(o,op)()
                if not isinstance(o, PeriodIndex):
                    expected = getattr(o.values, op)()
                else:
                    expected = pd.Period(ordinal=getattr(o.values, op)(), freq=o.freq)
                try:
                    self.assertEqual(result, expected)
                except ValueError:
                    # comparing tz-aware series with np.array results in ValueError
                    expected = expected.astype('M8[ns]').astype('int64')
                    self.assertEqual(result.value, expected)

    def test_nanops(self):
        # GH 7261
        for op in ['max','min']:
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
        obj = Index(np.arange(5,dtype='int64'))
        self.assertEqual(obj.argmin(),0)
        self.assertEqual(obj.argmax(),4)

        obj = Index([np.nan, 1, np.nan, 2])
        self.assertEqual(obj.argmin(),1)
        self.assertEqual(obj.argmax(),3)

        obj = Index([np.nan])
        self.assertEqual(obj.argmin(),-1)
        self.assertEqual(obj.argmax(),-1)

        obj = Index([pd.NaT, datetime(2011, 11, 1), datetime(2011,11,2),pd.NaT])
        self.assertEqual(obj.argmin(),1)
        self.assertEqual(obj.argmax(),2)

        obj = Index([pd.NaT])
        self.assertEqual(obj.argmin(),-1)
        self.assertEqual(obj.argmax(),-1)

    def test_value_counts_unique_nunique(self):
        for o in self.objs:
            klass = type(o)
            values = o.values

            # create repeated values, 'n'th element is repeated by n+1 times
            if isinstance(o, PeriodIndex):
                # freq must be specified because repeat makes freq ambiguous
                expected_index = o[::-1]
                o = klass(np.repeat(values, range(1, len(o) + 1)), freq=o.freq)
            elif isinstance(o, Index):
                expected_index = values[::-1]
                o = klass(np.repeat(values, range(1, len(o) + 1)))
            else:
                expected_index = values[::-1]
                idx = np.repeat(o.index.values, range(1, len(o) + 1))
                o = klass(np.repeat(values, range(1, len(o) + 1)), index=idx)

            expected_s = Series(range(10, 0, -1), index=expected_index, dtype='int64')
            tm.assert_series_equal(o.value_counts(), expected_s)

            result = o.unique()
            if isinstance(o, (DatetimeIndex, PeriodIndex)):
                self.assertTrue(isinstance(result, o.__class__))
                self.assertEqual(result.name, o.name)
                self.assertEqual(result.freq, o.freq)

            self.assert_numpy_array_equal(result, values)

            self.assertEqual(o.nunique(), len(np.unique(o.values)))

        for null_obj in [np.nan, None]:
            for o in self.objs:
                klass = type(o)
                values = o.values

                if ((isinstance(o, Int64Index) and not isinstance(o,
                    (DatetimeIndex, PeriodIndex)))):
                    # skips int64 because it doesn't allow to include nan or None
                    continue

                # special assign to the numpy array
                if o.values.dtype == 'datetime64[ns]' or isinstance(o, PeriodIndex):
                    values[0:2] = pd.tslib.iNaT
                else:
                    values[0:2] = null_obj

                # create repeated values, 'n'th element is repeated by n+1 times
                if isinstance(o, PeriodIndex):
                    # freq must be specified because repeat makes freq ambiguous
                    expected_index = o
                    o = klass(np.repeat(values, range(1, len(o) + 1)), freq=o.freq)
                elif isinstance(o, Index):
                    expected_index = values
                    o = klass(np.repeat(values, range(1, len(o) + 1)))
                else:
                    expected_index = values
                    idx = np.repeat(o.index.values, range(1, len(o) + 1))
                    o = klass(np.repeat(values, range(1, len(o) + 1)), index=idx)

                expected_s_na = Series(list(range(10, 2, -1)) +[3], index=expected_index[9:0:-1], dtype='int64')
                expected_s = Series(list(range(10, 2, -1)), index=expected_index[9:1:-1], dtype='int64')

                tm.assert_series_equal(o.value_counts(dropna=False), expected_s_na)
                tm.assert_series_equal(o.value_counts(), expected_s)

                # numpy_array_equal cannot compare arrays includes nan
                result = o.unique()
                self.assert_numpy_array_equal(result[1:], values[2:])

                if isinstance(o, (DatetimeIndex, PeriodIndex)):
                    self.assertTrue(result.asi8[0] == pd.tslib.iNaT)
                else:
                    self.assertTrue(pd.isnull(result[0]))

                self.assertEqual(o.nunique(), 8)
                self.assertEqual(o.nunique(dropna=False), 9)

    def test_value_counts_inferred(self):
        klasses = [Index, Series]
        for klass in klasses:
            s_values = ['a', 'b', 'b', 'b', 'b', 'c', 'd', 'd', 'a', 'a']
            s = klass(s_values)
            expected = Series([4, 3, 2, 1], index=['b', 'a', 'd', 'c'])
            tm.assert_series_equal(s.value_counts(), expected)

            self.assert_numpy_array_equal(s.unique(), np.unique(s_values))
            self.assertEqual(s.nunique(), 4)
            # don't sort, have to sort after the fact as not sorting is platform-dep
            hist = s.value_counts(sort=False)
            hist.sort()
            expected = Series([3, 1, 4, 2], index=list('acbd'))
            expected.sort()
            tm.assert_series_equal(hist, expected)

            # sort ascending
            hist = s.value_counts(ascending=True)
            expected = Series([1, 2, 3, 4], index=list('cdab'))
            tm.assert_series_equal(hist, expected)

            # relative histogram.
            hist = s.value_counts(normalize=True)
            expected = Series([.4, .3, .2, .1], index=['b', 'a', 'd', 'c'])
            tm.assert_series_equal(hist, expected)

            # bins
            self.assertRaises(TypeError, lambda bins: s.value_counts(bins=bins), 1)

            s1 = Series([1, 1, 2, 3])
            res1 = s1.value_counts(bins=1)
            exp1 = Series({0.998: 4})
            tm.assert_series_equal(res1, exp1)
            res1n = s1.value_counts(bins=1, normalize=True)
            exp1n = Series({0.998: 1.0})
            tm.assert_series_equal(res1n, exp1n)

            self.assert_numpy_array_equal(s1.unique(), np.array([1, 2, 3]))
            self.assertEqual(s1.nunique(), 3)

            res4 = s1.value_counts(bins=4)
            exp4 = Series({0.998: 2, 1.5: 1, 2.0: 0, 2.5: 1}, index=[0.998, 2.5, 1.5, 2.0])
            tm.assert_series_equal(res4, exp4)
            res4n = s1.value_counts(bins=4, normalize=True)
            exp4n = Series({0.998: 0.5, 1.5: 0.25, 2.0: 0.0, 2.5: 0.25}, index=[0.998, 2.5, 1.5, 2.0])
            tm.assert_series_equal(res4n, exp4n)

            # handle NA's properly
            s_values = ['a', 'b', 'b', 'b', np.nan, np.nan, 'd', 'd', 'a', 'a', 'b']
            s = klass(s_values)
            expected = Series([4, 3, 2], index=['b', 'a', 'd'])
            tm.assert_series_equal(s.value_counts(), expected)

            self.assert_numpy_array_equal(s.unique(), np.array(['a', 'b', np.nan, 'd'], dtype='O'))
            self.assertEqual(s.nunique(), 3)

            s = klass({})
            expected = Series([], dtype=np.int64)
            tm.assert_series_equal(s.value_counts(), expected)
            self.assert_numpy_array_equal(s.unique(), np.array([]))
            self.assertEqual(s.nunique(), 0)

            # GH 3002, datetime64[ns]
            txt = "\n".join(['xxyyzz20100101PIE', 'xxyyzz20100101GUM', 'xxyyzz20100101EGG',
                             'xxyyww20090101EGG', 'foofoo20080909PIE', 'foofoo20080909GUM'])
            f = StringIO(txt)
            df = pd.read_fwf(f, widths=[6, 8, 3], names=["person_id", "dt", "food"],
                             parse_dates=["dt"])

            s = klass(df['dt'].copy())

            idx = pd.to_datetime(['2010-01-01 00:00:00Z', '2008-09-09 00:00:00Z', '2009-01-01 00:00:00X'])
            expected_s = Series([3, 2, 1], index=idx)
            tm.assert_series_equal(s.value_counts(), expected_s)

            expected = np.array(['2010-01-01 00:00:00Z', '2009-01-01 00:00:00Z', '2008-09-09 00:00:00Z'],
                                dtype='datetime64[ns]')
            if isinstance(s, DatetimeIndex):
                expected = DatetimeIndex(expected)
                self.assertTrue(s.unique().equals(expected))
            else:
                self.assert_numpy_array_equal(s.unique(), expected)

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
            self.assert_numpy_array_equal(unique[:3], expected)
            self.assertTrue(unique[3] is pd.NaT or unique[3].astype('int64') == pd.tslib.iNaT)

            self.assertEqual(s.nunique(), 3)
            self.assertEqual(s.nunique(dropna=False), 4)

            # timedelta64[ns]
            td = df.dt - df.dt + timedelta(1)
            td = klass(td)

            result = td.value_counts()
            expected_s = Series([6], index=[86400000000000])
            self.assertEqual(result.index.dtype, 'int64')
            tm.assert_series_equal(result, expected_s)

            # get nanoseconds to compare
            expected = np.array([86400000000000])
            self.assert_numpy_array_equal(td.unique(), expected)
            self.assertEqual(td.nunique(), 1)

            td2 = timedelta(1) + (df.dt - df.dt)
            td2 = klass(td2)
            result2 = td2.value_counts()

            self.assertEqual(result2.index.dtype, 'int64')
            tm.assert_series_equal(result2, expected_s)

            self.assert_numpy_array_equal(td.unique(), expected)
            self.assertEqual(td.nunique(), 1)

    def test_factorize(self):
        for o in self.objs:
            exp_arr = np.array(range(len(o)))
            labels, uniques = o.factorize()

            self.assert_numpy_array_equal(labels, exp_arr)
            if isinstance(o, Series):
                expected = Index(o.values)
                self.assert_numpy_array_equal(uniques, expected)
            else:
                self.assertTrue(uniques.equals(o))

        for o in self.objs:
            # sort by value, and create duplicates
            if isinstance(o, Series):
                o.sort()
            else:
                indexer = o.argsort()
                o = o.take(indexer)
            n = o[5:].append(o)

            exp_arr = np.array([5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
            labels, uniques = n.factorize(sort=True)

            self.assert_numpy_array_equal(labels, exp_arr)
            if isinstance(o, Series):
                expected = Index(o.values)
                self.assert_numpy_array_equal(uniques, expected)
            else:
                self.assertTrue(uniques.equals(o))

            exp_arr = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4])
            labels, uniques = n.factorize(sort=False)
            self.assert_numpy_array_equal(labels, exp_arr)

            if isinstance(o, Series):
                expected = Index(np.concatenate([o.values[5:10], o.values[:5]]))
                self.assert_numpy_array_equal(uniques, expected)
            else:
                expected = o[5:].append(o[:5])
                self.assertTrue(uniques.equals(expected))

    def test_duplicated_drop_duplicates(self):
        # GH 4060
        for original in self.objs:

            if isinstance(original, Index):
                # original doesn't have duplicates
                expected = Index([False] * len(original))
                tm.assert_index_equal(original.duplicated(), expected)
                result = original.drop_duplicates()
                tm.assert_index_equal(result, original)
                self.assertFalse(result is original)

                # create repeated values, 3rd and 5th values are duplicated
                idx = original[list(range(len(original))) + [5, 3]]
                expected = Index([False] * len(original) + [True, True])
                tm.assert_index_equal(idx.duplicated(), expected)
                tm.assert_index_equal(idx.drop_duplicates(), original)

                last_base = [False] * len(idx)
                last_base[3] = True
                last_base[5] = True
                expected = Index(last_base)
                tm.assert_index_equal(idx.duplicated(take_last=True), expected)
                tm.assert_index_equal(idx.drop_duplicates(take_last=True),
                                      idx[~np.array(last_base)])

                with tm.assertRaisesRegexp(TypeError,
                                           "drop_duplicates\(\) got an unexpected keyword argument"):
                    idx.drop_duplicates(inplace=True)

            else:
                expected = Series([False] * len(original), index=original.index)
                tm.assert_series_equal(original.duplicated(), expected)
                result = original.drop_duplicates()
                tm.assert_series_equal(result, original)
                self.assertFalse(result is original)

                idx = original.index[list(range(len(original))) + [5, 3]]
                values = original.values[list(range(len(original))) + [5, 3]]
                s = Series(values, index=idx)

                expected = Series([False] * len(original) + [True, True], index=idx)
                tm.assert_series_equal(s.duplicated(), expected)
                tm.assert_series_equal(s.drop_duplicates(), original)

                last_base = [False] * len(idx)
                last_base[3] = True
                last_base[5] = True
                expected = Series(last_base, index=idx)
                expected
                tm.assert_series_equal(s.duplicated(take_last=True), expected)
                tm.assert_series_equal(s.drop_duplicates(take_last=True),
                                       s[~np.array(last_base)])

                s.drop_duplicates(inplace=True)
                tm.assert_series_equal(s, original)


class TestDatetimeIndexOps(Ops):
    tz = [None, 'UTC', 'Asia/Tokyo', 'US/Eastern',
          'dateutil/Asia/Singapore', 'dateutil/US/Pacific']

    def setUp(self):
        super(TestDatetimeIndexOps, self).setUp()
        mask = lambda x: isinstance(x, DatetimeIndex) or isinstance(x, PeriodIndex) or is_datetimelike(x)
        self.is_valid_objs  = [ o for o in self.objs if mask(o) ]
        self.not_valid_objs = [ o for o in self.objs if not mask(o) ]

    def test_ops_properties(self):
        self.check_ops_properties(['year','month','day','hour','minute','second','weekofyear','week','dayofweek','dayofyear','quarter'])
        self.check_ops_properties(['date','time','microsecond','nanosecond', 'is_month_start', 'is_month_end', 'is_quarter_start',
                                   'is_quarter_end', 'is_year_start', 'is_year_end'], lambda x: isinstance(x,DatetimeIndex))

    def test_ops_properties_basic(self):

        # sanity check that the behavior didn't change
        # GH7206
        for op in ['year','day','second','weekday']:
            self.assertRaises(TypeError, lambda x: getattr(self.dt_series,op))

        # attribute access should still work!
        s = Series(dict(year=2000,month=1,day=10))
        self.assertEquals(s.year,2000)
        self.assertEquals(s.month,1)
        self.assertEquals(s.day,10)
        self.assertRaises(AttributeError, lambda : s.weekday)

    def test_asobject_tolist(self):
        idx = pd.date_range(start='2013-01-01', periods=4, freq='M', name='idx')
        expected_list = [pd.Timestamp('2013-01-31'), pd.Timestamp('2013-02-28'),
                         pd.Timestamp('2013-03-31'), pd.Timestamp('2013-04-30')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = pd.date_range(start='2013-01-01', periods=4, freq='M', name='idx', tz='Asia/Tokyo')
        expected_list = [pd.Timestamp('2013-01-31', tz='Asia/Tokyo'),
                         pd.Timestamp('2013-02-28', tz='Asia/Tokyo'),
                         pd.Timestamp('2013-03-31', tz='Asia/Tokyo'),
                         pd.Timestamp('2013-04-30', tz='Asia/Tokyo')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = DatetimeIndex([datetime(2013, 1, 1), datetime(2013, 1, 2),
                             pd.NaT, datetime(2013, 1, 4)], name='idx')
        expected_list = [pd.Timestamp('2013-01-01'), pd.Timestamp('2013-01-02'),
                         pd.NaT, pd.Timestamp('2013-01-04')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

    def test_minmax(self):
        for tz in self.tz:
            # monotonic
            idx1 = pd.DatetimeIndex([pd.NaT, '2011-01-01', '2011-01-02',
                                     '2011-01-03'], tz=tz)
            self.assertTrue(idx1.is_monotonic)

            # non-monotonic
            idx2 = pd.DatetimeIndex(['2011-01-01', pd.NaT, '2011-01-03',
                                     '2011-01-02', pd.NaT], tz=tz)
            self.assertFalse(idx2.is_monotonic)

            for idx in [idx1, idx2]:
                self.assertEqual(idx.min(), pd.Timestamp('2011-01-01', tz=tz))
                self.assertEqual(idx.max(), pd.Timestamp('2011-01-03', tz=tz))

        for op in ['min', 'max']:
            # Return NaT
            obj = DatetimeIndex([])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = DatetimeIndex([pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = DatetimeIndex([pd.NaT, pd.NaT, pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

    def test_representation(self):
        idx1 = DatetimeIndex([], freq='D')
        idx2 = DatetimeIndex(['2011-01-01'], freq='D')
        idx3 = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'],
                             freq='H', tz='Asia/Tokyo')
        idx6 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT],
                             tz='US/Eastern')

        exp1 = """<class 'pandas.tseries.index.DatetimeIndex'>
Length: 0, Freq: D, Timezone: None"""
        exp2 = """<class 'pandas.tseries.index.DatetimeIndex'>
[2011-01-01]
Length: 1, Freq: D, Timezone: None"""
        exp3 = """<class 'pandas.tseries.index.DatetimeIndex'>
[2011-01-01, 2011-01-02]
Length: 2, Freq: D, Timezone: None"""
        exp4 = """<class 'pandas.tseries.index.DatetimeIndex'>
[2011-01-01, ..., 2011-01-03]
Length: 3, Freq: D, Timezone: None"""
        exp5 = """<class 'pandas.tseries.index.DatetimeIndex'>
[2011-01-01 09:00:00+09:00, ..., 2011-01-01 11:00:00+09:00]
Length: 3, Freq: H, Timezone: Asia/Tokyo"""
        exp6 = """<class 'pandas.tseries.index.DatetimeIndex'>
[2011-01-01 09:00:00-05:00, ..., NaT]
Length: 3, Freq: None, Timezone: US/Eastern"""

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5, idx6],
                                 [exp1, exp2, exp3, exp4, exp5, exp6]):
            for func in ['__repr__', '__unicode__', '__str__']:
                result = getattr(idx, func)()
                self.assertEqual(result, expected)

    def test_resolution(self):
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T', 'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day',
                                   'hour', 'minute', 'second', 'millisecond', 'microsecond']):
            for tz in [None, 'Asia/Tokyo', 'US/Eastern']:
                idx = pd.date_range(start='2013-04-01', periods=30, freq=freq, tz=tz)
                self.assertEqual(idx.resolution, expected)

    def test_add_iadd(self):
        for tz in self.tz:
            # union
            rng1 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other1 = pd.date_range('1/6/2000', freq='D', periods=5, tz=tz)
            expected1 = pd.date_range('1/1/2000', freq='D', periods=10, tz=tz)

            rng2 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other2 = pd.date_range('1/4/2000', freq='D', periods=5, tz=tz)
            expected2 = pd.date_range('1/1/2000', freq='D', periods=8, tz=tz)

            rng3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other3 = pd.DatetimeIndex([], tz=tz)
            expected3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)

            for rng, other, expected in [(rng1, other1, expected1), (rng2, other2, expected2),
                                         (rng3, other3, expected3)]:
                result_add = rng + other
                result_union = rng.union(other)

                tm.assert_index_equal(result_add, expected)
                tm.assert_index_equal(result_union, expected)
                rng += other
                tm.assert_index_equal(rng, expected)

            # offset
            offsets = [pd.offsets.Hour(2), timedelta(hours=2), np.timedelta64(2, 'h')]

            for delta in offsets:
                rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
                result = rng + delta
                expected = pd.date_range('2000-01-01 02:00', '2000-02-01 02:00', tz=tz)
                tm.assert_index_equal(result, expected)
                rng += delta
                tm.assert_index_equal(rng, expected)

            # int
            rng = pd.date_range('2000-01-01 09:00', freq='H', periods=10, tz=tz)
            result = rng + 1
            expected = pd.date_range('2000-01-01 10:00', freq='H', periods=10, tz=tz)
            tm.assert_index_equal(result, expected)
            rng += 1
            tm.assert_index_equal(rng, expected)

    def test_sub_isub(self):
        for tz in self.tz:
            # diff
            rng1 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other1 = pd.date_range('1/6/2000', freq='D', periods=5, tz=tz)
            expected1 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)

            rng2 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other2 = pd.date_range('1/4/2000', freq='D', periods=5, tz=tz)
            expected2 = pd.date_range('1/1/2000', freq='D', periods=3, tz=tz)

            rng3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other3 = pd.DatetimeIndex([], tz=tz)
            expected3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)

            for rng, other, expected in [(rng1, other1, expected1), (rng2, other2, expected2),
                                         (rng3, other3, expected3)]:
                result_add = rng - other
                result_union = rng.diff(other)

                tm.assert_index_equal(result_add, expected)
                tm.assert_index_equal(result_union, expected)
                rng -= other
                tm.assert_index_equal(rng, expected)

            # offset
            offsets = [pd.offsets.Hour(2), timedelta(hours=2), np.timedelta64(2, 'h')]

            for delta in offsets:
                rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
                result = rng - delta
                expected = pd.date_range('1999-12-31 22:00', '2000-01-31 22:00', tz=tz)
                tm.assert_index_equal(result, expected)
                rng -= delta
                tm.assert_index_equal(rng, expected)

            # int
            rng = pd.date_range('2000-01-01 09:00', freq='H', periods=10, tz=tz)
            result = rng - 1
            expected = pd.date_range('2000-01-01 08:00', freq='H', periods=10, tz=tz)
            tm.assert_index_equal(result, expected)
            rng -= 1
            tm.assert_index_equal(rng, expected)

    def test_value_counts_unique(self):
        # GH 7735
        for tz in [None, 'UTC', 'Asia/Tokyo', 'US/Eastern']:
            idx = pd.date_range('2011-01-01 09:00', freq='H', periods=10)
            # create repeated values, 'n'th element is repeated by n+1 times
            idx = DatetimeIndex(np.repeat(idx.values, range(1, len(idx) + 1)), tz=tz)

            exp_idx = pd.date_range('2011-01-01 18:00', freq='-1H', periods=10, tz=tz)
            expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')
            tm.assert_series_equal(idx.value_counts(), expected)

            expected = pd.date_range('2011-01-01 09:00', freq='H', periods=10, tz=tz)
            tm.assert_index_equal(idx.unique(), expected)

            idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 09:00', '2013-01-01 09:00',
                                 '2013-01-01 08:00', '2013-01-01 08:00', pd.NaT], tz=tz)

            exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00'], tz=tz)
            expected = Series([3, 2], index=exp_idx)
            tm.assert_series_equal(idx.value_counts(), expected)

            exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00', pd.NaT], tz=tz)
            expected = Series([3, 2, 1], index=exp_idx)
            tm.assert_series_equal(idx.value_counts(dropna=False), expected)

            tm.assert_index_equal(idx.unique(), exp_idx)


class TestPeriodIndexOps(Ops):

    def setUp(self):
        super(TestPeriodIndexOps, self).setUp()
        mask = lambda x: isinstance(x, DatetimeIndex) or isinstance(x, PeriodIndex) or is_datetimelike(x)
        self.is_valid_objs  = [ o for o in self.objs if mask(o) ]
        self.not_valid_objs = [ o for o in self.objs if not mask(o) ]

    def test_ops_properties(self):
        self.check_ops_properties(['year','month','day','hour','minute','second','weekofyear','week','dayofweek','dayofyear','quarter'])
        self.check_ops_properties(['qyear'], lambda x: isinstance(x,PeriodIndex))

    def test_asobject_tolist(self):
        idx = pd.period_range(start='2013-01-01', periods=4, freq='M', name='idx')
        expected_list = [pd.Period('2013-01-31', freq='M'), pd.Period('2013-02-28', freq='M'),
                         pd.Period('2013-03-31', freq='M'), pd.Period('2013-04-30', freq='M')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = PeriodIndex(['2013-01-01', '2013-01-02', 'NaT', '2013-01-04'], freq='D', name='idx')
        expected_list = [pd.Period('2013-01-01', freq='D'), pd.Period('2013-01-02', freq='D'),
                         pd.Period('NaT', freq='D'), pd.Period('2013-01-04', freq='D')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        for i in [0, 1, 3]:
            self.assertTrue(result[i], expected[i])
        self.assertTrue(result[2].ordinal, pd.tslib.iNaT)
        self.assertTrue(result[2].freq, 'D')
        self.assertEqual(result.name, expected.name)

        result_list = idx.tolist()
        for i in [0, 1, 3]:
            self.assertTrue(result_list[i], expected_list[i])
        self.assertTrue(result_list[2].ordinal, pd.tslib.iNaT)
        self.assertTrue(result_list[2].freq, 'D')

    def test_minmax(self):

        # monotonic
        idx1 = pd.PeriodIndex([pd.NaT, '2011-01-01', '2011-01-02',
                               '2011-01-03'], freq='D')
        self.assertTrue(idx1.is_monotonic)

        # non-monotonic
        idx2 = pd.PeriodIndex(['2011-01-01', pd.NaT, '2011-01-03',
                                '2011-01-02', pd.NaT], freq='D')
        self.assertFalse(idx2.is_monotonic)

        for idx in [idx1, idx2]:
            self.assertEqual(idx.min(), pd.Period('2011-01-01', freq='D'))
            self.assertEqual(idx.max(), pd.Period('2011-01-03', freq='D'))

        for op in ['min', 'max']:
            # Return NaT
            obj = PeriodIndex([], freq='M')
            result = getattr(obj, op)()
            self.assertEqual(result.ordinal, tslib.iNaT)
            self.assertEqual(result.freq, 'M')

            obj = PeriodIndex([pd.NaT], freq='M')
            result = getattr(obj, op)()
            self.assertEqual(result.ordinal, tslib.iNaT)
            self.assertEqual(result.freq, 'M')

            obj = PeriodIndex([pd.NaT, pd.NaT, pd.NaT], freq='M')
            result = getattr(obj, op)()
            self.assertEqual(result.ordinal, tslib.iNaT)
            self.assertEqual(result.freq, 'M')

    def test_representation(self):
        # GH 7601
        idx1 = PeriodIndex([], freq='D')
        idx2 = PeriodIndex(['2011-01-01'], freq='D')
        idx3 = PeriodIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = PeriodIndex(['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = PeriodIndex(['2011', '2012', '2013'], freq='A')
        idx6 = PeriodIndex(['2011-01-01 09:00', '2012-02-01 10:00', 'NaT'], freq='H')

        idx7 = pd.period_range('2013Q1', periods=1, freq="Q")
        idx8 = pd.period_range('2013Q1', periods=2, freq="Q")
        idx9 = pd.period_range('2013Q1', periods=3, freq="Q")

        exp1 = """<class 'pandas.tseries.period.PeriodIndex'>
Length: 0, Freq: D"""
        exp2 = """<class 'pandas.tseries.period.PeriodIndex'>
[2011-01-01]
Length: 1, Freq: D"""
        exp3 = """<class 'pandas.tseries.period.PeriodIndex'>
[2011-01-01, 2011-01-02]
Length: 2, Freq: D"""
        exp4 = """<class 'pandas.tseries.period.PeriodIndex'>
[2011-01-01, ..., 2011-01-03]
Length: 3, Freq: D"""
        exp5 = """<class 'pandas.tseries.period.PeriodIndex'>
[2011, ..., 2013]
Length: 3, Freq: A-DEC"""
        exp6 = """<class 'pandas.tseries.period.PeriodIndex'>
[2011-01-01 09:00, ..., NaT]
Length: 3, Freq: H"""
        exp7 = """<class 'pandas.tseries.period.PeriodIndex'>
[2013Q1]
Length: 1, Freq: Q-DEC"""
        exp8 = """<class 'pandas.tseries.period.PeriodIndex'>
[2013Q1, 2013Q2]
Length: 2, Freq: Q-DEC"""
        exp9 = """<class 'pandas.tseries.period.PeriodIndex'>
[2013Q1, ..., 2013Q3]
Length: 3, Freq: Q-DEC"""

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx9],
                                 [exp1, exp2, exp3, exp4, exp5, exp6, exp7, exp8, exp9]):
            for func in ['__repr__', '__unicode__', '__str__']:
                result = getattr(idx, func)()
                self.assertEqual(result, expected)

    def test_resolution(self):
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T', 'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day',
                                   'hour', 'minute', 'second', 'millisecond', 'microsecond']):

            idx = pd.period_range(start='2013-04-01', periods=30, freq=freq)
            self.assertEqual(idx.resolution, expected)

    def test_add_iadd(self):
        # union
        rng1 = pd.period_range('1/1/2000', freq='D', periods=5)
        other1 = pd.period_range('1/6/2000', freq='D', periods=5)
        expected1 = pd.period_range('1/1/2000', freq='D', periods=10)

        rng2 = pd.period_range('1/1/2000', freq='D', periods=5)
        other2 = pd.period_range('1/4/2000', freq='D', periods=5)
        expected2 = pd.period_range('1/1/2000', freq='D', periods=8)

        rng3 = pd.period_range('1/1/2000', freq='D', periods=5)
        other3 = pd.PeriodIndex([], freq='D')
        expected3 = pd.period_range('1/1/2000', freq='D', periods=5)

        rng4 = pd.period_range('2000-01-01 09:00', freq='H', periods=5)
        other4 = pd.period_range('2000-01-02 09:00', freq='H', periods=5)
        expected4 = pd.PeriodIndex(['2000-01-01 09:00', '2000-01-01 10:00',
                                    '2000-01-01 11:00', '2000-01-01 12:00',
                                    '2000-01-01 13:00', '2000-01-02 09:00',
                                    '2000-01-02 10:00', '2000-01-02 11:00',
                                    '2000-01-02 12:00', '2000-01-02 13:00'],
                                   freq='H')

        rng5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:03',
                               '2000-01-01 09:05'], freq='T')
        other5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:05'
                                 '2000-01-01 09:08'], freq='T')
        expected5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:03',
                                    '2000-01-01 09:05', '2000-01-01 09:08'],
                                   freq='T')

        rng6 = pd.period_range('2000-01-01', freq='M', periods=7)
        other6 = pd.period_range('2000-04-01', freq='M', periods=7)
        expected6 = pd.period_range('2000-01-01', freq='M', periods=10)

        rng7 = pd.period_range('2003-01-01', freq='A', periods=5)
        other7 = pd.period_range('1998-01-01', freq='A', periods=8)
        expected7 = pd.period_range('1998-01-01', freq='A', periods=10)

        for rng, other, expected in [(rng1, other1, expected1), (rng2, other2, expected2),
                                     (rng3, other3, expected3), (rng4, other4, expected4),
                                     (rng5, other5, expected5), (rng6, other6, expected6),
                                     (rng7, other7, expected7)]:

            result_add = rng + other
            result_union = rng.union(other)

            tm.assert_index_equal(result_add, expected)
            tm.assert_index_equal(result_union, expected)
            # GH 6527
            rng += other
            tm.assert_index_equal(rng, expected)

        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng + pd.offsets.YearEnd(5)
        expected = pd.period_range('2019', '2029', freq='A')
        tm.assert_index_equal(result, expected)
        rng += pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1), pd.offsets.Minute(),
                  np.timedelta64(365, 'D'), timedelta(365)]:
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng + o

        rng = pd.period_range('2014-01', '2016-12', freq='M')
        result = rng + pd.offsets.MonthEnd(5)
        expected = pd.period_range('2014-06', '2017-05', freq='M')
        tm.assert_index_equal(result, expected)
        rng += pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1), pd.offsets.Minute(),
                  np.timedelta64(365, 'D'), timedelta(365)]:
            rng = pd.period_range('2014-01', '2016-12', freq='M')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng + o

        # Tick
        offsets = [pd.offsets.Day(3), timedelta(days=3), np.timedelta64(3, 'D'),
                   pd.offsets.Hour(72), timedelta(minutes=60*24*3), np.timedelta64(72, 'h')]
        for delta in offsets:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            result = rng + delta
            expected = pd.period_range('2014-05-04', '2014-05-18', freq='D')
            tm.assert_index_equal(result, expected)
            rng += delta
            tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1), pd.offsets.Minute(),
                  np.timedelta64(4, 'h'), timedelta(hours=23)]:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng + o

        offsets = [pd.offsets.Hour(2), timedelta(hours=2), np.timedelta64(2, 'h'),
                   pd.offsets.Minute(120), timedelta(minutes=120), np.timedelta64(120, 'm')]
        for delta in offsets:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
            result = rng + delta
            expected = pd.period_range('2014-01-01 12:00', '2014-01-05 12:00', freq='H')
            tm.assert_index_equal(result, expected)
            rng += delta
            tm.assert_index_equal(rng, expected)

        for delta in [pd.offsets.YearBegin(2), timedelta(minutes=30), np.timedelta64(30, 's')]:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                result = rng + delta
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng += delta

        # int
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng + 1
        expected = pd.period_range('2000-01-01 10:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng += 1
        tm.assert_index_equal(rng, expected)

    def test_sub_isub(self):
        # diff
        rng1 = pd.period_range('1/1/2000', freq='D', periods=5)
        other1 = pd.period_range('1/6/2000', freq='D', periods=5)
        expected1 = pd.period_range('1/1/2000', freq='D', periods=5)

        rng2 = pd.period_range('1/1/2000', freq='D', periods=5)
        other2 = pd.period_range('1/4/2000', freq='D', periods=5)
        expected2 = pd.period_range('1/1/2000', freq='D', periods=3)

        rng3 = pd.period_range('1/1/2000', freq='D', periods=5)
        other3 = pd.PeriodIndex([], freq='D')
        expected3 = pd.period_range('1/1/2000', freq='D', periods=5)

        rng4 = pd.period_range('2000-01-01 09:00', freq='H', periods=5)
        other4 = pd.period_range('2000-01-02 09:00', freq='H', periods=5)
        expected4 = rng4

        rng5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:03',
                               '2000-01-01 09:05'], freq='T')
        other5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:05'], freq='T')
        expected5 = pd.PeriodIndex(['2000-01-01 09:03'], freq='T')

        rng6 = pd.period_range('2000-01-01', freq='M', periods=7)
        other6 = pd.period_range('2000-04-01', freq='M', periods=7)
        expected6 = pd.period_range('2000-01-01', freq='M', periods=3)

        rng7 = pd.period_range('2003-01-01', freq='A', periods=5)
        other7 = pd.period_range('1998-01-01', freq='A', periods=8)
        expected7 = pd.period_range('2006-01-01', freq='A', periods=2)

        for rng, other, expected in [(rng1, other1, expected1), (rng2, other2, expected2),
                                     (rng3, other3, expected3), (rng4, other4, expected4),
                                     (rng5, other5, expected5), (rng6, other6, expected6),
                                     (rng7, other7, expected7),]:
            result_add = rng - other
            result_union = rng.diff(other)

            tm.assert_index_equal(result_add, expected)
            tm.assert_index_equal(result_union, expected)
            rng -= other
            tm.assert_index_equal(rng, expected)

        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng - pd.offsets.YearEnd(5)
        expected = pd.period_range('2009', '2019', freq='A')
        tm.assert_index_equal(result, expected)
        rng -= pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1), pd.offsets.Minute(),
                  np.timedelta64(365, 'D'), timedelta(365)]:
            rng = pd.period_range('2014', '2024', freq='A')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng - o

        rng = pd.period_range('2014-01', '2016-12', freq='M')
        result = rng - pd.offsets.MonthEnd(5)
        expected = pd.period_range('2013-08', '2016-07', freq='M')
        tm.assert_index_equal(result, expected)
        rng -= pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1), pd.offsets.Minute(),
                  np.timedelta64(365, 'D'), timedelta(365)]:
            rng = pd.period_range('2014-01', '2016-12', freq='M')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng - o

        # Tick
        offsets = [pd.offsets.Day(3), timedelta(days=3), np.timedelta64(3, 'D'),
                   pd.offsets.Hour(72), timedelta(minutes=60*24*3), np.timedelta64(72, 'h')]
        for delta in offsets:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            result = rng - delta
            expected = pd.period_range('2014-04-28', '2014-05-12', freq='D')
            tm.assert_index_equal(result, expected)
            rng -= delta
            tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1), pd.offsets.Minute(),
                  np.timedelta64(4, 'h'), timedelta(hours=23)]:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng - o

        offsets = [pd.offsets.Hour(2), timedelta(hours=2), np.timedelta64(2, 'h'),
                   pd.offsets.Minute(120), timedelta(minutes=120), np.timedelta64(120, 'm')]
        for delta in offsets:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
            result = rng - delta
            expected = pd.period_range('2014-01-01 08:00', '2014-01-05 08:00', freq='H')
            tm.assert_index_equal(result, expected)
            rng -= delta
            tm.assert_index_equal(rng, expected)

        for delta in [pd.offsets.YearBegin(2), timedelta(minutes=30), np.timedelta64(30, 's')]:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                result = rng + delta
            with tm.assertRaisesRegexp(ValueError, 'Input has different freq from Period'):
                rng += delta

        # int
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng - 1
        expected = pd.period_range('2000-01-01 08:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng -= 1
        tm.assert_index_equal(rng, expected)

    def test_value_counts_unique(self):
        # GH 7735
        idx = pd.period_range('2011-01-01 09:00', freq='H', periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = PeriodIndex(np.repeat(idx.values, range(1, len(idx) + 1)), freq='H')

        exp_idx = PeriodIndex(['2011-01-01 18:00', '2011-01-01 17:00', '2011-01-01 16:00',
                               '2011-01-01 15:00', '2011-01-01 14:00', '2011-01-01 13:00',
                               '2011-01-01 12:00', '2011-01-01 11:00', '2011-01-01 10:00',
                               '2011-01-01 09:00'], freq='H')
        expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')
        tm.assert_series_equal(idx.value_counts(), expected)

        expected = pd.period_range('2011-01-01 09:00', freq='H', periods=10)
        tm.assert_index_equal(idx.unique(), expected)

        idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 09:00', '2013-01-01 09:00',
                           '2013-01-01 08:00', '2013-01-01 08:00', pd.NaT], freq='H')

        exp_idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 08:00'], freq='H')
        expected = Series([3, 2], index=exp_idx)
        tm.assert_series_equal(idx.value_counts(), expected)

        exp_idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 08:00', pd.NaT], freq='H')
        expected = Series([3, 2, 1], index=exp_idx)
        tm.assert_series_equal(idx.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)


if __name__ == '__main__':
    import nose

    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
