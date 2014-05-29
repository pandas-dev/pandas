import re
from datetime import datetime, timedelta
import numpy as np
import pandas.compat as compat
import pandas as pd
from pandas.compat import u, StringIO
from pandas.core.base import FrozenList, FrozenNDArray, DatetimeIndexOpsMixin
from pandas.util.testing import assertRaisesRegexp, assert_isinstance
from pandas import Series, Index, Int64Index, DatetimeIndex, PeriodIndex
from pandas import _np_version_under1p7
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

    def test_ops(self):
        if _np_version_under1p7:
            raise nose.SkipTest("test only valid in numpy >= 1.7")
        for op in ['max','min']:
            for o in self.objs:
                result = getattr(o,op)()
                expected = getattr(o.values,op)()
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

            # explicitly create DatetimeIndex
            obj = DatetimeIndex([])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = DatetimeIndex([pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = DatetimeIndex([pd.NaT, pd.NaT, pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))


    def test_value_counts_unique_nunique(self):
        for o in self.objs:
            klass = type(o)
            values = o.values

            # create repeated values, 'n'th element is repeated by n+1 times
            if isinstance(o, PeriodIndex):
                # freq must be specified because repeat makes freq ambiguous
                o = klass(np.repeat(values, range(1, len(o) + 1)), freq=o.freq)
            else:
                o = klass(np.repeat(values, range(1, len(o) + 1)))

            expected_s = Series(range(10, 0, -1), index=values[::-1], dtype='int64')
            tm.assert_series_equal(o.value_counts(), expected_s)

            if isinstance(o, DatetimeIndex):
                # DatetimeIndex.unique returns DatetimeIndex
                self.assertTrue(o.unique().equals(klass(values)))
            else:
                self.assert_numpy_array_equal(o.unique(), values)

            self.assertEqual(o.nunique(), len(np.unique(o.values)))

        for null_obj in [np.nan, None]:
            for o in self.objs:
                klass = type(o)
                values = o.values

                if o.values.dtype == 'int64':
                    # skips int64 because it doesn't allow to include nan or None
                    continue

                if o.values.dtype == 'datetime64[ns]' and _np_version_under1p7:
                    # Unable to assign None
                    continue

                # special assign to the numpy array
                if o.values.dtype == 'datetime64[ns]':
                    values[0:2] = pd.tslib.iNaT
                else:
                    values[0:2] = null_obj

                # create repeated values, 'n'th element is repeated by n+1 times
                if isinstance(o, PeriodIndex):
                    o = klass(np.repeat(values, range(1, len(o) + 1)), freq=o.freq)
                else:
                    o = klass(np.repeat(values, range(1, len(o) + 1)))

                if isinstance(o, DatetimeIndex):
                    # DatetimeIndex: nan is casted to Nat and included
                    expected_s = Series(list(range(10, 2, -1)) + [3], index=values[9:0:-1])
                else:
                    # nan is excluded
                    expected_s = Series(range(10, 2, -1), index=values[9:1:-1], dtype='int64')

                tm.assert_series_equal(o.value_counts(), expected_s)

                # numpy_array_equal cannot compare arrays includes nan
                result = o.unique()
                self.assert_numpy_array_equal(result[1:], values[2:])

                if isinstance(o, DatetimeIndex):
                    self.assertTrue(result[0] is pd.NaT)
                else:
                    self.assertTrue(pd.isnull(result[0]))

                if isinstance(o, DatetimeIndex):
                    self.assertEqual(o.nunique(), 9)
                else:
                    self.assertEqual(o.nunique(), 8)

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
            expected_s[pd.NaT] = 1
            tm.assert_series_equal(result, expected_s)

            unique = s.unique()
            self.assertEqual(unique.dtype, 'datetime64[ns]')
            # numpy_array_equal cannot compare pd.NaT
            self.assert_numpy_array_equal(unique[:3], expected)
            self.assertTrue(unique[3] is pd.NaT or unique[3].astype('int64') == pd.tslib.iNaT)

            self.assertEqual(s.nunique(), 4)

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


class TestDatetimeIndexOps(Ops):
    _allowed = '_allow_datetime_index_ops'

    def setUp(self):
        super(TestDatetimeIndexOps, self).setUp()
        mask = lambda x: x._allow_datetime_index_ops or x._allow_period_index_ops
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

class TestPeriodIndexOps(Ops):
    _allowed = '_allow_period_index_ops'

    def setUp(self):
        super(TestPeriodIndexOps, self).setUp()
        mask = lambda x: x._allow_datetime_index_ops or x._allow_period_index_ops
        self.is_valid_objs  = [ o for o in self.objs if mask(o) ]
        self.not_valid_objs = [ o for o in self.objs if not mask(o) ]

    def test_ops_properties(self):
        self.check_ops_properties(['year','month','day','hour','minute','second','weekofyear','week','dayofweek','dayofyear','quarter'])
        self.check_ops_properties(['qyear'], lambda x: isinstance(x,PeriodIndex))

if __name__ == '__main__':
    import nose

    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
