import re
import numpy as np
import pandas.compat as compat
import pandas as pd
from pandas.compat import u
from pandas.core.base import FrozenList, FrozenNDArray
from pandas.util.testing import assertRaisesRegexp, assert_isinstance
from pandas import Series, Index, DatetimeIndex, PeriodIndex
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
        self.period_index  = tm.makePeriodIndex(10)
        self.string_index  = tm.makeStringIndex(10)

        arr = np.random.randn(10)
        self.int_series    = Series(arr, index=self.int_index)
        self.float_series  = Series(arr, index=self.int_index)
        self.dt_series     = Series(arr, index=self.dt_index)
        self.period_series = Series(arr, index=self.period_index)
        self.string_series = Series(arr, index=self.string_index)

        self.objs = [ getattr(self,"{0}_{1}".format(t,f)) for t in ['int','float','dt','period','string'] for f in ['index','series'] ]

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
                    self.assertRaises(TypeError, lambda : getattr(o,op))

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
                self.assertEqual(result, expected)

class TestDatetimeIndexOps(Ops):
    _allowed = '_allow_datetime_index_ops'

    def setUp(self):
        super(TestDatetimeIndexOps, self).setUp()
        mask = lambda x: x._allow_datetime_index_ops or x._allow_period_index_ops
        self.is_valid_objs  = [ o for o in self.objs if mask(o) ]
        self.not_valid_objs = [ o for o in self.objs if not mask(o) ]

    def test_ops_properties(self):
        self.check_ops_properties(['year','month','day','hour','minute','second','weekofyear','week','dayofweek','dayofyear','quarter'])
        self.check_ops_properties(['date','time','microsecond','nanosecond'], lambda x: isinstance(x,DatetimeIndex))

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
