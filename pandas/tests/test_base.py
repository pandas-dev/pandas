import re
import numpy as np
import pandas.compat as compat
from pandas.compat import u
from pandas.core.base import FrozenList, FrozenNDArray
from pandas.util.testing import assertRaisesRegexp, assert_isinstance
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
        self.assert_(not isinstance(self.container.view(np.ndarray), FrozenNDArray))
        self.assert_(self.container.view() is not self.container)
        self.assert_(np.array_equal(self.container, original))
        # shallow copy should be the same too
        assert_isinstance(self.container._shallow_copy(), FrozenNDArray)
        # setting should not be allowed
        def testit(container): container[0] = 16

        self.check_mutable_error(testit, self.container)

    def test_values(self):
        original = self.container.view(np.ndarray).copy()
        n = original[0] + 15
        vals = self.container.values()
        self.assert_(np.array_equal(original, vals))
        self.assert_(original is not vals)
        vals[0] = n
        self.assert_(np.array_equal(self.container, original))
        self.assertEqual(vals[0], n)


if __name__ == '__main__':
    import nose

    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
