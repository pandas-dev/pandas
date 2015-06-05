# -*- coding: utf-8 -*-
import warnings

import nose

import sys
import pandas.util
from pandas.util.decorators import deprecate_kwarg
import pandas.util.testing as tm

class TestDecorators(tm.TestCase):
    def setUp(self):
        @deprecate_kwarg('old', 'new')
        def _f1(new=False):
            return new

        @deprecate_kwarg('old', 'new', {'yes': True, 'no': False})
        def _f2(new=False):
            return new

        @deprecate_kwarg('old', 'new', lambda x: x+1)
        def _f3(new=0):
            return new

        self.f1 = _f1
        self.f2 = _f2
        self.f3 = _f3

    def test_deprecate_kwarg(self):
        x = 78
        with tm.assert_produces_warning(FutureWarning):
            result = self.f1(old=x)
        self.assertIs(result, x)
        with tm.assert_produces_warning(None):
            self.f1(new=x)

    def test_dict_deprecate_kwarg(self):
        x = 'yes'
        with tm.assert_produces_warning(FutureWarning):
            result = self.f2(old=x)
        self.assertEqual(result, True)

    def test_missing_deprecate_kwarg(self):
        x = 'bogus'
        with tm.assert_produces_warning(FutureWarning):
            result = self.f2(old=x)
        self.assertEqual(result, 'bogus')

    def test_callable_deprecate_kwarg(self):
        x = 5
        with tm.assert_produces_warning(FutureWarning):
            result = self.f3(old=x)
        self.assertEqual(result, x+1)
        with tm.assertRaises(TypeError):
            self.f3(old='hello')

    def test_bad_deprecate_kwarg(self):
        with tm.assertRaises(TypeError):
            @deprecate_kwarg('old', 'new', 0)
            def f4(new=None):
                pass


class TestTesting(tm.TestCase):

    def test_warning(self):

        with tm.assert_produces_warning(FutureWarning):
            self.assertEquals(1, 1)

        with tm.assert_produces_warning(FutureWarning):
            self.assertNotEquals(1, 2)

        with tm.assert_produces_warning(FutureWarning):
            self.assert_(True)

        with tm.assert_produces_warning(FutureWarning):
            self.assertAlmostEquals(1.0, 1.0000000001)

        with tm.assert_produces_warning(FutureWarning):
            self.assertNotAlmostEquals(1, 2)

    def test_locale(self):
        if sys.platform == 'win32':
            raise nose.SkipTest("skipping on win platforms as locale not available")

        #GH9744
        locales = pandas.util.testing.get_locales()
        self.assertTrue(len(locales) >= 1)

def test_rands():
    r = tm.rands(10)
    assert(len(r) == 10)


def test_rands_array():
    arr = tm.rands_array(5, size=10)
    assert(arr.shape == (10,))
    assert(len(arr[0]) == 5)

    arr = tm.rands_array(7, size=(10, 10))
    assert(arr.shape == (10, 10))
    assert(len(arr[1, 1]) == 7)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
