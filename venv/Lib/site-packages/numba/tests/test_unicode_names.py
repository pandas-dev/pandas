# -*- coding: utf-8 -*-


from numba import njit, cfunc
from numba.tests.support import TestCase, unittest
from numba.core import cgutils

unicode_name1 = u"""
def unicode_name1(ಠ_ರೃ, ಠਊಠ):
    return (ಠ_ರೃ) + (ಠਊಠ)
"""

unicode_name2 = u"""
def Ծ_Ծ(ಠ_ರೃ, ಠਊಠ):
    return (ಠ_ರೃ) + (ಠਊಠ)
"""


class TestUnicodeNames(TestCase):
    def make_testcase(self, src, fname):
        glb = {}
        exec(src, glb)
        fn = glb[fname]
        return fn

    def test_unicode_name1(self):
        fn = self.make_testcase(unicode_name1, 'unicode_name1')
        cfn = njit(fn)
        self.assertEqual(cfn(1, 2), 3)

    def test_unicode_name2(self):
        fn = self.make_testcase(unicode_name2, 'Ծ_Ծ')
        cfn = njit(fn)
        self.assertEqual(cfn(1, 2), 3)

    def test_cfunc(self):
        fn = self.make_testcase(unicode_name2, 'Ծ_Ծ')
        cfn = cfunc("int32(int32, int32)")(fn)
        self.assertEqual(cfn.ctypes(1, 2), 3)


class TestUnicodeUtils(TestCase):
    def test_normalize_ir_text(self):
        # non-unicode input
        out = cgutils.normalize_ir_text('abc')
        # str returned
        self.assertIsInstance(out, str)
        # try encoding to latin
        out.encode('latin1')

    def test_normalize_ir_text_unicode(self):
        # unicode input
        out = cgutils.normalize_ir_text(unicode_name2)
        # str returned
        self.assertIsInstance(out, str)
        # try encoding to latin
        out.encode('latin1')


if __name__ == '__main__':
    unittest.main()
