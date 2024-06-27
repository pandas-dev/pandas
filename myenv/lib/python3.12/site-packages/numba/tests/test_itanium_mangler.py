# -*- coding: utf-8 -*-
from numba import int32, int64, uint32, uint64, float32, float64
from numba.core.types import range_iter32_type
from numba.core import itanium_mangler
import unittest


class TestItaniumManager(unittest.TestCase):
    def test_ident(self):
        got = itanium_mangler.mangle_identifier("apple")
        expect = "5apple"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_identifier("ap_ple")
        expect = "6ap_ple"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_identifier("apple213")
        expect = "8apple213"
        self.assertEqual(expect, got)

    def test_types(self):
        got = itanium_mangler.mangle_type(int32)
        expect = "i"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_type(int64)
        expect = "x"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_type(uint32)
        expect = "j"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_type(uint64)
        expect = "y"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_type(float32)
        expect = "f"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle_type(float64)
        expect = "d"
        self.assertEqual(expect, got)

    def test_function(self):
        got = itanium_mangler.mangle("what", [int32, float32])
        expect = "_Z4whatif"
        self.assertEqual(expect, got)

        got = itanium_mangler.mangle("a_little_brown_fox", [uint64,
                                                            uint32,
                                                            float64])
        expect = "_Z18a_little_brown_foxyjd"
        self.assertEqual(expect, got)

    def test_custom_type(self):
        got = itanium_mangler.mangle_type(range_iter32_type)
        name = str(range_iter32_type)
        expect = "{n}{name}".format(n=len(name), name=name)
        self.assertEqual(expect, got)

    def test_mangle_literal(self):
        # check int
        got = itanium_mangler.mangle_value(123)
        expect = "Li123E"
        self.assertEqual(expect, got)
        # check float (not handled using standard)
        got = itanium_mangler.mangle_value(12.3)
        self.assertRegex(got, r'^\d+_12_[0-9a-z][0-9a-z]3$')

    def test_mangle_unicode(self):
        name = u'f∂ƒ©z'
        got = itanium_mangler.mangle_identifier(name)
        self.assertRegex(got, r'^\d+f(_[a-z0-9][a-z0-9])+z$')


if __name__ == '__main__':
    unittest.main()
