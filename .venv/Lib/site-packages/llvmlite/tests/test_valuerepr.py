import math
import sys
import unittest

from llvmlite.ir import (
    Constant, FloatType, DoubleType, LiteralStructType, IntType,
    ArrayType, HalfType)
from llvmlite.tests import TestCase


int8 = IntType(8)
int16 = IntType(16)


PY36_OR_LATER = sys.version_info[:2] >= (3, 6)


class TestValueRepr(TestCase):

    def test_double_repr(self):
        def check_repr(val, expected):
            c = Constant(DoubleType(), val)
            self.assertEqual(str(c), expected)
        check_repr(math.pi, "double 0x400921fb54442d18")
        check_repr(float('inf'), "double 0x7ff0000000000000")
        check_repr(float('-inf'), "double 0xfff0000000000000")

    def test_float_repr(self):
        def check_repr(val, expected):
            c = Constant(FloatType(), val)
            self.assertEqual(str(c), expected)
        check_repr(math.pi, "float 0x400921fb60000000")
        check_repr(float('inf'), "float 0x7ff0000000000000")
        check_repr(float('-inf'), "float 0xfff0000000000000")

    @unittest.skipUnless(PY36_OR_LATER, 'py36+ only')
    def test_half_repr(self):
        def check_repr(val, expected):
            c = Constant(HalfType(), val)
            self.assertEqual(str(c), expected)
        check_repr(math.pi, "half 0x4009200000000000")
        check_repr(float('inf'), "half 0x7ff0000000000000")
        check_repr(float('-inf'), "half 0xfff0000000000000")

    def test_struct_repr(self):
        tp = LiteralStructType([int8, int16])
        c = Constant(tp, (Constant(int8, 100), Constant(int16, 1000)))
        self.assertEqual(str(c), "{i8, i16} {i8 100, i16 1000}")

    def test_array_repr(self):
        tp = ArrayType(int8, 3)
        values = [Constant(int8, x) for x in (5, 10, -15)]
        c = Constant(tp, values)
        self.assertEqual(str(c), "[3 x i8] [i8 5, i8 10, i8 -15]")
        c = Constant(tp, bytearray(b"\x01\x02\x03"))
        self.assertEqual(str(c), '[3 x i8] c"\\01\\02\\03"')


if __name__ == "__main__":
    unittest.main()
