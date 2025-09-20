import unittest
from numba import jit
from numba.core import types


force_pyobj_flags = {'forceobj': True}
no_pyobj_flags = {'nopython': True}


def return_int(a, b):
    return a + b


class TestPythonInt(unittest.TestCase):

    # Issue #474: ints should be returned rather than longs under Python 2,
    # as much as possible.

    def test_int_return_type(self, flags=force_pyobj_flags,
                             int_type=types.int64, operands=(3, 4)):
        pyfunc = return_int
        cfunc = jit((int_type, int_type), **flags)(pyfunc)
        expected = pyfunc(*operands)
        got = cfunc(*operands)
        self.assertIs(type(got), type(expected))
        self.assertEqual(got, expected)

    def test_int_return_type_npm(self):
        self.test_int_return_type(flags=no_pyobj_flags)

    def test_unsigned_int_return_type(self, flags=force_pyobj_flags):
        self.test_int_return_type(int_type=types.uint64, flags=flags)

    def test_unsigned_int_return_type_npm(self):
        self.test_unsigned_int_return_type(flags=no_pyobj_flags)

    def test_long_int_return_type(self, flags=force_pyobj_flags):
        # Same but returning a 64-bit integer.  The return type should be
        # `int` on 64-bit builds, `long` on 32-bit ones (or Windows).
        self.test_int_return_type(flags=flags, operands=(2**33, 2**40))

    def test_long_int_return_type_npm(self):
        self.test_long_int_return_type(flags=no_pyobj_flags)

    def test_longer_int_return_type(self, flags=force_pyobj_flags):
        # This won't be supported in nopython mode.
        self.test_int_return_type(flags=flags, operands=(2**70, 2**75))


if __name__ == '__main__':
    unittest.main()
