import numpy as np

import unittest
from numba import njit
from numba.core.errors import TypingError
from numba import jit, typeof
from numba.core import types
from numba.tests.support import TestCase


a0 = np.array(42)

s1 = np.int32(64)

a1 = np.arange(12)
a2 = a1[::2]
a3 = a1.reshape((3, 4)).T

dt = np.dtype([('x', np.int8), ('y', 'S3')])

a4 = np.arange(32, dtype=np.int8).view(dt)
a5 = a4[::-2]

# A recognizable data string
a6 = np.frombuffer(b"XXXX_array_contents_XXXX", dtype=np.float32)


myarray = np.array([1, ])


def getitem0(i):
    return a0[()]


def getitem1(i):
    return a1[i]


def getitem2(i):
    return a2[i]


def getitem3(i):
    return a3[i]


def getitem4(i):
    return a4[i]


def getitem5(i):
    return a5[i]


def getitem6(i):
    return a6[i]


def use_arrayscalar_const():
    return s1


def write_to_global_array():
    myarray[0] = 1


def bytes_as_const_array():
    return np.frombuffer(b'foo', dtype=np.uint8)


class TestConstantArray(TestCase):
    """
    Test array constants.
    """

    def check_array_const(self, pyfunc):
        cfunc = njit((types.int32,))(pyfunc)
        for i in [0, 1, 2]:
            np.testing.assert_array_equal(pyfunc(i), cfunc(i))

    def test_array_const_0d(self):
        self.check_array_const(getitem0)

    def test_array_const_1d_contig(self):
        self.check_array_const(getitem1)

    def test_array_const_1d_noncontig(self):
        self.check_array_const(getitem2)

    def test_array_const_2d(self):
        self.check_array_const(getitem3)

    def test_record_array_const_contig(self):
        self.check_array_const(getitem4)

    def test_record_array_const_noncontig(self):
        self.check_array_const(getitem5)

    def test_array_const_alignment(self):
        """
        Issue #1933: the array declaration in the LLVM IR must have
        the right alignment specified.
        """
        sig = (types.intp,)
        cfunc = jit(sig, nopython=True)(getitem6)
        ir = cfunc.inspect_llvm(sig)
        for line in ir.splitlines():
            if 'XXXX_array_contents_XXXX' in line:
                self.assertIn("constant [24 x i8]", line)  # sanity check
                # Should be the ABI-required alignment for float32
                # on most platforms...
                self.assertIn(", align 4", line)
                break
        else:
            self.fail("could not find array declaration in LLVM IR")

    def test_arrayscalar_const(self):
        pyfunc = use_arrayscalar_const
        cfunc = njit((),)(pyfunc)
        self.assertEqual(pyfunc(), cfunc())

    def test_write_to_global_array(self):
        pyfunc = write_to_global_array
        with self.assertRaises(TypingError):
            njit((),)(pyfunc)

    def test_issue_1850(self):
        """
        This issue is caused by an unresolved bug in numpy since version 1.6.
        See numpy GH issue #3147.
        """
        constarr = np.array([86])

        def pyfunc():
            return constarr[0]

        cfunc = njit((),)(pyfunc)
        out = cfunc()
        self.assertEqual(out, 86)

    @TestCase.run_test_in_subprocess # isolate MCJIT use
    def test_too_big_to_freeze(self):
        """
        Test issue https://github.com/numba/numba/issues/2188 where freezing
        a constant array into the code that's prohibitively long and consumes
        too much RAM.
        """
        def test(biggie):
            expect = np.copy(biggie)
            self.assertEqual(typeof(biggie), typeof(expect))

            def pyfunc():
                return biggie

            cfunc = njit((),)(pyfunc)
            # Check that the array is not frozen into the LLVM IR.
            # LLVM size must be less than the array size.
            self.assertLess(len(cfunc.inspect_llvm((),)), biggie.nbytes)
            # Run and test result
            out = cfunc()
            self.assertIs(biggie, out)
            # Remove all local references to biggie
            del out
            biggie = None  # del biggie is syntax error in py2
            # Run again and verify result
            out = cfunc()
            np.testing.assert_equal(expect, out)
            self.assertEqual(typeof(expect), typeof(out))

        nelem = 10**7   # 10 million items

        c_array = np.arange(nelem).reshape(nelem)
        f_array = np.asfortranarray(np.random.random((2, nelem // 2)))
        self.assertEqual(typeof(c_array).layout, 'C')
        self.assertEqual(typeof(f_array).layout, 'F')
        # Test C contig
        test(c_array)
        # Test F contig
        test(f_array)


class TestConstantBytes(TestCase):
    def test_constant_bytes(self):
        pyfunc = bytes_as_const_array
        cfunc = njit((),)(pyfunc)
        np.testing.assert_array_equal(pyfunc(), cfunc())


if __name__ == '__main__':
    unittest.main()
