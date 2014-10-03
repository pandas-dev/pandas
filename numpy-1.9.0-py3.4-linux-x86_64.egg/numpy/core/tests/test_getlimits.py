""" Test functions for limits module.

"""
from __future__ import division, absolute_import, print_function

from numpy.testing import *

from numpy.core import finfo, iinfo
from numpy import half, single, double, longdouble
import numpy as np

##################################################

class TestPythonFloat(TestCase):
    def test_singleton(self):
        ftype = finfo(float)
        ftype2 = finfo(float)
        assert_equal(id(ftype), id(ftype2))

class TestHalf(TestCase):
    def test_singleton(self):
        ftype = finfo(half)
        ftype2 = finfo(half)
        assert_equal(id(ftype), id(ftype2))

class TestSingle(TestCase):
    def test_singleton(self):
        ftype = finfo(single)
        ftype2 = finfo(single)
        assert_equal(id(ftype), id(ftype2))

class TestDouble(TestCase):
    def test_singleton(self):
        ftype = finfo(double)
        ftype2 = finfo(double)
        assert_equal(id(ftype), id(ftype2))

class TestLongdouble(TestCase):
    def test_singleton(self,level=2):
        ftype = finfo(longdouble)
        ftype2 = finfo(longdouble)
        assert_equal(id(ftype), id(ftype2))

class TestIinfo(TestCase):
    def test_basic(self):
        dts = list(zip(['i1', 'i2', 'i4', 'i8',
                   'u1', 'u2', 'u4', 'u8'],
                  [np.int8, np.int16, np.int32, np.int64,
                   np.uint8, np.uint16, np.uint32, np.uint64]))
        for dt1, dt2 in dts:
            assert_equal(iinfo(dt1).min, iinfo(dt2).min)
            assert_equal(iinfo(dt1).max, iinfo(dt2).max)
        self.assertRaises(ValueError, iinfo, 'f4')

    def test_unsigned_max(self):
        types = np.sctypes['uint']
        for T in types:
            assert_equal(iinfo(T).max, T(-1))

class TestRepr(TestCase):
    def test_iinfo_repr(self):
        expected = "iinfo(min=-32768, max=32767, dtype=int16)"
        assert_equal(repr(np.iinfo(np.int16)), expected)

    def test_finfo_repr(self):
        expected = "finfo(resolution=1e-06, min=-3.4028235e+38," + \
                   " max=3.4028235e+38, dtype=float32)"
        # Python 2.5 float formatting on Windows adds an extra 0 to the
        # exponent.  So test for both.  Once 2.5 compatibility is dropped, this
        # can simply use `assert_equal(repr(np.finfo(np.float32)), expected)`.
        expected_win25 = "finfo(resolution=1e-006, min=-3.4028235e+038," + \
                         " max=3.4028235e+038, dtype=float32)"

        actual = repr(np.finfo(np.float32))
        if not actual == expected:
            if not actual == expected_win25:
                msg = build_err_msg([actual, desired], verbose=True)
                raise AssertionError(msg)


def test_instances():
    iinfo(10)
    finfo(3.0)

if __name__ == "__main__":
    run_module_suite()
