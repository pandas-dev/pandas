import numpy as np

from numba import float32, jit, njit
from numba.np.ufunc import Vectorize
from numba.core.errors import TypingError
from numba.tests.support import TestCase
import unittest


dtype = np.float32
a = np.arange(80, dtype=dtype).reshape(8, 10)
b = a.copy()
c = a.copy(order='F')
d = np.arange(16 * 20, dtype=dtype).reshape(16, 20)[::2, ::2]


def add(a, b):
    return a + b


def add_multiple_args(a, b, c, d):
    return a + b + c + d


def gufunc_add(a, b):
    result = 0.0
    for i in range(a.shape[0]):
        result += a[i] * b[i]

    return result


def ufunc_reduce(ufunc, arg):
    for i in range(arg.ndim):
        arg = ufunc.reduce(arg)
    return arg


vectorizers = [
    Vectorize,
    # ParallelVectorize,
    # StreamVectorize,
    # CudaVectorize,
    # GUFuncVectorize,
]


class TestUFuncs(TestCase):

    def _test_ufunc_attributes(self, cls, a, b, *args):
        "Test ufunc attributes"
        vectorizer = cls(add, *args)
        vectorizer.add(float32(float32, float32))
        ufunc = vectorizer.build_ufunc()

        info = (cls, a.ndim)
        self.assertPreciseEqual(ufunc(a, b), a + b, msg=info)
        self.assertPreciseEqual(ufunc_reduce(ufunc, a), np.sum(a), msg=info)
        self.assertPreciseEqual(ufunc.accumulate(a), np.add.accumulate(a),
                                msg=info)
        self.assertPreciseEqual(ufunc.outer(a, b), np.add.outer(a, b), msg=info)

    def _test_broadcasting(self, cls, a, b, c, d):
        "Test multiple args"
        vectorizer = cls(add_multiple_args)
        vectorizer.add(float32(float32, float32, float32, float32))
        ufunc = vectorizer.build_ufunc()

        info = (cls, a.shape)
        self.assertPreciseEqual(ufunc(a, b, c, d), a + b + c + d, msg=info)

    def test_ufunc_attributes(self):
        for v in vectorizers: # 1D
            self._test_ufunc_attributes(v, a[0], b[0])
        for v in vectorizers: # 2D
            self._test_ufunc_attributes(v, a, b)
        for v in vectorizers: # 3D
            self._test_ufunc_attributes(v, a[:, np.newaxis, :],
                                        b[np.newaxis, :, :])

    def test_broadcasting(self):
        for v in vectorizers: # 1D
            self._test_broadcasting(v, a[0], b[0], c[0], d[0])
        for v in vectorizers: # 2D
            self._test_broadcasting(v, a, b, c, d)
        for v in vectorizers: # 3D
            self._test_broadcasting(v, a[:, np.newaxis, :], b[np.newaxis, :, :],
                                    c[:, np.newaxis, :], d[np.newaxis, :, :])

    def test_implicit_broadcasting(self):
        for v in vectorizers:
            vectorizer = v(add)
            vectorizer.add(float32(float32, float32))
            ufunc = vectorizer.build_ufunc()

            broadcasting_b = b[np.newaxis, :, np.newaxis, np.newaxis, :]
            self.assertPreciseEqual(ufunc(a, broadcasting_b),
                                    a + broadcasting_b)

    def test_ufunc_exception_on_write_to_readonly(self):
        z = np.ones(10)
        z.flags.writeable = False # flip write bit

        tests = []
        expect = "ufunc 'sin' called with an explicit output that is read-only"
        tests.append((jit(nopython=True), TypingError, expect))
        tests.append((jit(forceobj=True), ValueError,
                      "output array is read-only"))

        for dec, exc, msg in tests:
            def test(x):
                a = np.ones(x.shape, x.dtype) # do not copy RO attribute from x
                np.sin(a, x)

            with self.assertRaises(exc) as raises:
                dec(test)(z)

            self.assertIn(msg, str(raises.exception))

    def test_optional_type_handling(self):
        # Tests ufunc compilation with Optional type

        @njit
        def inner(x, y):
            if y > 2:
                z = None
            else:
                z = np.ones(4)
            return np.add(x, z)

        # This causes `z` to be np.ones(4) at runtime, success
        self.assertPreciseEqual(inner(np.arange(4), 1),
                                np.arange(1, 5).astype(np.float64))

        with self.assertRaises(TypeError) as raises:
            # This causes `z` to be None at runtime, TypeError raised on the
            # type cast of the Optional.
            inner(np.arange(4), 3)

        msg = "expected array(float64, 1d, C), got None"
        self.assertIn(msg, str(raises.exception))


class TestUFuncsMisc(TestCase):
    # Test for miscellaneous ufunc issues

    def test_exp2(self):
        # See issue #8898, and TargetLibraryInfo based fix in #9336
        @njit
        def foo(x):
            return np.exp2(x)

        for ty in (np.int8, np.uint16):
            x = ty(2)
            expected = foo.py_func(x)
            got = foo(x)
            self.assertPreciseEqual(expected, got)

    def test_log2(self):
        # See issue #8898, and TargetLibraryInfo based fix in #9336
        @njit
        def foo(x):
            return np.log2(x)

        for ty in (np.int8, np.uint16):
            x = ty(2)
            expected = foo.py_func(x)
            got = foo(x)
            self.assertPreciseEqual(expected, got)


if __name__ == '__main__':
    unittest.main()
