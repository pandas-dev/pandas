import gc
from itertools import product

import numpy as np
from numpy.polynomial import polynomial as poly
from numpy.polynomial import polyutils as pu

from numba import jit, njit
from numba.tests.support import (TestCase, needs_lapack,
                                 EnableNRTStatsMixin, MemoryLeakMixin)
from numba.core.errors import TypingError


def roots_fn(p):
    return np.roots(p)


def polyadd(c1,c2):
    return poly.polyadd(c1,c2)


def polysub(c1,c2):
    return poly.polysub(c1,c2)


def polymul(c1,c2):
    return poly.polymul(c1,c2)


def trimseq(seq):
    return pu.trimseq(seq)


def polyasseries1(a):
    res = pu.as_series(a)
    return res


def polyasseries2(a, trim):
    res = pu.as_series(a, trim)
    return res


def polydiv(c1, c2):
    res = poly.polydiv(c1, c2)
    return res


def polyval2(x, c):
    res = poly.polyval(x, c)
    return res


def polyval3T(x, c):
    res = poly.polyval(x, c, True)
    return res


def polyval3F(x, c):
    res = poly.polyval(x, c, False)
    return res


def polyint(c, m=1):
    res = poly.polyint(c, m)
    return res


class TestPolynomialBase(EnableNRTStatsMixin, TestCase):
    """
    Provides setUp and common data/error modes for testing polynomial functions.
    """

    # supported dtypes
    dtypes = (np.float64, np.float32, np.complex128, np.complex64)

    def setUp(self):
        # Collect leftovers from previous test cases before checking for leaks
        gc.collect()
        super(TestPolynomialBase, self).setUp()

    def assert_error(self, cfunc, args, msg, err=ValueError):
        with self.assertRaises(err) as raises:
            cfunc(*args)
        self.assertIn(msg, str(raises.exception))

    def assert_1d_input(self, cfunc, args):
        msg = "Input must be a 1d array."
        self.assert_error(cfunc, args, msg)


class TestPoly1D(TestPolynomialBase):

    def assert_no_domain_change(self, name, cfunc, args):
        msg = name + "() argument must not cause a domain change."
        self.assert_error(cfunc, args, msg)

    @needs_lapack
    def test_roots(self):

        cfunc = jit(nopython=True)(roots_fn)

        default_resolution = np.finfo(np.float64).resolution

        def check(a, **kwargs):
            expected = roots_fn(a, **kwargs)
            got = cfunc(a, **kwargs)

            # eigen decomposition used so type specific impl
            # will be used in numba whereas a wide type impl
            # will be used in numpy, so compare using a more
            # fuzzy comparator

            if a.dtype in self.dtypes:
                resolution = np.finfo(a.dtype).resolution
            else:
                # this is for integer types when roots() will cast to float64
                resolution = default_resolution

            np.testing.assert_allclose(
                expected,
                got,
                rtol=10 * resolution,
                atol=100 * resolution  # zeros tend to be fuzzy
            )

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # test vectors in real space
        # contrived examples to trip branches
        r_vectors = (
            np.array([1]),
            np.array([1, 3, 2]),
            np.array([0, 0, 0]),
            np.array([1, 6, 11, 6]),
            np.array([0, 0, 0, 1, 3, 2]),
            np.array([1, 1, 0, 0, 0]),
            np.array([0, 0, 1, 0, 0, 0])
        )

        # test loop real space
        for v, dtype in \
                product(r_vectors, [np.int32, np.int64] + list(self.dtypes)):
            a = v.astype(dtype)
            check(a)

        c_vectors = (
            np.array([1 + 1j]),
            np.array([1, 3 + 1j, 2]),
            np.array([0, 0 + 0j, 0]),
            np.array([1, 6 + 1j, 11, 6]),
            np.array([0, 0, 0, 1 + 1j, 3, 2]),
            np.array([1 + 1j, 1, 0, 0, 0]),
            np.array([0, 0, 1 + 1j, 0, 0, 0])
        )

        # test loop complex space
        for v, dtype in product(c_vectors, self.dtypes[2:]):
            a = v.astype(dtype)
            check(a)

        # check input with dimension > 1 raises
        self.assert_1d_input(cfunc, (np.arange(4.).reshape(2, 2),))

        # check real input with complex roots raises
        x = np.array([7., 2., 0., 1.])
        self.assert_no_domain_change("eigvals", cfunc, (x,))
        # but works fine if type conv to complex first
        cfunc(x.astype(np.complex128))


class TestPolynomial(MemoryLeakMixin, TestCase):

    #
    # tests for Polyutils functions
    #

    def test_trimseq_basic(self):
        pyfunc = trimseq
        cfunc = njit(trimseq)

        def inputs():
            for i in range(5):
                yield np.array([1] + [0] * i)

        for coefs in inputs():
            self.assertPreciseEqual(pyfunc(coefs), cfunc(coefs))

    def test_trimseq_exception(self):
        cfunc = njit(trimseq)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc")
        self.assertIn('The argument "seq" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as e:
            cfunc(np.arange(10).reshape(5, 2))
        self.assertIn('Coefficient array is not 1-d',
                      str(e.exception))

        with self.assertRaises(TypingError) as e:
            cfunc((1, 2, 3, 0))
        self.assertIn('Unsupported type UniTuple(int64, 4) for argument "seq"',
                      str(e.exception))

    def test_pu_as_series_basic(self):
        pyfunc1 = polyasseries1
        cfunc1 = njit(polyasseries1)
        pyfunc2 = polyasseries2
        cfunc2 = njit(polyasseries2)

        def inputs():
            yield np.arange(4)
            yield np.arange(6).reshape((2,3))
            yield (1, np.arange(3), np.arange(2, dtype=np.float32))
            yield ([1, 2, 3, 4, 0], [1, 2, 3])
            yield ((0, 0, 1e-3, 0, 1e-5, 0, 0), (1, 2, 3, 4, 5, 6, 7))
            yield ((0, 0, 1e-3, 0, 1e-5, 0, 0), (1j, 2, 3j, 4j, 5, 6j, 7))
            yield (2, [1.1, 0.])
            yield ([1, 2, 3, 0], )
            yield ((1, 2, 3, 0), )
            yield (np.array([1, 2, 3, 0]), )
            yield [np.array([1, 2, 3, 0]), np.array([1, 2, 3, 0])]
            yield [np.array([1,2,3]), ]

        for input in inputs():
            self.assertPreciseEqual(pyfunc1(input), cfunc1(input))
            self.assertPreciseEqual(pyfunc2(input, False), cfunc2(input, False))
            self.assertPreciseEqual(pyfunc2(input, True), cfunc2(input, True))

    def test_pu_as_series_exception(self):
        cfunc1 = njit(polyasseries1)
        cfunc2 = njit(polyasseries2)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc1("abc")
        self.assertIn('The argument "alist" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc2("abc", True)
        self.assertIn('The argument "alist" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc2(np.arange(4), "abc")
        self.assertIn('The argument "trim" must be boolean',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc1(([1, 2, 3], np.arange(16).reshape(4,4)))
        self.assertIn('Coefficient array is not 1-d',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc1(np.arange(8).reshape((2, 2, 2)))
        self.assertIn('Coefficient array is not 1-d',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc1([np.array([[1,2,3],[1,2,3]]), ])
        self.assertIn('Coefficient array is not 1-d',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc1(np.array([[]], dtype=np.float64))
        self.assertIn('Coefficient array is empty',
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc1(([1, 2, 3], np.array([], dtype=np.float64),
                    np.array([1,2,1])))
        self.assertIn('Coefficient array is empty',
                      str(raises.exception))

    def _test_polyarithm_basic(self, pyfunc, ignore_sign_on_zero=False):
        # test suite containing tests for polyadd, polysub, polymul, polydiv
        cfunc = njit(pyfunc)

        def inputs():
            # basic, taken from https://github.com/numpy/numpy/blob/48a8277855849be094a5979c48d9f5f1778ee4de/numpy/polynomial/tests/test_polynomial.py#L58-L123 # noqa: E501
            for i in range(5):
                for j in range(5):
                    p1 = np.array([0] * i + [1])
                    p2 = np.array([0] * j + [1])
                    yield p1, p2
            # test lists, tuples, scalars
            yield [1, 2, 3], [1, 2, 3]
            yield [1, 2, 3], (1, 2, 3)
            yield (1, 2, 3), [1, 2, 3]
            yield [1, 2, 3], 3
            yield 3, (1, 2, 3)
            # test different dtypes
            yield np.array([1, 2, 3]), np.array([1.0, 2.0, 3.0])
            yield np.array([1j, 2j, 3j]), np.array([1.0, 2.0, 3.0])
            yield np.array([1, 2, 3]), np.array([1j, 2j, 3j])
            yield (1, 2, 3), 3.0
            yield (1, 2, 3), 3j
            yield (1, 1e-3, 3), (1, 2, 3)

        for p1, p2 in inputs():
            self.assertPreciseEqual(pyfunc(p1,p2), cfunc(p1,p2),
                                    ignore_sign_on_zero=ignore_sign_on_zero)

    def _test_polyarithm_exception(self, pyfunc):
        # test suite containing tests for polyadd, polysub, polymul, polydiv
        cfunc = njit(pyfunc)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc", np.array([1,2,3]))
        self.assertIn('The argument "c1" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(np.array([1,2,3]), "abc")
        self.assertIn('The argument "c2" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as e:
            cfunc(np.arange(10).reshape(5, 2), np.array([1, 2, 3]))
        self.assertIn('Coefficient array is not 1-d',
                      str(e.exception))

        with self.assertRaises(TypingError) as e:
            cfunc(np.array([1, 2, 3]), np.arange(10).reshape(5, 2))
        self.assertIn('Coefficient array is not 1-d',
                      str(e.exception))

    def test_polyadd_basic(self):
        self._test_polyarithm_basic(polyadd)

    def test_polyadd_exception(self):
        self._test_polyarithm_exception(polyadd)

    def test_polysub_basic(self):
        self._test_polyarithm_basic(polysub, ignore_sign_on_zero=True)

    def test_polysub_exception(self):
        self._test_polyarithm_exception(polysub)

    def test_polymul_basic(self):
        self._test_polyarithm_basic(polymul)

    def test_polymul_exception(self):
        self._test_polyarithm_exception(polymul)

    def test_poly_polydiv_basic(self):
        pyfunc = polydiv
        cfunc = njit(polydiv)
        self._test_polyarithm_basic(polydiv)

        def inputs():
            # Based on https://github.com/numpy/numpy/blob/160c16f055d4d2fce072004e286d8075b31955cd/numpy/polynomial/tests/test_polynomial.py#L99-L114 # noqa: E501
            # check scalar division
            yield [2], [2]
            yield [2, 2], [2]
            # check rest.
            for i in range(5):
                for j in range(5):
                    ci = [0] * i + [1, 2]
                    cj = [0] * j + [1, 2]
                    tgt = poly.polyadd(ci, cj)
                    yield tgt, ci
            yield np.array([1,0,0,0,0,0,-1]), np.array([1,0,0,-1])

        for c1, c2 in inputs():
            self.assertPreciseEqual(pyfunc(c1, c2), cfunc(c1, c2))

    def test_poly_polydiv_exception(self):
        self._test_polyarithm_exception(polydiv)
        cfunc = njit(polydiv)
        # Based on https://github.com/numpy/numpy/blob/160c16f055d4d2fce072004e286d8075b31955cd/numpy/polynomial/tests/test_polynomial.py#L97 # noqa: E501
        # check zero division
        with self.assertRaises(ZeroDivisionError) as _:
            cfunc([1], [0])

    def test_poly_polyval_basic(self):
        pyfunc2 = polyval2
        cfunc2 = njit(polyval2)
        pyfunc3T = polyval3T
        cfunc3T = njit(polyval3T)
        pyfunc3F = polyval3F
        cfunc3F = njit(polyval3F)

        def inputs():
            # Based on https://github.com/numpy/numpy/blob/160c16f055d4d2fce072004e286d8075b31955cd/numpy/polynomial/tests/test_polynomial.py#L137-L157 # noqa: E501
            # check empty input
            yield np.array([], dtype=np.float64), [1]
            yield 1, [1,2,3]
            yield np.arange(4).reshape(2,2), [1,2,3]
            # check normal input
            for i in range(5):
                yield np.linspace(-1, 1), [0] * i + [1]
            yield np.linspace(-1, 1), [0, -1, 0, 1]
            # check that shape is preserved
            for i in range(3):
                dims = [2] * i
                x = np.zeros(dims)
                yield x, [1]
                yield x, [1, 0]
                yield x, [1, 0, 0]
            # Check that behaviour corresponds to tensor = False
            yield np.array([1, 2]), np.arange(4).reshape(2,2)
            yield [1, 2], np.arange(4).reshape(2,2)

        for x, c in inputs():
            self.assertPreciseEqual(pyfunc2(x, c), cfunc2(x, c))
            # test tensor argument
            self.assertPreciseEqual(pyfunc3T(x, c), cfunc3T(x, c))
            self.assertPreciseEqual(pyfunc3F(x, c), cfunc3F(x, c))

    def test_poly_polyval_exception(self):
        cfunc2 = njit(polyval2)
        cfunc3T = njit(polyval3T)
        cfunc3F = njit(polyval3F)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc2(3, "abc")
        self.assertIn('The argument "c" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc2("abc", 3)
        self.assertIn('The argument "x" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc2("abc", "def")
        self.assertIn('The argument "x" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc3T(3, "abc")
        self.assertIn('The argument "c" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc3T("abc", 3)
        self.assertIn('The argument "x" must be array-like',
                      str(raises.exception))

        @njit
        def polyval3(x, c, tensor):
            res = poly.polyval(x, c, tensor)
            return res
        with self.assertRaises(TypingError) as raises:
            polyval3(3, 3, "abc")
        self.assertIn('The argument "tensor" must be boolean',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc3F("abc", "def")
        self.assertIn('The argument "x" must be array-like',
                      str(raises.exception))

    def test_poly_polyint_basic(self):
        pyfunc = polyint
        cfunc = njit(polyint)
        # basic
        self.assertPreciseEqual(pyfunc([1,2,3]), cfunc([1,2,3]))
        # Based on https://github.com/numpy/numpy/blob/160c16f055d4d2fce072004e286d8075b31955cd/numpy/polynomial/tests/test_polynomial.py#L314-L381 # noqa: E501
        # test integration of zero polynomial
        for i in range(2, 5):
            self.assertPreciseEqual(pyfunc([0], m=i), cfunc([0], m=i))

        # check single integration with integration constant
        for i in range(5):
            pol = [0] * i + [1]
            self.assertPreciseEqual(pyfunc(pol, m=1), pyfunc(pol, m=1))

        # check multiple integrations with default k
        for i in range(5):
            for j in range(2, 5):
                pol = [0] * i + [1]
                self.assertPreciseEqual(pyfunc(pol, m=j), cfunc(pol, m=j))

        # test multidimensional arrays
        c2 = np.array([[0,1], [0,2]])
        self.assertPreciseEqual(pyfunc(c2), cfunc(c2))
        c3 = np.arange(8).reshape((2,2,2))
        self.assertPreciseEqual(pyfunc(c3), cfunc(c3))

    def test_poly_polyint_exception(self):
        cfunc = njit(polyint)

        self.disable_leak_check()

        with self.assertRaises(TypingError) as raises:
            cfunc("abc")
        self.assertIn('The argument "c" must be array-like',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(np.array([1,2,3]), "abc")
        self.assertIn('The argument "m" must be an integer',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(['a', 'b', 'c'], 1)
        self.assertIn('Input dtype must be scalar.',
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc(('a', 'b', 'c'), 1)
        self.assertIn('Input dtype must be scalar.',
                      str(raises.exception))

    #
    # tests for Polynomial class
    #

    def test_Polynomial_constructor(self):

        def pyfunc3(c, dom, win):
            p = poly.Polynomial(c, dom, win)
            return p
        cfunc3 = njit(pyfunc3)

        def pyfunc1(c):
            p = poly.Polynomial(c)
            return p
        cfunc1 = njit(pyfunc1)
        list1 = (np.array([0, 1]), np.array([0., 1.]))
        list2 = (np.array([0, 1]), np.array([0., 1.]))
        list3 = (np.array([0, 1]), np.array([0., 1.]))
        for c in list1:
            for dom in list2:
                for win in list3:
                    p1 = pyfunc3(c, dom, win)
                    p2 = cfunc3(c, dom, win)
                    q1 = pyfunc1(c)
                    q2 = cfunc1(c)
                    self.assertPreciseEqual(p1, p2)
                    self.assertPreciseEqual(p1.coef, p2.coef)
                    self.assertPreciseEqual(p1.domain, p2.domain)
                    self.assertPreciseEqual(p1.window, p2.window)
                    self.assertPreciseEqual(q1.coef, q2.coef)
                    self.assertPreciseEqual(q1.domain, q2.domain)
                    self.assertPreciseEqual(q1.window, q2.window)

    def test_Polynomial_exeption(self):
        def pyfunc3(c, dom, win):
            p = poly.Polynomial(c, dom, win)
            return p
        cfunc3 = njit(pyfunc3)

        self.disable_leak_check()

        input2 = np.array([1, 2])
        input3 = np.array([1, 2, 3])
        input2D = np.arange(4).reshape((2, 2))

        with self.assertRaises(ValueError) as raises:
            cfunc3(input2, input3, input2)
        self.assertIn("Domain has wrong number of elements.",
                      str(raises.exception))

        with self.assertRaises(ValueError) as raises:
            cfunc3(input2, input2, input3)
        self.assertIn("Window has wrong number of elements.",
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            cfunc3(input2D, input2, input2)
        self.assertIn("Coefficient array is not 1-d",
                      str(raises.exception))
