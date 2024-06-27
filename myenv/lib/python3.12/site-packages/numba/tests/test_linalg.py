import contextlib
import gc
from itertools import product, cycle
import sys
import warnings
from numbers import Number, Integral
import platform

import numpy as np

from numba import jit, njit, typeof
from numba.core import errors
from numba.tests.support import (TestCase, tag, needs_lapack, needs_blas,
                                 _is_armv7l, EnableNRTStatsMixin)
from .matmul_usecase import matmul_usecase
import unittest


def dot2(a, b):
    return np.dot(a, b)


def dot3(a, b, out):
    return np.dot(a, b, out=out)


def vdot(a, b):
    return np.vdot(a, b)


class TestProduct(EnableNRTStatsMixin, TestCase):
    """
    Tests for dot products.
    """

    dtypes = (np.float64, np.float32, np.complex128, np.complex64)

    def setUp(self):
        # Collect leftovers from previous test cases before checking for leaks
        gc.collect()
        super(TestProduct, self).setUp()

    def sample_vector(self, n, dtype):
        # Be careful to generate only exactly representable float values,
        # to avoid rounding discrepancies between Numpy and Numba
        base = np.arange(n)
        if issubclass(dtype, np.complexfloating):
            return (base * (1 - 0.5j) + 2j).astype(dtype)
        else:
            return (base * 0.5 - 1).astype(dtype)

    def sample_matrix(self, m, n, dtype):
        return self.sample_vector(m * n, dtype).reshape((m, n))

    @contextlib.contextmanager
    def check_contiguity_warning(self, pyfunc):
        """
        Check performance warning(s) for non-contiguity.
        """
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', errors.NumbaPerformanceWarning)
            yield
        self.assertGreaterEqual(len(w), 1)
        self.assertIs(w[0].category, errors.NumbaPerformanceWarning)
        self.assertIn("faster on contiguous arrays", str(w[0].message))
        self.assertEqual(w[0].filename, pyfunc.__code__.co_filename)
        # This works because our functions are one-liners
        self.assertEqual(w[0].lineno, pyfunc.__code__.co_firstlineno + 1)

    def check_func(self, pyfunc, cfunc, args):
        with self.assertNoNRTLeak():
            expected = pyfunc(*args)
            got = cfunc(*args)
            self.assertPreciseEqual(got, expected, ignore_sign_on_zero=True)
            del got, expected


    def _aligned_copy(self, arr):
        # This exists for armv7l because NumPy wants aligned arrays for the
        # `out` arg of functions, but np.empty/np.copy doesn't seem to always
        # produce them, in particular for complex dtypes
        size = (arr.size + 1) * arr.itemsize + 1
        datasize = arr.size * arr.itemsize
        tmp = np.empty(size, dtype=np.uint8)
        for i in range(arr.itemsize + 1):
            new = tmp[i : i + datasize].view(dtype=arr.dtype)
            if new.flags.aligned:
                break
        else:
            raise Exception("Could not obtain aligned array")
        if arr.flags.c_contiguous:
            new = np.reshape(new, arr.shape, order='C')
        else:
            new = np.reshape(new, arr.shape, order='F')
        new[:] = arr[:]
        assert new.flags.aligned
        return new

    def check_func_out(self, pyfunc, cfunc, args, out):
        copier = self._aligned_copy if _is_armv7l else np.copy
        with self.assertNoNRTLeak():
            expected = copier(out)
            got = copier(out)
            self.assertIs(pyfunc(*args, out=expected), expected)
            self.assertIs(cfunc(*args, out=got), got)
            self.assertPreciseEqual(got, expected, ignore_sign_on_zero=True)
            del got, expected

    def assert_mismatching_sizes(self, cfunc, args, is_out=False):
        with self.assertRaises(ValueError) as raises:
            cfunc(*args)
        msg = ("incompatible output array size" if is_out else
               "incompatible array sizes")
        self.assertIn(msg, str(raises.exception))

    def assert_mismatching_dtypes(self, cfunc, args, func_name="np.dot()"):
        with self.assertRaises(errors.TypingError) as raises:
            cfunc(*args)
        self.assertIn("%s arguments must all have the same dtype"
                      % (func_name,),
                      str(raises.exception))

    def check_dot_vv(self, pyfunc, func_name):
        n = 3
        cfunc = jit(nopython=True)(pyfunc)
        for dtype in self.dtypes:
            a = self.sample_vector(n, dtype)
            b = self.sample_vector(n, dtype)
            self.check_func(pyfunc, cfunc, (a, b))
            # Non-contiguous
            self.check_func(pyfunc, cfunc, (a[::-1], b[::-1]))

        # Mismatching sizes
        a = self.sample_vector(n - 1, np.float64)
        b = self.sample_vector(n, np.float64)
        self.assert_mismatching_sizes(cfunc, (a, b))
        # Mismatching dtypes
        a = self.sample_vector(n, np.float32)
        b = self.sample_vector(n, np.float64)
        self.assert_mismatching_dtypes(cfunc, (a, b), func_name=func_name)

    @needs_blas
    def test_dot_vv(self):
        """
        Test vector * vector np.dot()
        """
        self.check_dot_vv(dot2, "np.dot()")

    @needs_blas
    def test_vdot(self):
        """
        Test np.vdot()
        """
        self.check_dot_vv(vdot, "np.vdot()")

    def check_dot_vm(self, pyfunc2, pyfunc3, func_name):

        def samples(m, n):
            for order in 'CF':
                a = self.sample_matrix(m, n, np.float64).copy(order=order)
                b = self.sample_vector(n, np.float64)
                yield a, b
            for dtype in self.dtypes:
                a = self.sample_matrix(m, n, dtype)
                b = self.sample_vector(n, dtype)
                yield a, b
            # Non-contiguous
            yield a[::-1], b[::-1]

        cfunc2 = jit(nopython=True)(pyfunc2)
        if pyfunc3 is not None:
            cfunc3 = jit(nopython=True)(pyfunc3)

        for m, n in [(2, 3),
                     (3, 0),
                     (0, 3)
                     ]:
            for a, b in samples(m, n):
                self.check_func(pyfunc2, cfunc2, (a, b))
                self.check_func(pyfunc2, cfunc2, (b, a.T))
            if pyfunc3 is not None:
                for a, b in samples(m, n):
                    out = np.empty(m, dtype=a.dtype)
                    self.check_func_out(pyfunc3, cfunc3, (a, b), out)
                    self.check_func_out(pyfunc3, cfunc3, (b, a.T), out)

        # Mismatching sizes
        m, n = 2, 3
        a = self.sample_matrix(m, n - 1, np.float64)
        b = self.sample_vector(n, np.float64)
        self.assert_mismatching_sizes(cfunc2, (a, b))
        self.assert_mismatching_sizes(cfunc2, (b, a.T))
        if pyfunc3 is not None:
            out = np.empty(m, np.float64)
            self.assert_mismatching_sizes(cfunc3, (a, b, out))
            self.assert_mismatching_sizes(cfunc3, (b, a.T, out))
            a = self.sample_matrix(m, m, np.float64)
            b = self.sample_vector(m, np.float64)
            out = np.empty(m - 1, np.float64)
            self.assert_mismatching_sizes(cfunc3, (a, b, out), is_out=True)
            self.assert_mismatching_sizes(cfunc3, (b, a.T, out), is_out=True)
        # Mismatching dtypes
        a = self.sample_matrix(m, n, np.float32)
        b = self.sample_vector(n, np.float64)
        self.assert_mismatching_dtypes(cfunc2, (a, b), func_name)
        if pyfunc3 is not None:
            a = self.sample_matrix(m, n, np.float64)
            b = self.sample_vector(n, np.float64)
            out = np.empty(m, np.float32)
            self.assert_mismatching_dtypes(cfunc3, (a, b, out), func_name)

    @needs_blas
    def test_dot_vm(self):
        """
        Test vector * matrix and matrix * vector np.dot()
        """
        self.check_dot_vm(dot2, dot3, "np.dot()")

    def check_dot_mm(self, pyfunc2, pyfunc3, func_name):

        def samples(m, n, k):
            for order_a, order_b in product('CF', 'CF'):
                a = self.sample_matrix(m, k, np.float64).copy(order=order_a)
                b = self.sample_matrix(k, n, np.float64).copy(order=order_b)
                yield a, b
            for dtype in self.dtypes:
                a = self.sample_matrix(m, k, dtype)
                b = self.sample_matrix(k, n, dtype)
                yield a, b
            # Non-contiguous
            yield a[::-1], b[::-1]

        cfunc2 = jit(nopython=True)(pyfunc2)
        if pyfunc3 is not None:
            cfunc3 = jit(nopython=True)(pyfunc3)

        # Test generic matrix * matrix as well as "degenerate" cases
        # where one of the outer dimensions is 1 (i.e. really represents
        # a vector, which may select a different implementation),
        # one of the matrices is empty, or both matrices are empty.
        for m, n, k in [(2, 3, 4),  # Generic matrix * matrix
                        (1, 3, 4),  # 2d vector * matrix
                        (1, 1, 4),  # 2d vector * 2d vector
                        (0, 3, 2),  # Empty matrix * matrix, empty output
                        (3, 0, 2),  # Matrix * empty matrix, empty output
                        (0, 0, 3),  # Both arguments empty, empty output
                        (3, 2, 0),  # Both arguments empty, nonempty output
                        ]:
            for a, b in samples(m, n, k):
                self.check_func(pyfunc2, cfunc2, (a, b))
                self.check_func(pyfunc2, cfunc2, (b.T, a.T))
            if pyfunc3 is not None:
                for a, b in samples(m, n, k):
                    out = np.empty((m, n), dtype=a.dtype)
                    self.check_func_out(pyfunc3, cfunc3, (a, b), out)
                    out = np.empty((n, m), dtype=a.dtype)
                    self.check_func_out(pyfunc3, cfunc3, (b.T, a.T), out)

        # Mismatching sizes
        m, n, k = 2, 3, 4
        a = self.sample_matrix(m, k - 1, np.float64)
        b = self.sample_matrix(k, n, np.float64)
        self.assert_mismatching_sizes(cfunc2, (a, b))
        if pyfunc3 is not None:
            out = np.empty((m, n), np.float64)
            self.assert_mismatching_sizes(cfunc3, (a, b, out))
            a = self.sample_matrix(m, k, np.float64)
            b = self.sample_matrix(k, n, np.float64)
            out = np.empty((m, n - 1), np.float64)
            self.assert_mismatching_sizes(cfunc3, (a, b, out), is_out=True)
        # Mismatching dtypes
        a = self.sample_matrix(m, k, np.float32)
        b = self.sample_matrix(k, n, np.float64)
        self.assert_mismatching_dtypes(cfunc2, (a, b), func_name)
        if pyfunc3 is not None:
            a = self.sample_matrix(m, k, np.float64)
            b = self.sample_matrix(k, n, np.float64)
            out = np.empty((m, n), np.float32)
            self.assert_mismatching_dtypes(cfunc3, (a, b, out), func_name)

    @needs_blas
    def test_dot_mm(self):
        """
        Test matrix * matrix np.dot()
        """
        self.check_dot_mm(dot2, dot3, "np.dot()")

    @needs_blas
    def test_matmul_vv(self):
        """
        Test vector @ vector
        """
        self.check_dot_vv(matmul_usecase, "'@'")

    @needs_blas
    def test_matmul_vm(self):
        """
        Test vector @ matrix and matrix @ vector
        """
        self.check_dot_vm(matmul_usecase, None, "'@'")

    @needs_blas
    def test_matmul_mm(self):
        """
        Test matrix @ matrix
        """
        self.check_dot_mm(matmul_usecase, None, "'@'")

    @needs_blas
    def test_contiguity_warnings(self):
        m, k, n = 2, 3, 4
        dtype = np.float64
        a = self.sample_matrix(m, k, dtype)[::-1]
        b = self.sample_matrix(k, n, dtype)[::-1]
        out = np.empty((m, n), dtype)

        cfunc = jit(nopython=True)(dot2)
        with self.check_contiguity_warning(cfunc.py_func):
            cfunc(a, b)
        cfunc = jit(nopython=True)(dot3)
        with self.check_contiguity_warning(cfunc.py_func):
            cfunc(a, b, out)

        a = self.sample_vector(n, dtype)[::-1]
        b = self.sample_vector(n, dtype)[::-1]

        cfunc = jit(nopython=True)(vdot)
        with self.check_contiguity_warning(cfunc.py_func):
            cfunc(a, b)


# Implementation definitions for the purpose of jitting.

def invert_matrix(a):
    return np.linalg.inv(a)


def cholesky_matrix(a):
    return np.linalg.cholesky(a)


def eig_matrix(a):
    return np.linalg.eig(a)


def eigvals_matrix(a):
    return np.linalg.eigvals(a)


def eigh_matrix(a):
    return np.linalg.eigh(a)


def eigvalsh_matrix(a):
    return np.linalg.eigvalsh(a)


def svd_matrix(a, full_matrices=1):
    return np.linalg.svd(a, full_matrices)


def qr_matrix(a):
    return np.linalg.qr(a)


def lstsq_system(A, B, rcond=-1):
    return np.linalg.lstsq(A, B, rcond)


def solve_system(A, B):
    return np.linalg.solve(A, B)


def pinv_matrix(A, rcond=1e-15):  # 1e-15 from numpy impl
    return np.linalg.pinv(A)


def slogdet_matrix(a):
    return np.linalg.slogdet(a)


def det_matrix(a):
    return np.linalg.det(a)


def norm_matrix(a, ord=None):
    return np.linalg.norm(a, ord)


def cond_matrix(a, p=None):
    return np.linalg.cond(a, p)


def matrix_rank_matrix(a, tol=None):
    return np.linalg.matrix_rank(a, tol)


def matrix_power_matrix(a, n):
    return np.linalg.matrix_power(a, n)


def trace_matrix(a, offset=0):
    return np.trace(a, offset)


def trace_matrix_no_offset(a):
    return np.trace(a)


def outer_matrix(a, b, out=None):
    return np.outer(a, b, out=out)


def kron_matrix(a, b):
    return np.kron(a, b)


class TestLinalgBase(EnableNRTStatsMixin, TestCase):
    """
    Provides setUp and common data/error modes for testing np.linalg functions.
    """

    # supported dtypes
    dtypes = (np.float64, np.float32, np.complex128, np.complex64)

    def setUp(self):
        # Collect leftovers from previous test cases before checking for leaks
        gc.collect()
        super(TestLinalgBase, self).setUp()

    def sample_vector(self, n, dtype):
        # Be careful to generate only exactly representable float values,
        # to avoid rounding discrepancies between Numpy and Numba
        base = np.arange(n)
        if issubclass(dtype, np.complexfloating):
            return (base * (1 - 0.5j) + 2j).astype(dtype)
        else:
            return (base * 0.5 + 1).astype(dtype)

    def specific_sample_matrix(
            self, size, dtype, order, rank=None, condition=None):
        """
        Provides a sample matrix with an optionally specified rank or condition
        number.

        size: (rows, columns), the dimensions of the returned matrix.
        dtype: the dtype for the returned matrix.
        order: the memory layout for the returned matrix, 'F' or 'C'.
        rank: the rank of the matrix, an integer value, defaults to full rank.
        condition: the condition number of the matrix (defaults to 1.)

        NOTE: Only one of rank or condition may be set.
        """

        # default condition
        d_cond = 1.

        if len(size) != 2:
            raise ValueError("size must be a length 2 tuple.")

        if order not in ['F', 'C']:
            raise ValueError("order must be one of 'F' or 'C'.")

        if dtype not in [np.float32, np.float64, np.complex64, np.complex128]:
            raise ValueError("dtype must be a numpy floating point type.")

        if rank is not None and condition is not None:
            raise ValueError("Only one of rank or condition can be specified.")

        if condition is None:
            condition = d_cond

        if condition < 1:
            raise ValueError("Condition number must be >=1.")

        np.random.seed(0)  # repeatable seed
        m, n = size

        if m < 0 or n < 0:
            raise ValueError("Negative dimensions given for matrix shape.")

        minmn = min(m, n)
        if rank is None:
            rv = minmn
        else:
            if rank <= 0:
                raise ValueError("Rank must be greater than zero.")
            if not isinstance(rank, Integral):
                raise ValueError("Rank must an integer.")
            rv = rank
            if rank > minmn:
                raise ValueError("Rank given greater than full rank.")

        if m == 1 or n == 1:
            # vector, must be rank 1 (enforced above)
            # condition of vector is also 1
            if condition != d_cond:
                raise ValueError(
                    "Condition number was specified for a vector (always 1.).")
            maxmn = max(m, n)
            Q = self.sample_vector(maxmn, dtype).reshape(m, n)
        else:
            # Build a sample matrix via combining SVD like inputs.

            # Create matrices of left and right singular vectors.
            # This could use Modified Gram-Schmidt and perhaps be quicker,
            # at present it uses QR decompositions to obtain orthonormal
            # matrices.
            tmp = self.sample_vector(m * m, dtype).reshape(m, m)
            U, _ = np.linalg.qr(tmp)
            # flip the second array, else for m==n the identity matrix appears
            tmp = self.sample_vector(n * n, dtype)[::-1].reshape(n, n)
            V, _ = np.linalg.qr(tmp)
            # create singular values.
            sv = np.linspace(d_cond, condition, rv)
            S = np.zeros((m, n))
            idx = np.nonzero(np.eye(m, n))
            S[idx[0][:rv], idx[1][:rv]] = sv
            Q = np.dot(np.dot(U, S), V.T)  # construct
            Q = np.array(Q, dtype=dtype, order=order)  # sort out order/type

        return Q

    def assert_error(self, cfunc, args, msg, err=ValueError):
        with self.assertRaises(err) as raises:
            cfunc(*args)
        self.assertIn(msg, str(raises.exception))

    def assert_non_square(self, cfunc, args):
        msg = "Last 2 dimensions of the array must be square."
        self.assert_error(cfunc, args, msg, np.linalg.LinAlgError)

    def assert_wrong_dtype(self, name, cfunc, args):
        msg = "np.linalg.%s() only supported on float and complex arrays" % name
        self.assert_error(cfunc, args, msg, errors.TypingError)

    def assert_wrong_dimensions(self, name, cfunc, args, la_prefix=True):
        prefix = "np.linalg" if la_prefix else "np"
        msg = "%s.%s() only supported on 2-D arrays" % (prefix, name)
        self.assert_error(cfunc, args, msg, errors.TypingError)

    def assert_no_nan_or_inf(self, cfunc, args):
        msg = "Array must not contain infs or NaNs."
        self.assert_error(cfunc, args, msg, np.linalg.LinAlgError)

    def assert_contig_sanity(self, got, expected_contig):
        """
        This checks that in a computed result from numba (array, possibly tuple
        of arrays) all the arrays are contiguous in memory and that they are
        all at least one of "C_CONTIGUOUS" or "F_CONTIGUOUS". The computed
        result of the contiguousness is then compared against a hardcoded
        expected result.

        got: is the computed results from numba
        expected_contig: is "C" or "F" and is the expected type of
                        contiguousness across all input values
                        (and therefore tests).
        """

        if isinstance(got, tuple):
            # tuple present, check all results
            for a in got:
                self.assert_contig_sanity(a, expected_contig)
        else:
            if not isinstance(got, Number):
                # else a single array is present
                c_contig = got.flags.c_contiguous
                f_contig = got.flags.f_contiguous

                # check that the result (possible set of) is at least one of
                # C or F contiguous.
                msg = "Results are not at least one of all C or F contiguous."
                self.assertTrue(c_contig | f_contig, msg)

                msg = "Computed contiguousness does not match expected."
                if expected_contig == "C":
                    self.assertTrue(c_contig, msg)
                elif expected_contig == "F":
                    self.assertTrue(f_contig, msg)
                else:
                    raise ValueError("Unknown contig")

    def assert_raise_on_singular(self, cfunc, args):
        msg = "Matrix is singular to machine precision."
        self.assert_error(cfunc, args, msg, err=np.linalg.LinAlgError)

    def assert_is_identity_matrix(self, got, rtol=None, atol=None):
        """
        Checks if a matrix is equal to the identity matrix.
        """
        # check it is square
        self.assertEqual(got.shape[-1], got.shape[-2])
        # create identity matrix
        eye = np.eye(got.shape[-1], dtype=got.dtype)
        resolution = 5 * np.finfo(got.dtype).resolution
        if rtol is None:
            rtol = 10 * resolution
        if atol is None:
            atol = 100 * resolution  # zeros tend to be fuzzy
        # check it matches
        np.testing.assert_allclose(got, eye, rtol, atol)

    def assert_invalid_norm_kind(self, cfunc, args):
        """
        For use in norm() and cond() tests.
        """
        msg = "Invalid norm order for matrices."
        self.assert_error(cfunc, args, msg, ValueError)

    def assert_raise_on_empty(self, cfunc, args):
        msg = 'Arrays cannot be empty'
        self.assert_error(cfunc, args, msg, np.linalg.LinAlgError)


class TestTestLinalgBase(TestCase):
    """
    The sample matrix code TestLinalgBase.specific_sample_matrix()
    is a bit involved, this class tests it works as intended.
    """

    def test_specific_sample_matrix(self):

        # add a default test to the ctor, it never runs so doesn't matter
        inst = TestLinalgBase('specific_sample_matrix')

        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]

        # test loop
        for size, dtype, order in product(sizes, inst.dtypes, 'FC'):

            m, n = size
            minmn = min(m, n)

            # test default full rank
            A = inst.specific_sample_matrix(size, dtype, order)
            self.assertEqual(A.shape, size)
            self.assertEqual(np.linalg.matrix_rank(A), minmn)

            # test reduced rank if a reduction is possible
            if minmn > 1:
                rank = minmn - 1
                A = inst.specific_sample_matrix(size, dtype, order, rank=rank)
                self.assertEqual(A.shape, size)
                self.assertEqual(np.linalg.matrix_rank(A), rank)

            resolution = 5 * np.finfo(dtype).resolution

            # test default condition
            A = inst.specific_sample_matrix(size, dtype, order)
            self.assertEqual(A.shape, size)
            np.testing.assert_allclose(np.linalg.cond(A),
                                       1.,
                                       rtol=resolution,
                                       atol=resolution)

            # test specified condition if matrix is > 1D
            if minmn > 1:
                condition = 10.
                A = inst.specific_sample_matrix(
                    size, dtype, order, condition=condition)
                self.assertEqual(A.shape, size)
                np.testing.assert_allclose(np.linalg.cond(A),
                                           10.,
                                           rtol=resolution,
                                           atol=resolution)

        # check errors are raised appropriately
        def check_error(args, msg, err=ValueError):
            with self.assertRaises(err) as raises:
                inst.specific_sample_matrix(*args)
            self.assertIn(msg, str(raises.exception))

        # check the checker runs ok
        with self.assertRaises(AssertionError) as raises:
            msg = "blank"
            check_error(((2, 3), np.float64, 'F'), msg, err=ValueError)

        # check invalid inputs...

        # bad size
        msg = "size must be a length 2 tuple."
        check_error(((1,), np.float64, 'F'), msg, err=ValueError)

        # bad order
        msg = "order must be one of 'F' or 'C'."
        check_error(((2, 3), np.float64, 'z'), msg, err=ValueError)

        # bad type
        msg = "dtype must be a numpy floating point type."
        check_error(((2, 3), np.int32, 'F'), msg, err=ValueError)

        # specifying both rank and condition
        msg = "Only one of rank or condition can be specified."
        check_error(((2, 3), np.float64, 'F', 1, 1), msg, err=ValueError)

        # specifying negative condition
        msg = "Condition number must be >=1."
        check_error(((2, 3), np.float64, 'F', None, -1), msg, err=ValueError)

        # specifying negative matrix dimension
        msg = "Negative dimensions given for matrix shape."
        check_error(((2, -3), np.float64, 'F'), msg, err=ValueError)

        # specifying negative rank
        msg = "Rank must be greater than zero."
        check_error(((2, 3), np.float64, 'F', -1), msg, err=ValueError)

        # specifying a rank greater than maximum rank
        msg = "Rank given greater than full rank."
        check_error(((2, 3), np.float64, 'F', 4), msg, err=ValueError)

        # specifying a condition number for a vector
        msg = "Condition number was specified for a vector (always 1.)."
        check_error(((1, 3), np.float64, 'F', None, 10), msg, err=ValueError)

        # specifying a non integer rank
        msg = "Rank must an integer."
        check_error(((2, 3), np.float64, 'F', 1.5), msg, err=ValueError)


class TestLinalgInv(TestLinalgBase):
    """
    Tests for np.linalg.inv.
    """

    @needs_lapack
    def test_linalg_inv(self):
        """
        Test np.linalg.inv
        """
        n = 10
        cfunc = jit(nopython=True)(invert_matrix)

        def check(a, **kwargs):
            expected = invert_matrix(a)
            got = cfunc(a)
            self.assert_contig_sanity(got, "F")

            use_reconstruction = False

            # try strict
            try:
                np.testing.assert_array_almost_equal_nulp(got, expected,
                                                          nulp=10)
            except AssertionError:
                # fall back to reconstruction
                use_reconstruction = True

            if use_reconstruction:
                rec = np.dot(got, a)
                self.assert_is_identity_matrix(rec)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a)

        for dtype, order in product(self.dtypes, 'CF'):
            a = self.specific_sample_matrix((n, n), dtype, order)
            check(a)

        # 0 dimensioned matrix
        check(np.empty((0, 0)))

        # Non square matrix
        self.assert_non_square(cfunc, (np.ones((2, 3)),))

        # Wrong dtype
        self.assert_wrong_dtype("inv", cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions("inv", cfunc, (np.ones(10),))

        # Singular matrix
        self.assert_raise_on_singular(cfunc, (np.zeros((2, 2)),))

    @needs_lapack
    def test_no_input_mutation(self):
        X = np.array([[1., 3, 2, 7,],
                      [-5, 4, 2, 3,],
                      [9, -3, 1, 1,],
                      [2, -2, 2, 8,]], order='F')

        X_orig = np.copy(X)

        @jit(nopython=True)
        def ainv(X, test):
            if test:
                # not executed, but necessary to trigger A ordering in X
                X = X[1:2, :]
            return np.linalg.inv(X)

        expected = ainv.py_func(X, False)
        np.testing.assert_allclose(X, X_orig)

        got = ainv(X, False)
        np.testing.assert_allclose(X, X_orig)

        np.testing.assert_allclose(expected, got)


class TestLinalgCholesky(TestLinalgBase):
    """
    Tests for np.linalg.cholesky.
    """

    def sample_matrix(self, m, dtype, order):
        # pd. (positive definite) matrix has eigenvalues in Z+
        np.random.seed(0)  # repeatable seed
        A = np.random.rand(m, m)
        # orthonormal q needed to form up q^{-1}*D*q
        # no "orth()" in numpy
        q, _ = np.linalg.qr(A)
        L = np.arange(1, m + 1)  # some positive eigenvalues
        Q = np.dot(np.dot(q.T, np.diag(L)), q)  # construct
        Q = np.array(Q, dtype=dtype, order=order)  # sort out order/type
        return Q

    def assert_not_pd(self, cfunc, args):
        msg = "Matrix is not positive definite."
        self.assert_error(cfunc, args, msg, np.linalg.LinAlgError)

    @needs_lapack
    def test_linalg_cholesky(self):
        """
        Test np.linalg.cholesky
        """
        n = 10
        cfunc = jit(nopython=True)(cholesky_matrix)

        def check(a):
            expected = cholesky_matrix(a)
            got = cfunc(a)
            use_reconstruction = False
            # check that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "C")

            # try strict
            try:
                np.testing.assert_array_almost_equal_nulp(got, expected,
                                                          nulp=10)
            except AssertionError:
                # fall back to reconstruction
                use_reconstruction = True

            # try via reconstruction
            if use_reconstruction:
                rec = np.dot(got, np.conj(got.T))
                resolution = 5 * np.finfo(a.dtype).resolution
                np.testing.assert_allclose(
                    a,
                    rec,
                    rtol=resolution,
                    atol=resolution
                )

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a)

        for dtype, order in product(self.dtypes, 'FC'):
            a = self.sample_matrix(n, dtype, order)
            check(a)

        # 0 dimensioned matrix
        check(np.empty((0, 0)))

        rn = "cholesky"
        # Non square matrices
        self.assert_non_square(cfunc, (np.ones((2, 3), dtype=np.float64),))

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64),))

        # not pd
        self.assert_not_pd(cfunc,
                           (np.ones(4, dtype=np.float64).reshape(2, 2),))


class TestLinalgEigenSystems(TestLinalgBase):
    """
    Tests for np.linalg.eig/eigvals.
    """

    def sample_matrix(self, m, dtype, order):
        # This is a tridiag with the same but skewed values on the diagonals
        v = self.sample_vector(m, dtype)
        Q = np.diag(v)
        idx = np.nonzero(np.eye(Q.shape[0], Q.shape[1], 1))
        Q[idx] = v[1:]
        idx = np.nonzero(np.eye(Q.shape[0], Q.shape[1], -1))
        Q[idx] = v[:-1]
        Q = np.array(Q, dtype=dtype, order=order)
        return Q

    def assert_no_domain_change(self, name, cfunc, args):
        msg = name + "() argument must not cause a domain change."
        self.assert_error(cfunc, args, msg)

    def _check_worker(self, cfunc, name, expected_res_len,
                      check_for_domain_change):
        def check(*args):
            expected = cfunc.py_func(*args)
            got = cfunc(*args)
            a = args[0]
            # check that the returned tuple is same length
            self.assertEqual(len(expected), len(got))
            # and that dimension is correct
            res_is_tuple = False
            if isinstance(got, tuple):
                res_is_tuple = True
                self.assertEqual(len(got), expected_res_len)
            else:  # its an array
                self.assertEqual(got.ndim, expected_res_len)

            # and that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "F")

            use_reconstruction = False
            # try plain match of each array to np first
            for k in range(len(expected)):
                try:
                    np.testing.assert_array_almost_equal_nulp(
                        got[k], expected[k], nulp=10)
                except AssertionError:
                    # plain match failed, test by reconstruction
                    use_reconstruction = True

            # If plain match fails then reconstruction is used.
            # this checks that A*V ~== V*diag(W)
            # i.e. eigensystem ties out
            # this is required as numpy uses only double precision lapack
            # routines and computation of eigenvectors is numerically
            # sensitive, numba uses the type specific routines therefore
            # sometimes comes out with a different (but entirely
            # valid) answer (eigenvectors are not unique etc.).
            # This is only applicable if eigenvectors are computed
            # along with eigenvalues i.e. result is a tuple.
            resolution = 5 * np.finfo(a.dtype).resolution
            if use_reconstruction:
                if res_is_tuple:
                    w, v = got
                    # modify 'a' if hermitian eigensystem functionality is
                    # being tested. 'L' for use lower part is default and
                    # the only thing used at present so we conjugate transpose
                    # the lower part into the upper for use in the
                    # reconstruction. By construction the sample matrix is
                    # tridiag so this is just a question of copying the lower
                    # diagonal into the upper and conjugating on the way.
                    if name[-1] == 'h':
                        idxl = np.nonzero(np.eye(a.shape[0], a.shape[1], -1))
                        idxu = np.nonzero(np.eye(a.shape[0], a.shape[1], 1))
                        cfunc(*args)
                        # upper idx must match lower for default uplo="L"
                        # if complex, conjugate
                        a[idxu] = np.conj(a[idxl])
                        # also, only the real part of the diagonals is
                        # considered in the calculation so the imag is zeroed
                        # out for the purposes of use in reconstruction.
                        a[np.diag_indices(a.shape[0])] = np.real(np.diag(a))

                    lhs = np.dot(a, v)
                    rhs = np.dot(v, np.diag(w))

                    np.testing.assert_allclose(
                        lhs.real,
                        rhs.real,
                        rtol=resolution,
                        atol=resolution
                    )
                    if np.iscomplexobj(v):
                        np.testing.assert_allclose(
                            lhs.imag,
                            rhs.imag,
                            rtol=resolution,
                            atol=resolution
                        )
                else:
                    # This isn't technically reconstruction but is here to
                    # deal with that the order of the returned eigenvalues
                    # may differ in the case of routines just returning
                    # eigenvalues and there's no true reconstruction
                    # available with which to perform a check.
                    np.testing.assert_allclose(
                        np.sort(expected),
                        np.sort(got),
                        rtol=resolution,
                        atol=resolution
                    )

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(*args)
        return check

    def checker_for_linalg_eig(
            self, name, func, expected_res_len, check_for_domain_change=None):
        """
        Test np.linalg.eig
        """
        n = 10
        cfunc = jit(nopython=True)(func)
        check = self._check_worker(cfunc, name, expected_res_len,
                                   check_for_domain_change)


        # The main test loop
        for dtype, order in product(self.dtypes, 'FC'):
            a = self.sample_matrix(n, dtype, order)
            check(a)

        # Test both a real and complex type as the impls are different
        for ty in [np.float32, np.complex64]:

            # 0 dimensioned matrix
            check(np.empty((0, 0), dtype=ty))

            # Non square matrices
            self.assert_non_square(cfunc, (np.ones((2, 3), dtype=ty),))

            # Wrong dtype
            self.assert_wrong_dtype(name, cfunc,
                                    (np.ones((2, 2), dtype=np.int32),))

            # Dimension issue
            self.assert_wrong_dimensions(name, cfunc, (np.ones(10, dtype=ty),))

            # no nans or infs
            self.assert_no_nan_or_inf(cfunc,
                                      (np.array([[1., 2., ], [np.inf, np.nan]],
                                                dtype=ty),))

        if check_for_domain_change:
            # By design numba does not support dynamic return types, numpy does
            # and uses this in the case of returning eigenvalues/vectors of
            # a real matrix. The return type of np.linalg.eig(), when
            # operating on a matrix in real space depends on the values present
            # in the matrix itself (recalling that eigenvalues are the roots of the
            # characteristic polynomial of the system matrix, which will by
            # construction depend on the values present in the system matrix).
            # This test asserts that if a domain change is required on the return
            # type, i.e. complex eigenvalues from a real input, an error is raised.
            # For complex types, regardless of the value of the imaginary part of
            # the returned eigenvalues, a complex type will be returned, this
            # follows numpy and fits in with numba.

            # First check that the computation is valid (i.e. in complex space)
            A = np.array([[1, -2], [2, 1]])
            check(A.astype(np.complex128))
            # and that the imaginary part is nonzero
            l, _ = func(A)
            self.assertTrue(np.any(l.imag))

            # Now check that the computation fails in real space
            for ty in [np.float32, np.float64]:
                self.assert_no_domain_change(name, cfunc, (A.astype(ty),))

    @needs_lapack
    def test_linalg_eig(self):
        self.checker_for_linalg_eig("eig", eig_matrix, 2, True)

    @needs_lapack
    def test_linalg_eigvals(self):
        self.checker_for_linalg_eig("eigvals", eigvals_matrix, 1, True)

    @needs_lapack
    def test_linalg_eigh(self):
        self.checker_for_linalg_eig("eigh", eigh_matrix, 2, False)

    @needs_lapack
    def test_linalg_eigvalsh(self):
        self.checker_for_linalg_eig("eigvalsh", eigvalsh_matrix, 1, False)

    @needs_lapack
    def test_no_input_mutation(self):
        # checks inputs are not mutated

        for c in (('eig', 2, True),
                  ('eigvals', 1, True),
                  ('eigh', 2, False),
                  ('eigvalsh', 1, False)):

            m, nout, domain_change = c

            meth = getattr(np.linalg, m)

            @jit(nopython=True)
            def func(X, test):
                if test:
                    # not executed, but necessary to trigger A ordering in X
                    X = X[1:2, :]
                return meth(X)

            check = self._check_worker(func, m, nout, domain_change)

            for dtype in (np.float64, np.complex128):
                with self.subTest(meth=meth, dtype=dtype):
                    # trivial system, doesn't matter, just checking if it gets
                    # mutated
                    X = np.array([[10., 1, 0, 1],
                                [1, 9, 0, 0],
                                [0, 0, 8, 0],
                                [1, 0, 0, 7],
                                ], order='F', dtype=dtype)

                    X_orig = np.copy(X)

                    expected = func.py_func(X, False)
                    np.testing.assert_allclose(X, X_orig)

                    got = func(X, False)
                    np.testing.assert_allclose(X, X_orig)

                    check(X, False)


class TestLinalgSvd(TestLinalgBase):
    """
    Tests for np.linalg.svd.
    """

    # This checks that A ~= U*S*V**H, i.e. SV decomposition ties out.  This is
    # required as NumPy uses only double precision LAPACK routines and
    # computation of SVD is numerically sensitive. Numba uses type-specific
    # routines and therefore sometimes comes out with a different answer to
    # NumPy (orthonormal bases are not unique, etc.).

    def check_reconstruction(self, a, got, expected):
        u, sv, vt = got

        # Check they are dimensionally correct
        for k in range(len(expected)):
            self.assertEqual(got[k].shape, expected[k].shape)

        # Columns in u and rows in vt dictates the working size of s
        s = np.zeros((u.shape[1], vt.shape[0]))
        np.fill_diagonal(s, sv)

        rec = np.dot(np.dot(u, s), vt)
        resolution = np.finfo(a.dtype).resolution
        np.testing.assert_allclose(
            a,
            rec,
            rtol=10 * resolution,
            atol=100 * resolution  # zeros tend to be fuzzy
        )

    @needs_lapack
    def test_linalg_svd(self):
        """
        Test np.linalg.svd
        """
        cfunc = jit(nopython=True)(svd_matrix)

        def check(a, **kwargs):
            expected = svd_matrix(a, **kwargs)
            got = cfunc(a, **kwargs)
            # check that the returned tuple is same length
            self.assertEqual(len(expected), len(got))
            # and that length is 3
            self.assertEqual(len(got), 3)
            # and that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "F")

            use_reconstruction = False
            # try plain match of each array to np first
            for k in range(len(expected)):

                try:
                    np.testing.assert_array_almost_equal_nulp(
                        got[k], expected[k], nulp=10)
                except AssertionError:
                    # plain match failed, test by reconstruction
                    use_reconstruction = True

            if use_reconstruction:
                self.check_reconstruction(a, got, expected)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # test: column vector, tall, wide, square, row vector
        # prime sizes
        sizes = [(7, 1), (7, 5), (5, 7), (3, 3), (1, 7)]

        # flip on reduced or full matrices
        full_matrices = (True, False)

        # test loop
        for size, dtype, fmat, order in \
                product(sizes, self.dtypes, full_matrices, 'FC'):

            a = self.specific_sample_matrix(size, dtype, order)
            check(a, full_matrices=fmat)

        rn = "svd"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64),))

        # no nans or infs
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2., ], [np.inf, np.nan]],
                                            dtype=np.float64),))
        # empty
        for sz in [(0, 1), (1, 0), (0, 0)]:
            args = (np.empty(sz), True)
            self.assert_raise_on_empty(cfunc, args)

    @needs_lapack
    def test_no_input_mutation(self):
        X = np.array([[1., 3, 2, 7,],
                      [-5, 4, 2, 3,],
                      [9, -3, 1, 1,],
                      [2, -2, 2, 8,]], order='F')

        X_orig = np.copy(X)

        @jit(nopython=True)
        def func(X, test):
            if test:
                # not executed, but necessary to trigger A ordering in X
                X = X[1:2, :]
            return np.linalg.svd(X)

        expected = func.py_func(X, False)
        np.testing.assert_allclose(X, X_orig)

        got = func(X, False)
        np.testing.assert_allclose(X, X_orig)

        try:
            for e_a, g_a in zip(expected, got):
                np.testing.assert_allclose(e_a, g_a)
        except AssertionError:
            self.check_reconstruction(X, got, expected)


class TestLinalgQr(TestLinalgBase):
    """
    Tests for np.linalg.qr.
    """

    @needs_lapack
    def test_linalg_qr(self):
        """
        Test np.linalg.qr
        """
        cfunc = jit(nopython=True)(qr_matrix)

        def check(a, **kwargs):
            expected = qr_matrix(a, **kwargs)
            got = cfunc(a, **kwargs)

            # check that the returned tuple is same length
            self.assertEqual(len(expected), len(got))
            # and that length is 2
            self.assertEqual(len(got), 2)
            # and that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "F")

            use_reconstruction = False
            # try plain match of each array to np first
            for k in range(len(expected)):
                try:
                    np.testing.assert_array_almost_equal_nulp(
                        got[k], expected[k], nulp=10)
                except AssertionError:
                    # plain match failed, test by reconstruction
                    use_reconstruction = True

            # if plain match fails then reconstruction is used.
            # this checks that A ~= Q*R and that (Q^H)*Q = I
            # i.e. QR decomposition ties out
            # this is required as numpy uses only double precision lapack
            # routines and computation of qr is numerically
            # sensitive, numba using the type specific routines therefore
            # sometimes comes out with a different answer (orthonormal bases
            # are not unique etc.).
            if use_reconstruction:
                q, r = got

                # check they are dimensionally correct
                for k in range(len(expected)):
                    self.assertEqual(got[k].shape, expected[k].shape)

                # check A=q*r
                rec = np.dot(q, r)
                resolution = np.finfo(a.dtype).resolution
                np.testing.assert_allclose(
                    a,
                    rec,
                    rtol=10 * resolution,
                    atol=100 * resolution  # zeros tend to be fuzzy
                )

                # check q is orthonormal
                self.assert_is_identity_matrix(np.dot(np.conjugate(q.T), q))

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # test: column vector, tall, wide, square, row vector
        # prime sizes
        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]

        # test loop
        for size, dtype, order in \
                product(sizes, self.dtypes, 'FC'):
            a = self.specific_sample_matrix(size, dtype, order)
            check(a)

        rn = "qr"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64),))

        # no nans or infs
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2., ], [np.inf, np.nan]],
                                            dtype=np.float64),))

        # empty
        for sz in [(0, 1), (1, 0), (0, 0)]:
            self.assert_raise_on_empty(cfunc, (np.empty(sz),))

    @needs_lapack
    def test_no_input_mutation(self):
        X = np.array([[1., 3, 2, 7,],
                      [-5, 4, 2, 3,],
                      [9, -3, 1, 1,],
                      [2, -2, 2, 8,]], order='F')

        X_orig = np.copy(X)

        @jit(nopython=True)
        def func(X, test):
            if test:
                # not executed, but necessary to trigger A ordering in X
                X = X[1:2, :]
            return np.linalg.qr(X)

        expected = func.py_func(X, False)
        np.testing.assert_allclose(X, X_orig)

        got = func(X, False)
        np.testing.assert_allclose(X, X_orig)

        for e_a, g_a in zip(expected, got):
            np.testing.assert_allclose(e_a, g_a)


class TestLinalgSystems(TestLinalgBase):
    """
    Base class for testing "system" solvers from np.linalg.
    Namely np.linalg.solve() and np.linalg.lstsq().
    """

    # check for RHS with dimension > 2 raises
    def assert_wrong_dimensions_1D(self, name, cfunc, args, la_prefix=True):
        prefix = "np.linalg" if la_prefix else "np"
        msg = "%s.%s() only supported on 1 and 2-D arrays" % (prefix, name)
        self.assert_error(cfunc, args, msg, errors.TypingError)

    # check that a dimensionally invalid system raises
    def assert_dimensionally_invalid(self, cfunc, args):
        msg = "Incompatible array sizes, system is not dimensionally valid."
        self.assert_error(cfunc, args, msg, np.linalg.LinAlgError)

    # check that args with differing dtypes raise
    def assert_homogeneous_dtypes(self, name, cfunc, args):
        msg = "np.linalg.%s() only supports inputs that have homogeneous dtypes." % name
        self.assert_error(cfunc, args, msg, errors.TypingError)


class TestLinalgLstsq(TestLinalgSystems):
    """
    Tests for np.linalg.lstsq.
    """

    # NOTE: The testing of this routine is hard as it has to handle numpy
    # using double precision routines on single precision input, this has
    # a knock on effect especially in rank deficient cases and cases where
    # conditioning is generally poor. As a result computed ranks can differ
    # and consequently the calculated residual can differ.
    # The tests try and deal with this as best as they can through the use
    # of reconstruction and measures like residual norms.
    # Suggestions for improvements are welcomed!

    @needs_lapack
    def test_linalg_lstsq(self):
        """
        Test np.linalg.lstsq
        """
        cfunc = jit(nopython=True)(lstsq_system)

        def check(A, B, **kwargs):
            expected = lstsq_system(A, B, **kwargs)
            got = cfunc(A, B, **kwargs)

            # check that the returned tuple is same length
            self.assertEqual(len(expected), len(got))
            # and that length is 4
            self.assertEqual(len(got), 4)
            # and that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "C")

            use_reconstruction = False

            # check the ranks are the same and continue to a standard
            # match if that is the case (if ranks differ, then output
            # in e.g. residual array is of different size!).
            try:
                self.assertEqual(got[2], expected[2])
                # try plain match of each array to np first
                for k in range(len(expected)):
                    try:
                        np.testing.assert_array_almost_equal_nulp(
                            got[k], expected[k], nulp=10)
                    except AssertionError:
                        # plain match failed, test by reconstruction
                        use_reconstruction = True
            except AssertionError:
                use_reconstruction = True

            if use_reconstruction:
                x, res, rank, s = got

                # indicies in the output which are ndarrays
                out_array_idx = [0, 1, 3]

                try:
                    # check the ranks are the same
                    self.assertEqual(rank, expected[2])
                    # check they are dimensionally correct, skip [2] = rank.
                    for k in out_array_idx:
                        if isinstance(expected[k], np.ndarray):
                            self.assertEqual(got[k].shape, expected[k].shape)
                except AssertionError:
                    # check the rank differs by 1. (numerical fuzz)
                    self.assertTrue(abs(rank - expected[2]) < 2)

                # check if A*X = B
                resolution = np.finfo(A.dtype).resolution
                try:
                    # this will work so long as the conditioning is
                    # ok and the rank is full
                    rec = np.dot(A, x)
                    np.testing.assert_allclose(
                        B,
                        rec,
                        rtol=10 * resolution,
                        atol=10 * resolution
                    )
                except AssertionError:
                    # system is probably under/over determined and/or
                    # poorly conditioned. Check slackened equality
                    # and that the residual norm is the same.
                    for k in out_array_idx:
                        try:
                            np.testing.assert_allclose(
                                expected[k],
                                got[k],
                                rtol=100 * resolution,
                                atol=100 * resolution
                            )
                        except AssertionError:
                            # check the fail is likely due to bad conditioning
                            c = np.linalg.cond(A)
                            self.assertGreater(10 * c, (1. / resolution))

                        # make sure the residual 2-norm is ok
                        # if this fails its probably due to numpy using double
                        # precision LAPACK routines for singles.
                        res_expected = np.linalg.norm(
                            B - np.dot(A, expected[0]))
                        res_got = np.linalg.norm(B - np.dot(A, x))
                        # rtol = 10. as all the systems are products of orthonormals
                        # and on the small side (rows, cols) < 100.
                        np.testing.assert_allclose(
                            res_expected, res_got, rtol=10.)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(A, B, **kwargs)

        # test: column vector, tall, wide, square, row vector
        # prime sizes, the A's
        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]
        # compatible B's for Ax=B must have same number of rows and 1 or more
        # columns

        # This test takes ages! So combinations are trimmed via cycling

        # gets a dtype
        cycle_dt = cycle(self.dtypes)

        orders = ['F', 'C']
        # gets a memory order flag
        cycle_order = cycle(orders)

        # a specific condition number to use in the following tests
        # there is nothing special about it other than it is not magic
        specific_cond = 10.

        # inner test loop, extracted as there's additional logic etc required
        # that'd end up with this being repeated a lot
        def inner_test_loop_fn(A, dt, **kwargs):
            # test solve Ax=B for (column, matrix) B, same dtype as A
            b_sizes = (1, 13)

            for b_size in b_sizes:

                # check 2D B
                b_order = next(cycle_order)
                B = self.specific_sample_matrix(
                    (A.shape[0], b_size), dt, b_order)
                check(A, B, **kwargs)

                # check 1D B
                b_order = next(cycle_order)
                tmp = B[:, 0].copy(order=b_order)
                check(A, tmp, **kwargs)

        # test loop
        for a_size in sizes:

            dt = next(cycle_dt)
            a_order = next(cycle_order)

            # A full rank, well conditioned system
            A = self.specific_sample_matrix(a_size, dt, a_order)

            # run the test loop
            inner_test_loop_fn(A, dt)

            m, n = a_size
            minmn = min(m, n)

            # operations that only make sense with a 2D matrix system
            if m != 1 and n != 1:

                # Test a rank deficient system
                r = minmn - 1
                A = self.specific_sample_matrix(
                    a_size, dt, a_order, rank=r)
                # run the test loop
                inner_test_loop_fn(A, dt)

                # Test a system with a given condition number for use in
                # testing the rcond parameter.
                # This works because the singular values in the
                # specific_sample_matrix code are linspace (1, cond, [0... if
                # rank deficient])
                A = self.specific_sample_matrix(
                    a_size, dt, a_order, condition=specific_cond)
                # run the test loop
                rcond = 1. / specific_cond
                approx_half_rank_rcond = minmn * rcond
                inner_test_loop_fn(A, dt,
                                   rcond=approx_half_rank_rcond)

        # check empty arrays
        empties = [
        [(0, 1), (1,)], # empty A, valid b
        [(1, 0), (1,)], # empty A, valid b
        [(1, 1), (0,)], # valid A, empty 1D b
        [(1, 1), (1, 0)],  # valid A, empty 2D b
        ]

        for A, b in empties:
            args = (np.empty(A), np.empty(b))
            self.assert_raise_on_empty(cfunc, args)

        # Test input validation
        ok = np.array([[1., 2.], [3., 4.]], dtype=np.float64)

        # check ok input is ok
        cfunc, (ok, ok)

        # check bad inputs
        rn = "lstsq"

        # Wrong dtype
        bad = np.array([[1, 2], [3, 4]], dtype=np.int32)
        self.assert_wrong_dtype(rn, cfunc, (ok, bad))
        self.assert_wrong_dtype(rn, cfunc, (bad, ok))

        # different dtypes
        bad = np.array([[1, 2], [3, 4]], dtype=np.float32)
        self.assert_homogeneous_dtypes(rn, cfunc, (ok, bad))
        self.assert_homogeneous_dtypes(rn, cfunc, (bad, ok))

        # Dimension issue
        bad = np.array([1, 2], dtype=np.float64)
        self.assert_wrong_dimensions(rn, cfunc, (bad, ok))

        # no nans or infs
        bad = np.array([[1., 2., ], [np.inf, np.nan]], dtype=np.float64)
        self.assert_no_nan_or_inf(cfunc, (ok, bad))
        self.assert_no_nan_or_inf(cfunc, (bad, ok))

        # check 1D is accepted for B (2D is done previously)
        # and then that anything of higher dimension raises
        oneD = np.array([1., 2.], dtype=np.float64)
        cfunc, (ok, oneD)
        bad = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], dtype=np.float64)
        self.assert_wrong_dimensions_1D(rn, cfunc, (ok, bad))

        # check a dimensionally invalid system raises (1D and 2D cases
        # checked)
        bad1D = np.array([1.], dtype=np.float64)
        bad2D = np.array([[1.], [2.], [3.]], dtype=np.float64)
        self.assert_dimensionally_invalid(cfunc, (ok, bad1D))
        self.assert_dimensionally_invalid(cfunc, (ok, bad2D))

    @needs_lapack
    def test_issue3368(self):
        X = np.array([[1., 7.54, 6.52],
                      [1., 2.70, 4.00],
                      [1., 2.50, 3.80],
                      [1., 1.15, 5.64],
                      [1., 4.22, 3.27],
                      [1., 1.41, 5.70],], order='F')

        X_orig = np.copy(X)
        y = np.array([1., 2., 3., 4., 5., 6.])

        @jit(nopython=True)
        def f2(X, y, test):
            if test:
                # never executed, but necessary to trigger the bug
                X = X[1:2, :]
            return np.linalg.lstsq(X, y)

        f2(X, y, False)
        np.testing.assert_allclose(X, X_orig)


class TestLinalgSolve(TestLinalgSystems):
    """
    Tests for np.linalg.solve.
    """

    @needs_lapack
    def test_linalg_solve(self):
        """
        Test np.linalg.solve
        """
        cfunc = jit(nopython=True)(solve_system)

        def check(a, b, **kwargs):
            expected = solve_system(a, b, **kwargs)
            got = cfunc(a, b, **kwargs)

            # check that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "F")

            use_reconstruction = False
            # try plain match of the result first
            try:
                np.testing.assert_array_almost_equal_nulp(
                    got, expected, nulp=10)
            except AssertionError:
                # plain match failed, test by reconstruction
                use_reconstruction = True

            # If plain match fails then reconstruction is used,
            # this checks that AX ~= B.
            # Plain match can fail due to numerical fuzziness associated
            # with system size and conditioning, or more simply from
            # numpy using double precision routines for computation that
            # could be done in single precision (which is what numba does).
            # Therefore minor differences in results can appear due to
            # e.g. numerical roundoff being different between two precisions.
            if use_reconstruction:
                # check they are dimensionally correct
                self.assertEqual(got.shape, expected.shape)

                # check AX=B
                rec = np.dot(a, got)
                resolution = np.finfo(a.dtype).resolution
                np.testing.assert_allclose(
                    b,
                    rec,
                    rtol=10 * resolution,
                    atol=100 * resolution  # zeros tend to be fuzzy
                )

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, b, **kwargs)

        # test: prime size squares
        sizes = [(1, 1), (3, 3), (7, 7)]

        # test loop
        for size, dtype, order in \
                product(sizes, self.dtypes, 'FC'):
            A = self.specific_sample_matrix(size, dtype, order)

            b_sizes = (1, 13)

            for b_size, b_order in product(b_sizes, 'FC'):
                # check 2D B
                B = self.specific_sample_matrix(
                    (A.shape[0], b_size), dtype, b_order)
                check(A, B)

                # check 1D B
                tmp = B[:, 0].copy(order=b_order)
                check(A, tmp)

        # check empty
        cfunc(np.empty((0, 0)), np.empty((0,)))

        # Test input validation
        ok = np.array([[1., 0.], [0., 1.]], dtype=np.float64)

        # check ok input is ok
        cfunc(ok, ok)

        # check bad inputs
        rn = "solve"

        # Wrong dtype
        bad = np.array([[1, 0], [0, 1]], dtype=np.int32)
        self.assert_wrong_dtype(rn, cfunc, (ok, bad))
        self.assert_wrong_dtype(rn, cfunc, (bad, ok))

        # different dtypes
        bad = np.array([[1, 2], [3, 4]], dtype=np.float32)
        self.assert_homogeneous_dtypes(rn, cfunc, (ok, bad))
        self.assert_homogeneous_dtypes(rn, cfunc, (bad, ok))

        # Dimension issue
        bad = np.array([1, 0], dtype=np.float64)
        self.assert_wrong_dimensions(rn, cfunc, (bad, ok))

        # no nans or infs
        bad = np.array([[1., 0., ], [np.inf, np.nan]], dtype=np.float64)
        self.assert_no_nan_or_inf(cfunc, (ok, bad))
        self.assert_no_nan_or_inf(cfunc, (bad, ok))

        # check 1D is accepted for B (2D is done previously)
        # and then that anything of higher dimension raises
        ok_oneD = np.array([1., 2.], dtype=np.float64)
        cfunc(ok, ok_oneD)
        bad = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], dtype=np.float64)
        self.assert_wrong_dimensions_1D(rn, cfunc, (ok, bad))

        # check an invalid system raises (1D and 2D cases checked)
        bad1D = np.array([1.], dtype=np.float64)
        bad2D = np.array([[1.], [2.], [3.]], dtype=np.float64)
        self.assert_dimensionally_invalid(cfunc, (ok, bad1D))
        self.assert_dimensionally_invalid(cfunc, (ok, bad2D))

        # check that a singular system raises
        bad2D = self.specific_sample_matrix((2, 2), np.float64, 'C', rank=1)
        self.assert_raise_on_singular(cfunc, (bad2D, ok))

    @needs_lapack
    def test_no_input_mutation(self):
        X = np.array([[1., 1, 1, 1],
                      [0., 1, 1, 1],
                      [0., 0, 1, 1],
                      [1., 0, 0, 1],], order='F')

        X_orig = np.copy(X)
        y = np.array([1., 2., 3., 4])
        y_orig = np.copy(y)

        @jit(nopython=True)
        def func(X, y, test):
            if test:
                # not executed, triggers A order in X
                X = X[1:2, :]
            return np.linalg.solve(X, y)

        expected = func.py_func(X, y, False)
        np.testing.assert_allclose(X, X_orig)
        np.testing.assert_allclose(y, y_orig)

        got = func(X, y, False)
        np.testing.assert_allclose(X, X_orig)
        np.testing.assert_allclose(y, y_orig)

        np.testing.assert_allclose(expected, got)


class TestLinalgPinv(TestLinalgBase):
    """
    Tests for np.linalg.pinv.
    """

    @needs_lapack
    def test_linalg_pinv(self):
        """
        Test np.linalg.pinv
        """
        cfunc = jit(nopython=True)(pinv_matrix)

        def check(a, **kwargs):
            expected = pinv_matrix(a, **kwargs)
            got = cfunc(a, **kwargs)

            # check that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "F")

            use_reconstruction = False
            # try plain match of each array to np first

            try:
                np.testing.assert_array_almost_equal_nulp(
                    got, expected, nulp=10)
            except AssertionError:
                # plain match failed, test by reconstruction
                use_reconstruction = True

            # If plain match fails then reconstruction is used.
            # This can occur due to numpy using double precision
            # LAPACK when single can be used, this creates round off
            # problems. Also, if the matrix has machine precision level
            # zeros in its singular values then the singular vectors are
            # likely to vary depending on round off.
            if use_reconstruction:

                # check they are dimensionally correct
                self.assertEqual(got.shape, expected.shape)

                # check pinv(A)*A~=eye
                # if the problem is numerical fuzz then this will probably
                # work, if the problem is rank deficiency then it won't!
                rec = np.dot(got, a)
                try:
                    self.assert_is_identity_matrix(rec)
                except AssertionError:
                    # check A=pinv(pinv(A))
                    resolution = 5 * np.finfo(a.dtype).resolution
                    rec = cfunc(got)
                    np.testing.assert_allclose(
                        rec,
                        a,
                        rtol=10 * resolution,
                        atol=100 * resolution  # zeros tend to be fuzzy
                    )
                    if a.shape[0] >= a.shape[1]:
                        # if it is overdetermined or fully determined
                        # use numba lstsq function (which is type specific) to
                        # compute the inverse and check against that.
                        lstsq = jit(nopython=True)(lstsq_system)
                        lstsq_pinv = lstsq(
                            a, np.eye(
                                a.shape[0]).astype(
                                a.dtype), **kwargs)[0]
                        np.testing.assert_allclose(
                            got,
                            lstsq_pinv,
                            rtol=10 * resolution,
                            atol=100 * resolution  # zeros tend to be fuzzy
                        )
                    # check the 2 norm of the difference is small
                    self.assertLess(np.linalg.norm(got - expected), resolution)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # test: column vector, tall, wide, square, row vector
        # prime sizes
        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]

        # When required, a specified condition number
        specific_cond = 10.

        # test loop
        for size, dtype, order in \
                product(sizes, self.dtypes, 'FC'):
            # check a full rank matrix
            a = self.specific_sample_matrix(size, dtype, order)
            check(a)

            m, n = size
            if m != 1 and n != 1:
                # check a rank deficient matrix
                minmn = min(m, n)
                a = self.specific_sample_matrix(size, dtype, order,
                                                condition=specific_cond)
                rcond = 1. / specific_cond
                approx_half_rank_rcond = minmn * rcond
                check(a, rcond=approx_half_rank_rcond)

        # check empty
        for sz in [(0, 1), (1, 0)]:
            check(np.empty(sz))

        rn = "pinv"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64),))

        # no nans or infs
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2., ], [np.inf, np.nan]],
                                            dtype=np.float64),))

    @needs_lapack
    def test_issue5870(self):
        # testing for mutation of input matrix
        @jit(nopython=True)
        def some_fn(v):
            return np.linalg.pinv(v[0])

        v_data = np.array([[1., 3, 2, 7,],
                           [-5, 4, 2, 3,],
                           [9, -3, 1, 1,],
                           [2, -2, 2, 8,]], order='F')

        v_orig = np.copy(v_data)
        reshaped_v = v_data.reshape((1, 4, 4))

        expected = some_fn.py_func(reshaped_v)
        np.testing.assert_allclose(v_data, v_orig)

        got = some_fn(reshaped_v)
        np.testing.assert_allclose(v_data, v_orig)

        np.testing.assert_allclose(expected, got)


class TestLinalgDetAndSlogdet(TestLinalgBase):
    """
    Tests for np.linalg.det. and np.linalg.slogdet.
    Exactly the same inputs are used for both tests as
    det() is a trivial function of slogdet(), the tests
    are therefore combined.
    """

    def check_det(self, cfunc, a, **kwargs):
        expected = det_matrix(a, **kwargs)
        got = cfunc(a, **kwargs)

        resolution = 5 * np.finfo(a.dtype).resolution

        # check the determinants are the same
        np.testing.assert_allclose(got, expected, rtol=resolution)

        # Ensure proper resource management
        with self.assertNoNRTLeak():
            cfunc(a, **kwargs)

    def check_slogdet(self, cfunc, a, **kwargs):
        expected = slogdet_matrix(a, **kwargs)
        got = cfunc(a, **kwargs)

        # As numba returns python floats types and numpy returns
        # numpy float types, some more adjustment and different
        # types of comparison than those used with array based
        # results is required.

        # check that the returned tuple is same length
        self.assertEqual(len(expected), len(got))
        # and that length is 2
        self.assertEqual(len(got), 2)

        # check that the domain of the results match
        for k in range(2):
            self.assertEqual(
                np.iscomplexobj(got[k]),
                np.iscomplexobj(expected[k]))

        # turn got[0] into the same dtype as `a`
        # this is so checking with nulp will work
        got_conv = a.dtype.type(got[0])
        np.testing.assert_array_almost_equal_nulp(
            got_conv, expected[0], nulp=10)
        # compare log determinant magnitude with a more fuzzy value
        # as numpy values come from higher precision lapack routines
        resolution = 5 * np.finfo(a.dtype).resolution
        np.testing.assert_allclose(
            got[1], expected[1], rtol=resolution, atol=resolution)

        # Ensure proper resource management
        with self.assertNoNRTLeak():
            cfunc(a, **kwargs)

    def do_test(self, rn, check, cfunc):

        # test: 1x1 as it is unusual, 4x4 as it is even and 7x7 as is it odd!
        sizes = [(1, 1), (4, 4), (7, 7)]

        # test loop
        for size, dtype, order in \
                product(sizes, self.dtypes, 'FC'):
            # check a full rank matrix
            a = self.specific_sample_matrix(size, dtype, order)
            check(cfunc, a)

        # use a matrix of zeros to trip xgetrf U(i,i)=0 singular test
        for dtype, order in product(self.dtypes, 'FC'):
            a = np.zeros((3, 3), dtype=dtype)
            check(cfunc, a)

        # check empty
        check(cfunc, np.empty((0, 0)))

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64),))

        # no nans or infs
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2., ], [np.inf, np.nan]],
                                            dtype=np.float64),))

    @needs_lapack
    def test_linalg_det(self):
        cfunc = jit(nopython=True)(det_matrix)
        self.do_test("det", self.check_det, cfunc)

    @needs_lapack
    def test_linalg_slogdet(self):
        cfunc = jit(nopython=True)(slogdet_matrix)
        self.do_test("slogdet", self.check_slogdet, cfunc)

    @needs_lapack
    def test_no_input_mutation(self):
        X = np.array([[1., 3, 2, 7,],
                      [-5, 4, 2, 3,],
                      [9, -3, 1, 1,],
                      [2, -2, 2, 8,]], order='F')

        X_orig = np.copy(X)

        @jit(nopython=True)
        def func(X, test):
            if test:
                # not executed, but necessary to trigger A ordering in X
                X = X[1:2, :]
            return np.linalg.slogdet(X)

        expected = func.py_func(X, False)
        np.testing.assert_allclose(X, X_orig)

        got = func(X, False)
        np.testing.assert_allclose(X, X_orig)

        np.testing.assert_allclose(expected, got)

# Use TestLinalgSystems as a base to get access to additional
# testing for 1 and 2D inputs.


class TestLinalgNorm(TestLinalgSystems):
    """
    Tests for np.linalg.norm.
    """

    @needs_lapack
    def test_linalg_norm(self):
        """
        Test np.linalg.norm
        """
        cfunc = jit(nopython=True)(norm_matrix)

        def check(a, **kwargs):
            expected = norm_matrix(a, **kwargs)
            got = cfunc(a, **kwargs)

            # All results should be in the real domain
            self.assertTrue(not np.iscomplexobj(got))

            resolution = 5 * np.finfo(a.dtype).resolution

            # check the norms are the same to the arg `a` precision
            np.testing.assert_allclose(got, expected, rtol=resolution)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # Check 1D inputs
        sizes = [1, 4, 7]
        nrm_types = [None, np.inf, -np.inf, 0, 1, -1, 2, -2, 5, 6.7, -4.3]

        # standard 1D input
        for size, dtype, nrm_type in \
                product(sizes, self.dtypes, nrm_types):
            a = self.sample_vector(size, dtype)
            check(a, ord=nrm_type)

        # sliced 1D input
        for dtype, nrm_type in \
                product(self.dtypes, nrm_types):
            a = self.sample_vector(10, dtype)[::3]
            check(a, ord=nrm_type)

        # Check 2D inputs:
        # test: column vector, tall, wide, square, row vector
        # prime sizes
        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]
        nrm_types = [None, np.inf, -np.inf, 1, -1, 2, -2]

        # standard 2D input
        for size, dtype, order, nrm_type in \
                product(sizes, self.dtypes, 'FC', nrm_types):
            # check a full rank matrix
            a = self.specific_sample_matrix(size, dtype, order)
            check(a, ord=nrm_type)

        # check 2D slices work for the case where xnrm2 is called from
        # BLAS (ord=None) to make sure it is working ok.
        nrm_types = [None]
        for dtype, nrm_type, order in \
                product(self.dtypes, nrm_types, 'FC'):
            a = self.specific_sample_matrix((17, 13), dtype, order)
            # contig for C order
            check(a[:3], ord=nrm_type)

            # contig for Fortran order
            check(a[:, 3:], ord=nrm_type)

            # contig for neither order
            check(a[1, 4::3], ord=nrm_type)

        # check that numba returns zero for empty arrays. Numpy returns zero
        # for most norm types and raises ValueError for +/-np.inf.
        # there is not a great deal of consistency in Numpy's response so
        # it is not being emulated in Numba
        for dtype, nrm_type, order in \
                product(self.dtypes, nrm_types, 'FC'):
            a = np.empty((0,), dtype=dtype, order=order)
            self.assertEqual(cfunc(a, nrm_type), 0.0)
            a = np.empty((0, 0), dtype=dtype, order=order)
            self.assertEqual(cfunc(a, nrm_type), 0.0)

        rn = "norm"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue, reuse the test from the TestLinalgSystems class
        self.assert_wrong_dimensions_1D(
            rn, cfunc, (np.ones(
                12, dtype=np.float64).reshape(
                2, 2, 3),))

        # no nans or infs for 2d case when SVD is used (e.g 2-norm)
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2.], [np.inf, np.nan]],
                                            dtype=np.float64), 2))

        # assert 2D input raises for an invalid norm kind kwarg
        self.assert_invalid_norm_kind(cfunc, (np.array([[1., 2.], [3., 4.]],
                                                       dtype=np.float64), 6))


class TestLinalgCond(TestLinalgBase):
    """
    Tests for np.linalg.cond.
    """

    @needs_lapack
    def test_linalg_cond(self):
        """
        Test np.linalg.cond
        """

        cfunc = jit(nopython=True)(cond_matrix)

        def check(a, **kwargs):
            expected = cond_matrix(a, **kwargs)
            got = cfunc(a, **kwargs)

            # All results should be in the real domain
            self.assertTrue(not np.iscomplexobj(got))

            resolution = 5 * np.finfo(a.dtype).resolution

            # check the cond is the same to the arg `a` precision
            np.testing.assert_allclose(got, expected, rtol=resolution)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # valid p values (used to indicate norm type)
        ps = [None, np.inf, -np.inf, 1, -1, 2, -2]
        sizes = [(3, 3), (7, 7)]

        for size, dtype, order, p in \
                product(sizes, self.dtypes, 'FC', ps):
            a = self.specific_sample_matrix(size, dtype, order)
            check(a, p=p)

        # When p=None non-square matrices are accepted.
        sizes = [(7, 1), (11, 5), (5, 11), (1, 7)]
        for size, dtype, order in \
                product(sizes, self.dtypes, 'FC'):
            a = self.specific_sample_matrix(size, dtype, order)
            check(a)

        # empty
        for sz in [(0, 1), (1, 0), (0, 0)]:
            self.assert_raise_on_empty(cfunc, (np.empty(sz),))

        # singular systems to trip divide-by-zero
        x = np.array([[1, 0], [0, 0]], dtype=np.float64)
        check(x)
        check(x, p=2)
        x = np.array([[0, 0], [0, 0]], dtype=np.float64)
        check(x, p=-2)

        # try an ill-conditioned system with 2-norm, make sure np raises an
        # overflow warning as the result is `+inf` and that the result from
        # numba matches.
        with warnings.catch_warnings():
            a = np.array([[1.e308, 0], [0, 0.1]], dtype=np.float64)
            warnings.simplefilter("ignore", RuntimeWarning)
            check(a)

        rn = "cond"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64),))

        # no nans or infs when p="None" (default for kwarg).
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2., ], [np.inf, np.nan]],
                                            dtype=np.float64),))

        # assert raises for an invalid norm kind kwarg
        self.assert_invalid_norm_kind(cfunc, (np.array([[1., 2.], [3., 4.]],
                                                       dtype=np.float64), 6))


class TestLinalgMatrixRank(TestLinalgSystems):
    """
    Tests for np.linalg.matrix_rank.
    """

    @needs_lapack
    def test_linalg_matrix_rank(self):
        """
        Test np.linalg.matrix_rank
        """

        cfunc = jit(nopython=True)(matrix_rank_matrix)

        def check(a, **kwargs):
            expected = matrix_rank_matrix(a, **kwargs)
            got = cfunc(a, **kwargs)

            # Ranks are integral so comparison should be trivial.
            # check the rank is the same
            np.testing.assert_allclose(got, expected)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]

        for size, dtype, order in \
                product(sizes, self.dtypes, 'FC'):
            # check full rank system
            a = self.specific_sample_matrix(size, dtype, order)
            check(a)

            # If the system is a matrix, check rank deficiency is reported
            # correctly. Check all ranks from 0 to (full rank - 1).
            tol = 1e-13
            # first check 1 to (full rank - 1)
            for k in range(1, min(size) - 1):
                # check rank k
                a = self.specific_sample_matrix(size, dtype, order, rank=k)
                self.assertEqual(cfunc(a), k)
                check(a)
                # check provision of a tolerance works as expected
                # create a (m x n) diagonal matrix with a singular value
                # guaranteed below the tolerance 1e-13
                m, n = a.shape
                a[:, :] = 0.  # reuse `a`'s memory
                idx = np.nonzero(np.eye(m, n))
                if np.iscomplexobj(a):
                    b = 1. + np.random.rand(k) + 1.j +\
                        1.j * np.random.rand(k)
                    # min singular value is sqrt(2)*1e-14
                    b[0] = 1e-14 + 1e-14j
                else:
                    b = 1. + np.random.rand(k)
                    b[0] = 1e-14  # min singular value is 1e-14
                a[idx[0][:k], idx[1][:k]] = b.astype(dtype)
                # rank should be k-1 (as tol is present)
                self.assertEqual(cfunc(a, tol), k - 1)
                check(a, tol=tol)
            # then check zero rank
            a[:, :] = 0.
            self.assertEqual(cfunc(a), 0)
            check(a)
            # add in a singular value that is small
            if np.iscomplexobj(a):
                a[-1, -1] = 1e-14 + 1e-14j
            else:
                a[-1, -1] = 1e-14
            # check the system has zero rank to a given tolerance
            self.assertEqual(cfunc(a, tol), 0)
            check(a, tol=tol)

        # check the zero vector returns rank 0 and a nonzero vector
        # returns rank 1.
        for dt in self.dtypes:
            a = np.zeros((5), dtype=dt)
            self.assertEqual(cfunc(a), 0)
            check(a)
            # make it a nonzero vector
            a[0] = 1.
            self.assertEqual(cfunc(a), 1)
            check(a)

        # empty
        for sz in [(0, 1), (1, 0), (0, 0)]:
            for tol in [None, 1e-13]:
                self.assert_raise_on_empty(cfunc, (np.empty(sz), tol))

        rn = "matrix_rank"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32),))

        # Dimension issue
        self.assert_wrong_dimensions_1D(
            rn, cfunc, (np.ones(
                12, dtype=np.float64).reshape(
                2, 2, 3),))

        # no nans or infs for 2D case
        self.assert_no_nan_or_inf(cfunc,
                                  (np.array([[1., 2., ], [np.inf, np.nan]],
                                            dtype=np.float64),))

    @needs_lapack
    def test_no_input_mutation(self):
        # this is here to test no input mutation by
        # numba.np.linalg._compute_singular_values
        # which is the workhorse for norm with 2d input, rank and cond.

        X = np.array([[1., 3, 2, 7,],
                      [-5, 4, 2, 3,],
                      [9, -3, 1, 1,],
                      [2, -2, 2, 8,]], order='F')

        X_orig = np.copy(X)

        @jit(nopython=True)
        def func(X, test):
            if test:
                # not executed, but necessary to trigger A ordering in X
                X = X[1:2, :]
            return np.linalg.matrix_rank(X)

        expected = func.py_func(X, False)
        np.testing.assert_allclose(X, X_orig)

        got = func(X, False)
        np.testing.assert_allclose(X, X_orig)

        np.testing.assert_allclose(expected, got)


class TestLinalgMatrixPower(TestLinalgBase):
    """
    Tests for np.linalg.matrix_power.
    """

    def assert_int_exponenent(self, cfunc, args):
        # validate first arg is ok
        cfunc(args[0], 1)
        # pass in both args and assert fail
        with self.assertRaises(errors.TypingError):
            cfunc(*args)

    @needs_lapack
    def test_linalg_matrix_power(self):
        cfunc = jit(nopython=True)(matrix_power_matrix)

        def check(a, pwr):
            expected = matrix_power_matrix(a, pwr)
            got = cfunc(a, pwr)

            # check that the computed results are contig and in the same way
            self.assert_contig_sanity(got, "C")

            res = 7 * np.finfo(a.dtype).resolution
            np.testing.assert_allclose(got, expected, rtol=res, atol=res)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, pwr)

        sizes = [(1, 1), (5, 5), (7, 7)]
        powers = [-33, -17] + list(range(-10, 10)) + [17, 33]

        for size, pwr, dtype, order in \
                product(sizes, powers, self.dtypes, 'FC'):
            a = self.specific_sample_matrix(size, dtype, order)
            check(a, pwr)
            a = np.empty((0, 0), dtype=dtype, order=order)
            check(a, pwr)

        rn = "matrix_power"

        # Wrong dtype
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32), 1))

        # not an int power
        self.assert_wrong_dtype(rn, cfunc,
                                (np.ones((2, 2), dtype=np.int32), 1))

        # non square system
        args = (np.ones((3, 5)), 1)
        msg = 'input must be a square array'
        self.assert_error(cfunc, args, msg)

        # Dimension issue
        self.assert_wrong_dimensions(rn, cfunc,
                                     (np.ones(10, dtype=np.float64), 1))

        # non-integer supplied as exponent
        self.assert_int_exponenent(cfunc, (np.ones((2, 2)), 1.2))

        # singular matrix is not invertible
        self.assert_raise_on_singular(cfunc, (np.array([[0., 0], [1, 1]]), -1))


class TestTrace(TestLinalgBase):
    """
    Tests for np.trace.
    """

    def setUp(self):
        super(TestTrace, self).setUp()
        # compile two versions, one with and one without the offset kwarg
        self.cfunc_w_offset = jit(nopython=True)(trace_matrix)
        self.cfunc_no_offset = jit(nopython=True)(trace_matrix_no_offset)

    def assert_int_offset(self, cfunc, a, **kwargs):
        # validate first arg is ok
        cfunc(a)
        # pass in kwarg and assert fail
        with self.assertRaises(errors.TypingError):
            cfunc(a, **kwargs)

    def test_trace(self):

        def check(a, **kwargs):
            if 'offset' in kwargs:
                expected = trace_matrix(a, **kwargs)
                cfunc = self.cfunc_w_offset
            else:
                expected = trace_matrix_no_offset(a, **kwargs)
                cfunc = self.cfunc_no_offset

            got = cfunc(a, **kwargs)

            res = 5 * np.finfo(a.dtype).resolution
            np.testing.assert_allclose(got, expected, rtol=res, atol=res)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, **kwargs)

        # test: column vector, tall, wide, square, row vector
        # prime sizes
        sizes = [(7, 1), (11, 5), (5, 11), (3, 3), (1, 7)]

        # offsets to cover the range of the matrix sizes above
        offsets = [-13, -12, -11] + list(range(-10, 10)) + [11, 12, 13]

        for size, offset, dtype, order in \
                product(sizes, offsets, self.dtypes, 'FC'):
            a = self.specific_sample_matrix(size, dtype, order)
            check(a, offset=offset)
            if offset == 0:
                check(a)
            a = np.empty((0, 0), dtype=dtype, order=order)
            check(a, offset=offset)
            if offset == 0:
                check(a)

        rn = "trace"

        # Dimension issue
        self.assert_wrong_dimensions(rn, self.cfunc_w_offset,
                                     (np.ones(10, dtype=np.float64), 1), False)
        self.assert_wrong_dimensions(rn, self.cfunc_no_offset,
                                     (np.ones(10, dtype=np.float64),), False)

        # non-integer supplied as exponent
        self.assert_int_offset(
            self.cfunc_w_offset, np.ones(
                (2, 2)), offset=1.2)

    def test_trace_w_optional_input(self):
        "Issue 2314"
        @jit("(optional(float64[:,:]),)", nopython=True)
        def tested(a):
            return np.trace(a)

        a = np.ones((5, 5), dtype=np.float64)
        tested(a)

        with self.assertRaises(TypeError) as raises:
            tested(None)

        errmsg = str(raises.exception)
        self.assertEqual('expected array(float64, 2d, A), got None', errmsg)


class TestBasics(TestLinalgSystems):  # TestLinalgSystems for 1d test

    order1 = cycle(['F', 'C', 'C', 'F'])
    order2 = cycle(['C', 'F', 'C', 'F'])

    # test: column vector, matrix, row vector, 1d sizes
    # (7, 1, 3) and two scalars
    sizes = [(7, 1), (3, 3), (1, 7), (7,), (1,), (3,), 3., 5.]

    def _assert_wrong_dim(self, rn, cfunc):
        # Dimension issue
        self.assert_wrong_dimensions_1D(
            rn, cfunc, (np.array([[[1]]], dtype=np.float64), np.ones(1)), False)
        self.assert_wrong_dimensions_1D(
            rn, cfunc, (np.ones(1), np.array([[[1]]], dtype=np.float64)), False)

    def _gen_input(self, size, dtype, order):
        if not isinstance(size, tuple):
            return size
        else:
            if len(size) == 1:
                return self.sample_vector(size[0], dtype)
            else:
                return self.sample_vector(
                    size[0] * size[1],
                    dtype).reshape(
                    size, order=order)

    def _get_input(self, size1, size2, dtype):
        a = self._gen_input(size1, dtype, next(self.order1))
        b = self._gen_input(size2, dtype, next(self.order2))
        # force domain consistency as underlying ufuncs require it
        if np.iscomplexobj(a):
            b = b + 1j
        if np.iscomplexobj(b):
            a = a + 1j
        return (a, b)

    def test_outer(self):
        cfunc = jit(nopython=True)(outer_matrix)

        def check(a, b, **kwargs):

            # check without kwargs
            expected = outer_matrix(a, b)
            got = cfunc(a, b)

            res = 5 * np.finfo(np.asarray(a).dtype).resolution
            np.testing.assert_allclose(got, expected, rtol=res, atol=res)

            # if kwargs present check with them too
            if 'out' in kwargs:
                got = cfunc(a, b, **kwargs)
                np.testing.assert_allclose(got, expected, rtol=res,
                                           atol=res)
                np.testing.assert_allclose(kwargs['out'], expected,
                                           rtol=res, atol=res)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, b, **kwargs)

        dts = cycle(self.dtypes)
        for size1, size2 in product(self.sizes, self.sizes):
            dtype = next(dts)
            (a, b) = self._get_input(size1, size2, dtype)
            check(a, b)
            c = np.empty((np.asarray(a).size, np.asarray(b).size),
                            dtype=np.asarray(a).dtype)
            check(a, b, out=c)

        self._assert_wrong_dim("outer", cfunc)

    def test_kron(self):
        cfunc = jit(nopython=True)(kron_matrix)

        def check(a, b, **kwargs):

            expected = kron_matrix(a, b)
            got = cfunc(a, b)

            res = 5 * np.finfo(np.asarray(a).dtype).resolution
            np.testing.assert_allclose(got, expected, rtol=res, atol=res)

            # Ensure proper resource management
            with self.assertNoNRTLeak():
                cfunc(a, b)

        for size1, size2, dtype in \
                product(self.sizes, self.sizes, self.dtypes):
            (a, b) = self._get_input(size1, size2, dtype)
            check(a, b)

        self._assert_wrong_dim("kron", cfunc)

        args = (np.empty(10)[::2], np.empty(10)[::2])
        msg = "only supports 'C' or 'F' layout"
        self.assert_error(cfunc, args, msg, err=errors.TypingError)


class TestHelpers(TestCase):
    def test_copy_to_fortran_order(self):
        from numba.np.linalg import _copy_to_fortran_order

        def check(udt, expectfn, shapes, dtypes, orders):
            for shape, dtype, order in product(shapes, dtypes, orders):
                a = np.arange(np.prod(shape)).reshape(shape, order=order)

                r = udt(a)
                # check correct operation
                self.assertPreciseEqual(expectfn(a), r)
                # check new copy has made
                self.assertNotEqual(a.ctypes.data, r.ctypes.data)

        @njit
        def direct_call(a):
            return _copy_to_fortran_order(a)

        shapes = [(3, 4), (3, 2, 5)]
        dtypes = [np.intp]
        orders = ['C', 'F']
        check(direct_call, np.asfortranarray, shapes, dtypes, orders)


        @njit
        def slice_to_any(a):
            # make a 'any' layout slice
            sliced = a[::2][0]
            return _copy_to_fortran_order(sliced)

        shapes = [(3, 3, 4), (3, 3, 2, 5)]
        dtypes = [np.intp]
        orders = ['C', 'F']

        def expected_slice_to_any(a):
            # make a 'any' layout slice
            sliced = a[::2][0]
            return np.asfortranarray(sliced)

        check(slice_to_any, expected_slice_to_any, shapes, dtypes, orders)

if __name__ == '__main__':
    unittest.main()
