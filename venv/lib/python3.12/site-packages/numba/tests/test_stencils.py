#
# Copyright (c) 2017 Intel Corporation
# SPDX-License-Identifier: BSD-2-Clause
#

import numpy as np
from contextlib import contextmanager

import numba
from numba import njit, stencil
from numba.core import types, registry
from numba.core.compiler import compile_extra, Flags
from numba.core.cpu import ParallelOptions
from numba.tests.support import skip_parfors_unsupported, _32bit
from numba.core.errors import LoweringError, TypingError, NumbaValueError
import unittest


skip_unsupported = skip_parfors_unsupported


@stencil
def stencil1_kernel(a):
    return 0.25 * (a[0, 1] + a[1, 0] + a[0, -1] + a[-1, 0])


@stencil(neighborhood=((-5, 0), ))
def stencil2_kernel(a):
    cum = a[-5]
    for i in range(-4, 1):
        cum += a[i]
    return 0.3 * cum


@stencil(cval=1.0)
def stencil3_kernel(a):
    return 0.25 * a[-2, 2]


@stencil
def stencil_multiple_input_kernel(a, b):
    return 0.25 * (a[0, 1] + a[1, 0] + a[0, -1] + a[-1, 0] +
                   b[0, 1] + b[1, 0] + b[0, -1] + b[-1, 0])


@stencil
def stencil_multiple_input_kernel_var(a, b, w):
    return w * (a[0, 1] + a[1, 0] + a[0, -1] + a[-1, 0] +
                b[0, 1] + b[1, 0] + b[0, -1] + b[-1, 0])


@stencil
def stencil_multiple_input_mixed_types_2d(a, b, f):
    return a[0, 0] if f[0, 0] else b[0, 0]


@stencil(standard_indexing=("b",))
def stencil_with_standard_indexing_1d(a, b):
    return a[-1] * b[0] + a[0] * b[1]


@stencil(standard_indexing=("b",))
def stencil_with_standard_indexing_2d(a, b):
    return (a[0, 1] * b[0, 1] + a[1, 0] * b[1, 0]
            + a[0, -1] * b[0, -1] + a[-1, 0] * b[-1, 0])


@njit
def addone_njit(a):
    return a + 1


if not _32bit: # prevent compilation on unsupported 32bit targets
    @njit(parallel=True)
    def addone_pjit(a):
        return a + 1


class TestStencilBase(unittest.TestCase):

    _numba_parallel_test_ = False

    def __init__(self, *args):
        # flags for njit()
        self.cflags = Flags()
        self.cflags.nrt = True

        super(TestStencilBase, self).__init__(*args)

    def _compile_this(self, func, sig, flags):
        return compile_extra(registry.cpu_target.typing_context,
                             registry.cpu_target.target_context, func, sig,
                             None, flags, {})

    def compile_parallel(self, func, sig, **kws):
        flags = Flags()
        flags.nrt = True
        options = True if not kws else kws
        flags.auto_parallel = ParallelOptions(options)
        return self._compile_this(func, sig, flags)

    def compile_njit(self, func, sig):
        return self._compile_this(func, sig, flags=self.cflags)

    def compile_all(self, pyfunc, *args, **kwargs):
        sig = tuple([numba.typeof(x) for x in args])
        # compile with parallel=True
        cpfunc = self.compile_parallel(pyfunc, sig)
        # compile a standard njit of the original function
        cfunc = self.compile_njit(pyfunc, sig)
        return cfunc, cpfunc

    def check(self, no_stencil_func, pyfunc, *args):
        cfunc, cpfunc = self.compile_all(pyfunc, *args)
        # results without stencil macro
        expected = no_stencil_func(*args)
        # python result
        py_output = pyfunc(*args)

        # njit result
        njit_output = cfunc.entry_point(*args)

        # parfor result
        parfor_output = cpfunc.entry_point(*args)

        np.testing.assert_almost_equal(py_output, expected, decimal=3)
        np.testing.assert_almost_equal(njit_output, expected, decimal=3)
        np.testing.assert_almost_equal(parfor_output, expected, decimal=3)

        # make sure parfor set up scheduling
        self.assertIn('@do_scheduling', cpfunc.library.get_llvm_str())


class TestStencil(TestStencilBase):

    def __init__(self, *args, **kwargs):
        super(TestStencil, self).__init__(*args, **kwargs)

    @skip_unsupported
    def test_stencil1(self):
        """Tests whether the optional out argument to stencil calls works.
        """
        def test_with_out(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.zeros(n**2).reshape((n, n))
            B = stencil1_kernel(A, out=B)
            return B

        def test_without_out(n):
            A = np.arange(n**2).reshape((n, n))
            B = stencil1_kernel(A)
            return B

        def test_impl_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.zeros(n**2).reshape((n, n))
            for i in range(1, n - 1):
                for j in range(1, n - 1):
                    B[i, j] = 0.25 * (A[i, j + 1] +
                                      A[i + 1, j] + A[i, j - 1] + A[i - 1, j])
            return B

        n = 100
        self.check(test_impl_seq, test_with_out, n)
        self.check(test_impl_seq, test_without_out, n)

    @skip_unsupported
    def test_stencil2(self):
        """Tests whether the optional neighborhood argument to the stencil
        decorate works.
        """
        def test_seq(n):
            A = np.arange(n)
            B = stencil2_kernel(A)
            return B

        def test_impl_seq(n):
            A = np.arange(n)
            B = np.zeros(n)
            for i in range(5, len(A)):
                B[i] = 0.3 * sum(A[i - 5:i + 1])
            return B

        n = 100
        self.check(test_impl_seq, test_seq, n)
        # variable length neighborhood in numba.stencil call
        # only supported in parallel path

        def test_seq(n, w):
            A = np.arange(n)

            def stencil2_kernel(a, w):
                cum = a[-w]
                for i in range(-w + 1, w + 1):
                    cum += a[i]
                return 0.3 * cum
            B = numba.stencil(stencil2_kernel, neighborhood=((-w, w), ))(A, w)
            return B

        def test_impl_seq(n, w):
            A = np.arange(n)
            B = np.zeros(n)
            for i in range(w, len(A) - w):
                B[i] = 0.3 * sum(A[i - w:i + w + 1])
            return B
        n = 100
        w = 5
        cpfunc = self.compile_parallel(test_seq, (types.intp, types.intp))
        expected = test_impl_seq(n, w)
        # parfor result
        parfor_output = cpfunc.entry_point(n, w)
        np.testing.assert_almost_equal(parfor_output, expected, decimal=3)
        self.assertIn('@do_scheduling', cpfunc.library.get_llvm_str())
        # test index_offsets

        def test_seq(n, w, offset):
            A = np.arange(n)

            def stencil2_kernel(a, w):
                cum = a[-w + 1]
                for i in range(-w + 1, w + 1):
                    cum += a[i + 1]
                return 0.3 * cum
            B = numba.stencil(stencil2_kernel, neighborhood=((-w, w), ),
                              index_offsets=(-offset, ))(A, w)
            return B

        offset = 1
        cpfunc = self.compile_parallel(test_seq, (types.intp, types.intp,
                                                  types.intp))
        parfor_output = cpfunc.entry_point(n, w, offset)
        np.testing.assert_almost_equal(parfor_output, expected, decimal=3)
        self.assertIn('@do_scheduling', cpfunc.library.get_llvm_str())
        # test slice in kernel

        def test_seq(n, w, offset):
            A = np.arange(n)

            def stencil2_kernel(a, w):
                return 0.3 * np.sum(a[-w + 1:w + 2])
            B = numba.stencil(stencil2_kernel, neighborhood=((-w, w), ),
                              index_offsets=(-offset, ))(A, w)
            return B

        offset = 1
        cpfunc = self.compile_parallel(test_seq, (types.intp, types.intp,
                                                  types.intp))
        parfor_output = cpfunc.entry_point(n, w, offset)
        np.testing.assert_almost_equal(parfor_output, expected, decimal=3)
        self.assertIn('@do_scheduling', cpfunc.library.get_llvm_str())

    @skip_unsupported
    def test_stencil3(self):
        """Tests whether a non-zero optional cval argument to the stencil
        decorator works.  Also tests integer result type.
        """
        def test_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = stencil3_kernel(A)
            return B

        test_njit = njit(test_seq)
        test_par = njit(test_seq, parallel=True)

        n = 5
        seq_res = test_seq(n)
        njit_res = test_njit(n)
        par_res = test_par(n)

        self.assertTrue(seq_res[0, 0] == 1.0 and seq_res[4, 4] == 1.0)
        self.assertTrue(njit_res[0, 0] == 1.0 and njit_res[4, 4] == 1.0)
        self.assertTrue(par_res[0, 0] == 1.0 and par_res[4, 4] == 1.0)

    @skip_unsupported
    def test_stencil_standard_indexing_1d(self):
        """Tests standard indexing with a 1d array.
        """
        def test_seq(n):
            A = np.arange(n)
            B = [3.0, 7.0]
            C = stencil_with_standard_indexing_1d(A, B)
            return C

        def test_impl_seq(n):
            A = np.arange(n)
            B = [3.0, 7.0]
            C = np.zeros(n)

            for i in range(1, n):
                C[i] = A[i - 1] * B[0] + A[i] * B[1]
            return C

        n = 100
        self.check(test_impl_seq, test_seq, n)

    @skip_unsupported
    def test_stencil_standard_indexing_2d(self):
        """Tests standard indexing with a 2d array and multiple stencil calls.
        """
        def test_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.ones((3, 3))
            C = stencil_with_standard_indexing_2d(A, B)
            D = stencil_with_standard_indexing_2d(C, B)
            return D

        def test_impl_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.ones((3, 3))
            C = np.zeros(n**2).reshape((n, n))
            D = np.zeros(n**2).reshape((n, n))

            for i in range(1, n - 1):
                for j in range(1, n - 1):
                    C[i, j] = (A[i, j + 1] * B[0, 1] + A[i + 1, j] * B[1, 0] +
                               A[i, j - 1] * B[0, -1] + A[i - 1, j] * B[-1, 0])
            for i in range(1, n - 1):
                for j in range(1, n - 1):
                    D[i, j] = (C[i, j + 1] * B[0, 1] + C[i + 1, j] * B[1, 0] +
                               C[i, j - 1] * B[0, -1] + C[i - 1, j] * B[-1, 0])
            return D

        n = 5
        self.check(test_impl_seq, test_seq, n)

    @skip_unsupported
    def test_stencil_multiple_inputs(self):
        """Tests whether multiple inputs of the same size work.
        """
        def test_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.arange(n**2).reshape((n, n))
            C = stencil_multiple_input_kernel(A, B)
            return C

        def test_impl_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.arange(n**2).reshape((n, n))
            C = np.zeros(n**2).reshape((n, n))
            for i in range(1, n - 1):
                for j in range(1, n - 1):
                    C[i, j] = 0.25 * \
                        (A[i, j + 1] + A[i + 1, j]
                         + A[i, j - 1] + A[i - 1, j]
                         + B[i, j + 1] + B[i + 1, j]
                         + B[i, j - 1] + B[i - 1, j])
            return C

        n = 3
        self.check(test_impl_seq, test_seq, n)
        # test stencil with a non-array input

        def test_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.arange(n**2).reshape((n, n))
            w = 0.25
            C = stencil_multiple_input_kernel_var(A, B, w)
            return C
        self.check(test_impl_seq, test_seq, n)

    @skip_unsupported
    def test_stencil_mixed_types(self):
        def test_impl_seq(n):
            A = np.arange(n ** 2).reshape((n, n))
            B = n ** 2 - np.arange(n ** 2).reshape((n, n))
            S = np.eye(n, dtype=np.bool_)
            O = np.zeros((n, n), dtype=A.dtype)
            for i in range(0, n):
                for j in range(0, n):
                    O[i, j] = A[i, j] if S[i, j] else B[i, j]
            return O

        def test_seq(n):
            A = np.arange(n ** 2).reshape((n, n))
            B = n ** 2 - np.arange(n ** 2).reshape((n, n))
            S = np.eye(n, dtype=np.bool_)
            O = stencil_multiple_input_mixed_types_2d(A, B, S)
            return O

        n = 3
        self.check(test_impl_seq, test_seq, n)

    @skip_unsupported
    def test_stencil_call(self):
        """Tests 2D numba.stencil calls.
        """
        def test_impl1(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.zeros(n**2).reshape((n, n))
            numba.stencil(lambda a: 0.25 * (a[0, 1] + a[1, 0] + a[0, -1]
                                            + a[-1, 0]))(A, out=B)
            return B

        def test_impl2(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.zeros(n**2).reshape((n, n))

            def sf(a):
                return 0.25 * (a[0, 1] + a[1, 0] + a[0, -1] + a[-1, 0])
            B = numba.stencil(sf)(A)
            return B

        def test_impl_seq(n):
            A = np.arange(n**2).reshape((n, n))
            B = np.zeros(n**2).reshape((n, n))
            for i in range(1, n - 1):
                for j in range(1, n - 1):
                    B[i, j] = 0.25 * (A[i, j + 1] + A[i + 1, j]
                                      + A[i, j - 1] + A[i - 1, j])
            return B

        n = 100
        self.check(test_impl_seq, test_impl1, n)
        self.check(test_impl_seq, test_impl2, n)

    @skip_unsupported
    def test_stencil_call_1D(self):
        """Tests 1D numba.stencil calls.
        """
        def test_impl(n):
            A = np.arange(n)
            B = np.zeros(n)
            numba.stencil(lambda a: 0.3 * (a[-1] + a[0] + a[1]))(A, out=B)
            return B

        def test_impl_seq(n):
            A = np.arange(n)
            B = np.zeros(n)
            for i in range(1, n - 1):
                B[i] = 0.3 * (A[i - 1] + A[i] + A[i + 1])
            return B

        n = 100
        self.check(test_impl_seq, test_impl, n)

    @skip_unsupported
    def test_stencil_call_const(self):
        """Tests numba.stencil call that has an index that can be inferred as
        constant from a unary expr. Otherwise, this would raise an error since
        neighborhood length is not specified.
        """
        def test_impl1(n):
            A = np.arange(n)
            B = np.zeros(n)
            c = 1
            numba.stencil(lambda a,c : 0.3 * (a[-c] + a[0] + a[c]))(A, c, out=B)
            return B

        def test_impl2(n):
            A = np.arange(n)
            B = np.zeros(n)
            c = 2
            numba.stencil(
                lambda a,c : 0.3 * (a[1 - c] + a[0] + a[c - 1]))(A, c, out=B)
            return B

        # recursive expr case
        def test_impl3(n):
            A = np.arange(n)
            B = np.zeros(n)
            c = 2
            numba.stencil(
                lambda a,c : 0.3 * (a[-c + 1] + a[0] + a[c - 1]))(A, c, out=B)
            return B

        # multi-constant case
        def test_impl4(n):
            A = np.arange(n)
            B = np.zeros(n)
            d = 1
            c = 2
            numba.stencil(
                lambda a,c,d : 0.3 * (a[-c + d] + a[0] + a[c - d]))(A, c, d,
                                                                    out=B)
            return B

        def test_impl_seq(n):
            A = np.arange(n)
            B = np.zeros(n)
            c = 1
            for i in range(1, n - 1):
                B[i] = 0.3 * (A[i - c] + A[i] + A[i + c])
            return B

        n = 100
        # constant inference is only possible in parallel path
        cpfunc1 = self.compile_parallel(test_impl1, (types.intp,))
        cpfunc2 = self.compile_parallel(test_impl2, (types.intp,))
        cpfunc3 = self.compile_parallel(test_impl3, (types.intp,))
        cpfunc4 = self.compile_parallel(test_impl4, (types.intp,))
        expected = test_impl_seq(n)
        # parfor result
        parfor_output1 = cpfunc1.entry_point(n)
        parfor_output2 = cpfunc2.entry_point(n)
        parfor_output3 = cpfunc3.entry_point(n)
        parfor_output4 = cpfunc4.entry_point(n)
        np.testing.assert_almost_equal(parfor_output1, expected, decimal=3)
        np.testing.assert_almost_equal(parfor_output2, expected, decimal=3)
        np.testing.assert_almost_equal(parfor_output3, expected, decimal=3)
        np.testing.assert_almost_equal(parfor_output4, expected, decimal=3)

        # check error in regular Python path
        with self.assertRaises(NumbaValueError) as e:
            test_impl4(4)

        self.assertIn("stencil kernel index is not constant, "
                      "'neighborhood' option required", str(e.exception))
        # check error in njit path
        # TODO: ValueError should be thrown instead of LoweringError
        with self.assertRaises((LoweringError, NumbaValueError)) as e:
            njit(test_impl4)(4)

        self.assertIn("stencil kernel index is not constant, "
                      "'neighborhood' option required", str(e.exception))

    @skip_unsupported
    def test_stencil_parallel_off(self):
        """Tests 1D numba.stencil calls without parallel translation
           turned off.
        """
        def test_impl(A):
            return numba.stencil(lambda a: 0.3 * (a[-1] + a[0] + a[1]))(A)

        cpfunc = self.compile_parallel(test_impl, (numba.float64[:],),
                                       stencil=False)
        self.assertNotIn('@do_scheduling', cpfunc.library.get_llvm_str())

    @skip_unsupported
    def test_stencil_nested1(self):
        """Tests whether nested stencil decorator works.
        """
        @njit(parallel=True)
        def test_impl(n):
            @stencil
            def fun(a):
                c = 2
                return a[-c + 1]
            B = fun(n)
            return B

        def test_impl_seq(n):
            B = np.zeros(len(n), dtype=int)
            for i in range(1, len(n)):
                B[i] = n[i - 1]
            return B

        n = np.arange(10)
        np.testing.assert_equal(test_impl(n), test_impl_seq(n))

    @skip_unsupported
    def test_out_kwarg_w_cval(self):
        """ Issue #3518, out kwarg did not work with cval."""
        # test const value that matches the arg dtype, and one that can be cast
        const_vals = [7, 7.0]

        def kernel(a):
            return (a[0, 0] - a[1, 0])

        for const_val in const_vals:
            stencil_fn = numba.stencil(kernel, cval=const_val)

            def wrapped():
                A = np.arange(12).reshape((3, 4))
                ret = np.ones_like(A)
                stencil_fn(A, out=ret)
                return ret

            # stencil function case
            A = np.arange(12).reshape((3, 4))
            expected = np.full_like(A, -4)
            expected[-1, :] = const_val
            ret = np.ones_like(A)
            stencil_fn(A, out=ret)
            np.testing.assert_almost_equal(ret, expected)

            # wrapped function case, check njit, then njit(parallel=True)
            impls = self.compile_all(wrapped,)
            for impl in impls:
                got = impl.entry_point()
                np.testing.assert_almost_equal(got, expected)

        # now check exceptions for cval dtype mismatch with out kwarg dtype
        stencil_fn = numba.stencil(kernel, cval=1j)

        def wrapped():
            A = np.arange(12).reshape((3, 4))
            ret = np.ones_like(A)
            stencil_fn(A, out=ret)
            return ret

        A = np.arange(12).reshape((3, 4))
        ret = np.ones_like(A)
        with self.assertRaises(NumbaValueError) as e:
            stencil_fn(A, out=ret)
        msg = "cval type does not match stencil return type."
        self.assertIn(msg, str(e.exception))

        for compiler in [self.compile_njit, self.compile_parallel]:
            try:
                compiler(wrapped,())
            except (NumbaValueError, LoweringError) as e:
                self.assertIn(msg, str(e))
            else:
                raise AssertionError("Expected error was not raised")

    @skip_unsupported
    def test_out_kwarg_w_cval_np_attr(self):
        """ Test issue #7286 where the cval is a np attr/string-based numerical
        constant"""
        for cval in (np.nan, np.inf, -np.inf, float('inf'), -float('inf')):
            def kernel(a):
                return (a[0, 0] - a[1, 0])

            stencil_fn = numba.stencil(kernel, cval=cval)

            def wrapped():
                A = np.arange(12.).reshape((3, 4))
                ret = np.ones_like(A)
                stencil_fn(A, out=ret)
                return ret

            # stencil function case
            A = np.arange(12.).reshape((3, 4))
            expected = np.full_like(A, -4)
            expected[-1, :] = cval
            ret = np.ones_like(A)
            stencil_fn(A, out=ret)
            np.testing.assert_almost_equal(ret, expected)

            # wrapped function case, check njit, then njit(parallel=True)
            impls = self.compile_all(wrapped,)
            for impl in impls:
                got = impl.entry_point()
                np.testing.assert_almost_equal(got, expected)


@skip_unsupported
class TestManyStencils(TestStencilBase):
    # NOTE: the original implementation of this test used manipulations of the
    # Python AST repr of a kernel to create another implementation of the
    # stencil being tested so to act as another reference point when
    # comparing the various forms of @stencil calls. This implementation was
    # based on the cPython 3.7 version of the AST and proved too much effort to
    # continuously port to newer python versions. Ahead of dropping Python 3.7
    # support, all the kernel invocations were translated via the ``astor``
    # package ``astor.to_source()`` function to pure python source and this
    # source was hardcoded into the tests themselves. In the following tests,
    # regions demarked with dashed lines (----) and with the header
    # "Autogenerated kernel" correspond to these translations.

    def __init__(self, *args, **kwargs):
        super(TestManyStencils, self).__init__(*args, **kwargs)

    def check_against_expected(self, pyfunc, expected, *args, **kwargs):
        """
        For a given kernel:

        The expected result is available from argument `expected`.

        The following results are then computed:
        * from a pure @stencil decoration of the kernel.
        * from the njit of a trivial wrapper function around the pure @stencil
          decorated function.
        * from the njit(parallel=True) of a trivial wrapper function around
           the pure @stencil decorated function.

        The results are then compared.
        """

        options = kwargs.get('options', dict())
        expected_exception = kwargs.get('expected_exception')

        # DEBUG print output arrays
        DEBUG_OUTPUT = False

        # collect fails
        should_fail = []
        should_not_fail = []

        # runner that handles fails
        @contextmanager
        def errorhandler(exty=None, usecase=None):
            try:
                yield
            except Exception as e:
                if exty is not None:
                    lexty = exty if hasattr(exty, '__iter__') else [exty, ]
                    found = False
                    for ex in lexty:
                        found |= isinstance(e, ex)
                    if not found:
                        raise
                else:
                    should_not_fail.append(
                        (usecase, "%s: %s" %
                         (type(e), str(e))))
            else:
                if exty is not None:
                    should_fail.append(usecase)

        if isinstance(expected_exception, dict):
            stencil_ex = expected_exception['stencil']
            njit_ex = expected_exception['njit']
            parfor_ex = expected_exception['parfor']
        else:
            stencil_ex = expected_exception
            njit_ex = expected_exception
            parfor_ex = expected_exception

        stencil_args = {'func_or_mode': pyfunc}
        stencil_args.update(options)

        stencilfunc_output = None
        with errorhandler(stencil_ex, "@stencil"):
            stencil_func_impl = stencil(**stencil_args)
            # stencil result
            stencilfunc_output = stencil_func_impl(*args)

        # wrapped stencil impl, could this be generated?
        if len(args) == 1:
            def wrap_stencil(arg0):
                return stencil_func_impl(arg0)
        elif len(args) == 2:
            def wrap_stencil(arg0, arg1):
                return stencil_func_impl(arg0, arg1)
        elif len(args) == 3:
            def wrap_stencil(arg0, arg1, arg2):
                return stencil_func_impl(arg0, arg1, arg2)
        else:
            raise ValueError(
                "Up to 3 arguments can be provided, found %s" %
                len(args))

        sig = tuple([numba.typeof(x) for x in args])

        njit_output = None
        with errorhandler(njit_ex, "njit"):
            wrapped_cfunc = self.compile_njit(wrap_stencil, sig)
            # njit result
            njit_output = wrapped_cfunc.entry_point(*args)

        parfor_output = None
        with errorhandler(parfor_ex, "parfors"):
            wrapped_cpfunc = self.compile_parallel(wrap_stencil, sig)
            # parfor result
            parfor_output = wrapped_cpfunc.entry_point(*args)

        if DEBUG_OUTPUT:
            print("\n@stencil_output:\n", stencilfunc_output)
            print("\nnjit_output:\n", njit_output)
            print("\nparfor_output:\n", parfor_output)

        try:
            if not stencil_ex:
                np.testing.assert_almost_equal(
                    stencilfunc_output, expected, decimal=1)
                self.assertEqual(expected.dtype, stencilfunc_output.dtype)
        except Exception as e:
            should_not_fail.append(
                ('@stencil', "%s: %s" %
                    (type(e), str(e))))
            print("@stencil failed: %s" % str(e))

        try:
            if not njit_ex:
                np.testing.assert_almost_equal(
                    njit_output, expected, decimal=1)
                self.assertEqual(expected.dtype, njit_output.dtype)
        except Exception as e:
            should_not_fail.append(('njit', "%s: %s" % (type(e), str(e))))
            print("@njit failed: %s" % str(e))

        try:
            if not parfor_ex:
                np.testing.assert_almost_equal(
                    parfor_output, expected, decimal=1)
                self.assertEqual(expected.dtype, parfor_output.dtype)
                try:
                    self.assertIn(
                        '@do_scheduling',
                        wrapped_cpfunc.library.get_llvm_str())
                except AssertionError:
                    msg = 'Could not find `@do_scheduling` in LLVM IR'
                    raise AssertionError(msg)
        except Exception as e:
            should_not_fail.append(
                ('parfors', "%s: %s" %
                    (type(e), str(e))))
            print("@njit(parallel=True) failed: %s" % str(e))

        if DEBUG_OUTPUT:
            print("\n\n")

        if should_fail:
            msg = ["%s" % x for x in should_fail]
            raise RuntimeError(("The following implementations should have "
                                "raised an exception but did not:\n%s") % msg)

        if should_not_fail:
            impls = ["%s" % x[0] for x in should_not_fail]
            errs = ''.join(["%s: Message: %s\n\n" %
                            x for x in should_not_fail])
            str1 = ("The following implementations should not have raised an "
                    "exception but did:\n%s\n" % impls)
            str2 = "Errors were:\n\n%s" % errs
            raise RuntimeError(str1 + str2)

    def check_exceptions(self, pyfunc, *args, **kwargs):
        """
        For a given kernel:

        The expected result is computed from a pyStencil version of the
        stencil.

        The following results are then computed:
        * from a pure @stencil decoration of the kernel.
        * from the njit of a trivial wrapper function around the pure @stencil
          decorated function.
        * from the njit(parallel=True) of a trivial wrapper function around
           the pure @stencil decorated function.

        The results are then compared.
        """
        options = kwargs.get('options', dict())
        expected_exception = kwargs.get('expected_exception')

        # collect fails
        should_fail = []
        should_not_fail = []

        # runner that handles fails
        @contextmanager
        def errorhandler(exty=None, usecase=None):
            try:
                yield
            except Exception as e:
                if exty is not None:
                    lexty = exty if hasattr(exty, '__iter__') else [exty, ]
                    found = False
                    for ex in lexty:
                        found |= isinstance(e, ex)
                    if not found:
                        raise
                else:
                    should_not_fail.append(
                        (usecase, "%s: %s" %
                         (type(e), str(e))))
            else:
                if exty is not None:
                    should_fail.append(usecase)

        if isinstance(expected_exception, dict):
            stencil_ex = expected_exception['stencil']
            njit_ex = expected_exception['njit']
            parfor_ex = expected_exception['parfor']
        else:
            stencil_ex = expected_exception
            njit_ex = expected_exception
            parfor_ex = expected_exception

        stencil_args = {'func_or_mode': pyfunc}
        stencil_args.update(options)

        with errorhandler(stencil_ex, "@stencil"):
            stencil_func_impl = stencil(**stencil_args)
            # stencil result
            stencil_func_impl(*args)

        # wrapped stencil impl, could this be generated?
        if len(args) == 1:
            def wrap_stencil(arg0):
                return stencil_func_impl(arg0)
        elif len(args) == 2:
            def wrap_stencil(arg0, arg1):
                return stencil_func_impl(arg0, arg1)
        elif len(args) == 3:
            def wrap_stencil(arg0, arg1, arg2):
                return stencil_func_impl(arg0, arg1, arg2)
        else:
            raise ValueError(
                "Up to 3 arguments can be provided, found %s" %
                len(args))

        sig = tuple([numba.typeof(x) for x in args])

        with errorhandler(njit_ex, "njit"):
            wrapped_cfunc = self.compile_njit(wrap_stencil, sig)
            # njit result
            wrapped_cfunc.entry_point(*args)

        with errorhandler(parfor_ex, "parfors"):
            wrapped_cpfunc = self.compile_parallel(wrap_stencil, sig)
            # parfor result
            wrapped_cpfunc.entry_point(*args)

        if should_fail:
            msg = ["%s" % x for x in should_fail]
            raise RuntimeError(("The following implementations should have "
                                "raised an exception but did not:\n%s") % msg)

        if should_not_fail:
            impls = ["%s" % x[0] for x in should_not_fail]
            errs = ''.join(["%s: Message: %s\n\n" %
                            x for x in should_not_fail])
            str1 = ("The following implementations should not have raised an "
                    "exception but did:\n%s\n" % impls)
            str2 = "Errors were:\n\n%s" % errs
            raise RuntimeError(str1 + str2)

    def exception_dict(self, **kwargs):
        d = dict()
        d['pyStencil'] = None
        d['stencil'] = None
        d['njit'] = None
        d['parfor'] = None
        for k, v in kwargs.items():
            d[k] = v
        return d

    def check_stencil_arrays(self, *args, **kwargs):
        neighborhood = kwargs.get('neighborhood')
        init_shape = args[0].shape
        if neighborhood is not None:
            if len(init_shape) != len(neighborhood):
                raise ValueError('Invalid neighborhood supplied')
        for x in args[1:]:
            if hasattr(x, 'shape'):
                if init_shape != x.shape:
                    raise ValueError('Input stencil arrays do not commute')

    def test_basic00(self):
        """rel index"""
        def kernel(a):
            return a[0, 0]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic01(self):
        """rel index add const"""
        def kernel(a):
            return a[0, 1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic02(self):
        """rel index add const"""
        def kernel(a):
            return a[0, -1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + -1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic03(self):
        """rel index add const"""
        def kernel(a):
            return a[1, 0]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + 1, __b + 0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic04(self):
        """rel index add const"""
        def kernel(a):
            return a[-1, 0]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(1, a.shape[0]):
                    __b0[__a, __b] = a[__a + -1, __b + 0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic05(self):
        """rel index add const"""
        def kernel(a):
            return a[-1, 1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(1, a.shape[0]):
                    __b0[__a, __b] = a[__a + -1, __b + 1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic06(self):
        """rel index add const"""
        def kernel(a):
            return a[1, -1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1]):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + 1, __b + -1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic07(self):
        """rel index add const"""
        def kernel(a):
            return a[1, 1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + 1, __b + 1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic08(self):
        """rel index add const"""
        def kernel(a):
            return a[-1, -1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1]):
                for __a in range(1, a.shape[0]):
                    __b0[__a, __b] = a[__a + -1, __b + -1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic09(self):
        """rel index add const"""
        def kernel(a):
            return a[-2, 2]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 2):
                for __a in range(2, a.shape[0]):
                    __b0[__a, __b] = a[__a + -2, __b + 2]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic10(self):
        """rel index add const"""
        def kernel(a):
            return a[0, 0] + a[1, 0]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + 0, __b + 0] + a[__a + 1, __b + 0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic11(self):
        """rel index add const"""
        def kernel(a):
            return a[-1, 0] + a[1, 0]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(1, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + -1, __b + 0] + a[__a + 1, __b + 0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic12(self):
        """rel index add const"""
        def kernel(a):
            return a[-1, 1] + a[1, -1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(1, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + -1, __b + 1] + a[__a + 1, __b + -1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic13(self):
        """rel index add const"""
        def kernel(a):
            return a[-1, -1] + a[1, 1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(1, a.shape[0] - 1):
                    __b0[__a, __b] = a[__a + -1, __b + -1] + a[__a + 1, __b + 1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic14(self):
        """rel index add domain change const"""
        def kernel(a):
            return a[0, 0] + 1j
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 0] + 1.0j
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic14b(self):
        """rel index add domain change const"""
        def kernel(a):
            t = 1.j
            return a[0, 0] + t
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    t = 1.0j
                    __b0[__a, __b] = a[__a + 0, __b + 0] + t
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic15(self):
        """two rel index, add const"""
        def kernel(a):
            return a[0, 0] + a[1, 0] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      a[__a + 1, __b + 0] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic17(self):
        """two rel index boundary test, add const"""
        def kernel(a):
            return a[0, 0] + a[2, 0] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0] - 2):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      a[__a + 2, __b + 0] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic18(self):
        """two rel index boundary test, add const"""
        def kernel(a):
            return a[0, 0] + a[-2, 0] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(2, a.shape[0]):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      a[__a + -2, __b + 0] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic19(self):
        """two rel index boundary test, add const"""
        def kernel(a):
            return a[0, 0] + a[0, 3] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 3):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      a[__a + 0, __b + 3] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic20(self):
        """two rel index boundary test, add const"""
        def kernel(a):
            return a[0, 0] + a[0, -3] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(3, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      a[__a + 0, __b + -3] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic21(self):
        """same rel, add const"""
        def kernel(a):
            return a[0, 0] + a[0, 0] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      a[__a + 0, __b + 0] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic22(self):
        """rel idx const expr folding, add const"""
        def kernel(a):
            return a[1 + 0, 0] + a[0, 0] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 1, __b + 0] +
                                      a[__a + 0, __b + 0] + 1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic23(self):
        """rel idx, work in body"""
        def kernel(a):
            x = np.sin(10 + a[2, 1])
            return a[1 + 0, 0] + a[0, 0] + x
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0] - 2):
                    x = np.sin(10 + a[__a + 2, __b + 1])
                    __b0[__a, __b] = (a[__a + 1, __b + 0] +
                                      a[__a + 0, __b + 0] + x)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic23a(self):
        """rel idx, dead code should not impact rel idx"""
        def kernel(a):
            x = np.sin(10 + a[2, 1]) # noqa: F841 # dead code expected
            return a[1 + 0, 0] + a[0, 0]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0] - 2):
                    x = np.sin(10 + a[__a + 2, __b + 1]) # noqa: F841
                    __b0[__a, __b] = a[__a + 1, __b + 0] + a[__a + 0, __b + 0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic24(self):
        """1d idx on 2d arr"""
        a = np.arange(12).reshape(3, 4)

        def kernel(a):
            return a[0] + 1.

        self.check_exceptions(kernel, a, expected_exception=[TypingError,])

    def test_basic25(self):
        """no idx on 2d arr"""
        a = np.arange(12).reshape(3, 4)

        def kernel(a):
            return 1.
        self.check_exceptions(kernel, a, expected_exception=[ValueError,
                                                             NumbaValueError,])

    def test_basic26(self):
        """3d arr"""

        def kernel(a):
            return a[0, 0, 0] - a[0, 1, 0] + 1.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __c in range(0, a.shape[2]):
                for __b in range(0, a.shape[1] - 1):
                    for __a in range(0, a.shape[0]):
                        __b0[__a, __b, __c] = (a[__a + 0, __b + 0, __c + 0] -
                                               a[__a + 0, __b + 1, __c + 0] +
                                               1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(64).reshape(4, 8, 2)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic27(self):
        """4d arr"""
        def kernel(a):
            return a[0, 0, 0, 0] - a[0, 1, 0, -1] + 1.

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __d in range(1, a.shape[3]):
                for __c in range(0, a.shape[2]):
                    for __b in range(0, a.shape[1] - 1):
                        for __a in range(0, a.shape[0]):
                            __b0[__a, __b, __c, __d] = (a[__a + 0, __b + 0,
                                                          __c + 0, __d + 0] -
                                                        a[__a + 0, __b + 1,
                                                          __c + 0, __d + -1] +
                                                        1.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(128).reshape(4, 8, 2, 2)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic28(self):
        """type widen """
        def kernel(a):
            return a[0, 0] + np.float64(10.)

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 0] + np.float64(10.0)
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4).astype(np.float32)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic29(self):
        """const index from func """
        a = np.arange(12.).reshape(3, 4)

        def kernel(a):
            return a[0, int(np.cos(0))]
        self.check_exceptions(kernel, a, expected_exception=[ValueError,
                                                             NumbaValueError,
                                                             LoweringError])

    def test_basic30(self):
        """signed zeros"""
        def kernel(a):
            return a[-0, -0]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + -0, __b + -0]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12).reshape(3, 4).astype(np.float32)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic31(self):
        """does a const propagate? 2D"""
        def kernel(a):
            t = 1
            return a[t, 0]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0] - 1):
                    t = 1
                    __b0[__a, __b] = a[__a + t, __b + 0]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12).reshape(3, 4).astype(np.float32)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    @unittest.skip("constant folding not implemented")
    def test_basic31b(self):
        """does a const propagate?"""
        a = np.arange(12.).reshape(3, 4) # noqa: F841

        def kernel(a):
            s = 1
            t = 1 - s
            return a[t, 0]

        #TODO: add check should this be implemented

    def test_basic31c(self):
        """does a const propagate? 1D"""
        def kernel(a):
            t = 1
            return a[t]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __a in range(0, a.shape[0] - 1):
                t = 1
                __b0[__a,] = a[__a + t]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic32(self):
        """typed int index"""
        a = np.arange(12.).reshape(3, 4)

        def kernel(a):
            return a[np.int8(1), 0]
        self.check_exceptions(kernel, a, expected_exception=[ValueError,
                                                             NumbaValueError,
                                                             LoweringError])

    def test_basic33(self):
        """add 0d array"""
        def kernel(a):
            return a[0, 0] + np.array(1)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 0] + np.array(1)
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic34(self):
        """More complex rel index with dependency on addition rel index"""
        def kernel(a):
            g = 4. + a[0, 1]
            return g + (a[0, 1] + a[1, 0] + a[0, -1] + np.sin(a[-2, 0]))
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(2, a.shape[0] - 1):
                    g = 4.0 + a[__a + 0, __b + 1]
                    __b0[__a, __b] = g + (a[__a + 0, __b + 1] +
                                          a[__a + 1, __b + 0] +
                                          a[__a + 0, __b + -1] +
                                          np.sin(a[__a + -2, __b + 0]))
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(144).reshape(12, 12)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic35(self):
        """simple cval where cval is int but castable to dtype of float"""
        def kernel(a):
            return a[0, 1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 5, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 1]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a, options={'cval': 5})

    def test_basic36(self):
        """more complex with cval"""
        def kernel(a):
            return a[0, 1] + a[0, -1] + a[1, -1] + a[1, -1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 5.0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 1] +
                                      a[__a + 0, __b + -1] +
                                      a[__a + 1, __b + -1] +
                                      a[__a + 1, __b + -1])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a, options={'cval': 5})

    def test_basic37(self):
        """cval is expr"""
        def kernel(a):
            return a[0, 1] + a[0, -1] + a[1, -1] + a[1, -1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 68.0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 1] +
                                      a[__a + 0, __b + -1] +
                                      a[__a + 1, __b + -1] +
                                      a[__a + 1, __b + -1])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a,
                                    options={'cval': 5 + 63.})

    def test_basic38(self):
        """cval is complex"""
        def kernel(a):
            return a[0, 1] + a[0, -1] + a[1, -1] + a[1, -1]
        a = np.arange(12.).reshape(3, 4)
        ex = self.exception_dict(
            stencil=NumbaValueError,
            parfor=NumbaValueError,
            njit=NumbaValueError)
        self.check_exceptions(kernel, a, options={'cval': 1.j},
                              expected_exception=ex)

    def test_basic39(self):
        """cval is func expr"""
        def kernel(a):
            return a[0, 1] + a[0, -1] + a[1, -1] + a[1, -1]

        cval = np.sin(3.) + np.cos(2)

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, cval, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 1] +
                                      a[__a + 0, __b + -1] +
                                      a[__a + 1, __b + -1] +
                                      a[__a + 1, __b + -1])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a,
                                    options={'cval': cval})

    def test_basic40(self):
        """2 args!"""
        def kernel(a, b):
            return a[0, 1] + b[0, -2]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(2, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[__a + 0, __b + -2]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b)

    def test_basic41(self):
        """2 args! rel arrays wildly not same size!"""
        def kernel(a, b):
            return a[0, 1] + b[0, -2]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(1.).reshape(1, 1)
        self.check_exceptions(kernel, a, b, expected_exception=[ValueError,
                                                                AssertionError])

    def test_basic42(self):
        """2 args! rel arrays very close in size"""
        def kernel(a, b):
            return a[0, 1] + b[0, -2]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(9.).reshape(3, 3)
        self.check_exceptions(kernel, a, b, expected_exception=[ValueError,
                                                                AssertionError])

    def test_basic43(self):
        """2 args more complexity"""
        def kernel(a, b):
            return a[0, 1] + a[1, 2] + b[-2, 0] + b[0, -1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 2):
                for __a in range(2, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 1] +
                                      a[__a + 1, __b + 2] +
                                      b[__a + -2, __b + 0] +
                                      b[__a + 0, __b + -1])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(30.).reshape(5, 6)
        b = np.arange(30.).reshape(5, 6)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b)

    def test_basic44(self):
        """2 args, has assignment before use"""
        def kernel(a, b):
            a[0, 1] = 12
            return a[0, 1]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        self.check_exceptions(kernel, a, b, expected_exception=[NumbaValueError,
                                                                LoweringError])

    def test_basic45(self):
        """2 args, has assignment and then cross dependency"""
        def kernel(a, b):
            a[0, 1] = 12
            return a[0, 1] + a[1, 0]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        self.check_exceptions(kernel, a, b, expected_exception=[NumbaValueError,
                                                                LoweringError])

    def test_basic46(self):
        """2 args, has cross relidx assignment"""
        def kernel(a, b):
            a[0, 1] = b[1, 2]
            return a[0, 1] + a[1, 0]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        self.check_exceptions(kernel, a, b, expected_exception=[NumbaValueError,
                                                                LoweringError])

    def test_basic47(self):
        """3 args"""
        def kernel(a, b, c):
            return a[0, 1] + b[1, 0] + c[-1, 0]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, c, neighborhood):
            self.check_stencil_arrays(a, b, c, neighborhood=neighborhood)
            __retdtype = kernel(a, b, c)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(1, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 1] +
                                      b[__a + 1, __b + 0] +
                                      c[__a + -1, __b + 0])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        c = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, c, None)
        self.check_against_expected(kernel, expected, a, b, c)

    # matches pyStencil, but all ought to fail
    # probably hard to detect?
    def test_basic48(self):
        """2 args, has assignment before use via memory alias"""
        def kernel(a):
            c = a.T
            c[:, :] = 10
            return a[0, 1]

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a,neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    c = a.T
                    c[:, :] = 10
                    __b0[__a, __b] = a[__a + 0, __b + 1]
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic49(self):
        """2 args, standard_indexing on second"""
        def kernel(a, b):
            return a[0, 1] + b[0, 3]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[0, 3]
            return __b0

        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b'})

    @unittest.skip("dynamic range checking not implemented")
    def test_basic50(self):
        """2 args, standard_indexing OOB"""
        def kernel(a, b):
            return a[0, 1] + b[0, 15]
        #TODO: add check should this be implemented

    def test_basic51(self):
        """2 args, standard_indexing, no relidx"""
        def kernel(a, b):
            return a[0, 1] + b[0, 2]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        self.check_exceptions(kernel, a, b,
                              options={'standard_indexing': ['a', 'b']},
                              expected_exception=[ValueError, NumbaValueError])

    def test_basic52(self):
        """3 args, standard_indexing on middle arg """
        def kernel(a, b, c):
            return a[0, 1] + b[0, 1] + c[1, 2]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, c, neighborhood):
            self.check_stencil_arrays(a, c, neighborhood=neighborhood)
            __retdtype = kernel(a, b, c)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 2):
                for __a in range(0, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + 0, __b + 1] + b[0, 1] +
                                      c[__a + 1, __b + 2])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(4.).reshape(2, 2)
        c = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, c, None)
        self.check_against_expected(kernel, expected, a, b, c,
                                    options={'standard_indexing': 'b'})

    def test_basic53(self):
        """2 args, standard_indexing on variable that does not exist"""
        def kernel(a, b):
            return a[0, 1] + b[0, 2]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        ex = self.exception_dict(
            stencil=Exception,
            parfor=NumbaValueError,
            njit=Exception)
        self.check_exceptions(kernel, a, b, options={'standard_indexing': 'c'},
                              expected_exception=ex)

    def test_basic54(self):
        """2 args, standard_indexing, index from var"""
        def kernel(a, b):
            t = 2
            return a[0, 1] + b[0, t]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    t = 2
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[0, t]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b'})

    def test_basic55(self):
        """2 args, standard_indexing, index from more complex var"""
        def kernel(a, b):
            s = 1
            t = 2 - s
            return a[0, 1] + b[0, t]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    s = 1
                    t = 2 - s
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[0, t]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b'})

    def test_basic56(self):
        """2 args, standard_indexing, added complexity """
        def kernel(a, b):
            s = 1
            acc = 0
            for k in b[0, :]:
                acc += k
            t = 2 - s - 1
            return a[0, 1] + b[0, t] + acc
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    s = 1
                    acc = 0
                    for k in b[(0), :]:
                        acc += k
                    t = 2 - s - 1
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[0, t] + acc
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b'})

    def test_basic57(self):
        """2 args, standard_indexing, split index operation """
        def kernel(a, b):
            c = b[0]
            return a[0, 1] + c[1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    c = b[0]
                    __b0[__a, __b] = a[__a + 0, __b + 1] + c[1]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b'})

    def test_basic58(self):
        """2 args, standard_indexing, split index with broadcast mutation """
        def kernel(a, b):
            c = b[0] + 1
            return a[0, 1] + c[1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    c = b[0] + 1
                    __b0[__a, __b] = a[__a + 0, __b + 1] + c[1]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b'})

    def test_basic59(self):
        """3 args, mix of array, relative and standard indexing and const"""
        def kernel(a, b, c):
            return a[0, 1] + b[1, 1] + c
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, c, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b, c)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[1, 1] + c
            return __b0
        # ----------------------------------------------------------------------

        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        c = 10
        expected = __kernel(a, b, c, None)
        self.check_against_expected(kernel, expected, a, b, c,
                                    options={'standard_indexing': ['b', 'c']})

    def test_basic60(self):
        """3 args, mix of array, relative and standard indexing,
        tuple pass through"""
        def kernel(a, b, c):
            return a[0, 1] + b[1, 1] + c[0]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        c = (10,)
        # parfors does not support tuple args for stencil kernels
        ex = self.exception_dict(parfor=NumbaValueError)
        self.check_exceptions(kernel, a, b, c,
                              options={'standard_indexing': ['b', 'c']},
                              expected_exception=ex)

    def test_basic61(self):
        """2 args, standard_indexing on first"""
        def kernel(a, b):
            return a[0, 1] + b[1, 1]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        self.check_exceptions(kernel, a, b,
                              options={'standard_indexing': 'a'},
                              expected_exception=Exception)

    def test_basic62(self):
        """2 args, standard_indexing and cval"""
        def kernel(a, b):
            return a[0, 1] + b[1, 1]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 10.0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = a[__a + 0, __b + 1] + b[1, 1]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12.).reshape(3, 4)
        expected = __kernel(a, b, None)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'standard_indexing': 'b',
                                             'cval': 10.})

    def test_basic63(self):
        """2 args, standard_indexing applied to relative, should fail,
        non-const idx"""
        def kernel(a, b):
            return a[0, b[0, 1]]
        a = np.arange(12.).reshape(3, 4)
        b = np.arange(12).reshape(3, 4)
        ex = self.exception_dict(
            stencil=NumbaValueError,
            parfor=NumbaValueError,
            njit=NumbaValueError)
        self.check_exceptions(kernel, a, b, options={'standard_indexing': 'b'},
                              expected_exception=ex)

    # stencil, njit, parfors all fail. Does this make sense?
    def test_basic64(self):
        """1 arg that uses standard_indexing"""
        def kernel(a):
            return a[0, 0]
        a = np.arange(12.).reshape(3, 4)
        self.check_exceptions(kernel, a, options={'standard_indexing': 'a'},
                              expected_exception=[ValueError, NumbaValueError])

    def test_basic65(self):
        """basic induced neighborhood test"""
        def kernel(a):
            cumul = 0
            for i in range(-29, 1):
                cumul += a[i]
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(29, a.shape[0]):
                cumul = 0
                for i in range(-29, 1):
                    cumul += a[__an + i]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-29, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    # Should this work? a[0] is out of neighborhood?
    def test_basic66(self):
        """basic const neighborhood test"""
        def kernel(a):
            cumul = 0
            for i in range(-29, 1):
                cumul += a[0]
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(29, a.shape[0]):
                cumul = 0
                for i in range(-29, 1):
                    cumul += a[__an + 0]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-29, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic67(self):
        """basic 2d induced neighborhood test"""
        def kernel(a):
            cumul = 0
            for i in range(-5, 1):
                for j in range(-10, 1):
                    cumul += a[i, j]
            return cumul / (10 * 5)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(10, a.shape[1]):
                for __an in range(5, a.shape[0]):
                    cumul = 0
                    for i in range(-5, 1):
                        for j in range(-10, 1):
                            cumul += a[__an + i, __bn + j]
                    __b0[__an, __bn] = cumul / 50
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        nh = ((-5, 0), (-10, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic67b(self):
        """basic 2d induced 1D neighborhood"""
        def kernel(a):
            cumul = 0
            for j in range(-10, 1):
                cumul += a[0, j]
            return cumul / (10 * 5)
        a = np.arange(10. * 20.).reshape(10, 20)
        self.check_exceptions(kernel, a, options={'neighborhood': ((-10, 0),)},
                              expected_exception=[TypingError, ValueError])

    # Should this work or is it UB? a[i, 0] is out of neighborhood?
    def test_basic68(self):
        """basic 2d one induced, one cost neighborhood test"""
        def kernel(a):
            cumul = 0
            for i in range(-5, 1):
                for j in range(-10, 1):
                    cumul += a[i, 0]
            return cumul / (10 * 5)

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(10, a.shape[1]):
                for __an in range(5, a.shape[0]):
                    cumul = 0
                    for i in range(-5, 1):
                        for j in range(-10, 1):
                            cumul += a[__an + i, __bn + 0]
                    __b0[__an, __bn] = cumul / 50
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        nh = ((-5, 0), (-10, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    # Should this work or is it UB? a[0, 0] is out of neighborhood?
    def test_basic69(self):
        """basic 2d two cost neighborhood test"""
        def kernel(a):
            cumul = 0
            for i in range(-5, 1):
                for j in range(-10, 1):
                    cumul += a[0, 0]
            return cumul / (10 * 5)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(10, a.shape[1]):
                for __an in range(5, a.shape[0]):
                    cumul = 0
                    for i in range(-5, 1):
                        for j in range(-10, 1):
                            cumul += a[__an + 0, __bn + 0]
                    __b0[__an, __bn] = cumul / 50
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        nh = ((-5, 0), (-10, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic70(self):
        """neighborhood adding complexity"""
        def kernel(a):
            cumul = 0
            zz = 12.
            for i in range(-5, 1):
                t = zz + i
                for j in range(-10, 1):
                    cumul += a[i, j] + t
            return cumul / (10 * 5)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(10, a.shape[1]):
                for __an in range(5, a.shape[0]):
                    cumul = 0
                    zz = 12.0
                    for i in range(-5, 1):
                        t = zz + i
                        for j in range(-10, 1):
                            cumul += a[__an + i, __bn + j] + t
                    __b0[__an, __bn] = cumul / 50
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        nh = ((-5, 0), (-10, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic71(self):
        """neighborhood, type change"""
        def kernel(a):
            cumul = 0
            for i in range(-29, 1):
                k = 0.
                if i > -15:
                    k = 1j
                cumul += a[i] + k
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(29, a.shape[0]):
                cumul = 0
                for i in range(-29, 1):
                    k = 0.0
                    if i > -15:
                        k = 1.0j
                    cumul += a[__an + i] + k
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-29, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic72(self):
        """neighborhood, narrower range than specified"""
        def kernel(a):
            cumul = 0
            for i in range(-19, -3):
                cumul += a[i]
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(29, a.shape[0]):
                cumul = 0
                for i in range(-19, -3):
                    cumul += a[__an + i]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-29, 0),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic73(self):
        """neighborhood, +ve range"""
        def kernel(a):
            cumul = 0
            for i in range(5, 11):
                cumul += a[i]
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(0, a.shape[0] - 10):
                cumul = 0
                for i in range(5, 11):
                    cumul += a[__an + i]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((5, 10),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic73b(self):
        """neighborhood, -ve range"""
        def kernel(a):
            cumul = 0
            for i in range(-10, -4):
                cumul += a[i]
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(10, a.shape[0]):
                cumul = 0
                for i in range(-10, -4):
                    cumul += a[__an + i]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-10, -5),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic74(self):
        """neighborhood, -ve->+ve range span"""
        def kernel(a):
            cumul = 0
            for i in range(-5, 11):
                cumul += a[i]
            return cumul / 30
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(5, a.shape[0] - 10):
                cumul = 0
                for i in range(-5, 11):
                    cumul += a[__an + i]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-5, 10),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic75(self):
        """neighborhood, -ve->-ve range span"""
        def kernel(a):
            cumul = 0
            for i in range(-10, -1):
                cumul += a[i]
            return cumul / 30

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(10, a.shape[0]):
                cumul = 0
                for i in range(-10, -1):
                    cumul += a[__an + i]
                __b0[__an,] = cumul / 30
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(60.)
        nh = ((-10, -2),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic76(self):
        """neighborhood, mixed range span"""
        def kernel(a):
            cumul = 0
            zz = 12.
            for i in range(-3, 0):
                t = zz + i
                for j in range(-3, 4):
                    cumul += a[i, j] + t
            return cumul / (10 * 5)

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(3, a.shape[1] - 3):
                for __an in range(3, a.shape[0]):
                    cumul = 0
                    zz = 12.0
                    for i in range(-3, 0):
                        t = zz + i
                        for j in range(-3, 4):
                            cumul += a[__an + i, __bn + j] + t
                    __b0[__an, __bn] = cumul / 50
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        nh = ((-3, -1), (-3, 3),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    def test_basic77(self):
        """ neighborhood, two args """
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b[i, j]
            return cumul / (9.)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(3, a.shape[1]):
                for __an in range(3, a.shape[0]):
                    cumul = 0
                    for i in range(-3, 1):
                        for j in range(-3, 1):
                            cumul += (a[__an + i, __bn + j] +
                                      b[__an + i, __bn + j])
                    __b0[__an, __bn] = cumul / 9.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        b = np.arange(10. * 20.).reshape(10, 20)
        nh = ((-3, 0), (-3, 0),)
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh})

    def test_basic78(self):
        """ neighborhood, two args, -ve range, -ve range """
        def kernel(a, b):
            cumul = 0
            for i in range(-6, -2):
                for j in range(-7, -1):
                    cumul += a[i, j] + b[i, j]
            return cumul / (9.)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(7, a.shape[1]):
                for __an in range(6, a.shape[0]):
                    cumul = 0
                    for i in range(-6, -2):
                        for j in range(-7, -1):
                            cumul += (a[__an + i, __bn + j] +
                                      b[__an + i, __bn + j])
                    __b0[__an, __bn] = cumul / 9.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(15. * 20.).reshape(15, 20)
        b = np.arange(15. * 20.).reshape(15, 20)
        nh = ((-6, -3), (-7, -2),)
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh})

    def test_basic78b(self):
        """ neighborhood, two args, -ve range, +ve range """
        def kernel(a, b):
            cumul = 0
            for i in range(-6, -2):
                for j in range(2, 10):
                    cumul += a[i, j] + b[i, j]
            return cumul / (9.)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(0, a.shape[1] - 9):
                for __an in range(6, a.shape[0]):
                    cumul = 0
                    for i in range(-6, -2):
                        for j in range(2, 10):
                            cumul += (a[__an + i, __bn + j] +
                                      b[__an + i, __bn + j])
                    __b0[__an, __bn] = cumul / 9.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(15. * 20.).reshape(15, 20)
        b = np.arange(15. * 20.).reshape(15, 20)
        nh = ((-6, -3), (2, 9),)
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh})

    def test_basic79(self):
        """ neighborhood, two incompatible args """
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b[i, j]
            return cumul / (9.)
        a = np.arange(10. * 20.).reshape(10, 20)
        b = np.arange(10. * 20.).reshape(10, 10, 2)
        ex = self.exception_dict(
            stencil=TypingError,
            parfor=TypingError,
            njit=TypingError)
        self.check_exceptions(kernel, a, b, options={'neighborhood':
                                                     ((-3, 0), (-3, 0),)},
                              expected_exception=ex)

    def test_basic80(self):
        """ neighborhood, type change """
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b
            return cumul / (9.)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(3, a.shape[1]):
                for __an in range(3, a.shape[0]):
                    cumul = 0
                    for i in range(-3, 1):
                        for j in range(-3, 1):
                            cumul += a[__an + i, __bn + j] + b
                    __b0[__an, __bn] = cumul / 9.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        b = 12.j
        nh = ((-3, 0), (-3, 0))
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh})

    def test_basic81(self):
        """ neighborhood, dimensionally incompatible arrays """
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b[i]
            return cumul / (9.)
        a = np.arange(10. * 20.).reshape(10, 20)
        b = a[0].copy()
        ex = self.exception_dict(
            stencil=TypingError,
            parfor=AssertionError,
            njit=TypingError)
        self.check_exceptions(kernel, a, b,
                              options={'neighborhood': ((-3, 0), (-3, 0))},
                              expected_exception=ex)

    def test_basic82(self):
        """ neighborhood, with standard_indexing"""
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b[1, 3]
            return cumul / (9.)
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(3, a.shape[1]):
                for __an in range(3, a.shape[0]):
                    cumul = 0
                    for i in range(-3, 1):
                        for j in range(-3, 1):
                            cumul += a[__an + i, __bn + j] + b[1, 3]
                    __b0[__an, __bn] = cumul / 9.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        b = a.copy()
        nh = ((-3, 0), (-3, 0))
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh,
                                             'standard_indexing': 'b'})

    def test_basic83(self):
        """ neighborhood, with standard_indexing and cval"""
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b[1, 3]
            return cumul / (9.)
        a = np.arange(10. * 20.).reshape(10, 20)
        b = a.copy()

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 1.5, dtype=type(__retdtype))
            for __bn in range(3, a.shape[1]):
                for __an in range(3, a.shape[0]):
                    cumul = 0
                    for i in range(-3, 1):
                        for j in range(-3, 1):
                            cumul += a[__an + i, __bn + j] + b[1, 3]
                    __b0[__an, __bn] = cumul / 9.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        b = a.copy()
        nh = ((-3, 0), (-3, 0))
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh,
                                             'standard_indexing': 'b',
                                             'cval': 1.5,})

    def test_basic84(self):
        """ kernel calls njit """
        def kernel(a):
            return a[0, 0] + addone_njit(a[0, 1])
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      addone_njit.py_func(a[__a + 0, __b + 1]))
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic85(self):
        """ kernel calls njit(parallel=True)"""
        def kernel(a):
            return a[0, 0] + addone_pjit(a[0, 1])

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 1):
                for __a in range(0, a.shape[0]):
                    __b0[__a, __b] = (a[__a + 0, __b + 0] +
                                      addone_pjit.py_func(a[__a + 0, __b + 1]))
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    # njit/parfors fail correctly, but the error message isn't very informative
    def test_basic86(self):
        """ bad kwarg """
        def kernel(a):
            return a[0, 0]

        a = np.arange(10. * 20.).reshape(10, 20)
        self.check_exceptions(kernel, a, options={'bad': 10},
                              expected_exception=[ValueError, TypingError])

    def test_basic87(self):
        """ reserved arg name in use """
        def kernel(__sentinel__):
            return __sentinel__[0, 0]
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(__sentinel__, neighborhood):
            self.check_stencil_arrays(__sentinel__, neighborhood=neighborhood)
            __retdtype = kernel(__sentinel__)
            __b0 = np.full(__sentinel__.shape, 0, dtype=type(__retdtype))
            for __b in range(0, __sentinel__.shape[1]):
                for __a in range(0, __sentinel__.shape[0]):
                    __b0[__a, __b] = __sentinel__[__a + 0, __b + 0]
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic88(self):
        """ use of reserved word """
        def kernel(a, out):
            return out * a[0, 1]
        a = np.arange(12.).reshape(3, 4)
        ex = self.exception_dict(
            stencil=NumbaValueError,
            parfor=NumbaValueError,
            njit=NumbaValueError)
        self.check_exceptions(kernel, a, 1.0, options={}, expected_exception=ex)

    def test_basic89(self):
        """ basic multiple return"""
        def kernel(a):
            if a[0, 1] > 10:
                return 10.
            elif a[0, 3] < 8:
                return a[0, 0]
            else:
                return 7.
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1] - 3):
                for __a in range(0, a.shape[0]):
                    if a[__a + 0, __b + 1] > 10:
                        __b0[__a, __b] = 10.0
                    elif a[__a + 0, __b + 3] < 8:
                        __b0[__a, __b] = a[__a + 0, __b + 0]
                    else:
                        __b0[__a, __b] = 7.0
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic90(self):
        """ neighborhood, with standard_indexing and cval, multiple returns"""
        def kernel(a, b):
            cumul = 0
            for i in range(-3, 1):
                for j in range(-3, 1):
                    cumul += a[i, j] + b[1, 3]
            res = cumul / (9.)
            if res > 200.0:
                return res + 1.0
            else:
                return res
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, b, neighborhood):
            self.check_stencil_arrays(a, b, neighborhood=neighborhood)
            __retdtype = kernel(a, b)
            __b0 = np.full(a.shape, 1.5, dtype=type(__retdtype))
            for __bn in range(3, a.shape[1]):
                for __an in range(3, a.shape[0]):
                    cumul = 0
                    for i in range(-3, 1):
                        for j in range(-3, 1):
                            cumul += a[__an + i, __bn + j] + b[1, 3]
                    res = cumul / 9.0
                    if res > 200.0:
                        __b0[__an, __bn] = res + 1.0
                    else:
                        __b0[__an, __bn] = res
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        b = a.copy()
        nh = ((-3, 0), (-3, 0))
        expected = __kernel(a, b, nh)
        self.check_against_expected(kernel, expected, a, b,
                                    options={'neighborhood': nh,
                                             'standard_indexing': 'b',
                                             'cval': 1.5,})

    def test_basic91(self):
        """ Issue #3454, const(int) == const(int) evaluating incorrectly. """
        def kernel(a):
            b = 0
            if (2 == 0):
                b = 2
            return a[0, 0] + b
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(0, a.shape[1]):
                for __a in range(0, a.shape[0]):
                    b = 0
                    if 2 == 0:
                        b = 2
                    __b0[__a, __b] = a[__a + 0, __b + 0] + b
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(10. * 20.).reshape(10, 20)
        expected = __kernel(a, None)
        self.check_against_expected(kernel, expected, a)

    def test_basic92(self):
        """ Issue #3497, bool return type evaluating incorrectly. """
        def kernel(a):
            return (a[-1, -1] ^ a[-1, 0] ^ a[-1, 1] ^
                    a[0, -1] ^ a[0, 0] ^ a[0, 1] ^
                    a[1, -1] ^ a[1, 0] ^ a[1, 1])

        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(1, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + -1, __b + -1] ^
                                      a[__a + -1, __b + 0] ^
                                      a[__a + -1, __b + 1] ^
                                      a[__a + 0, __b + -1] ^
                                      a[__a + 0, __b + 0] ^
                                      a[__a + 0, __b + 1] ^
                                      a[__a + 1, __b + -1] ^
                                      a[__a + 1, __b + 0] ^
                                      a[__a + 1, __b + 1])
            return __b0
        # ----------------------------------------------------------------------
        A = np.array(np.arange(20) % 2).reshape(4, 5).astype(np.bool_)
        expected = __kernel(A, None)
        self.check_against_expected(kernel, expected, A)

    def test_basic93(self):
        """ Issue #3497, bool return type evaluating incorrectly. """
        def kernel(a):
            return (a[-1, -1] ^ a[-1, 0] ^ a[-1, 1] ^
                    a[0, -1] ^ a[0, 0] ^ a[0, 1] ^
                    a[1, -1] ^ a[1, 0] ^ a[1, 1])
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 1, dtype=type(__retdtype))
            for __b in range(1, a.shape[1] - 1):
                for __a in range(1, a.shape[0] - 1):
                    __b0[__a, __b] = (a[__a + -1, __b + -1] ^
                                      a[__a + -1, __b + 0] ^
                                      a[__a + -1, __b + 1] ^
                                      a[__a + 0, __b + -1] ^
                                      a[__a + 0, __b + 0] ^
                                      a[__a + 0, __b + 1] ^
                                      a[__a + 1, __b + -1] ^
                                      a[__a + 1, __b + 0] ^
                                      a[__a + 1, __b + 1])
            return __b0
        # ----------------------------------------------------------------------
        A = np.array(np.arange(20) % 2).reshape(4, 5).astype(np.bool_)
        expected = __kernel(A, None)
        self.check_against_expected(kernel, expected, A, options={'cval': True})

    def test_basic94(self):
        """ Issue #3528. Support for slices. """
        def kernel(a):
            return np.median(a[-1:2, -1:2])
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __bn in range(1, a.shape[1] - 1):
                for __an in range(1, a.shape[0] - 1):
                    __b0[__an, __bn] = np.median(a[__an + -1:__an + 2,
                                                   __bn + -1:__bn + 2])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(20, dtype=np.uint32).reshape(4, 5)
        nh = ((-1, 1), (-1, 1),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    @unittest.skip("not yet supported")
    def test_basic95(self):
        """ Slice, calculate neighborhood. """
        def kernel(a):
            return np.median(a[-1:2, -3:4])
        #TODO: add check should this be implemented

    def test_basic96(self):
        """ 1D slice. """
        def kernel(a):
            return np.median(a[-1:2])
        # ----------------------------------------------------------------------
        # Autogenerated kernel

        def __kernel(a, neighborhood):
            self.check_stencil_arrays(a, neighborhood=neighborhood)
            __retdtype = kernel(a)
            __b0 = np.full(a.shape, 0, dtype=type(__retdtype))
            for __an in range(1, a.shape[0] - 1):
                __b0[__an,] = np.median(a[__an + -1:__an + 2])
            return __b0
        # ----------------------------------------------------------------------
        a = np.arange(20, dtype=np.uint32)
        nh = ((-1, 1),)
        expected = __kernel(a, nh)
        self.check_against_expected(kernel, expected, a,
                                    options={'neighborhood': nh})

    @unittest.skip("not yet supported")
    def test_basic97(self):
        """ 2D slice and index. """
        def kernel(a):
            return np.median(a[-1:2, 3])
        #TODO: add check should this be implemented

    def test_basic98(self):
        """ Test issue #7286 where the cval is a np attr/string-based numerical
        constant"""
        for cval in (np.nan, np.inf, -np.inf, float('inf'), -float('inf')):
            def kernel(a):
                return a[0, 0]
            ## -----------------------------------------------------------------
            ## Autogenerated kernel

            def __kernel(a, neighborhood):
                self.check_stencil_arrays(a, neighborhood=neighborhood)
                __retdtype = kernel(a)
                __b0 = np.full(a.shape, cval, dtype=type(__retdtype))
                for __bn in range(1, a.shape[1] - 1):
                    for __an in range(1, a.shape[0] - 1):
                        __b0[__an, __bn] = a[__an + 0, __bn + 0]
                return __b0

            ## -----------------------------------------------------------------
            a = np.arange(6.).reshape((2, 3))
            nh = ((-1, 1), (-1, 1),)
            expected = __kernel(a, nh)
            self.check_against_expected(kernel, expected, a,
                                        options={'neighborhood': nh,
                                                 'cval':cval})


if __name__ == "__main__":
    unittest.main()
