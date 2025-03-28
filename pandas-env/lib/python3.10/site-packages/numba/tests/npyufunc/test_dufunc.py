import itertools
import pickle
import textwrap
import warnings

import numpy as np

from numba import njit, vectorize
from numba.tests.support import MemoryLeakMixin, TestCase
from numba.core.errors import (TypingError, NumbaNotImplementedError,
                               NumbaExperimentalFeatureWarning)
import unittest
from numba.np.ufunc import dufunc
from numba.np.numpy_support import from_dtype


def pyuadd(a0, a1):
    return a0 + a1


def pysub(a0, a1):
    return a0 - a1


def pymult(a0, a1):
    return a0 * a1


def pydiv(a0, a1):
    return a0 // a1


def pymin(a0, a1):
    return a0 if a0 < a1 else a1


class TestDUFunc(MemoryLeakMixin, unittest.TestCase):

    def nopython_dufunc(self, pyfunc):
        return dufunc.DUFunc(pyfunc, targetoptions=dict(nopython=True))

    def test_frozen(self):
        duadd = self.nopython_dufunc(pyuadd)
        self.assertFalse(duadd._frozen)
        duadd._frozen = True
        self.assertTrue(duadd._frozen)
        with self.assertRaises(ValueError):
            duadd._frozen = False
        with self.assertRaises(TypeError):
            duadd(np.linspace(0,1,10), np.linspace(1,2,10))

    def test_scalar(self):
        duadd = self.nopython_dufunc(pyuadd)
        self.assertEqual(pyuadd(1,2), duadd(1,2))

    def test_npm_call(self):
        duadd = self.nopython_dufunc(pyuadd)

        @njit
        def npmadd(a0, a1, o0):
            duadd(a0, a1, o0)
        X = np.linspace(0,1.9,20)
        X0 = X[:10]
        X1 = X[10:]
        out0 = np.zeros(10)
        npmadd(X0, X1, out0)
        np.testing.assert_array_equal(X0 + X1, out0)
        Y0 = X0.reshape((2,5))
        Y1 = X1.reshape((2,5))
        out1 = np.zeros((2,5))
        npmadd(Y0, Y1, out1)
        np.testing.assert_array_equal(Y0 + Y1, out1)
        Y2 = X1[:5]
        out2 = np.zeros((2,5))
        npmadd(Y0, Y2, out2)
        np.testing.assert_array_equal(Y0 + Y2, out2)

    def test_npm_call_implicit_output(self):
        duadd = self.nopython_dufunc(pyuadd)

        @njit
        def npmadd(a0, a1):
            return duadd(a0, a1)
        X = np.linspace(0,1.9,20)
        X0 = X[:10]
        X1 = X[10:]
        out0 = npmadd(X0, X1)
        np.testing.assert_array_equal(X0 + X1, out0)
        Y0 = X0.reshape((2,5))
        Y1 = X1.reshape((2,5))
        out1 = npmadd(Y0, Y1)
        np.testing.assert_array_equal(Y0 + Y1, out1)
        Y2 = X1[:5]
        out2 = npmadd(Y0, Y2)
        np.testing.assert_array_equal(Y0 + Y2, out2)
        out3 = npmadd(1.,2.)
        self.assertEqual(out3, 3.)

    def test_ufunc_props(self):
        duadd = self.nopython_dufunc(pyuadd)
        self.assertEqual(duadd.nin, 2)
        self.assertEqual(duadd.nout, 1)
        self.assertEqual(duadd.nargs, duadd.nin + duadd.nout)
        self.assertEqual(duadd.ntypes, 0)
        self.assertEqual(duadd.types, [])
        self.assertEqual(duadd.identity, None)
        duadd(1, 2)
        self.assertEqual(duadd.ntypes, 1)
        self.assertEqual(duadd.ntypes, len(duadd.types))
        self.assertIsNone(duadd.signature)

    def test_ufunc_props_jit(self):
        duadd = self.nopython_dufunc(pyuadd)
        duadd(1, 2)  # initialize types attribute

        attributes = {'nin': duadd.nin,
                      'nout': duadd.nout,
                      'nargs': duadd.nargs,
                      #'ntypes': duadd.ntypes,
                      #'types': duadd.types,
                      'identity': duadd.identity,
                      'signature': duadd.signature}

        def get_attr_fn(attr):
            fn = f'''
                def impl():
                    return duadd.{attr}
            '''
            l = {}
            exec(textwrap.dedent(fn), {'duadd': duadd}, l)
            return l['impl']

        for attr, val in attributes.items():
            cfunc = njit(get_attr_fn(attr))
            self.assertEqual(val, cfunc(),
                             f'Attribute differs from original: {attr}')

        # We don't expose [n]types attributes as they are dynamic attributes
        # and can change as the user calls the ufunc
        # cfunc = njit(get_attr_fn('ntypes'))
        # self.assertEqual(cfunc(), 1)
        # duadd(1.1, 2.2)
        # self.assertEqual(cfunc(), 2)


class TestDUFuncAt(TestCase):
    def _compare_output(self, fn, ufunc, a, *args):
        expected = a.copy()
        got = a.copy()
        ufunc.at(expected, *args)
        fn(got, *args)
        self.assertPreciseEqual(expected, got)

    def _generate_jit(self, ufunc):
        if ufunc.nin == 2:
            vec = vectorize()(lambda a, b: ufunc(a, b))
        else:
            vec = vectorize()(lambda a: ufunc(a))

        @njit
        def fn(*args):
            return vec.at(*args)
        return fn

    def test_numpy_ufunc_at_basic(self):
        # tests taken from: https://github.com/numpy/numpy/blob/27d8c43eb958b4ecee59b4d66908750759a9afc2/numpy/core/tests/test_ufunc.py#L1974-L2003  # noqa: E501
        # NumPy also test this function with a Rational array dtype. We skip
        # this test as Numba doesn't support Rational
        a = np.arange(10, dtype=int)

        add_at = self._generate_jit(np.add)
        negative_at = self._generate_jit(np.negative)

        negative_vec = vectorize()(lambda a: np.negative(a))

        @njit
        def negative_jit_2(a, indices, b):
            return negative_vec.at(a, indices, b)

        # basic testing
        self._compare_output(add_at, np.add, a, [2, 5, 2], 1)

        # missing second operand
        err_msg = 'second operand needed for ufunc'
        with self.assertRaisesRegex(TypingError, err_msg):
            add_at(a.copy(), [2, 5, 3], None)

        self._compare_output(negative_at, np.negative, a.copy(), [2, 5, 3])

        b = np.array([100, 100, 100])
        self._compare_output(add_at, np.add, a.copy(), [2, 5, 2], b)

        # extraneous second operand
        err_msg = 'second operand provided when ufunc is unary'
        with self.assertRaisesRegex(TypingError, err_msg):
            negative_jit_2(a.copy(), [2, 5, 3], [1, 2, 3])

        with self.assertRaises(TypingError):
            add_at(a.copy(), [2, 5, 3], [[1, 2], 1])

    def test_ufunc_at_inner_loop(self):
        typecodes = np.typecodes['Complex']
        ufuncs = (np.add, np.subtract, np.multiply)
        for typecode in typecodes:

            try:
                from_dtype(np.dtype(typecode))
            except NumbaNotImplementedError:
                continue

            for ufunc in ufuncs:
                a = np.ones(10, dtype=typecode)
                indx = np.concatenate([np.ones(6, dtype=np.intp),
                                       np.full(18, 4, dtype=np.intp)])
                value = a.dtype.type(1j)
                ufunc_at = self._generate_jit(ufunc)
                ufunc_at(a, indx, value)
                expected = np.ones_like(a)
                if ufunc is np.multiply:
                    expected[1] = expected[4] = -1
                else:
                    expected[1] += 6 * (value if ufunc is np.add else -value)
                    expected[4] += 18 * (value if ufunc is np.add else -value)

        self.assertPreciseEqual(a, expected)

    def test_ufunc_at_ellipsis(self):
        # Make sure the indexed loop check does not choke on iters
        # with subspaces
        arr = np.zeros(5, dtype=int)
        add_at = self._generate_jit(np.add)
        self._compare_output(add_at, np.add, arr, slice(None),
                             np.ones(5, dtype=int))

    def test_ufunc_at_negative(self):
        arr = np.ones(5, dtype=np.int32)
        indx = np.arange(5)
        at = self._generate_jit(np.negative)
        at(arr, indx)
        assert np.all(arr == [-1, -1, -1, -1, -1])

    def test_ufunc_at_large(self):
        # NumPy issue gh-23457
        indices = np.zeros(8195, dtype=np.int16)
        b = np.zeros(8195, dtype=float)
        b[0] = 10
        b[1] = 5
        b[8192:] = 100
        a = np.zeros(1, dtype=float)
        add_at = self._generate_jit(np.add)
        add_at(a, indices, b)
        assert a[0] == b.sum()

    def test_cast_index_fastpath(self):
        arr = np.zeros(10)
        values = np.ones(100000)
        add_at = self._generate_jit(np.add)
        # index must be cast, which may be buffered in chunks:
        index = np.zeros(len(values), dtype=np.uint8)
        add_at(arr, index, values)
        assert arr[0] == len(values)

    def test_ufunc_at_scalar_value_fastpath(self):
        values = (np.ones(1), np.ones(()), np.float64(1.), 1.)
        for value in values:
            arr = np.zeros(1000)
            # index must be cast, which may be buffered in chunks:
            index = np.repeat(np.arange(1000), 2)
            add_at = self._generate_jit(np.add)
            add_at(arr, index, value)
            np.testing.assert_array_equal(arr, np.full_like(arr, 2 * value))

    def test_ufunc_at_multiD(self):
        a = np.arange(9).reshape(3, 3)
        b = np.array([[100, 100, 100], [200, 200, 200], [300, 300, 300]])
        add_at = self._generate_jit(np.add)
        add_at(a, (slice(None), np.asarray([1, 2, 1])), b)
        self.assertPreciseEqual(a, np.array(
            [[0, 201, 102], [3, 404, 205], [6, 607, 308]]))

        a = np.arange(27).reshape(3, 3, 3)
        b = np.array([100, 200, 300])
        add_at(a, (slice(None), slice(None), np.asarray([1, 2, 1])), b)
        self.assertPreciseEqual(a, np.array(
                                [[[0, 401, 202],
                                  [3, 404, 205],
                                  [6, 407, 208]],

                                 [[9, 410, 211],
                                  [12, 413, 214],
                                  [15, 416, 217]],

                                 [[18, 419, 220],
                                  [21, 422, 223],
                                  [24, 425, 226]]]))

        a = np.arange(9).reshape(3, 3)
        b = np.array([[100, 100, 100], [200, 200, 200], [300, 300, 300]])
        add_at(a, (np.asarray([1, 2, 1]), slice(None)), b)
        self.assertPreciseEqual(a, np.asarray(
            [[0, 1, 2], [403, 404, 405], [206, 207, 208]]))

        a = np.arange(27).reshape(3, 3, 3)
        b = np.array([100, 200, 300])
        add_at(a, (slice(None), np.asarray([1, 2, 1]), slice(None)), b)
        self.assertPreciseEqual(a, np.asarray(
                                [[[0,  1,  2],
                                  [203, 404, 605],
                                  [106, 207, 308]],

                                 [[9,  10, 11],
                                  [212, 413, 614],
                                  [115, 216, 317]],

                                 [[18, 19, 20],
                                  [221, 422, 623],
                                  [124, 225, 326]]]))

        a = np.arange(9).reshape(3, 3)
        b = np.array([100, 200, 300])
        add_at(a, (0, np.asarray([1, 2, 1])), b)
        self.assertPreciseEqual(a, np.asarray(
            [[0, 401, 202], [3, 4, 5], [6, 7, 8]]))

        a = np.arange(27).reshape(3, 3, 3)
        b = np.array([100, 200, 300])
        add_at(a, (np.asarray([1, 2, 1]), 0, slice(None)), b)
        self.assertPreciseEqual(a, np.asarray(
                                [[[0,  1,  2],
                                  [3,  4,  5],
                                  [6,  7,  8]],

                                 [[209, 410, 611],
                                  [12,  13, 14],
                                  [15,  16, 17]],

                                 [[118, 219, 320],
                                  [21,  22, 23],
                                  [24,  25, 26]]]))

        a = np.arange(27).reshape(3, 3, 3)
        b = np.array([100, 200, 300])
        add_at = self._generate_jit(np.add)
        add_at(a, (slice(None), slice(None), slice(None)), b)
        self.assertPreciseEqual(a, np.asarray(
                                [[[100, 201, 302],
                                  [103, 204, 305],
                                  [106, 207, 308]],

                                 [[109, 210, 311],
                                  [112, 213, 314],
                                  [115, 216, 317]],

                                 [[118, 219, 320],
                                  [121, 222, 323],
                                  [124, 225, 326]]]))

    def test_ufunc_at_0D(self):
        a = np.array(0)
        add_at = self._generate_jit(np.add)
        add_at(a, (), 1)
        self.assertPreciseEqual(a, np.array(1))

        with self.assertRaises(TypingError):
            add_at(a, 0, 1)

        b = np.arange(3)
        add_at(b, 0, 1)
        self.assertPreciseEqual(b, np.array([1, 1, 2]))

        # NumPy checks for IndexError but we can't call a jit function with an
        # empty list as Numba raises "can't compute fingerprint of empty list"
        with self.assertRaises(ValueError):
            add_at(a, [], 1)

    def test_ufunc_at_dtypes(self):
        # Test mixed dtypes
        a = np.arange(10)
        power_at = self._generate_jit(np.power)
        power_at(a, [1, 2, 3, 2], 3.5)
        self.assertPreciseEqual(a, np.array([0, 1, 4414, 46, 4, 5, 6, 7, 8, 9]))

    def test_ufunc_at_boolean(self):
        # Test boolean indexing and boolean ufuncs
        a = np.arange(10)
        index = a % 2 == 0
        equal_at = self._generate_jit(np.equal)
        # boolean indexing not supported
        equal_at(a, index, [0, 2, 4, 6, 8])
        self.assertPreciseEqual(a, np.array([1, 1, 1, 3, 1, 5, 1, 7, 1, 9]))

    def test_ufunc_at_boolean2(self):
        # Test unary operator
        a = np.arange(10, dtype='u4')
        invert_at = self._generate_jit(np.invert)
        invert_at(a, [2, 5, 2])
        self.assertPreciseEqual(a, np.array([0, 1, 2, 3, 4, 5 ^ 0xffffffff, 6,
                                             7, 8, 9], dtype=np.uint32))

    def test_ufunc_at_advanced(self):
        # Test empty subspace
        orig = np.arange(4)
        a = orig[:, None][:, 0:0]
        add_at = self._generate_jit(np.add)
        add_at(a, [0, 1], 3)
        self.assertPreciseEqual(orig, np.arange(4))

    @unittest.expectedFailure
    def test_ufunc_at_advanced_2(self):
        # Test with swapped byte order
        index = np.array([1, 2, 1], np.dtype('i').newbyteorder())
        values = np.array([1, 2, 3, 4], np.dtype('f').newbyteorder())
        add_at = self._generate_jit(np.add)
        add_at(values, index, 3)
        self.assertPreciseEqual(values, [1, 8, 6, 4])

    def test_ufunc_at_advanced_3(self):
        # Test exception thrown
        values = np.array(['a', 1], dtype=object)
        add_at = self._generate_jit(np.add)
        with self.assertRaises(TypingError):
            add_at(values, [0, 1], 1)
        self.assertPreciseEqual(values, np.array(['a', 1], dtype=object))

    def test_ufunc_at_advanced_4(self):
        # Test multiple output ufuncs raise error, NumPy gh-5665
        modf_at = self._generate_jit(np.modf)
        # NumPy raises ValueError as modf returns multiple outputs
        with self.assertRaises(TypingError):
            modf_at(np.arange(10), [1])

    def test_ufunc_at_advanced_5(self):
        # Test maximum
        maximum_at = self._generate_jit(np.maximum)
        a = np.array([1, 2, 3])
        maximum_at(a, [0], 0)
        self.assertPreciseEqual(a, np.array([1, 2, 3]))

    def test_ufunc_at_negative_indexes(self):
        dtypes = np.typecodes['AllInteger'] + np.typecodes['Float']
        ufuncs = (np.add, np.subtract, np.divide, np.minimum, np.maximum)

        for dtype in dtypes:

            if dtype in ('e',):  # skip float16 as we don't have an impl. for it
                continue

            try:
                from_dtype(np.dtype(dtype))
            except NumbaNotImplementedError:
                continue

            for ufunc in ufuncs:
                a = np.arange(0, 10).astype(dtype)
                indxs = np.array([-1, 1, -1, 2]).astype(np.intp)
                vals = np.array([1, 5, 2, 10], dtype=a.dtype)

                expected = a.copy()
                for i, v in zip(indxs, vals):
                    expected[i] = ufunc(expected[i], v)

                ufunc_at = self._generate_jit(ufunc)
                ufunc_at(a, indxs, vals)
                np.testing.assert_array_equal(a, expected)
                assert np.all(indxs == [-1, 1, -1, 2])

    @unittest.expectedFailure
    def test_ufunc_at_not_none_signature(self):
        # Test ufuncs with non-trivial signature raise a TypeError
        a = np.ones((2, 2, 2))
        b = np.ones((1, 2, 2))
        # matmul is a gufunc, thus, this will fail atm
        matmul_at = self._generate_jit(np.matmul)
        err_msg = 'does not support ufunc with non-trivial signature'
        with self.assertRaisesRegex(TypingError, err_msg):
            matmul_at(a, [0], b)

        # a = np.array([[[1, 2], [3, 4]]])
        # assert_raises(TypeError, np.linalg._umath_linalg.det.at, a, [0])

    def test_ufunc_at_no_loop_for_op(self):
        # str dtype does not have a ufunc loop for np.add
        arr = np.ones(10, dtype=str)
        add_at = self._generate_jit(np.add)
        # NumPy raises `np.core._exceptions._UFuncNoLoopError`
        with self.assertRaises(TypingError):
            add_at(arr, [0, 1], [0, 1])

    def test_ufunc_at_output_casting(self):
        arr = np.array([-1])
        equal_at = self._generate_jit(np.equal)
        equal_at(arr, [0], [0])
        assert arr[0] == 0

    def test_ufunc_at_broadcast_failure(self):
        arr = np.arange(5)
        add_at = self._generate_jit(np.add)

        # NumPy raises ValueError('array is not broadcastable to correct shape')
        msg = 'operands could not be broadcast together with remapped shapes'
        with self.assertRaisesRegex(ValueError, msg):
            add_at(arr, [0, 1], [1, 2, 3])

    def test_ufunc_at_dynamic(self):
        arr = np.arange(5)

        @vectorize
        def inc(x):
            return x + 1

        self.assertEqual(len(inc.types), 0)

        # trying to call inc.at should trigger compilation
        inc.at(arr, [1, 3])

        self.assertGreater(len(inc.types), 0)

    def test_ufunc_at_experimental_warning(self):
        arr = np.arange(5)
        add_at = self._generate_jit(np.add)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', NumbaExperimentalFeatureWarning)

            add_at(arr, [0, 3], 10)

        self.assertGreater(len(w), 0)
        self.assertIn('ufunc.at feature is experimental', str(w[0].message))


class TestDUFuncReduceNumPyTests(TestCase):
    # Tests taken from
    # https://github.com/numpy/numpy/blob/51ee17b6bd4ccec60a5483ee8bff94ad0c0e8585/numpy/_core/tests/test_ufunc.py  # noqa: E501

    def _generate_jit(self, ufunc, identity=None):
        if ufunc.nin == 2:
            vec = vectorize(identity=identity)(lambda a, b: ufunc(a, b))
        else:
            vec = vectorize(identity=identity)(lambda a: ufunc(a))

        @njit
        def fn(array, axis=0, initial=None):
            return vec.reduce(array, axis=axis, initial=initial)
        return fn

    @unittest.expectedFailure
    def test_numpy_scalar_reduction(self):
        # scalar reduction is not supported
        power_reduce = self._generate_jit(np.power)
        expected = np.power.reduce(3)
        got = power_reduce(3)
        self.assertPreciseEqual(expected, got)

    def check_identityless_reduction(self, a):
        def compare_output(a, b):
            # We don't use self.assertPreciseEqual as the dtype differs
            # between the value from the reduction and the hardcoded output
            np.testing.assert_equal(a, b)
        # test taken from:
        # https://github.com/numpy/numpy/blob/51ee17b6bd4ccec60a5483ee8bff94ad0c0e8585/numpy/_core/tests/test_ufunc.py#L1591  # noqa: E501

        minimum_reduce = self._generate_jit(np.minimum, identity='reorderable')

        # np.minimum.reduce is an identityless reduction

        # Verify that it sees the zero at various positions
        a[...] = 1
        a[1, 0, 0] = 0
        compare_output(minimum_reduce(a, axis=None), 0)
        compare_output(minimum_reduce(a, axis=(0, 1)), [0, 1, 1, 1])
        compare_output(minimum_reduce(a, axis=(0, 2)), [0, 1, 1])
        compare_output(minimum_reduce(a, axis=(1, 2)), [1, 0])
        compare_output(minimum_reduce(a, axis=0),
                       [[0, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
        compare_output(minimum_reduce(a, axis=1),
                       [[1, 1, 1, 1], [0, 1, 1, 1]])
        compare_output(minimum_reduce(a, axis=2),
                       [[1, 1, 1], [0, 1, 1]])
        compare_output(minimum_reduce(a, axis=()), a)

        a[...] = 1
        a[0, 1, 0] = 0
        compare_output(minimum_reduce(a, axis=None), 0)
        compare_output(minimum_reduce(a, axis=(0, 1)), [0, 1, 1, 1])
        compare_output(minimum_reduce(a, axis=(0, 2)), [1, 0, 1])
        compare_output(minimum_reduce(a, axis=(1, 2)), [0, 1])
        compare_output(minimum_reduce(a, axis=0),
                       [[1, 1, 1, 1], [0, 1, 1, 1], [1, 1, 1, 1]])
        compare_output(minimum_reduce(a, axis=1),
                       [[0, 1, 1, 1], [1, 1, 1, 1]])
        compare_output(minimum_reduce(a, axis=2),
                       [[1, 0, 1], [1, 1, 1]])
        compare_output(minimum_reduce(a, axis=()), a)

        a[...] = 1
        a[0, 0, 1] = 0
        compare_output(minimum_reduce(a, axis=None), 0)
        compare_output(minimum_reduce(a, axis=(0, 1)), [1, 0, 1, 1])
        compare_output(minimum_reduce(a, axis=(0, 2)), [0, 1, 1])
        compare_output(minimum_reduce(a, axis=(1, 2)), [0, 1])
        compare_output(minimum_reduce(a, axis=0),
                       [[1, 0, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
        compare_output(minimum_reduce(a, axis=1),
                       [[1, 0, 1, 1], [1, 1, 1, 1]])
        compare_output(minimum_reduce(a, axis=2),
                       [[0, 1, 1], [1, 1, 1]])
        compare_output(minimum_reduce(a, axis=()), a)

    def test_numpy_identityless_reduction_corder(self):
        a = np.empty((2, 3, 4), order='C')
        self.check_identityless_reduction(a)

    def test_numpy_identityless_reduction_forder(self):
        a = np.empty((2, 3, 4), order='F')
        self.check_identityless_reduction(a)

    def test_numpy_identityless_reduction_otherorder(self):
        a = np.empty((2, 4, 3), order='C').swapaxes(1, 2)
        self.check_identityless_reduction(a)

    def test_numpy_identityless_reduction_noncontig(self):
        a = np.empty((3, 5, 4), order='C').swapaxes(1, 2)
        a = a[1:, 1:, 1:]
        self.check_identityless_reduction(a)

    def test_numpy_identityless_reduction_noncontig_unaligned(self):
        a = np.empty((3 * 4 * 5 * 8 + 1,), dtype='i1')
        a = a[1:].view(dtype='f8')
        a.shape = (3, 4, 5)
        a = a[1:, 1:, 1:]
        self.check_identityless_reduction(a)

    def test_numpy_initial_reduction(self):
        # np.minimum.reduce is an identityless reduction
        add_reduce = self._generate_jit(np.add)
        min_reduce = self._generate_jit(np.minimum)
        max_reduce = self._generate_jit(np.maximum)

        # For cases like np.maximum(np.abs(...), initial=0)
        # More generally, a supremum over non-negative numbers.
        self.assertPreciseEqual(max_reduce(np.asarray([]), initial=0), 0.0)

        # For cases like reduction of an empty array over the reals.
        self.assertPreciseEqual(min_reduce(np.asarray([]), initial=np.inf),
                                np.inf)
        self.assertPreciseEqual(max_reduce(np.asarray([]), initial=-np.inf),
                                -np.inf)

        # Random tests
        self.assertPreciseEqual(min_reduce(np.asarray([5]), initial=4), 4)
        self.assertPreciseEqual(max_reduce(np.asarray([4]), initial=5), 5)
        self.assertPreciseEqual(max_reduce(np.asarray([5]), initial=4), 5)
        self.assertPreciseEqual(min_reduce(np.asarray([4]), initial=5), 4)

        # Check initial=None raises ValueError for both types of ufunc
        # reductions
        msg = 'zero-size array to reduction operation'
        for func in (add_reduce, min_reduce):
            with self.assertRaisesRegex(ValueError, msg):
                func(np.asarray([]), initial=None)

    def test_numpy_empty_reduction_and_identity(self):
        arr = np.zeros((0, 5))
        true_divide_reduce = self._generate_jit(np.true_divide)

        # OK, since the reduction itself is *not* empty, the result is
        expected = np.true_divide.reduce(arr, axis=1)
        got = true_divide_reduce(arr, axis=1)
        self.assertPreciseEqual(expected, got)
        self.assertPreciseEqual(got.shape, (0,))

        # Not OK, the reduction itself is empty and we have no identity
        msg = 'zero-size array to reduction operation'
        with self.assertRaisesRegex(ValueError, msg):
            true_divide_reduce(arr, axis=0)

        # Test that an empty reduction fails also if the result is empty
        arr = np.zeros((0, 0, 5))
        with self.assertRaisesRegex(ValueError, msg):
            true_divide_reduce(arr, axis=1)

        # Division reduction makes sense with `initial=1` (empty or not):
        expected = np.true_divide.reduce(arr, axis=1, initial=1)
        got = true_divide_reduce(arr, axis=1, initial=1)
        self.assertPreciseEqual(expected, got)

    def test_identityless_reduction_nonreorderable(self):
        a = np.array([[8.0, 2.0, 2.0], [1.0, 0.5, 0.25]])

        divide_reduce = self._generate_jit(np.divide)
        res = divide_reduce(a, axis=0)
        self.assertPreciseEqual(res, np.asarray([8.0, 4.0, 8.0]))

        res = divide_reduce(a, axis=1)
        self.assertPreciseEqual(res, np.asarray([2.0, 8.0]))

        res = divide_reduce(a, axis=())
        self.assertPreciseEqual(res, a)

        # will not raise as per Numba issue #9283
        # assert_raises(ValueError, np.divide.reduce, a, axis=(0, 1))

    def test_reduce_zero_axis(self):
        # If we have a n x m array and do a reduction with axis=1, then we are
        # doing n reductions, and each reduction takes an m-element array. For
        # a reduction operation without an identity, then:
        #   n > 0, m > 0: fine
        #   n = 0, m > 0: fine, doing 0 reductions of m-element arrays
        #   n > 0, m = 0: can't reduce a 0-element array, ValueError
        #   n = 0, m = 0: can't reduce a 0-element array, ValueError (for
        #     consistency with the above case)
        # This test doesn't actually look at return values, it just checks to
        # make sure that error we get an error in exactly those cases where we
        # expect one, and assumes the calculations themselves are done
        # correctly.

        def ok(f, *args, **kwargs):
            f(*args, **kwargs)

        def err(f, *args, **kwargs):
            with self.assertRaises(ValueError):
                f(*args, **kwargs)

        def t(expect, func, n, m):
            expect(func, np.zeros((n, m)), axis=1)
            expect(func, np.zeros((m, n)), axis=0)
            expect(func, np.zeros((n // 2, n // 2, m)), axis=2)
            expect(func, np.zeros((n // 2, m, n // 2)), axis=1)
            expect(func, np.zeros((n, m // 2, m // 2)), axis=(1, 2))
            expect(func, np.zeros((m // 2, n, m // 2)), axis=(0, 2))
            expect(func, np.zeros((m // 3, m // 3, m // 3,
                                   n // 2, n // 2)), axis=(0, 1, 2))
            # Check what happens if the inner (resp. outer) dimensions are a
            # mix of zero and non-zero:
            expect(func, np.zeros((10, m, n)), axis=(0, 1))
            expect(func, np.zeros((10, n, m)), axis=(0, 2))
            expect(func, np.zeros((m, 10, n)), axis=0)
            expect(func, np.zeros((10, m, n)), axis=1)
            expect(func, np.zeros((10, n, m)), axis=2)

        # np.maximum is just an arbitrary ufunc with no reduction identity
        maximum_reduce = self._generate_jit(np.maximum, identity='reorderable')
        self.assertEqual(np.maximum.identity, None)
        t(ok, maximum_reduce, 30, 30)
        t(ok, maximum_reduce, 0, 30)
        t(err, maximum_reduce, 30, 0)
        t(err, maximum_reduce, 0, 0)
        err(maximum_reduce, [])
        maximum_reduce(np.zeros((0, 0)), axis=())

        # all of the combinations are fine for a reduction that has an
        # identity
        add_reduce = self._generate_jit(np.add, identity=0)
        t(ok, add_reduce, 30, 30)
        t(ok, add_reduce, 0, 30)
        t(ok, add_reduce, 30, 0)
        t(ok, add_reduce, 0, 0)
        add_reduce(np.array([], dtype=np.int64))
        add_reduce(np.zeros((0, 0)), axis=())


class TestDUFuncReduce(TestCase):
    def _check_reduce(self, ufunc, dtype=None, initial=None):

        @njit
        def foo(a, axis, dtype, initial):
            return ufunc.reduce(a,
                                axis=axis,
                                dtype=dtype,
                                initial=initial)

        inputs = [
            np.arange(5),
            np.arange(4).reshape(2, 2),
            np.arange(40).reshape(5, 4, 2),
        ]
        for array in inputs:
            for axis in range(array.ndim):
                expected = foo.py_func(array, axis, dtype, initial)
                got = foo(array, axis, dtype, initial)
                self.assertPreciseEqual(expected, got)

    def _check_reduce_axis(self, ufunc, dtype, initial=None):

        @njit
        def foo(a, axis):
            return ufunc.reduce(a, axis=axis, initial=initial)

        def _check(*args):
            try:
                expected = foo.py_func(array, axis)
            except ValueError as e:
                self.assertEqual(e.args[0], exc_msg)
                with self.assertRaisesRegex(TypingError, exc_msg):
                    got = foo(array, axis)
            else:
                got = foo(array, axis)
                self.assertPreciseEqual(expected, got)

        exc_msg = (f"reduction operation '{ufunc.__name__}' is not "
                   "reorderable, so at most one axis may be specified")
        inputs = [
            np.arange(40, dtype=dtype).reshape(5, 4, 2),
            np.arange(10, dtype=dtype),
        ]
        for array in inputs:
            for i in range(1, array.ndim + 1):
                for axis in itertools.combinations(range(array.ndim), r=i):
                    _check(array, axis)

            # corner cases: Reduce over axis=() and axis=None
            for axis in ((), None):
                _check(array, axis)

    def test_add_reduce(self):
        duadd = vectorize('int64(int64, int64)', identity=0)(pyuadd)
        self._check_reduce(duadd)
        self._check_reduce_axis(duadd, dtype=np.int64)

    def test_mul_reduce(self):
        dumul = vectorize('int64(int64, int64)', identity=1)(pymult)
        self._check_reduce(dumul)

    def test_non_associative_reduce(self):
        dusub = vectorize('int64(int64, int64)', identity=None)(pysub)
        dudiv = vectorize('int64(int64, int64)', identity=None)(pydiv)
        self._check_reduce(dusub)
        self._check_reduce_axis(dusub, dtype=np.int64)
        self._check_reduce(dudiv)
        self._check_reduce_axis(dudiv, dtype=np.int64)

    def test_reduce_dtype(self):
        duadd = vectorize('float64(float64, int64)', identity=0)(pyuadd)
        self._check_reduce(duadd, dtype=np.float64)

    def test_min_reduce(self):
        dumin = vectorize('int64(int64, int64)', identity='reorderable')(pymin)
        self._check_reduce(dumin, initial=10)
        self._check_reduce_axis(dumin, dtype=np.int64)

    def test_add_reduce_initial(self):
        # Initial should be used as a start
        duadd = vectorize('int64(int64, int64)', identity=0)(pyuadd)
        self._check_reduce(duadd, dtype=np.int64, initial=100)

    def test_add_reduce_no_initial_or_identity(self):
        # don't provide an initial or identity value
        duadd = vectorize('int64(int64, int64)')(pyuadd)
        self._check_reduce(duadd, dtype=np.int64)

    def test_invalid_input(self):
        duadd = vectorize('float64(float64, int64)', identity=0)(pyuadd)

        @njit
        def foo(a):
            return duadd.reduce(a)

        exc_msg = 'The first argument "array" must be array-like'
        with self.assertRaisesRegex(TypingError, exc_msg):
            foo('a')

    def test_dufunc_negative_axis(self):
        duadd = vectorize('int64(int64, int64)', identity=0)(pyuadd)

        @njit
        def foo(a, axis):
            return duadd.reduce(a, axis=axis)

        a = np.arange(40).reshape(5, 4, 2)
        cases = (0, -1, (0, -1), (-1, -2), (1, -1), -3)
        for axis in cases:
            expected = duadd.reduce(a, axis)
            got = foo(a, axis)
            self.assertPreciseEqual(expected, got)

    def test_dufunc_invalid_axis(self):
        duadd = vectorize('int64(int64, int64)', identity=0)(pyuadd)

        @njit
        def foo(a, axis):
            return duadd.reduce(a, axis=axis)

        a = np.arange(40).reshape(5, 4, 2)
        cases = ((0, 0), (0, 1, 0), (0, -3), (-1, -1), (-1, 2))
        for axis in cases:
            msg = "duplicate value in 'axis'"
            with self.assertRaisesRegex(ValueError, msg):
                foo(a, axis)

        cases = (-4, 3, (0, -4),)
        for axis in cases:
            with self.assertRaisesRegex(ValueError, "Invalid axis"):
                foo(a, axis)


class TestDUFuncPickling(MemoryLeakMixin, unittest.TestCase):
    def check(self, ident, result_type):
        buf = pickle.dumps(ident)
        rebuilt = pickle.loads(buf)

        # Check reconstructed dufunc
        r = rebuilt(123)
        self.assertEqual(123, r)
        self.assertIsInstance(r, result_type)

        # Try to use reconstructed dufunc in @jit
        @njit
        def foo(x):
            return rebuilt(x)

        r = foo(321)
        self.assertEqual(321, r)
        self.assertIsInstance(r, result_type)

    def test_unrestricted(self):
        @vectorize
        def ident(x1):
            return x1

        self.check(ident, result_type=(int, np.integer))

    def test_restricted(self):
        @vectorize(["float64(float64)"])
        def ident(x1):
            return x1

        self.check(ident, result_type=float)


if __name__ == "__main__":
    unittest.main()
