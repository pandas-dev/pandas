import random
import numpy as np

from numba.tests.support import TestCase, captured_stdout
from numba import njit, literally
from numba.core import types
from numba.cpython.unsafe.tuple import tuple_setitem, build_full_slice_tuple
from numba.np.unsafe.ndarray import to_fixed_tuple, empty_inferred
from numba.core.unsafe.bytes import memcpy_region
from numba.core.unsafe.refcount import dump_refcount
from numba.cpython.unsafe.numbers import trailing_zeros, leading_zeros
from numba.core.errors import TypingError


class TestTupleIntrinsic(TestCase):
    """Tests for numba.unsafe.tuple
    """
    def test_tuple_setitem(self):
        @njit
        def foo(tup, idxs, vals):
            out_tup = tup
            for i, v in zip(idxs, vals):
                out_tup = tuple_setitem(out_tup, i, v)
            return tup, out_tup

        random.seed(123)
        for _ in range(20):
            # Random data
            n = random.randint(1, 10)
            tup = tuple([random.randint(0, n) for i in range(n)])
            vals = tuple([random.randint(10, 20) for i in range(n)])
            idxs = list(range(len(vals)))
            random.shuffle(idxs)
            idxs = tuple(idxs)
            # Expect
            expect_tup = tuple(tup)
            expect_out = np.asarray(expect_tup)
            expect_out[np.asarray(idxs)] = vals
            # Got
            got_tup, got_out = foo(tup, idxs, vals)
            # Check
            self.assertEqual(got_tup, expect_tup)
            self.assertEqual(got_out, tuple(expect_out))

    def test_slice_tuple(self):
        @njit
        def full_slice_array(a, n):
            # Since numba slices can't be boxed at the moment
            return a[build_full_slice_tuple(literally(n))]

        for n in range(1, 3):
            a = np.random.random(np.arange(n) + 1)
            for i in range(1, n + 1):
                np.testing.assert_array_equal(a, full_slice_array(a, i))
            with self.assertRaises(TypingError):
                # numpy would throw an IndexError here
                full_slice_array(a, n + 1)


class TestNdarrayIntrinsic(TestCase):
    """Tests for numba.unsafe.ndarray
    """
    def test_to_fixed_tuple(self):
        const = 3

        @njit
        def foo(array):
            a = to_fixed_tuple(array, length=1)
            b = to_fixed_tuple(array, 2)
            c = to_fixed_tuple(array, const)
            d = to_fixed_tuple(array, 0)
            return a, b, c, d

        np.random.seed(123)
        for _ in range(10):
            # Random data
            arr = np.random.random(3)
            # Run
            a, b, c, d = foo(arr)
            # Check
            self.assertEqual(a, tuple(arr[:1]))
            self.assertEqual(b, tuple(arr[:2]))
            self.assertEqual(c, tuple(arr[:3]))
            self.assertEqual(d, ())

        # Check error with ndim!=1
        with self.assertRaises(TypingError) as raises:
            foo(np.random.random((1, 2)))
        self.assertIn("Not supported on array.ndim=2",
                      str(raises.exception))

        # Check error with non-constant length
        @njit
        def tuple_with_length(array, length):
            return to_fixed_tuple(array, length)

        with self.assertRaises(TypingError) as raises:
            tuple_with_length(np.random.random(3), 1)
        expectmsg = "*length* argument must be a constant"
        self.assertIn(expectmsg, str(raises.exception))

    def test_issue_3586_variant1(self):
        @njit
        def func():
            S = empty_inferred((10,))
            a = 1.1
            for i in range(len(S)):
                S[i] = a + 2
            return S

        got = func()
        expect = np.asarray([3.1] * 10)
        np.testing.assert_array_equal(got, expect)

    def test_issue_3586_variant2(self):
        @njit
        def func():
            S = empty_inferred((10,))
            a = 1.1
            for i in range(S.size):
                S[i] = a + 2
            return S

        got = func()
        expect = np.asarray([3.1] * 10)
        np.testing.assert_array_equal(got, expect)


class TestBytesIntrinsic(TestCase):
    """Tests for numba.unsafe.bytes
    """
    def test_memcpy_region(self):
        @njit
        def foo(dst, dst_index, src, src_index, nbytes):
            # last arg is assume 1 byte alignment
            memcpy_region(dst.ctypes.data, dst_index,
                          src.ctypes.data, src_index, nbytes, 1)

        d = np.zeros(10, dtype=np.int8)
        s = np.arange(10, dtype=np.int8)

        # copy s[1:6] to d[4:9]
        foo(d, 4, s, 1, 5)

        expected = [0, 0, 0, 0, 1, 2, 3, 4, 5, 0]
        np.testing.assert_array_equal(d, expected)


class TestRefCount(TestCase):
    def test_dump_refcount(self):
        @njit
        def use_dump_refcount():
            a = np.ones(10)
            b = (a, a)
            dump_refcount(a)
            dump_refcount(b)

        # Capture output to sys.stdout
        with captured_stdout() as stream:
            use_dump_refcount()

        output = stream.getvalue()
        # Check that it printed
        pat = "dump refct of {}"
        aryty = types.float64[::1]
        tupty = types.Tuple.from_types([aryty] * 2)
        self.assertIn(pat.format(aryty), output)
        self.assertIn(pat.format(tupty), output)


class TestZeroCounts(TestCase):
    def test_zero_count(self):
        lz = njit(lambda x: leading_zeros(x))
        tz = njit(lambda x: trailing_zeros(x))

        evens = [2, 42, 126, 128]

        for T in types.unsigned_domain:
            self.assertTrue(tz(T(0)) == lz(T(0)) == T.bitwidth)
            for i in range(T.bitwidth):
                val = T(2 ** i)
                self.assertEqual(lz(val) + tz(val) + 1, T.bitwidth)
            for n in evens:
                self.assertGreater(tz(T(n)), 0)
                self.assertEqual(tz(T(n + 1)), 0)

        for T in types.signed_domain:
            self.assertTrue(tz(T(0)) == lz(T(0)) == T.bitwidth)
            for i in range(T.bitwidth - 1):
                val = T(2 ** i)
                self.assertEqual(lz(val) + tz(val) + 1, T.bitwidth)
                self.assertEqual(lz(-val), 0)
                self.assertEqual(tz(val), tz(-val))
            for n in evens:
                if not T.minval <= n <= T.maxval:
                    continue
                self.assertGreater(tz(T(n)), 0)
                self.assertEqual(tz(T(n + 1)), 0)

    def check_error_msg(self, func):
        cfunc = njit(lambda *x: func(*x))
        func_name = func._name

        unsupported_types = filter(
            lambda x: not isinstance(x, types.Integer), types.number_domain
        )
        for typ in sorted(unsupported_types, key=str):
            with self.assertRaises(TypingError) as e:
                cfunc(typ(2))
            self.assertIn(
                "{} is only defined for integers, but value passed was '{}'."
                .format(func_name, typ),
                str(e.exception),
            )

        # Testing w/ too many/few arguments
        def check(args, string):
            with self.assertRaises((TypingError, TypeError)) as e:
                cfunc(*args)
            self.assertIn(
                "{}() ".format(func_name),
                str(e.exception)
            )

        check((1, 2), "takes 2 positional arguments but 3 were given")
        check((), "missing 1 required positional argument")

    def test_trailing_zeros_error(self):
        self.check_error_msg(trailing_zeros)

    def test_leading_zeros_error(self):
        self.check_error_msg(leading_zeros)
