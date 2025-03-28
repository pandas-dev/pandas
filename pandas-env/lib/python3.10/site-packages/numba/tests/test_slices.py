from functools import partial
import itertools
from itertools import chain, product, starmap
import sys

import numpy as np

from numba import jit, literally, njit, typeof, TypingError
from numba.core import utils, types
from numba.tests.support import TestCase, MemoryLeakMixin
from numba.core.types.functions import _header_lead
import unittest


def slice_passing(sl):
    return sl.start, sl.stop, sl.step

def slice_constructor(*args):
    sl = slice(*args)
    return sl.start, sl.stop, sl.step

def slice_construct_and_use(args, l):
    sl = slice(*args)
    return l[sl]

def slice_indices(s, *indargs):
    return s.indices(*indargs)

class TestSlices(MemoryLeakMixin, TestCase):

    def test_slice_passing(self):
        """
        Check passing a slice object to a Numba function.
        """
        # NOTE this also checks slice attributes
        def check(a, b, c, d, e, f):
            sl = slice(a, b, c)
            got = cfunc(sl)
            self.assertPreciseEqual(got, (d, e, f))

        maxposint = sys.maxsize
        maxnegint = -maxposint - 1
        cfunc = jit(nopython=True)(slice_passing)

        # Positive steps
        start_cases = [(None, 0), (42, 42), (-1, -1)]
        stop_cases = [(None, maxposint), (9, 9), (-11, -11)]
        step_cases = [(None, 1), (12, 12)]
        for (a, d), (b, e), (c, f) in itertools.product(start_cases,
                                                        stop_cases,
                                                        step_cases):
            check(a, b, c, d, e, f)

        # Negative steps
        start_cases = [(None, maxposint), (42, 42), (-1, -1)]
        stop_cases = [(None, maxnegint), (9, 9), (-11, -11)]
        step_cases = [(-1, -1), (-12, -12)]
        for (a, d), (b, e), (c, f) in itertools.product(start_cases,
                                                        stop_cases,
                                                        step_cases):
            check(a, b, c, d, e, f)

        # Some member is neither integer nor None
        with self.assertRaises(TypeError):
            cfunc(slice(1.5, 1, 1))

    def test_slice_constructor(self):
        """
        Test the 'happy path' for slice() constructor in nopython mode.
        """
        maxposint = sys.maxsize
        maxnegint = -maxposint - 1
        a = np.arange(10)
        cfunc = jit(nopython=True)(slice_constructor)
        cfunc_use = jit(nopython=True)(slice_construct_and_use)
        for args, expected in [
            ((None,), (0, maxposint, 1)),
            ((5,), (0, 5, 1)),
            ((None, None), (0, maxposint, 1)),
            ((1, None), (1, maxposint, 1)),
            ((None, 2), (0, 2, 1)),
            ((1, 2), (1, 2, 1)),
            ((None, None, 3), (0, maxposint, 3)),
            ((None, 2, 3), (0, 2, 3)),
            ((1, None, 3), (1, maxposint, 3)),
            ((1, 2, 3), (1, 2, 3)),
            ((None, None, -1), (maxposint, maxnegint, -1)),
            ((10, None, -1), (10, maxnegint, -1)),
            ((None, 5, -1), (maxposint, 5, -1)),
            ((10, 5, -1), (10, 5, -1)),
        ]:
            got = cfunc(*args)
            self.assertPreciseEqual(got, expected)
            usage = slice_construct_and_use(args, a)
            cusage = cfunc_use(args, a)
            self.assertPreciseEqual(usage, cusage)

    def test_slice_constructor_cases(self):
        """
        Test that slice constructor behaves same in python and compiled code.
        """
        options = (None, -1, 0, 1)
        arg_cases = chain.from_iterable(
            product(options, repeat=n) for n in range(5)
        )
        array = np.arange(10)

        cfunc = jit(nopython=True)(slice_construct_and_use)

        self.disable_leak_check()
        for args in arg_cases:
            try:
                expected = slice_construct_and_use(args, array)
            except TypeError as py_type_e:
                # Catch cases of 0, or more than 3 arguments.
                # This becomes a typing error in numba
                n_args = len(args)
                self.assertRegex(
                    str(py_type_e),
                    r"slice expected at (most|least) (3|1) arguments?, got {}"
                    .format(n_args)
                )
                with self.assertRaises(TypingError) as numba_e:
                    cfunc(args, array)
                self.assertIn(
                    _header_lead,
                    str(numba_e.exception)
                )
                self.assertIn(
                    ", ".join(str(typeof(arg)) for arg in args),
                    str(numba_e.exception)
                )
            except Exception as py_e:
                with self.assertRaises(type(py_e)) as numba_e:
                    cfunc(args, array)
                self.assertIn(
                    str(py_e),
                    str(numba_e.exception)
                )
            else:
                self.assertPreciseEqual(expected, cfunc(args, array))

    def test_slice_indices(self):
        """Test that a numba slice returns same result for .indices as a python one."""
        slices = starmap(
            slice,
            product(
                chain(range(-5, 5), (None,)),
                chain(range(-5, 5), (None,)),
                chain(range(-5, 5), (None,))
            )
        )
        lengths = range(-2, 3)

        cfunc = jit(nopython=True)(slice_indices)

        for s, l in product(slices, lengths):
            try:
                expected = slice_indices(s, l)
            except Exception as py_e:
                with self.assertRaises(type(py_e)) as numba_e:
                    cfunc(s, l)
                self.assertIn(
                    str(py_e),
                    str(numba_e.exception)
                )
            else:
                self.assertPreciseEqual(expected, cfunc(s, l))

    def test_slice_indices_examples(self):
        """Tests for specific error cases."""
        cslice_indices = jit(nopython=True)(slice_indices)

        with self.assertRaises(TypingError) as e:
            cslice_indices(slice(None), 1, 2, 3)
        self.assertIn(
             "indices() takes exactly one argument (3 given)",
             str(e.exception)
        )

        with self.assertRaises(TypingError) as e:
            cslice_indices(slice(None, None, 0), 1.2)
        self.assertIn(
            "'%s' object cannot be interpreted as an integer" % typeof(1.2),
            str(e.exception)
        )

    def test_slice_from_constant(self):
        test_tuple = (1, 2, 3, 4)

        for ts in itertools.product(
            [None, 1, 2, 3], [None, 1, 2, 3], [None, 1, 2, -1, -2]
        ):
            ts = slice(*ts)

            @jit(nopython=True)
            def test_fn():
                return test_tuple[ts]

            self.assertEqual(test_fn(), test_fn.py_func())

    def test_literal_slice_distinct(self):
        sl1 = types.misc.SliceLiteral(slice(1, None, None))
        sl2 = types.misc.SliceLiteral(slice(None, None, None))
        sl3 = types.misc.SliceLiteral(slice(1, None, None))

        self.assertNotEqual(sl1, sl2)
        self.assertEqual(sl1, sl3)

    def test_literal_slice_boxing(self):
        # Tests that a literal slice can be
        # returned from a JIT function.
        @njit
        def f(x):
            return literally(x)

        slices = (
            slice(1, 4, 2),
            slice(1, 2),
            slice(1),
            slice(None, 1, 1),
            slice(1, None, 1),
            slice(None, None, 1),
            slice(None),
            slice(None, None, None)
        )
        for sl in slices:
            self.assertEqual(sl, f(sl))


    def test_literal_slice_freevar(self):
        # Tests passing a literal slice as a freevar
        # in a closure.
        z = slice(1, 2, 3)
        @njit
        def foo():
            return z

        self.assertEqual(z, foo())

    def test_literal_slice_maxint(self):
        # Tests that passing a slice with an integer
        # that exceeds the maxint size throws a reasonable
        # error message.
        @njit()
        def foo(z):
            return literally(z)

        maxval = int(2**63)
        with self.assertRaises(ValueError) as e:
            foo(slice(None, None, -maxval-1))
        self.assertIn(
            "Int value is too large",
            str(e.exception)
        )


if __name__ == '__main__':
    unittest.main()
