import unittest
from numba import jit, njit
from numba.core import types
from numba.tests.support import TestCase


force_pyobj_jit_opt = {'forceobj': True}
no_pyobj_jit_opt = {'nopython': True}


def assignments(a):
    b = c = str(a)
    return b + c


def assignments2(a):
    b = c = d = str(a)
    return b + c + d


# Use cases for issue #503

def var_propagate1(a, b):
    c = (a if a > b else b) + 5
    return c


def var_propagate2(a, b):
    c = 5 + (a if a > b else b + 12) / 2.0
    return c


def var_propagate3(a, b):
    c = 5 + (a > b and a or b)
    return c


def var_propagate4(a, b):
    c = 5 + (a - 1 and b + 1) or (a + 1 and b - 1)
    return c


# Issue #480
def chained_compare(a):
    return 1 < a < 3


# Issue #591
def stack_effect_error(x):
    i = 2
    c = 1
    if i == x:
        for i in range(3):
            c = i
    return i + c

# Some more issues with stack effect and blocks
def for_break(n, x):
    for i in range(n):
        n = 0
        if i == x:
            break
    else:
        n = i
    return i, n

# Issue #571
def var_swapping(a, b, c, d, e):
    a, b = b, a
    c, d, e = e, c, d
    a, b, c, d = b, c, d, a
    return a + b + c + d +e

class TestDataFlow(TestCase):

    def test_assignments(self, flags=force_pyobj_jit_opt):
        pyfunc = assignments
        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-1, 0, 1]:
            self.assertPreciseEqual(pyfunc(x), cfunc(x))

    def test_assignments2(self, flags=force_pyobj_jit_opt):
        pyfunc = assignments2
        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-1, 0, 1]:
            self.assertPreciseEqual(pyfunc(x), cfunc(x))

        if flags is force_pyobj_jit_opt:
            cfunc("a")

    # The dataflow analysis must be good enough for native mode
    # compilation to succeed, hence the use of njit in the following tests.

    def run_propagate_func(self, func, args):
        self.assertPreciseEqual(func(*args), func.py_func(*args))

    def test_var_propagate1(self):
        cfunc = njit((types.intp, types.intp))(var_propagate1)
        self.run_propagate_func(cfunc, (2, 3))
        self.run_propagate_func(cfunc, (3, 2))

    def test_var_propagate2(self):
        cfunc = njit((types.intp, types.intp))(var_propagate2)
        self.run_propagate_func(cfunc, (2, 3))
        self.run_propagate_func(cfunc, (3, 2))

    def test_var_propagate3(self):
        cfunc = njit((types.intp, types.intp))(var_propagate3)
        self.run_propagate_func(cfunc, (2, 3))
        self.run_propagate_func(cfunc, (3, 2))
        self.run_propagate_func(cfunc, (2, 0))
        self.run_propagate_func(cfunc, (-1, 0))
        self.run_propagate_func(cfunc, (0, 2))
        self.run_propagate_func(cfunc, (0, -1))

    def test_var_propagate4(self):
        cfunc = njit((types.intp, types.intp))(var_propagate4)
        self.run_propagate_func(cfunc, (1, 1))
        self.run_propagate_func(cfunc, (1, 0))
        self.run_propagate_func(cfunc, (1, -1))
        self.run_propagate_func(cfunc, (0, 1))
        self.run_propagate_func(cfunc, (0, 0))
        self.run_propagate_func(cfunc, (0, -1))
        self.run_propagate_func(cfunc, (-1, 1))
        self.run_propagate_func(cfunc, (-1, 0))
        self.run_propagate_func(cfunc, (-1, -1))

    def test_chained_compare(self, flags=force_pyobj_jit_opt):
        pyfunc = chained_compare
        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [0, 1, 2, 3, 4]:
            self.assertPreciseEqual(pyfunc(x), cfunc(x))

    def test_chained_compare_npm(self):
        self.test_chained_compare(no_pyobj_jit_opt)

    def test_stack_effect_error(self, flags=force_pyobj_jit_opt):
        # Issue #591: POP_BLOCK must undo all stack pushes done inside
        # the block.
        pyfunc = stack_effect_error
        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in (0, 1, 2, 3):
            self.assertPreciseEqual(pyfunc(x), cfunc(x))

    def test_stack_effect_error_npm(self):
        self.test_stack_effect_error(no_pyobj_jit_opt)

    def test_var_swapping(self, flags=force_pyobj_jit_opt):
        pyfunc = var_swapping
        cfunc = jit((types.int32,) * 5, **flags)(pyfunc)
        args = tuple(range(0, 10, 2))
        self.assertPreciseEqual(pyfunc(*args), cfunc(*args))

    def test_var_swapping_npm(self):
        self.test_var_swapping(no_pyobj_jit_opt)

    def test_for_break(self, flags=force_pyobj_jit_opt):
        # BREAK_LOOP must unwind the current inner syntax block.
        pyfunc = for_break
        cfunc = jit((types.intp, types.intp), **flags)(pyfunc)
        for (n, x) in [(4, 2), (4, 6)]:
            self.assertPreciseEqual(pyfunc(n, x), cfunc(n, x))

    def test_for_break_npm(self):
        self.test_for_break(no_pyobj_jit_opt)


if __name__ == '__main__':
    unittest.main()

