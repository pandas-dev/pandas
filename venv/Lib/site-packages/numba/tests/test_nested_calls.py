"""
Test problems in nested calls.
Usually due to invalid type conversion between function boundaries.
"""


from numba import int32, int64
from numba import jit, njit
from numba.core import types
from numba.extending import overload
from numba.tests.support import TestCase, tag
import unittest


@jit(nopython=True)
def f_inner(a, b, c):
    return a, b, c

def f(x, y, z):
    return f_inner(x, c=y, b=z)

@jit(nopython=True)
def g_inner(a, b=2, c=3):
    return a, b, c

def g(x, y, z):
    return g_inner(x, b=y), g_inner(a=z, c=x)

@jit(nopython=True)
def star_inner(a=5, *b):
    return a, b

def star(x, y, z):
    return star_inner(a=x), star_inner(x, y, z)

def star_call(x, y, z):
    return star_inner(x, *y), star_inner(*z)

@jit(nopython=True)
def argcast_inner(a, b):
    if b:
        # Here `a` is unified to int64 (from int32 originally)
        a = int64(0)
    return a

def argcast(a, b):
    return argcast_inner(int32(a), b)


def generated_inner(x, y=5, z=6):
    assert 0, "unreachable"


@overload(generated_inner)
def ol_generated_inner(x, y=5, z=6):
    if isinstance(x, types.Complex):
        def impl(x, y=5, z=6):
            return x + y, z
    else:
        def impl(x, y=5, z=6):
            return x - y, z
    return impl


def call_generated(a, b):
    return generated_inner(a, z=b)


def nested_annotated() -> int:
    def inner(inner_arg: int) -> int:
        return inner_arg
    return 1


class TestNestedCall(TestCase):

    def compile_func(self, pyfunc, objmode=False):
        def check(*args, **kwargs):
            expected = pyfunc(*args, **kwargs)
            result = f(*args, **kwargs)
            self.assertPreciseEqual(result, expected)
        flags = dict(forceobj=True) if objmode else dict(nopython=True)
        f = jit(**flags)(pyfunc)
        return f, check

    def test_boolean_return(self):
        @jit(nopython=True)
        def inner(x):
            return not x

        @jit(nopython=True)
        def outer(x):
            if inner(x):
                return True
            else:
                return False

        self.assertFalse(outer(True))
        self.assertTrue(outer(False))

    def test_named_args(self, objmode=False):
        """
        Test a nested function call with named (keyword) arguments.
        """
        cfunc, check = self.compile_func(f, objmode)
        check(1, 2, 3)
        check(1, y=2, z=3)

    def test_named_args_objmode(self):
        self.test_named_args(objmode=True)

    def test_default_args(self, objmode=False):
        """
        Test a nested function call using default argument values.
        """
        cfunc, check = self.compile_func(g, objmode)
        check(1, 2, 3)
        check(1, y=2, z=3)

    def test_default_args_objmode(self):
        self.test_default_args(objmode=True)

    def test_star_args(self):
        """
        Test a nested function call to a function with *args in its signature.
        """
        cfunc, check = self.compile_func(star)
        check(1, 2, 3)

    def test_star_call(self, objmode=False):
        """
        Test a function call with a *args.
        """
        cfunc, check = self.compile_func(star_call, objmode)
        check(1, (2,), (3,))

    def test_star_call_objmode(self):
        self.test_star_call(objmode=True)

    def test_argcast(self):
        """
        Issue #1488: implicitly casting an argument variable should not
        break nested calls.
        """
        cfunc, check = self.compile_func(argcast)
        check(1, 0)
        check(1, 1)

    def test_call_generated(self):
        """
        Test a nested function call to a generated jit function.
        """
        cfunc = jit(nopython=True)(call_generated)
        self.assertPreciseEqual(cfunc(1, 2), (-4, 2))
        self.assertPreciseEqual(cfunc(1j, 2), (1j + 5, 2))

    def test_nested_annotated(self):
        """
        Tested a nested function with annotations.
        """
        cfunc = njit()(nested_annotated)
        # should not fail
        cfunc()


if __name__ == '__main__':
    unittest.main()
