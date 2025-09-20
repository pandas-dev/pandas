"""
Tests for enum support.
"""


import numpy as np
import unittest
from numba import jit, vectorize, int8, int16, int32

from numba.tests.support import TestCase
from numba.tests.enum_usecases import (Color, Shape, Shake,
                                       Planet, RequestError,
                                       IntEnumWithNegatives)


def compare_usecase(a, b):
    return a == b, a != b, a is b, a is not b


def getattr_usecase(a):
    # Lookup of a enum member on its class
    return a is Color.red


def getitem_usecase(a):
    """Lookup enum member by string name"""
    return a is Color['red']


def identity_usecase(a, b, c):
    return (a is Shake.mint,
            b is Shape.circle,
            c is RequestError.internal_error,
            )


def make_constant_usecase(const):
    def constant_usecase(a):
        return a is const
    return constant_usecase


def return_usecase(a, b, pred):
    return a if pred else b


def int_coerce_usecase(x):
    # Implicit coercion of intenums to ints
    if x > RequestError.internal_error:
        return x - RequestError.not_found
    else:
        return x + Shape.circle

def int_cast_usecase(x):
    # Explicit coercion of intenums to ints
    if x > int16(RequestError.internal_error):
        return x - int32(RequestError.not_found)
    else:
        return x + int16(Shape.circle)


def vectorize_usecase(x):
    if x != RequestError.not_found:
        return RequestError['internal_error']
    else:
        return RequestError.dummy


class BaseEnumTest(object):

    def test_compare(self):
        pyfunc = compare_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for args in self.pairs:
            self.assertPreciseEqual(pyfunc(*args), cfunc(*args))

    def test_return(self):
        """
        Passing and returning enum members.
        """
        pyfunc = return_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for pair in self.pairs:
            for pred in (True, False):
                args = pair + (pred,)
                self.assertIs(pyfunc(*args), cfunc(*args))

    def check_constant_usecase(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)

        for arg in self.values:
            self.assertPreciseEqual(pyfunc(arg), cfunc(arg))

    def test_constant(self):
        self.check_constant_usecase(getattr_usecase)
        self.check_constant_usecase(getitem_usecase)
        self.check_constant_usecase(make_constant_usecase(self.values[0]))


class TestEnum(BaseEnumTest, TestCase):
    """
    Tests for Enum classes and members.
    """
    values = [Color.red, Color.green]

    pairs = [
        (Color.red, Color.red),
        (Color.red, Color.green),
        (Shake.mint, Shake.vanilla),
        (Planet.VENUS, Planet.MARS),
        (Planet.EARTH, Planet.EARTH),
        ]

    def test_identity(self):
        """
        Enum with equal values should not compare identical
        """
        pyfunc = identity_usecase
        cfunc = jit(nopython=True)(pyfunc)
        args = (Color.blue, Color.green, Shape.square)
        self.assertPreciseEqual(pyfunc(*args), cfunc(*args))


class TestIntEnum(BaseEnumTest, TestCase):
    """
    Tests for IntEnum classes and members.
    """
    values = [Shape.circle, Shape.square]

    pairs = [
        (Shape.circle, Shape.circle),
        (Shape.circle, Shape.square),
        (RequestError.not_found, RequestError.not_found),
        (RequestError.internal_error, RequestError.not_found),
        ]

    def test_int_coerce(self):
        pyfunc = int_coerce_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for arg in [300, 450, 550]:
            self.assertPreciseEqual(pyfunc(arg), cfunc(arg))

    def test_int_cast(self):
        pyfunc = int_cast_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for arg in [300, 450, 550]:
            self.assertPreciseEqual(pyfunc(arg), cfunc(arg))

    def test_vectorize(self):
        cfunc = vectorize(nopython=True)(vectorize_usecase)
        arg = np.array([2, 404, 500, 404])
        sol = np.array([vectorize_usecase(i) for i in arg], dtype=arg.dtype)
        self.assertPreciseEqual(sol, cfunc(arg))

    def test_hash(self):
        def pyfun(x):
            return hash(x)
        cfunc = jit(nopython=True)(pyfun)
        for member in IntEnumWithNegatives:
            self.assertPreciseEqual(pyfun(member), cfunc(member))

    def test_int_shape_cast(self):
        def pyfun_empty(x):
            return np.empty((x, x), dtype='int64').fill(-1)
        def pyfun_zeros(x):
            return np.zeros((x, x), dtype='int64')
        def pyfun_ones(x):
            return np.ones((x, x), dtype='int64')
        for pyfun in [pyfun_empty, pyfun_zeros, pyfun_ones]:
            cfunc = jit(nopython=True)(pyfun)
            for member in IntEnumWithNegatives:
                if member >= 0:
                    self.assertPreciseEqual(pyfun(member), cfunc(member))


if __name__ == '__main__':
    unittest.main()
