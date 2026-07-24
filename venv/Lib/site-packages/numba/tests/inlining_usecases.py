""" Test cases for inlining IR from another module """
from numba import jit, njit
from numba.core import types
from numba.core.extending import overload

_GLOBAL1 = 100


@njit(inline='always')
def bar():
    return _GLOBAL1 + 10


def baz_factory(a):
    b = 17 + a

    @njit(inline='always')
    def baz():
        return _GLOBAL1 + a - b
    return baz


def baz():
    return _GLOBAL1 + 10


@overload(baz, inline='always')
def baz_ol():
    def impl():
        return _GLOBAL1 + 10
    return impl


def bop_factory(a):
    b = 17 + a

    def bop():
        return _GLOBAL1 + a - b

    @overload(bop, inline='always')
    def baz():
        def impl():
            return _GLOBAL1 + a - b
        return impl

    return bop


@jit((types.int32,), nopython=True)
def inner(a):
    return a + 1


@jit((types.int32,), nopython=True)
def more(a):
    return inner(inner(a))


def outer_simple(a):
    return inner(a) * 2


def outer_multiple(a):
    return inner(a) * more(a)


@njit
def __dummy__():
    return
