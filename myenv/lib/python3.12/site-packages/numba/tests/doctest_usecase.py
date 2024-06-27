"""
Test that all docstrings are the same:

>>> len({f.__doc__ for f in (a, b, c, d)})
1
"""
from numba import guvectorize, int64, njit, vectorize


def a():
    """>>> x = 1"""
    return 1


@njit
def b():
    """>>> x = 1"""
    return 1


@guvectorize([(int64[:], int64, int64[:])], "(n),()->(n)")
def c(x, y, res):
    """>>> x = 1"""
    for i in range(x.shape[0]):
        res[i] = x[i] + y


@vectorize([int64(int64, int64)])
def d(x, y):
    """>>> x = 1"""
    return x + y
