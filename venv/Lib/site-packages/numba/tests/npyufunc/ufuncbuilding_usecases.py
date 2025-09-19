from numba import vectorize


def add(a, b):
    """An addition"""
    return a + b


def equals(a, b):
    return a == b


def mul(a, b):
    """A multiplication"""
    return a * b


def guadd(a, b, c):
    """A generalized addition"""
    x, y = c.shape
    for i in range(x):
        for j in range(y):
            c[i, j] = a[i, j] + b[i, j]


@vectorize(nopython=True)
def inner(a, b):
    return a + b


@vectorize(["int64(int64, int64)"], nopython=True)
def inner_explicit(a, b):
    return a + b


def outer(a, b):
    return inner(a, b)


def outer_explicit(a, b):
    return inner_explicit(a, b)


class Dummy:
    pass


def guadd_obj(a, b, c):
    Dummy()  # to force object mode
    x, y = c.shape
    for i in range(x):
        for j in range(y):
            c[i, j] = a[i, j] + b[i, j]


def guadd_scalar_obj(a, b, c):
    Dummy()  # to force object mode
    x, y = c.shape
    for i in range(x):
        for j in range(y):
            c[i, j] = a[i, j] + b


class MyException(Exception):
    pass


def guerror(a, b, c):
    raise MyException
