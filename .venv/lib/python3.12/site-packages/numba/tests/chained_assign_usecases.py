from numba import jit
import numpy as np


@jit
def inc(a):
    for i in range(len(a)):
        a[i] += 1
    return a


@jit
def inc1(a):
    a[0] += 1
    return a[0]


@jit
def inc2(a):
    a[0] += 1
    return a[0], a[0] + 1


def chain1(a):
    x = y = z = inc(a)
    return x + y + z


def chain2(v):
    a = np.zeros(2)
    a[0] = x = a[1] = v
    return a[0] + a[1] + (x / 2)


def unpack1(x, y):
    a, b = x, y
    return a + b / 2


def unpack2(x, y):
    a, b = c, d = inc1(x), inc1(y)
    return a + c / 2, b + d / 2


def chain3(x, y):
    a = (b, c) = (inc1(x), inc1(y))
    (d, e) = f = (inc1(x), inc1(y))
    return (a[0] + b / 2 + d + f[0]), (a[1] + c + e / 2 + f[1])


def unpack3(x):
    a, b = inc2(x)
    return a + b / 2


def unpack4(x):
    a, b = c, d = inc2(x)
    return a + c / 2, b + d / 2


def unpack5(x):
    a = b, c = inc2(x)
    d, e = f = inc2(x)
    return (a[0] + b / 2 + d + f[0]), (a[1] + c + e / 2 + f[1])


def unpack6(x, y):
    (a, b), (c, d) = (x, y), (y + 1, x + 1)
    return a + c / 2, b / 2 + d
