from pandas import *
from pandas.util.testing import rands
from pandas.compat import range

N = 1000
K = 50


def _random_index(howmany):
    return Index([rands(10) for _ in range(howmany)])

df = DataFrame(np.random.randn(N, K), index=_random_index(N),
               columns=_random_index(K))


def get1():
    for col in df.columns:
        for row in df.index:
            _ = df[col][row]


def get2():
    for col in df.columns:
        for row in df.index:
            _ = df.get_value(row, col)


def put1():
    for col in df.columns:
        for row in df.index:
            df[col][row] = 0


def put2():
    for col in df.columns:
        for row in df.index:
            df.set_value(row, col, 0)


def resize1():
    buf = DataFrame()
    for col in df.columns:
        for row in df.index:
            buf = buf.set_value(row, col, 5.)
    return buf


def resize2():
    from collections import defaultdict

    buf = defaultdict(dict)
    for col in df.columns:
        for row in df.index:
            buf[col][row] = 5.

    return DataFrame(buf)
