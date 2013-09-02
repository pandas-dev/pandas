#!/usr/bin/env python

"""
Microbenchmarks for comparison with R's "with" and "subset" functions
"""

from __future__ import print_function
from timeit import timeit


def bench_with(n=1e7, times=10, repeat=3):
    setup = "from pandas import DataFrame\n"
    setup += "from numpy.random import randn\n"
    setup += "df = DataFrame(randn(%d, 3), columns=list('abc'))\n" % n
    setup += "s = 'a + b * (c ** 2 + b ** 2 - a) / (a * c) ** 3'"
    print('DataFrame.eval:')
    print(timeit('df.eval(s)', setup=setup, repeat=repeat, number=times))


def bench_subset(n=1e7, times=10, repeat=3):
    setup = "from pandas import DataFrame\n"
    setup += "from numpy.random import randn\n"
    setup += "df = DataFrame(randn(%d, 3), columns=list('abc'))\n" % n
    setup += "s = 'a <= b <= (c ** 2 + b ** 2 - a) and b > c'"
    print('DataFrame.query:')
    print(timeit('df.query(s)', setup=setup, repeat=repeat, number=times))
    print('DataFrame.__getitem__:')
    print(timeit('df[s]', setup=setup, repeat=repeat, number=times))


def bench():
    bench_with()
    bench_subset()


if __name__ == '__main__':
    bench()
