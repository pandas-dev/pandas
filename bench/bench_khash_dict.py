"""
Some comparisons of khash.h to Python dict
"""
from __future__ import print_function

import numpy as np
import os

from vbench.api import Benchmark
from pandas.util.testing import rands
from pandas.compat import range
import pandas._tseries as lib
import pandas._sandbox as sbx
import time

import psutil

pid = os.getpid()
proc = psutil.Process(pid)


def object_test_data(n):
    pass


def string_test_data(n):
    return np.array([rands(10) for _ in range(n)], dtype='O')


def int_test_data(n):
    return np.arange(n, dtype='i8')

N = 1000000

#----------------------------------------------------------------------
# Benchmark 1: map_locations


def map_locations_python_object():
    arr = string_test_data(N)
    return _timeit(lambda: lib.map_indices_object(arr))


def map_locations_khash_object():
    arr = string_test_data(N)

    def f():
        table = sbx.PyObjectHashTable(len(arr))
        table.map_locations(arr)
    return _timeit(f)


def _timeit(f, iterations=10):
    start = time.time()
    for _ in range(iterations):
        foo = f()
    elapsed = time.time() - start
    return elapsed

#----------------------------------------------------------------------
# Benchmark 2: lookup_locations


def lookup_python(values):
    table = lib.map_indices_object(values)
    return _timeit(lambda: lib.merge_indexer_object(values, table))


def lookup_khash(values):
    table = sbx.PyObjectHashTable(len(values))
    table.map_locations(values)
    locs = table.lookup_locations(values)
    # elapsed = _timeit(lambda: table.lookup_locations2(values))
    return table


def leak(values):
    for _ in range(100):
        print(proc.get_memory_info())
        table = lookup_khash(values)
        # table.destroy()

arr = string_test_data(N)

#----------------------------------------------------------------------
# Benchmark 3: unique

#----------------------------------------------------------------------
# Benchmark 4: factorize
