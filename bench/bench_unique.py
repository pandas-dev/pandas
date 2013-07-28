from __future__ import print_function
from pandas import *
from pandas.util.testing import rands
from pandas.compat import range, zip
import pandas._tseries as lib
import numpy as np
import matplotlib.pyplot as plt

N = 50000
K = 10000

groups = np.array([rands(10) for _ in range(K)], dtype='O')
groups2 = np.array([rands(10) for _ in range(K)], dtype='O')

labels = np.tile(groups, N // K)
labels2 = np.tile(groups2, N // K)
data = np.random.randn(N)


def timeit(f, niter):
    import gc
    import time
    gc.disable()
    start = time.time()
    for _ in range(niter):
        f()
    elapsed = (time.time() - start) / niter
    gc.enable()
    return elapsed


def algo1():
    unique_labels = np.unique(labels)
    result = np.empty(len(unique_labels))
    for i, label in enumerate(unique_labels):
        result[i] = data[labels == label].sum()


def algo2():
    unique_labels = np.unique(labels)
    indices = lib.groupby_indices(labels)
    result = np.empty(len(unique_labels))

    for i, label in enumerate(unique_labels):
        result[i] = data.take(indices[label]).sum()


def algo3_nosort():
    rizer = lib.DictFactorizer()
    labs, counts = rizer.factorize(labels, sort=False)
    k = len(rizer.uniques)
    out = np.empty(k)
    lib.group_add(out, counts, data, labs)


def algo3_sort():
    rizer = lib.DictFactorizer()
    labs, counts = rizer.factorize(labels, sort=True)
    k = len(rizer.uniques)
    out = np.empty(k)
    lib.group_add(out, counts, data, labs)

import numpy as np
import random


# dict to hold results
counts = {}

# a hack to generate random key, value pairs.
# 5k keys, 100k values
x = np.tile(np.arange(5000, dtype='O'), 20)
random.shuffle(x)
xarr = x
x = [int(y) for y in x]
data = np.random.uniform(0, 1, 100000)


def f():
    # groupby sum
    for k, v in zip(x, data):
        try:
            counts[k] += v
        except KeyError:
            counts[k] = v


def f2():
    rizer = lib.DictFactorizer()
    labs, counts = rizer.factorize(xarr, sort=False)
    k = len(rizer.uniques)
    out = np.empty(k)
    lib.group_add(out, counts, data, labs)


def algo4():
    rizer = lib.DictFactorizer()
    labs1, _ = rizer.factorize(labels, sort=False)
    k1 = len(rizer.uniques)

    rizer = lib.DictFactorizer()
    labs2, _ = rizer.factorize(labels2, sort=False)
    k2 = len(rizer.uniques)

    group_id = labs1 * k2 + labs2
    max_group = k1 * k2

    if max_group > 1e6:
        rizer = lib.Int64Factorizer(len(group_id))
        group_id, _ = rizer.factorize(group_id.astype('i8'), sort=True)
        max_group = len(rizer.uniques)

    out = np.empty(max_group)
    counts = np.zeros(max_group, dtype='i4')
    lib.group_add(out, counts, data, group_id)

# cumtime  percall filename:lineno(function)
#   0.592    0.592 <string>:1(<module>)
  # 0.584    0.006 groupby_ex.py:37(algo3_nosort)
  # 0.535    0.005 {method 'factorize' of DictFactorizer' objects}
  # 0.047    0.000 {pandas._tseries.group_add}
  # 0.002    0.000 numeric.py:65(zeros_like)
  # 0.001    0.000 {method 'fill' of 'numpy.ndarray' objects}
  # 0.000    0.000 {numpy.core.multiarray.empty_like}
  # 0.000    0.000 {numpy.core.multiarray.empty}

# UNIQUE timings

# N = 10000000
# K = 500000

# groups = np.array([rands(10) for _ in range(K)], dtype='O')

# labels = np.tile(groups, N // K)
data = np.random.randn(N)

data = np.random.randn(N)

Ks = [100, 1000, 5000, 10000, 25000, 50000, 100000]

# Ks = [500000, 1000000, 2500000, 5000000, 10000000]

import psutil
import os
import gc

pid = os.getpid()
proc = psutil.Process(pid)


def dict_unique(values, expected_K, sort=False, memory=False):
    if memory:
        gc.collect()
        before_mem = proc.get_memory_info().rss

    rizer = lib.DictFactorizer()
    result = rizer.unique_int64(values)

    if memory:
        result = proc.get_memory_info().rss - before_mem
        return result

    if sort:
        result.sort()
    assert(len(result) == expected_K)
    return result


def khash_unique(values, expected_K, size_hint=False, sort=False,
                 memory=False):
    if memory:
        gc.collect()
        before_mem = proc.get_memory_info().rss

    if size_hint:
        rizer = lib.Factorizer(len(values))
    else:
        rizer = lib.Factorizer(100)

    result = []
    result = rizer.unique(values)

    if memory:
        result = proc.get_memory_info().rss - before_mem
        return result

    if sort:
        result.sort()
    assert(len(result) == expected_K)


def khash_unique_str(values, expected_K, size_hint=False, sort=False,
                     memory=False):
    if memory:
        gc.collect()
        before_mem = proc.get_memory_info().rss

    if size_hint:
        rizer = lib.StringHashTable(len(values))
    else:
        rizer = lib.StringHashTable(100)

    result = []
    result = rizer.unique(values)

    if memory:
        result = proc.get_memory_info().rss - before_mem
        return result

    if sort:
        result.sort()
    assert(len(result) == expected_K)


def khash_unique_int64(values, expected_K, size_hint=False, sort=False):
    if size_hint:
        rizer = lib.Int64HashTable(len(values))
    else:
        rizer = lib.Int64HashTable(100)

    result = []
    result = rizer.unique(values)

    if sort:
        result.sort()
    assert(len(result) == expected_K)


def hash_bench():
    numpy = []
    dict_based = []
    dict_based_sort = []
    khash_hint = []
    khash_nohint = []
    for K in Ks:
        print(K)
        # groups = np.array([rands(10) for _ in range(K)])
        # labels = np.tile(groups, N // K).astype('O')

        groups = np.random.randint(0, long(100000000000), size=K)
        labels = np.tile(groups, N // K)
        dict_based.append(timeit(lambda: dict_unique(labels, K), 20))
        khash_nohint.append(timeit(lambda: khash_unique_int64(labels, K), 20))
        khash_hint.append(timeit(lambda: khash_unique_int64(labels, K,
                                                            size_hint=True), 20))

        # memory, hard to get
        # dict_based.append(np.mean([dict_unique(labels, K, memory=True)
        #                            for _ in range(10)]))
        # khash_nohint.append(np.mean([khash_unique(labels, K, memory=True)
        #                              for _ in range(10)]))
        # khash_hint.append(np.mean([khash_unique(labels, K, size_hint=True, memory=True)
        #                            for _ in range(10)]))

        # dict_based_sort.append(timeit(lambda: dict_unique(labels, K,
        #                                                   sort=True), 10))
        # numpy.append(timeit(lambda: np.unique(labels), 10))

    # unique_timings = DataFrame({'numpy.unique' : numpy,
    #                             'dict, no sort' : dict_based,
    #                             'dict, sort' : dict_based_sort},
    #                            columns=['dict, no sort',
    #                                     'dict, sort', 'numpy.unique'],
    #                            index=Ks)

    unique_timings = DataFrame({'dict': dict_based,
                                'khash, preallocate': khash_hint,
                                'khash': khash_nohint},
                               columns=['khash, preallocate', 'khash', 'dict'],
                               index=Ks)

    unique_timings.plot(kind='bar', legend=False)
    plt.legend(loc='best')
    plt.title('Unique on 100,000 values, int64')
    plt.xlabel('Number of unique labels')
    plt.ylabel('Mean execution time')

    plt.show()
