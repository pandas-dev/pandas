from __future__ import print_function
from pandas import DataFrame
from pandas.compat import range, zip
import timeit

setup = """
from pandas import Series
import pandas._tseries as _tseries
from pandas.compat import range
import random
import numpy as np

def better_unique(values):
    uniques = _tseries.fast_unique(values)
    id_map = _tseries.map_indices_buf(uniques)
    labels = _tseries.get_unique_labels(values, id_map)
    return uniques, labels

tot = 100000

def get_test_data(ngroups=100, n=tot):
    unique_groups = range(ngroups)
    random.shuffle(unique_groups)
    arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
                         dtype=object)

    return arr

arr = get_test_data(ngroups=%d)
"""

group_sizes = [10, 100, 1000, 10000,
               20000, 30000, 40000,
               50000, 60000, 70000,
               80000, 90000, 100000]

numbers = [100, 100, 50] + [10] * 10

numpy = []
wes = []

for sz, n in zip(group_sizes, numbers):
    # wes_timer =  timeit.Timer(stmt='better_unique(arr)',
    #                           setup=setup % sz)
    wes_timer = timeit.Timer(stmt='_tseries.fast_unique(arr)',
                             setup=setup % sz)

    numpy_timer = timeit.Timer(stmt='np.unique(arr)',
                               setup=setup % sz)

    print(n)
    numpy_result = numpy_timer.timeit(number=n) / n
    wes_result = wes_timer.timeit(number=n) / n

    print('Groups: %d, NumPy: %s, Wes: %s' % (sz, numpy_result, wes_result))

    wes.append(wes_result)
    numpy.append(numpy_result)

result = DataFrame({'wes': wes, 'numpy': numpy}, index=group_sizes)


def make_plot(numpy, wes):
    pass

# def get_test_data(ngroups=100, n=100000):
#     unique_groups = range(ngroups)
#     random.shuffle(unique_groups)
#     arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

#     if len(arr) < n:
#         arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
#                          dtype=object)

#     return arr

# arr = get_test_data(ngroups=1000)
