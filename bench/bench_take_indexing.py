from __future__ import print_function
import numpy as np

from pandas import *
import pandas._tseries as lib

from pandas import DataFrame
import timeit
from pandas.compat import zip

setup = """
from pandas import Series
import pandas._tseries as lib
import random
import numpy as np

import random
n = %d
k = %d
arr = np.random.randn(n, k)
indexer = np.arange(n, dtype=np.int32)
indexer = indexer[::-1]
"""

sizes = [100, 1000, 10000, 100000]
iters = [1000, 1000, 100, 1]

fancy_2d = []
take_2d = []
cython_2d = []

n = 1000


def _timeit(stmt, size, k=5, iters=1000):
    timer = timeit.Timer(stmt=stmt, setup=setup % (sz, k))
    return timer.timeit(n) / n

for sz, its in zip(sizes, iters):
    print(sz)
    fancy_2d.append(_timeit('arr[indexer]', sz, iters=its))
    take_2d.append(_timeit('arr.take(indexer, axis=0)', sz, iters=its))
    cython_2d.append(_timeit('lib.take_axis0(arr, indexer)', sz, iters=its))

df = DataFrame({'fancy': fancy_2d,
                'take': take_2d,
                'cython': cython_2d})

print(df)

from pandas.rpy.common import r
r('mat <- matrix(rnorm(50000), nrow=10000, ncol=5)')
r('set.seed(12345')
r('indexer <- sample(1:10000)')
r('mat[indexer,]')
