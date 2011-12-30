from pandas import *
from pandas.util.testing import rands
import random

N = 10000
ngroups = 10

def get_test_data(ngroups=100, n=N):
    unique_groups = range(ngroups)
    arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
                         dtype=object)

    random.shuffle(arr)
    return arr

# aggregate multiple columns
df = DataFrame({'key1' : get_test_data(ngroups=ngroups),
                'key2' : get_test_data(ngroups=ngroups),
                'data1' : np.random.randn(N),
                'data2' : np.random.randn(N)})

df2 = DataFrame({'key1'  : get_test_data(ngroups=ngroups, n=N//10),
                 'key2'  : get_test_data(ngroups=ngroups//2, n=N//10),
                 'value' : np.random.randn(N // 10)})


import pandas.tools.merge as merge
reload(merge)

result = merge.merge(df, df2, on='key2')

from pandas.util.testing import rands
N = 10000
indices = np.array([rands(10) for _ in xrange(N)], dtype='O')

key = np.tile(indices, 10)
key2 = key.copy()
random.shuffle(key2)
indices2 = indices.copy()
random.shuffle(indices2)


left = DataFrame({'key' : key, 'key2':key2,
                  'value' : np.random.randn(100000)})
right = DataFrame({'key': indices, 'key2':indices2,
                   'value2' : np.random.randn(10000)})
