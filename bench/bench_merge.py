from pandas import *
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
