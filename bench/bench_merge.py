from pandas import *
import random

N = 10000
ngroups = 3

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

df2 = DataFrame({'key1'  : [0, 1, 2, 0, 1, 2],
                 'key2'  : [0, 1, 2, 0, 1, 2],
                 'value' : list('abcdef')})


import pandas.tools.merge as merge
reload(merge)

left, right = merge._get_group_keys([df['key1'], df['key2']],
                                        [df2['key1'], df2['key2']])

left, right = merge._get_group_keys([df['key1']], [df2['key1']])

