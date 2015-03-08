from pandas import *
from pandas.util.testing import rands
from pandas.compat import range

import string
import random

k = 20000
n = 10

foo = np.tile(np.array([rands(10) for _ in range(k)], dtype='O'), n)
foo2 = list(foo)
random.shuffle(foo)
random.shuffle(foo2)

df = DataFrame({'A': foo,
                'B': foo2,
                'C': np.random.randn(n * k)})

import pandas._sandbox as sbx


def f():
    table = sbx.StringHashTable(len(df))
    ret = table.factorize(df['A'])
    return ret


def g():
    table = sbx.PyObjectHashTable(len(df))
    ret = table.factorize(df['A'])
    return ret

ret = f()

"""
import pandas._tseries as lib

f = np.std


grouped = df.groupby(['A', 'B'])

label_list = [ping.labels for ping in grouped.groupings]
shape = [len(ping.ids) for ping in grouped.groupings]

from pandas.core.groupby import get_group_index


group_index = get_group_index(label_list, shape,
                              sort=True, xnull=True).astype('i4')

ngroups = np.prod(shape)

indexer = lib.groupsort_indexer(group_index, ngroups)

values = df['C'].values.take(indexer)
group_index = group_index.take(indexer)

f = lambda x: x.std(ddof=1)

grouper = lib.Grouper(df['C'], np.ndarray.std, group_index, ngroups)
result = grouper.get_result()

expected = grouped.std()
"""
