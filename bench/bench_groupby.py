from pandas import *

import string
import random

k = 200
n = 1000


foo = np.tile(list(range(k)), n)
foo2 = list(foo)
random.shuffle(foo)
random.shuffle(foo2)

df = DataFrame({'A' : foo,
                'B' : foo2,
                'C' : np.random.randn(n * k)})


import pandas._tseries as lib

f = np.std


grouped = df.groupby(['A', 'B'])

label_list = [ping.labels for ping in grouped.groupings]
shape = [len(ping.ids) for ping in grouped.groupings]

from pandas.core.groupby import get_group_index


group_index = get_group_index(label_list, shape).astype('i4')

ngroups = np.prod(shape)

indexer = lib.groupsort_indexer(group_index, ngroups)

values = df['C'].values.take(indexer)
group_index = group_index.take(indexer)

f = lambda x: x.std(ddof=1)

grouper = lib.Grouper(df['C'], np.ndarray.std, group_index, ngroups)
result = grouper.get_result()

expected = grouped.std()
