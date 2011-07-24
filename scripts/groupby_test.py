from collections import defaultdict

import numpy as np

from pandas import *

import pandas._tseries as tseries
import pandas.core.groupby as gp
reload(gp)

k = 1000
values = np.random.randn(8 * k)
key1 = np.array(['foo', 'bar', 'baz', 'bar', 'foo', 'baz', 'bar', 'baz'] * k,
                dtype=object)
key2 = np.array(['b', 'b', 'b', 'b', 'a', 'a', 'a', 'a' ] * k,
                dtype=object)
shape, labels, idicts = gp.labelize(key1, key2)

print tseries.group_labels(key1)

# print shape
# print labels
# print idicts

result = tseries.group_aggregate(values, labels, shape)

print tseries.groupby_indices(key2)

df = DataFrame({'key1' : key1,
                'key2' : key2,
                'values' : values})
k1 = df['key1']
k2 = df['key2']

# del df['key1']
# del df['key2']

# r2 = gp.multi_groupby(df, np.sum, k1, k2)

# print result

gen = gp.generate_groups(df['values'], labels, shape, axis=1,
                         factory=DataFrame)

res = defaultdict(dict)
for a, gen1 in gen:
    for b, group in gen1:
        print a, b
        print group
        # res[b][a] = group['values'].sum()
        res[b][a] = group.sum()

res = DataFrame(res)

