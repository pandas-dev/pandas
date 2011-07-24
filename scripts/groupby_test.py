import numpy as np

from pandas import *

import pandas._tseries as tseries
import pandas.core.groupby as gp
reload(gp)

k = 10000
values = np.random.randn(8 * k)
key1 = np.array(['foo', 'bar', 'baz', 'bar', 'foo', 'baz', 'bar', 'baz'] * k,
                dtype=object)
key2 = np.array(['b', 'b', 'b', 'b', 'a', 'a', 'a', 'a' ] * k,
                dtype=object)
shape, labels, idicts = tseries.labelize(key1, key2)

print tseries.group_labels(key1)

# print shape
# print labels
# print idicts

result = tseries.group_aggregate(values, labels, shape)

print tseries.groupby_indices(key2)

df = DataFrame({'key1' : key1,
                'key2' : key2,
                'values' : values})
k1 = df.pop('key1')
k2 = df.pop('key2')

r2 = gp.multi_groupby(df, np.sum, k1, k2)

# print result

