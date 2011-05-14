import numpy as np

from pandas import *
import pandas.core.sparse as spm
reload(spm)
from pandas.core.sparse import *

N = 10000.

arr1 = np.arange(N)
index = Index(np.arange(N))

off = N//10
arr1[off : 2 * off] = np.NaN
arr1[4*off: 5 * off] = np.NaN
arr1[8*off: 9 * off] = np.NaN

arr2 = np.arange(N)
arr2[3 * off // 2: 2 * off  + off // 2] = np.NaN
arr2[8 * off + off // 2: 9 * off + off // 2] = np.NaN

s1 = SparseSeries(arr1, index=index)
s2 = SparseSeries(arr2, index=index)

is1 = SparseSeries(arr1, kind='integer', index=index)
is2 = SparseSeries(arr2, kind='integer', index=index)

s1_dense = s1.to_dense()
s2_dense = s2.to_dense()

dm = DataMatrix.load('/Users/wesm/code/pandas/example')

data = {}
for col, ser in dm.iteritems():
    data[col] = SparseSeries(ser)

dwp = WidePanel.fromDict({'foo' : dm})

sdf = SparseDataFrame(data)

lp = stack_sparse_frame(sdf)

swp = SparseWidePanel({'A' : sdf,
                       'B' : sdf})
