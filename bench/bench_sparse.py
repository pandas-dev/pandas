import sys
import numpy as np

from pandas import *
import pandas.core.sparse as spm
import pandas.compat as compat
reload(spm)
from pandas.core.sparse import *

N = 10000.

arr1 = np.arange(N)
index = Index(np.arange(N))

off = N // 10
arr1[off: 2 * off] = np.NaN
arr1[4 * off: 5 * off] = np.NaN
arr1[8 * off: 9 * off] = np.NaN

arr2 = np.arange(N)
arr2[3 * off // 2: 2 * off + off // 2] = np.NaN
arr2[8 * off + off // 2: 9 * off + off // 2] = np.NaN

s1 = SparseSeries(arr1, index=index)
s2 = SparseSeries(arr2, index=index)

is1 = SparseSeries(arr1, kind='integer', index=index)
is2 = SparseSeries(arr2, kind='integer', index=index)

s1_dense = s1.to_dense()
s2_dense = s2.to_dense()

if 'linux' in sys.platform:
    pth = '/home/wesm/code/pandas/example'
else:
    pth = '/Users/wesm/code/pandas/example'

dm = DataFrame.load(pth)

sdf = dm.to_sparse()


def new_data_like(sdf):
    new_data = {}
    for col, series in compat.iteritems(sdf):
        new_data[col] = SparseSeries(np.random.randn(len(series.sp_values)),
                                     index=sdf.index,
                                     sparse_index=series.sp_index,
                                     fill_value=series.fill_value)

    return SparseDataFrame(new_data)

# data = {}
# for col, ser in dm.iteritems():
#     data[col] = SparseSeries(ser)

dwp = Panel.fromDict({'foo': dm})
# sdf = SparseDataFrame(data)


lp = stack_sparse_frame(sdf)


swp = SparsePanel({'A': sdf})
swp = SparsePanel({'A': sdf,
                   'B': sdf,
                   'C': sdf,
                   'D': sdf})

y = sdf
x = SparsePanel({'x1': sdf + new_data_like(sdf) / 10,
                 'x2': sdf + new_data_like(sdf) / 10})

dense_y = sdf
dense_x = x.to_dense()

# import hotshot, hotshot.stats
# prof = hotshot.Profile('test.prof')

# benchtime, stones = prof.runcall(ols, y=y, x=x)

# prof.close()

# stats = hotshot.stats.load('test.prof')

dense_model = ols(y=dense_y, x=dense_x)

import pandas.stats.plm as plm
import pandas.stats.interface as face
reload(plm)
reload(face)

# model = face.ols(y=y, x=x)
