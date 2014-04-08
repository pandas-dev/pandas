from pandas import *

K = 100
N = 100000
rng = DatetimeIndex('1/1/2000', periods=N, offset=datetools.Minute())

rng2 = np.asarray(rng).astype('M8[us]').astype('i8')

series = {}
for i in range(1, K + 1):
    data = np.random.randn(N)[:-i]
    this_rng = rng2[:-i]
    data[100:] = np.nan
    series[i] = SparseSeries(data, index=this_rng)
