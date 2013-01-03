from pandas import *
from pandas.util.testing import rands

n = 1000000
# indices = Index([rands(10) for _ in xrange(n)])


def sample(values, k):
    sampler = np.random.permutation(len(values))
    return values.take(sampler[:k])
sz = 500000
rng = np.arange(0, 10000000000000, 10000000)
stamps = np.datetime64(datetime.now()).view('i8') + rng
idx1 = np.sort(sample(stamps, sz))
idx2 = np.sort(sample(stamps, sz))
ts1 = Series(np.random.randn(sz), idx1)
ts2 = Series(np.random.randn(sz), idx2)


# subsample_size = 90000

# x = Series(np.random.randn(100000), indices)
# y = Series(np.random.randn(subsample_size),
#            index=sample(indices, subsample_size))


# lx = larry(np.random.randn(100000), [list(indices)])
# ly = larry(np.random.randn(subsample_size), [list(y.index)])

# Benchmark 1: Two 1-million length time series (int64-based index) with
# randomly chosen timestamps

# Benchmark 2: Join two 5-variate time series DataFrames (outer and inner join)

# df1 = DataFrame(np.random.randn(1000000, 5), idx1, columns=range(5))
# df2 = DataFrame(np.random.randn(1000000, 5), idx2, columns=range(5, 10))
