from pandas import *
from pandas.compat import range
try:
    import pandas.core.internals as internals
    reload(internals)
    import pandas.core.frame as frame
    reload(frame)
    from pandas.core.frame import DataFrame as DataMatrix
except ImportError:
    pass

N = 1000
K = 500


def horribly_unconsolidated():
    index = np.arange(N)

    df = DataMatrix(index=index)

    for i in range(K):
        df[i] = float(K)

    return df


def bench_reindex_index(df, it=100):
    new_idx = np.arange(0, N, 2)
    for i in range(it):
        df.reindex(new_idx)


def bench_reindex_columns(df, it=100):
    new_cols = np.arange(0, K, 2)
    for i in range(it):
        df.reindex(columns=new_cols)


def bench_join_index(df, it=10):
    left = df.reindex(index=np.arange(0, N, 2),
                      columns=np.arange(K // 2))
    right = df.reindex(columns=np.arange(K // 2 + 1, K))
    for i in range(it):
        joined = left.join(right)

if __name__ == '__main__':
    df = horribly_unconsolidated()
    left = df.reindex(index=np.arange(0, N, 2),
                      columns=np.arange(K // 2))
    right = df.reindex(columns=np.arange(K // 2 + 1, K))
    bench_join_index(df)
