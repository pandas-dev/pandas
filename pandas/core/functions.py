from pandas.core.common import isnull
import numpy as np

#-------------------------------------------------------------------------------
# NaN-friendly reductions and such

def reduce_mean(values, index, buckets, inclusive=False):
    def _reduceat_mean(values, mask, locs):
        the_sum = np.add.reduceat(values, locs)
        the_count = np.add.reduceat(-mask, locs)
        return the_sum / the_count
    return _reduce_generic(values, index, buckets, _reduceat_mean,
                           inclusive=inclusive, na_fill=0)


def _reduceat_var(values, mask, locs):
    XX = np.add.reduceat(values ** 2, locs)
    X = np.add.reduceat(values, locs)
    nobs = np.add.reduceat(-mask, locs)
    return (XX - X * X) / (nobs - 1)

def reduce_std(values, index, buckets, inclusive=False):
    result = _reduce_generic(values, index, buckets, _reduceat_var,
                             inclusive=inclusive, na_fill=0)
    return np.sqrt(result)

def reduce_prod(values, index, buckets, inclusive=False):
    def _reduceat_prod(values, mask, locs):
        return np.multiply.reduceat(values, locs)
    return _reduce_generic(values, index, buckets, _reduceat_prod,
                           inclusive=inclusive, na_fill=1)

def reduce_min(values, index, buckets, inclusive=False):
    def _reduceat_min(values, mask, locs):
        return np.minimum.reduceat(values, locs)
    return _reduce_generic(values, index, buckets, _reduceat_min,
                           inclusive=inclusive, na_fill=np.inf)

def reduce_max(values, index, buckets, inclusive=False):
    def _reduceat_max(values, mask, locs):
        return np.maximum.reduceat(values, locs)
    return _reduce_generic(values, index, buckets, _reduceat_max,
                           inclusive=inclusive, na_fill=-np.inf)

def _reduce_generic(values, index, buckets, freduce, inclusive=False,
                    na_fill=None):
    """

    """
    locs = _bucket_locs(index, buckets, inclusive=inclusive)

    values = np.asarray(values)
    mask = isnull(values)

    if na_fill is not None:
        values = values.copy()
        np.putmask(values, mask, na_fill)

    return freduce(values, mask, locs)

def _reduceat_count(values, mask, locs):
    return np.add.reduceat(-mask, locs)

def _bucket_locs(index, buckets, inclusive=False):
    if inclusive:
        locs = index.searchsorted(buckets, side='left')
    else:
        locs = index.searchsorted(buckets, side='right')

    return locs

def get_bucket(date, bucks):
    if date in bucks:
        idx = bucks.indexMap[date] + 1
    else:
        idx = bucks.searchsorted(date)
    return bucks[idx]

def dumb_way(series, buckets):
    sampled2 = hfseries.groupby(lambda x: get_bucket(x, buckets)).mean()
    sampled2 = sampled2.reindex(buckets)
    return sampled2

def ts_upsample(dates, buckets, values, aggfunc, inclusive=True):
    '''
    put something here
    '''
    nbuckets = len(buckets)
    nvalues = len(dates)
    output = np.empty(nbuckets, dtype=float)

    if inclusive:
        _check = lambda x, y: x < y
    else:
        _check = lambda x, y: x <= y

    j = 0
    for i, bound in enumerate(buckets):
        next_bound = buckets[i + 1]
        jstart = j

        while _check(dates[j], next_bound) and j < nvalues:
            j += 1

        output[i] = aggfunc(values[jstart:j])

    return Series(output, index=buckets)

if __name__ == '__main__':
    N = 1000000
    K = 1000

    values = np.random.randn(N)
    index = np.arange(N).astype(object)
    buckets = np.arange(0, N, N // K).astype(object)

    result = reduce_mean(values, index, buckets)

    import pandas._tseries as tseries
    tseries.ts_upsample_mean(index, buckets, values)
