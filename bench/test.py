import numpy as np
import itertools
import collections
import scipy.ndimage as ndi
from pandas.compat import zip, range

N = 10000

lat = np.random.randint(0, 360, N)
lon = np.random.randint(0, 360, N)
data = np.random.randn(N)


def groupby1(lat, lon, data):
    indexer = np.lexsort((lon, lat))
    lat = lat.take(indexer)
    lon = lon.take(indexer)
    sorted_data = data.take(indexer)

    keys = 1000. * lat + lon
    unique_keys = np.unique(keys)
    bounds = keys.searchsorted(unique_keys)

    result = group_agg(sorted_data, bounds, lambda x: x.mean())

    decoder = keys.searchsorted(unique_keys)

    return dict(zip(zip(lat.take(decoder), lon.take(decoder)), result))


def group_mean(lat, lon, data):
    indexer = np.lexsort((lon, lat))
    lat = lat.take(indexer)
    lon = lon.take(indexer)
    sorted_data = data.take(indexer)

    keys = 1000 * lat + lon
    unique_keys = np.unique(keys)

    result = ndi.mean(sorted_data, labels=keys, index=unique_keys)
    decoder = keys.searchsorted(unique_keys)

    return dict(zip(zip(lat.take(decoder), lon.take(decoder)), result))


def group_mean_naive(lat, lon, data):
    grouped = collections.defaultdict(list)
    for lt, ln, da in zip(lat, lon, data):
        grouped[(lt, ln)].append(da)

    averaged = dict((ltln, np.mean(da)) for ltln, da in grouped.items())

    return averaged


def group_agg(values, bounds, f):
    N = len(values)
    result = np.empty(len(bounds), dtype=float)
    for i, left_bound in enumerate(bounds):
        if i == len(bounds) - 1:
            right_bound = N
        else:
            right_bound = bounds[i + 1]

        result[i] = f(values[left_bound: right_bound])

    return result

# for i in range(10):
#     groupby1(lat, lon, data)
