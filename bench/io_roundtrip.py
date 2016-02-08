from __future__ import print_function
import time
import os
import numpy as np

import la
import pandas
from pandas.compat import range
from pandas import datetools, DatetimeIndex


def timeit(f, iterations):
    start = time.clock()

    for i in range(iterations):
        f()

    return time.clock() - start


def rountrip_archive(N, K=50, iterations=10):
    # Create data
    arr = np.random.randn(N, K)
    # lar = la.larry(arr)
    dma = pandas.DataFrame(arr,
                           DatetimeIndex('1/1/2000', periods=N,
                                     offset=datetools.Minute()))
    dma[201] = 'bar'

    # filenames
    filename_numpy = '/Users/wesm/tmp/numpy.npz'
    filename_larry = '/Users/wesm/tmp/archive.hdf5'
    filename_pandas = '/Users/wesm/tmp/pandas_tmp'

    # Delete old files
    try:
        os.unlink(filename_numpy)
    except:
        pass
    try:
        os.unlink(filename_larry)
    except:
        pass

    try:
        os.unlink(filename_pandas)
    except:
        pass

    # Time a round trip save and load
    # numpy_f = lambda: numpy_roundtrip(filename_numpy, arr, arr)
    # numpy_time = timeit(numpy_f, iterations) / iterations

    # larry_f = lambda: larry_roundtrip(filename_larry, lar, lar)
    # larry_time = timeit(larry_f, iterations) / iterations

    pandas_f = lambda: pandas_roundtrip(filename_pandas, dma, dma)
    pandas_time = timeit(pandas_f, iterations) / iterations
    print('pandas (HDF5) %7.4f seconds' % pandas_time)

    pickle_f = lambda: pandas_roundtrip(filename_pandas, dma, dma)
    pickle_time = timeit(pickle_f, iterations) / iterations
    print('pandas (pickle) %7.4f seconds' % pickle_time)

    # print('Numpy (npz)   %7.4f seconds' % numpy_time)
    # print('larry (HDF5)  %7.4f seconds' % larry_time)

    # Delete old files
    try:
        os.unlink(filename_numpy)
    except:
        pass
    try:
        os.unlink(filename_larry)
    except:
        pass

    try:
        os.unlink(filename_pandas)
    except:
        pass


def numpy_roundtrip(filename, arr1, arr2):
    np.savez(filename, arr1=arr1, arr2=arr2)
    npz = np.load(filename)
    arr1 = npz['arr1']
    arr2 = npz['arr2']


def larry_roundtrip(filename, lar1, lar2):
    io = la.IO(filename)
    io['lar1'] = lar1
    io['lar2'] = lar2
    lar1 = io['lar1']
    lar2 = io['lar2']


def pandas_roundtrip(filename, dma1, dma2):
    # What's the best way to code this?
    from pandas.io.pytables import HDFStore
    store = HDFStore(filename)
    store['dma1'] = dma1
    store['dma2'] = dma2
    dma1 = store['dma1']
    dma2 = store['dma2']


def pandas_roundtrip_pickle(filename, dma1, dma2):
    dma1.save(filename)
    dma1 = pandas.DataFrame.load(filename)
    dma2.save(filename)
    dma2 = pandas.DataFrame.load(filename)

if __name__ == '__main__':
    rountrip_archive(10000, K=200)
