from __future__ import print_function
from pandas.compat import range, lrange
import time
import os
import numpy as np

import la
import pandas


def timeit(f, iterations):
    start = time.clock()

    for i in range(iterations):
        f()

    return time.clock() - start


def roundtrip_archive(N, iterations=10):

    # Create data
    arr = np.random.randn(N, N)
    lar = la.larry(arr)
    dma = pandas.DataFrame(arr, lrange(N), lrange(N))

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
    numpy_f = lambda: numpy_roundtrip(filename_numpy, arr, arr)
    numpy_time = timeit(numpy_f, iterations) / iterations

    larry_f = lambda: larry_roundtrip(filename_larry, lar, lar)
    larry_time = timeit(larry_f, iterations) / iterations

    pandas_f = lambda: pandas_roundtrip(filename_pandas, dma, dma)
    pandas_time = timeit(pandas_f, iterations) / iterations

    print('Numpy (npz)   %7.4f seconds' % numpy_time)
    print('larry (HDF5)  %7.4f seconds' % larry_time)
    print('pandas (HDF5) %7.4f seconds' % pandas_time)


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
