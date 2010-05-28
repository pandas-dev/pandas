import time, os
import numpy as np

import la
import pandas

def timeit(f, iterations):
    start = time.clock()

    for i in xrange(iterations):
        f()

    return time.clock() - start

def rountrip_archive(N, iterations=10):

    # Create data
    arr = np.random.randn(N, N)
    lar = la.larry(arr)
    dma = pandas.DataMatrix(arr, range(N), range(N))

    # filenames
    filename_numpy = 'c:/temp/numpy.npz'
    filename_larry = 'c:/temp/archive.hdf5'
    filename_pandas = 'c:/temp/pandas_tmp'

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

    print 'Numpy (npz)   %7.4f seconds' % numpy_time
    print 'larry (HDF5)  %7.4f seconds' % larry_time
    print 'pandas (HDF5) %7.4f seconds' % pandas_time

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
    dma1 = pandas.DataMatrix.load(filename)
    dma2.save(filename)
    dma2 = pandas.DataMatrix.load(filename)


    In [65]: df1
    Out[65]:
                           A              B
    2000-01-03 00:00:00    -0.1174        -0.941
    2000-01-04 00:00:00    -0.6034        -0.008094
    2000-01-05 00:00:00    -0.3816        -0.9338
    2000-01-06 00:00:00    -0.3298        -0.9548
    2000-01-07 00:00:00    0.9576         0.4652
    2000-01-10 00:00:00    -0.7208        -1.131
    2000-01-11 00:00:00    1.568          0.8498
    2000-01-12 00:00:00    0.3717         -0.2323
    2000-01-13 00:00:00    -1.428         -1.997
    2000-01-14 00:00:00    -1.084         -0.271


    In [66]: df1.join?
    Type:           instancemethod
    Base Class:     <type 'instancemethod'>
    <bound method DataFrame.join of                        A              B
                2000-01-03  <...> 0:00:00    -1.428         -1.997
                2000-01-14 00:00:00    -1.084         -0.271
                >
    Namespace:      Interactive
    File:           h:\workspace\pandas\pandas\core\frame.py
    Definition:     df1.join(self, other, on=None, how=None)
    Docstring:
        Join columns with other DataFrame either on index or on a key
    column

    Parameters
    ----------
    other : DataFrame
        Index should be similar to one of the columns in this one
    on : string, default None
        Column name to use, otherwise join on index
    how : {'left', 'right', 'outer', 'inner'}
        default: 'left' for joining on index, None otherwise
        How to handle indexes of the two objects.
          * left: use calling frame's index
          * right: use input frame's index
          * outer: form union of indexes
          * inner: use intersection of indexes


    In [67]: df2
    Out[67]:
                           C              D
    2000-01-03 00:00:00    0.2833         -0.1937
    2000-01-05 00:00:00    1.868          1.207
    2000-01-07 00:00:00    -0.8586        -0.7367
    2000-01-11 00:00:00    2.121          0.9104
    2000-01-13 00:00:00    0.7856         0.9063


    In [68]: df1.join(df2)
    Out[68]:
                           A              B              C              D
    2000-01-03 00:00:00    -0.1174        -0.941         0.2833         -0.1937
    2000-01-04 00:00:00    -0.6034        -0.008094      NaN            NaN
    2000-01-05 00:00:00    -0.3816        -0.9338        1.868          1.207
    2000-01-06 00:00:00    -0.3298        -0.9548        NaN            NaN
    2000-01-07 00:00:00    0.9576         0.4652         -0.8586        -0.7367
    2000-01-10 00:00:00    -0.7208        -1.131         NaN            NaN
    2000-01-11 00:00:00    1.568          0.8498         2.121          0.9104
    2000-01-12 00:00:00    0.3717         -0.2323        NaN            NaN
    2000-01-13 00:00:00    -1.428         -1.997         0.7856         0.9063
    2000-01-14 00:00:00    -1.084         -0.271         NaN            NaN

    In [70]: df1.join(df2, how='inner')
    Out[70]:
                           A              B              C              D
    2000-01-03 00:00:00    -0.1174        -0.941         0.2833         -0.1937
    2000-01-05 00:00:00    -0.3816        -0.9338        1.868          1.207
    2000-01-07 00:00:00    0.9576         0.4652         -0.8586        -0.7367
    2000-01-11 00:00:00    1.568          0.8498         2.121          0.9104
    2000-01-13 00:00:00    -1.428         -1.997         0.7856         0.9063

    In [73]: df2
    Out[73]:
                           C              D              key
    2000-01-03 00:00:00    0.2833         -0.1937        0
    2000-01-05 00:00:00    1.868          1.207          1
    2000-01-07 00:00:00    -0.8586        -0.7367        0
    2000-01-11 00:00:00    2.121          0.9104         1
    2000-01-13 00:00:00    0.7856         0.9063         0


    In [74]: df3 = DataFrame({'code' : {0 : 'foo', 1 : 'bar'}})

    In [75]: df3
    Out[75]:
         code
    0    foo
    1    bar


    In [76]: df2.join(df3, on='key')
    Out[76]:
                           C              D              code           key
    2000-01-03 00:00:00    0.2833         -0.1937        foo            0
    2000-01-05 00:00:00    1.868          1.207          bar            1
    2000-01-07 00:00:00    -0.8586        -0.7367        foo            0
    2000-01-11 00:00:00    2.121          0.9104         bar            1
    2000-01-13 00:00:00    0.7856         0.9063         foo            0
