""" self-contained to write legacy pickle files """
from __future__ import print_function

# make sure we are < 0.13 compat (in py3)
try:
    from pandas.compat import zip, cPickle as pickle
except:
    import pickle

def _create_sp_series():

    import numpy as np
    from pandas import SparseSeries

    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    index = np.arange(15)
    arr[7:12] = nan
    arr[-1:] = nan

    bseries = SparseSeries(arr, kind='block')
    bseries.name = 'bseries'
    return bseries

def _create_sp_tsseries():

    import numpy as np
    from pandas import bdate_range, SparseTimeSeries

    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    index = np.arange(15)
    arr[7:12] = nan
    arr[-1:] = nan

    date_index = bdate_range('1/1/2011', periods=len(index))
    bseries = SparseTimeSeries(arr, index=date_index, kind='block')
    bseries.name = 'btsseries'
    return bseries

def _create_sp_frame():
    import numpy as np
    from pandas import bdate_range, SparseDataFrame

    nan = np.nan

    data = {'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
            'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
            'C': np.arange(10).astype(np.int64),
            'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}

    dates = bdate_range('1/1/2011', periods=10)
    return SparseDataFrame(data, index=dates)

def create_data():
    """ create the pickle data """

    import numpy as np
    import pandas
    from pandas import (Series,TimeSeries,DataFrame,Panel,
                        SparseSeries,SparseTimeSeries,SparseDataFrame,SparsePanel,
                        Index,MultiIndex,PeriodIndex,
                        date_range,bdate_range,Timestamp)
    nan = np.nan

    data = {
        'A': [0., 1., 2., 3., np.nan],
        'B': [0, 1, 0, 1, 0],
        'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D': date_range('1/1/2009', periods=5),
        'E' : [0., 1, Timestamp('20100101'),'foo',2.],
        }

    index = dict(int = Index(np.arange(10)),
                  date = date_range('20130101',periods=10))
    mi = dict(reg = MultiIndex.from_tuples(list(zip([['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                                                      ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']])),
                                                 names=['first', 'second']))
    series = dict(float = Series(data['A']),
                  int = Series(data['B']),
                  mixed = Series(data['E']),
                  ts = TimeSeries(np.arange(10).astype(np.int64),index=date_range('20130101',periods=10)))
    frame = dict(float = DataFrame(dict(A = series['float'], B = series['float'] + 1)),
                 int = DataFrame(dict(A = series['int']  , B = series['int']   + 1)),
                 mixed = DataFrame(dict([ (k,data[k]) for k in ['A','B','C','D']])))
    panel = dict(float = Panel(dict(ItemA = frame['float'], ItemB = frame['float']+1)))



    return dict( series = series,
                 frame = frame,
                 panel = panel,
                 index = index,
                 mi = mi,
                 sp_series = dict(float = _create_sp_series(),
                                  ts = _create_sp_tsseries()),
                 sp_frame = dict(float = _create_sp_frame())
                 )

def write_legacy_pickles():

    # force our cwd to be the first searched
    import sys
    sys.path.insert(0,'.')

    import os
    import numpy as np
    import pandas
    import pandas.util.testing as tm
    import platform as pl

    print("This script generates a pickle file for the current arch, system, and python version")

    version = pandas.__version__

    # construct a reasonable platform name
    f = '_'.join([ str(version), str(pl.machine()), str(pl.system().lower()), str(pl.python_version()) ])
    pth = '{0}.pickle'.format(f)

    fh = open(pth,'wb')
    pickle.dump(create_data(),fh,pickle.HIGHEST_PROTOCOL)
    fh.close()

    print("created pickle file: %s" % pth)

if __name__ == '__main__':
    write_legacy_pickles()
