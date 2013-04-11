""" self-contained to write legacy pickle files """

def _create_sp_series():

    import numpy as np
    from pandas import bdate_range, SparseSeries

    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=float)
    index = np.arange(15)
    arr[7:12] = nan
    arr[-1:] = nan

    date_index = bdate_range('1/1/2011', periods=len(index))
    bseries = SparseSeries(arr, index=index, kind='block')
    bseries.name = 'bseries'
    return bseries

def _create_sp_frame():
    import numpy as np
    from pandas import bdate_range, SparseDataFrame

    nan = np.nan

    data = {'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
            'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
            'C': np.arange(10),
            'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}
    
    dates = bdate_range('1/1/2011', periods=10)
    return SparseDataFrame(data, index=dates)

def create_data():
    """ create the pickle data """
    
    import numpy as np
    import pandas
    from pandas import (Series,DataFrame,Panel,
                        SparseSeries,SparseDataFrame,SparsePanel,
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
    
    index  = dict(int   = Index(np.arange(10)),
                  date  = date_range('20130101',periods=10))
    mi     = dict(reg   = MultiIndex.from_tuples(zip([['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                                                      ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]),
                                                 names=['first', 'second']))
    series = dict(float = Series(data['A']),
                  int   = Series(data['B']),
                  mixed = Series(data['E']))
    frame  = dict(float = DataFrame(dict(A = series['float'], B = series['float'] + 1)),
                  int   = DataFrame(dict(A = series['int']  , B = series['int']   + 1)),
                  mixed = DataFrame(dict([ (k,data[k]) for k in ['A','B','C','D']])))
    panel  = dict(float = Panel(dict(ItemA = frame['float'], ItemB = frame['float']+1)))

 

    return dict( series = series, 
                 frame  = frame, 
                 panel  = panel,
                 index  = index,
                 mi     = mi,
                 sp_series = dict(float = _create_sp_series()),
                 sp_frame  = dict(float = _create_sp_frame())
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
    import cPickle as pickle

    print("This script generates a pickle file for the current arch, system, and python version")

    base_dir, _ = os.path.split(os.path.abspath(__file__))
    base_dir = os.path.join(base_dir,'data/legacy_pickle')
    
    # could make this a parameter?
    version  = None


    if version is None:
        version = pandas.__version__
    pth = os.path.join(base_dir, str(version))
    try:
        os.mkdir(pth)
    except:
        pass

    # construct a reasonable platform name
    f = '_'.join([ str(pl.machine()), str(pl.system().lower()), str(pl.python_version()) ])
    pth = os.path.abspath(os.path.join(pth,'%s.pickle' % f))
    
    fh = open(pth,'wb')
    pickle.dump(create_data(),fh,pickle.HIGHEST_PROTOCOL)
    fh.close()
    
    print("created pickle file: %s" % pth)

if __name__ == '__main__':
    write_legacy_pickles()
