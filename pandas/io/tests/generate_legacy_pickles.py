""" self-contained to write legacy pickle files """

def create_data():
    """ create the pickle data """
    
    import numpy as np
    import pandas
    from pandas import (Series,DataFrame,Panel,
                        SparseSeries,SparseDataFrame,SparsePanel,
                        Index,MultiIndex,PeriodIndex,
                        date_range,Timestamp)

    data = {
        'A': [0., 1., 2., 3., np.nan],
        'B': [0, 1, 0, 1, 0],
        'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D': date_range('1/1/2009', periods=5),
        'E' : [0., 1, Timestamp('20100101'),'foo',2.],
        }
    
    series = dict(float = Series(data['A']),
                  int   = Series(data['B']),
                  mixed = Series(data['E']))
    frame  = dict(float = DataFrame(dict(A = series['float'], B = series['float'] + 1)),
                  int   = DataFrame(dict(A = series['int']  , B = series['int']   + 1)),
                  mixed = DataFrame(dict([ (k,data[k]) for k in ['A','B','C','D']])))
    panel  = dict(float = Panel(dict(ItemA = frame['float'], ItemB = frame['float']+1)))
    
    return dict( series = series, 
                 frame  = frame, 
                 panel  = panel )

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

    print "This script generates a pickle file for the current arch, system, and python version"

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
    
    print "created pickle file: %s" % pth

if __name__ == '__main__':
    write_legacy_pickles()
