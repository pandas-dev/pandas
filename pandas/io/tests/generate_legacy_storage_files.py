""" self-contained to write legacy storage (pickle/msgpack) files """
from __future__ import print_function
from distutils.version import LooseVersion
from pandas import (Series, DataFrame, Panel,
                    SparseSeries, SparseDataFrame, SparsePanel,
                    Index, MultiIndex, PeriodIndex, bdate_range, to_msgpack,
                    date_range, period_range, bdate_range, Timestamp, Categorical,
                    Period)
import os
import sys
import numpy as np
import pandas
import pandas.util.testing as tm
import platform as pl


def _create_sp_series():
    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    arr[7:12] = nan
    arr[-1:] = nan

    bseries = SparseSeries(arr, kind='block')
    bseries.name = 'bseries'
    return bseries


def _create_sp_tsseries():
    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    arr[7:12] = nan
    arr[-1:] = nan

    date_index = bdate_range('1/1/2011', periods=len(arr))
    bseries = SparseSeries(arr, index=date_index, kind='block')
    bseries.name = 'btsseries'
    return bseries


def _create_sp_frame():
    nan = np.nan

    data = {'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
            'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
            'C': np.arange(10).astype(np.int64),
            'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}

    dates = bdate_range('1/1/2011', periods=10)
    return SparseDataFrame(data, index=dates)


def create_data():
    """ create the pickle/msgpack data """

    data = {
        'A': [0., 1., 2., 3., np.nan],
        'B': [0, 1, 0, 1, 0],
        'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D': date_range('1/1/2009', periods=5),
        'E': [0., 1, Timestamp('20100101'), 'foo', 2.]
    }

    scalars = dict(timestamp=Timestamp('20130101'))
    if LooseVersion(pandas.__version__) >= '0.17.0':
        scalars['period'] = Period('2012','M')

    index = dict(int=Index(np.arange(10)),
                 date=date_range('20130101', periods=10),
                 period=period_range('2013-01-01', freq='M', periods=10))

    mi = dict(reg2=MultiIndex.from_tuples(tuple(zip(*[['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                                                      ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']])),
                                          names=['first', 'second']))
    series = dict(float=Series(data['A']),
                  int=Series(data['B']),
                  mixed=Series(data['E']),
                  ts=Series(np.arange(10).astype(np.int64), index=date_range('20130101',periods=10)),
                  mi=Series(np.arange(5).astype(np.float64),
                            index=MultiIndex.from_tuples(tuple(zip(*[[1, 1, 2, 2, 2], [3, 4, 3, 4, 5]])),
                                                         names=['one', 'two'])),
                  dup=Series(np.arange(5).astype(np.float64), index=['A', 'B', 'C', 'D', 'A']),
                  cat=Series(Categorical(['foo', 'bar', 'baz'])),
                  dt=Series(date_range('20130101',periods=5)),
                  dt_tz=Series(date_range('20130101',periods=5,tz='US/Eastern')))
    if LooseVersion(pandas.__version__) >= '0.17.0':
        series['period'] = Series([Period('2000Q1')] * 5)

    mixed_dup_df = DataFrame(data)
    mixed_dup_df.columns = list("ABCDA")
    frame = dict(float=DataFrame(dict(A=series['float'], B=series['float'] + 1)),
                 int=DataFrame(dict(A=series['int'], B=series['int'] + 1)),
                 mixed=DataFrame(dict([(k, data[k]) for k in ['A', 'B', 'C', 'D']])),
                 mi=DataFrame(dict(A=np.arange(5).astype(np.float64), B=np.arange(5).astype(np.int64)),
                              index=MultiIndex.from_tuples(tuple(zip(*[['bar', 'bar', 'baz', 'baz', 'baz'],
                                                                       ['one', 'two', 'one', 'two', 'three']])),
                                                           names=['first', 'second'])),
                 dup=DataFrame(np.arange(15).reshape(5, 3).astype(np.float64),
                               columns=['A', 'B', 'A']),
                 cat_onecol=DataFrame(dict(A=Categorical(['foo', 'bar']))),
                 cat_and_float=DataFrame(dict(A=Categorical(['foo', 'bar', 'baz']),
                                              B=np.arange(3).astype(np.int64))),
                 mixed_dup=mixed_dup_df,
                 dt_mixed_tzs=DataFrame(dict(A=Timestamp('20130102', tz='US/Eastern'), B=Timestamp('20130603', tz='CET')), index=range(5)),
                 )

    mixed_dup_panel = Panel(dict(ItemA=frame['float'], ItemB=frame['int']))
    mixed_dup_panel.items = ['ItemA', 'ItemA']
    panel = dict(float=Panel(dict(ItemA=frame['float'], ItemB=frame['float'] + 1)),
                 dup=Panel(np.arange(30).reshape(3, 5, 2).astype(np.float64),
                           items=['A', 'B', 'A']),
                 mixed_dup=mixed_dup_panel)

    return dict(series=series,
                frame=frame,
                panel=panel,
                index=index,
                scalars=scalars,
                mi=mi,
                sp_series=dict(float=_create_sp_series(),
                               ts=_create_sp_tsseries()),
                sp_frame=dict(float=_create_sp_frame()))


def create_pickle_data():
    data = create_data()

    # Pre-0.14.1 versions generated non-unpicklable mixed-type frames and
    # panels if their columns/items were non-unique.
    if LooseVersion(pandas.__version__) < '0.14.1':
        del data['frame']['mixed_dup']
        del data['panel']['mixed_dup']
    return data


def create_msgpack_data():
    data = create_data()
    if LooseVersion(pandas.__version__) < '0.17.0':
        del data['frame']['mixed_dup']
        del data['panel']['mixed_dup']
        del data['frame']['dup']
        del data['panel']['dup']
    # Not supported
    del data['sp_series']
    del data['sp_frame']
    del data['series']['cat']
    del data['frame']['cat_onecol']
    del data['frame']['cat_and_float']
    return data


def platform_name():
    return '_'.join([str(pandas.__version__), str(pl.machine()), str(pl.system().lower()), str(pl.python_version())])


def write_legacy_pickles(output_dir):

    # make sure we are < 0.13 compat (in py3)
    try:
        from pandas.compat import zip, cPickle as pickle
    except:
        import pickle

    version = pandas.__version__

    print("This script generates a storage file for the current arch, system, and python version")
    print("  pandas version: {0}".format(version))
    print("  output dir    : {0}".format(output_dir))
    print("  storage format: pickle")

    pth = '{0}.pickle'.format(platform_name())

    fh = open(os.path.join(output_dir, pth), 'wb')
    pickle.dump(create_pickle_data(), fh, pickle.HIGHEST_PROTOCOL)
    fh.close()

    print("created pickle file: %s" % pth)


def write_legacy_msgpack(output_dir):

    version = pandas.__version__

    print("This script generates a storage file for the current arch, system, and python version")
    print("  pandas version: {0}".format(version))
    print("  output dir    : {0}".format(output_dir))
    print("  storage format: msgpack")

    pth = '{0}.msgpack'.format(platform_name())
    to_msgpack(os.path.join(output_dir, pth), create_msgpack_data())

    print("created msgpack file: %s" % pth)


def write_legacy_file():
    # force our cwd to be the first searched
    sys.path.insert(0, '.')

    if len(sys.argv) != 3:
        exit("Specify output directory and storage type: generate_legacy_storage_files.py <output_dir> <storage_type>")

    output_dir = str(sys.argv[1])
    storage_type = str(sys.argv[2])

    if storage_type == 'pickle':
        write_legacy_pickles(output_dir=output_dir)
    elif storage_type == 'msgpack':
        write_legacy_msgpack(output_dir=output_dir)
    else:
        exit("storage_type must be one of {'pickle', 'msgpack'}")


if __name__ == '__main__':
    write_legacy_file()
