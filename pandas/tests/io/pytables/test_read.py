import pytest
from warnings import catch_warnings

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from pandas.tests.io.pytables.common import (ensure_clean_store, _maybe_remove,
                                             ensure_clean_path)
from pandas.util.testing import assert_frame_equal
from pandas.compat import (BytesIO, range, PY35, is_platform_windows)
import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas.io.pytables import read_hdf, HDFStore, TableIterator


def test_read_hdf_open_store():
    # GH10330
    # No check for non-string path_or-buf, and no test of open store
    df = DataFrame(np.random.rand(4, 5),
                   index=list('abcd'),
                   columns=list('ABCDE'))
    df.index.name = 'letters'
    df = df.set_index(keys='E', append=True)

    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', mode='w')
        direct = read_hdf(path, 'df')
        store = HDFStore(path, mode='r')
        indirect = read_hdf(store, 'df')
        tm.assert_frame_equal(direct, indirect)
        assert store.is_open
        store.close()


def test_read_hdf_iterator():
    df = DataFrame(np.random.rand(4, 5),
                   index=list('abcd'),
                   columns=list('ABCDE'))
    df.index.name = 'letters'
    df = df.set_index(keys='E', append=True)

    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', mode='w', format='t')
        direct = read_hdf(path, 'df')
        iterator = read_hdf(path, 'df', iterator=True)
        assert isinstance(iterator, TableIterator)
        indirect = next(iterator.__iter__())
        tm.assert_frame_equal(direct, indirect)
        iterator.store.close()


def test_read_hdf_errors():
    df = DataFrame(np.random.rand(4, 5),
                   index=list('abcd'),
                   columns=list('ABCDE'))

    with ensure_clean_path() as path:
        pytest.raises(IOError, read_hdf, path, 'key')
        df.to_hdf(path, 'df')
        store = HDFStore(path, mode='r')
        store.close()
        pytest.raises(IOError, read_hdf, store, 'df')


def test_read_hdf_generic_buffer_errors():
    pytest.raises(NotImplementedError, read_hdf, BytesIO(b''), 'df')


def test_read_nokey():
    df = DataFrame(np.random.rand(4, 5),
                   index=list('abcd'),
                   columns=list('ABCDE'))

    # Categorical dtype not supported for "fixed" format. So no need
    # to test with that dtype in the dataframe here.
    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', mode='a')
        reread = read_hdf(path)
        assert_frame_equal(df, reread)
        df.to_hdf(path, 'df2', mode='a')
        pytest.raises(ValueError, read_hdf, path)


def test_read_nokey_table():
    # GH13231
    df = DataFrame({'i': range(5),
                    'c': Series(list('abacd'), dtype='category')})

    with ensure_clean_path() as path:
        df.to_hdf(path, 'df', mode='a', format='table')
        reread = read_hdf(path)
        assert_frame_equal(df, reread)
        df.to_hdf(path, 'df2', mode='a', format='table')
        pytest.raises(ValueError, read_hdf, path)


def test_read_nokey_empty():
    with ensure_clean_path() as path:
        store = HDFStore(path)
        store.close()
        pytest.raises(ValueError, read_hdf, path)


@td.skip_if_no('pathlib')
def test_read_from_pathlib_path():
    # GH11773
    Path = pytest.importorskip('pathlib').Path
    expected = DataFrame(np.random.rand(4, 5),
                         index=list('abcd'),
                         columns=list('ABCDE'))
    with ensure_clean_path() as filename:
        path_obj = Path(filename)

        expected.to_hdf(path_obj, 'df', mode='a')
        actual = read_hdf(path_obj, 'df')

    tm.assert_frame_equal(expected, actual)


@td.skip_if_no('py.path')
def test_read_from_py_localpath():
    # GH11773
    LocalPath = pytest.importorskip('py.path').local
    expected = DataFrame(np.random.rand(4, 5),
                         index=list('abcd'),
                         columns=list('ABCDE'))
    with ensure_clean_path() as filename:
        path_obj = LocalPath(filename)

        expected.to_hdf(path_obj, 'df', mode='a')
        actual = read_hdf(path_obj, 'df')

    tm.assert_frame_equal(expected, actual)


@pytest.mark.parametrize('format', ['fixed', 'table'])
def test_read_hdf_series_mode_r(format):
    # GH 16583
    # Tests that reading a Series saved to an HDF file
    # still works if a mode='r' argument is supplied
    series = tm.makeFloatSeries()
    with ensure_clean_path() as path:
        series.to_hdf(path, key='data', format=format)
        result = pd.read_hdf(path, key='data', mode='r')
    tm.assert_series_equal(result, series)


def test_read_py2_hdf_file_in_py3():
    # GH 16781
    # tests reading a PeriodIndex DataFrame written in Python2 in Python3
    # the file was generated in Python 2.7 like so:
    # df = pd.DataFrame([1.,2,3], index=pd.PeriodIndex(
    #              ['2015-01-01', '2015-01-02', '2015-01-05'], freq='B'))
    # df.to_hdf('periodindex_0.20.1_x86_64_darwin_2.7.13.h5', 'p')
    expected = pd.DataFrame([1., 2, 3], index=pd.PeriodIndex(
        ['2015-01-01', '2015-01-02', '2015-01-05'], freq='B'))

    with ensure_clean_store(
            tm.get_data_path(
                'legacy_hdf/periodindex_0.20.1_x86_64_darwin_2.7.13.h5'),
            mode='r') as store:
        result = store['p']
        assert_frame_equal(result, expected)


def test_read_column():
    df = tm.makeTimeDataFrame()
    with ensure_clean_store() as store:
        _maybe_remove(store, 'df')
        store.append('df', df)

        # error
        pytest.raises(KeyError, store.select_column, 'df', 'foo')

        def f():
            store.select_column('df', 'index', where=['index>5'])
        pytest.raises(Exception, f)

        # valid
        result = store.select_column('df', 'index')
        tm.assert_almost_equal(result.values, Series(df.index).values)
        assert isinstance(result, Series)

        # not a data indexable column
        pytest.raises(
            ValueError, store.select_column, 'df', 'values_block_0')

        # a data column
        df2 = df.copy()
        df2['string'] = 'foo'
        store.append('df2', df2, data_columns=['string'])
        result = store.select_column('df2', 'string')
        tm.assert_almost_equal(result.values, df2['string'].values)

        # a data column with NaNs, result excludes the NaNs
        df3 = df.copy()
        df3['string'] = 'foo'
        df3.loc[4:6, 'string'] = np.nan
        store.append('df3', df3, data_columns=['string'])
        result = store.select_column('df3', 'string')
        tm.assert_almost_equal(result.values, df3['string'].values)

        # start/stop
        result = store.select_column('df3', 'string', start=2)
        tm.assert_almost_equal(result.values, df3['string'].values[2:])

        result = store.select_column('df3', 'string', start=-2)
        tm.assert_almost_equal(result.values, df3['string'].values[-2:])

        result = store.select_column('df3', 'string', stop=2)
        tm.assert_almost_equal(result.values, df3['string'].values[:2])

        result = store.select_column('df3', 'string', stop=-2)
        tm.assert_almost_equal(result.values, df3['string'].values[:-2])

        result = store.select_column('df3', 'string', start=2, stop=-2)
        tm.assert_almost_equal(result.values, df3['string'].values[2:-2])

        result = store.select_column('df3', 'string', start=-2, stop=2)
        tm.assert_almost_equal(result.values, df3['string'].values[-2:2])

        # GH 10392 - make sure column name is preserved
        df4 = DataFrame({'A': np.random.randn(10), 'B': 'foo'})
        store.append('df4', df4, data_columns=True)
        expected = df4['B']
        result = store.select_column('df4', 'B')
        tm.assert_series_equal(result, expected)


def test_pytables_native_read():
    with ensure_clean_store(
            tm.get_data_path('legacy_hdf/pytables_native.h5'),
            mode='r') as store:
        d2 = store['detector/readout']
        assert isinstance(d2, DataFrame)


@pytest.mark.skipif(PY35 and is_platform_windows(),
                    reason="native2 read fails oddly on windows / 3.5")
def test_pytables_native2_read():
    with ensure_clean_store(
            tm.get_data_path('legacy_hdf/pytables_native2.h5'),
            mode='r') as store:
        str(store)
        d1 = store['detector']
        assert isinstance(d1, DataFrame)


def test_legacy_table_read():
    # legacy table types
    with ensure_clean_store(
            tm.get_data_path('legacy_hdf/legacy_table.h5'),
            mode='r') as store:

        with catch_warnings(record=True):
            store.select('df1')
            store.select('df2')
            store.select('wp1')

            # force the frame
            store.select('df2', typ='legacy_frame')

            # old version warning
            pytest.raises(
                Exception, store.select, 'wp1', 'minor_axis=B')

            df2 = store.select('df2')
            result = store.select('df2', 'index>df2.index[2]')
            expected = df2[df2.index > df2.index[2]]
            assert_frame_equal(expected, result)
