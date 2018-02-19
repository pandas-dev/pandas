import os
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.tests.io.pytables.common import ensure_clean_path
from pandas.compat import PY36
from pandas.io import pytables as pytables
from pandas.io.pytables import HDFStore, ClosedFileError


def test_path_pathlib():
    df = tm.makeDataFrame()

    result = tm.round_trip_pathlib(
        lambda p: df.to_hdf(p, 'df'),
        lambda p: pd.read_hdf(p, 'df'))
    tm.assert_frame_equal(df, result)


def test_path_pathlib_hdfstore():
    df = tm.makeDataFrame()

    def writer(path):
        with pd.HDFStore(path) as store:
            df.to_hdf(store, 'df')

    def reader(path):
        with pd.HDFStore(path) as store:
            return pd.read_hdf(store, 'df')

    result = tm.round_trip_pathlib(writer, reader)
    tm.assert_frame_equal(df, result)


def test_pickle_path_localpath():
    df = tm.makeDataFrame()
    result = tm.round_trip_pathlib(
        lambda p: df.to_hdf(p, 'df'),
        lambda p: pd.read_hdf(p, 'df'))
    tm.assert_frame_equal(df, result)


def test_path_localpath_hdfstore():
    df = tm.makeDataFrame()

    def writer(path):
        with pd.HDFStore(path) as store:
            df.to_hdf(store, 'df')

    def reader(path):
        with pd.HDFStore(path) as store:
            return pd.read_hdf(store, 'df')

    result = tm.round_trip_localpath(writer, reader)
    tm.assert_frame_equal(df, result)


def test_multiple_open_close():
    # gh-4409: open & close multiple times
    with ensure_clean_path() as path:

        df = tm.makeDataFrame()
        df.to_hdf(path, 'df', mode='w', format='table')

        # single
        store = HDFStore(path)
        assert 'CLOSED' not in store.info()
        assert store.is_open

        store.close()
        assert 'CLOSED' in store.info()
        assert not store.is_open

    with ensure_clean_path() as path:
        if pytables._table_file_open_policy_is_strict:
            # multiples
            store1 = HDFStore(path)

            def f():
                HDFStore(path)
            pytest.raises(ValueError, f)
            store1.close()

        else:
            # multiples
            store1 = HDFStore(path)
            store2 = HDFStore(path)

            assert 'CLOSED' not in store1.info()
            assert 'CLOSED' not in store2.info()
            assert store1.is_open
            assert store2.is_open

            store1.close()
            assert 'CLOSED' in store1.info()
            assert not store1.is_open
            assert 'CLOSED' not in store2.info()
            assert store2.is_open

            store2.close()
            assert 'CLOSED' in store1.info()
            assert 'CLOSED' in store2.info()
            assert not store1.is_open
            assert not store2.is_open

            # nested close
            store = HDFStore(path, mode='w')
            store.append('df', df)

            store2 = HDFStore(path)
            store2.append('df2', df)
            store2.close()
            assert 'CLOSED' in store2.info()
            assert not store2.is_open

            store.close()
            assert 'CLOSED' in store.info()
            assert not store.is_open

            # double closing
            store = HDFStore(path, mode='w')
            store.append('df', df)

            store2 = HDFStore(path)
            store.close()
            assert 'CLOSED' in store.info()
            assert not store.is_open

            store2.close()
            assert 'CLOSED' in store2.info()
            assert not store2.is_open

    # ops on a closed store
    with ensure_clean_path() as path:
        df = tm.makeDataFrame()
        df.to_hdf(path, 'df', mode='w', format='table')

        store = HDFStore(path)
        store.close()

        pytest.raises(ClosedFileError, store.keys)
        pytest.raises(ClosedFileError, lambda: 'df' in store)
        pytest.raises(ClosedFileError, lambda: len(store))
        pytest.raises(ClosedFileError, lambda: store['df'])
        pytest.raises(AttributeError, lambda: store.df)
        pytest.raises(ClosedFileError, store.select, 'df')
        pytest.raises(ClosedFileError, store.get, 'df')
        pytest.raises(ClosedFileError, store.append, 'df2', df)
        pytest.raises(ClosedFileError, store.put, 'df3', df)
        pytest.raises(ClosedFileError, store.get_storer, 'df2')
        pytest.raises(ClosedFileError, store.remove, 'df2')

        def f():
            store.select('df')
        tm.assert_raises_regex(ClosedFileError, 'file is not open', f)


@pytest.mark.skipif(not PY36, reason="Need python 3.6")
def test_fspath():
    with tm.ensure_clean('foo.h5') as path:
        with pd.HDFStore(path) as store:
            assert os.fspath(store) == str(path)
