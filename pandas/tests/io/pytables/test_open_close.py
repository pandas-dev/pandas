import pytest
import os

import pandas.util.testing as tm
from pandas.io import pytables as pytables
from .common import ensure_clean_path
from pandas.io.pytables import (HDFStore, ClosedFileError,
                                PossibleDataLossError)


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


def test_reopen_handle():
    with ensure_clean_path() as path:
        store = HDFStore(path, mode='a')
        store['a'] = tm.makeTimeSeries()

        # invalid mode change
        pytest.raises(PossibleDataLossError, store.open, 'w')
        store.close()
        assert not store.is_open

        # truncation ok here
        store.open('w')
        assert store.is_open
        assert len(store) == 0
        store.close()
        assert not store.is_open

        store = HDFStore(path, mode='a')
        store['a'] = tm.makeTimeSeries()

        # reopen as read
        store.open('r')
        assert store.is_open
        assert len(store) == 1
        assert store._mode == 'r'
        store.close()
        assert not store.is_open

        # reopen as append
        store.open('a')
        assert store.is_open
        assert len(store) == 1
        assert store._mode == 'a'
        store.close()
        assert not store.is_open

        # reopen as append (again)
        store.open('a')
        assert store.is_open
        assert len(store) == 1
        assert store._mode == 'a'
        store.close()
        assert not store.is_open


def test_open_args():
    with ensure_clean_path() as path:
        df = tm.makeDataFrame()

        # create an in memory store
        store = HDFStore(path, mode='a', driver='H5FD_CORE',
                         driver_core_backing_store=0)
        store['df'] = df
        store.append('df2', df)

        tm.assert_frame_equal(store['df'], df)
        tm.assert_frame_equal(store['df2'], df)

        store.close()

        # the file should not have actually been written
        assert not os.path.exists(path)
