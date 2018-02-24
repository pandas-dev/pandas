import os
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.compat import PY36


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


@pytest.mark.skipif(not PY36, reason="Need python 3.6")
def test_fspath():
    with tm.ensure_clean('foo.h5') as path:
        with pd.HDFStore(path) as store:
            assert os.fspath(store) == str(path)
