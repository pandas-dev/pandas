import os
import tempfile
from contextlib import contextmanager

import pytest

import pandas as pd

from pandas.util.testing import assert_frame_equal


tables = pytest.importorskip('tables')


def safe_remove(path):
    if path is not None:
        try:
            os.remove(path)
        except OSError:
            pass


def create_tempfile(path):
    return os.path.join(tempfile.gettempdir(), path)


@contextmanager
def ensure_clean_path(path):
    try:
        if isinstance(path, list):
            filenames = [create_tempfile(p) for p in path]
            yield filenames
        else:
            filenames = [create_tempfile(path)]
            yield filenames[0]
    finally:
        for f in filenames:
            safe_remove(f)


@pytest.fixture
def pytables_hdf5_file():
    """Use PyTables to create a simple HDF5 file."""

    table_schema = {
        'c0': tables.Time64Col(pos=0),
        'c1': tables.StringCol(5, pos=1),
        'c2': tables.Int64Col(pos=2),
    }

    t0 = 1561105000.0

    testsamples = [
        {'c0': t0, 'c1': 'aaaaa', 'c2': 1},
        {'c0': t0 + 1, 'c1': 'bbbbb', 'c2': 2},
        {'c0': t0 + 2, 'c1': 'ccccc', 'c2': 10**5},
        {'c0': t0 + 3, 'c1': 'ddddd', 'c2': 4294967295},
    ]

    objname = 'pandas_test_timeseries'

    with ensure_clean_path('pytables_hdf5_file') as path:
        # The `ensure_clean_path` context mgr removes the temp file upon exit.
        with tables.open_file(path, mode='w') as f:
            t = f.create_table('/', name=objname, description=table_schema)
            for sample in testsamples:
                for key, value in sample.items():
                    t.row[key] = value
                t.row.append()

        yield path, objname, pd.DataFrame(testsamples)


class TestReadPyTablesHDF5:
    """
    A group of tests which covers reading HDF5 files written by plain PyTables
    (not written by pandas).
    """

    def test_read_complete(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        result = pd.read_hdf(path, key=objname)
        expected = df
        assert_frame_equal(result, expected)

    def test_read_with_start(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        # This is a regression test for pandas-dev/pandas/issues/11188
        result = pd.read_hdf(path, key=objname, start=1)
        expected = df[1:].reset_index(drop=True)
        assert_frame_equal(result, expected)

    def test_read_with_stop(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        # This is a regression test for pandas-dev/pandas/issues/11188
        result = pd.read_hdf(path, key=objname, stop=1)
        expected = df[:1].reset_index(drop=True)
        assert_frame_equal(result, expected)

    def test_read_with_startstop(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        # This is a regression test for pandas-dev/pandas/issues/11188
        result = pd.read_hdf(path, key=objname, start=1, stop=2)
        expected = df[1:2].reset_index(drop=True)
        assert_frame_equal(result, expected)
