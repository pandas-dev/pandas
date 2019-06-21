import os
import tempfile
from contextlib import contextmanager

import pytest

import pandas as pd

from pandas.util.testing import assert_frame_equal


tables = pytest.importorskip('tables')


@pytest.fixture
def pytables_hdf5_file(tmp_path):
    """Use PyTables to create a simple HDF5 file.

    There is no need for temporary file cleanup, pytest makes sure that there is
    no collision between tests and that storage needs to not grow indefinitely.
    """

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

    # The `tmp_path` fixture provides a temporary directory unique to the
    # individual test invocation. Create a file in there.
    h5path = tmp_path / 'written_with_pytables.h5'

    with tables.open_file(h5path, mode='w') as f:
        t = f.create_table('/', name=objname, description=table_schema)
        for sample in testsamples:
            for key, value in sample.items():
                t.row[key] = value
            t.row.append()

    return h5path, objname, pd.DataFrame(testsamples)


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
