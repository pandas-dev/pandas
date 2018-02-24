import pytest

import numpy as np
import pandas as pd
from pandas import DataFrame
import pandas.util.testing as tm
from .common import ensure_clean_path

tables = pytest.importorskip('tables')


def test_complibs_default_settings():
    # GH15943
    df = tm.makeDataFrame()

    # Set complevel and check if complib is automatically set to
    # default value
    with ensure_clean_path() as tmpfile:
        df.to_hdf(tmpfile, 'df', complevel=9)
        result = pd.read_hdf(tmpfile, 'df')
        tm.assert_frame_equal(result, df)

        with tables.open_file(tmpfile, mode='r') as h5file:
            for node in h5file.walk_nodes(where='/df', classname='Leaf'):
                assert node.filters.complevel == 9
                assert node.filters.complib == 'zlib'

    # Set complib and check to see if compression is disabled
    with ensure_clean_path() as tmpfile:
        df.to_hdf(tmpfile, 'df', complib='zlib')
        result = pd.read_hdf(tmpfile, 'df')
        tm.assert_frame_equal(result, df)

        with tables.open_file(tmpfile, mode='r') as h5file:
            for node in h5file.walk_nodes(where='/df', classname='Leaf'):
                assert node.filters.complevel == 0
                assert node.filters.complib is None

    # Check if not setting complib or complevel results in no compression
    with ensure_clean_path() as tmpfile:
        df.to_hdf(tmpfile, 'df')
        result = pd.read_hdf(tmpfile, 'df')
        tm.assert_frame_equal(result, df)

        with tables.open_file(tmpfile, mode='r') as h5file:
            for node in h5file.walk_nodes(where='/df', classname='Leaf'):
                assert node.filters.complevel == 0
                assert node.filters.complib is None

    # Check if file-defaults can be overridden on a per table basis
    with ensure_clean_path() as tmpfile:
        store = pd.HDFStore(tmpfile)
        store.append('dfc', df, complevel=9, complib='blosc')
        store.append('df', df)
        store.close()

        with tables.open_file(tmpfile, mode='r') as h5file:
            for node in h5file.walk_nodes(where='/df', classname='Leaf'):
                assert node.filters.complevel == 0
                assert node.filters.complib is None
            for node in h5file.walk_nodes(where='/dfc', classname='Leaf'):
                assert node.filters.complevel == 9
                assert node.filters.complib == 'blosc'


def test_complibs():
    # GH14478
    df = tm.makeDataFrame()

    # Building list of all complibs and complevels tuples
    all_complibs = tables.filters.all_complibs
    # Remove lzo if its not available on this platform
    if not tables.which_lib_version('lzo'):
        all_complibs.remove('lzo')
    # Remove bzip2 if its not available on this platform
    if not tables.which_lib_version("bzip2"):
        all_complibs.remove("bzip2")

    all_levels = range(0, 10)
    all_tests = [(lib, lvl) for lib in all_complibs for lvl in all_levels]

    for (lib, lvl) in all_tests:
        with ensure_clean_path() as tmpfile:
            gname = 'foo'

            # Write and read file to see if data is consistent
            df.to_hdf(tmpfile, gname, complib=lib, complevel=lvl)
            result = pd.read_hdf(tmpfile, gname)
            tm.assert_frame_equal(result, df)

            # Open file and check metadata
            # for correct amount of compression
            h5table = tables.open_file(tmpfile, mode='r')
            for node in h5table.walk_nodes(where='/' + gname,
                                           classname='Leaf'):
                assert node.filters.complevel == lvl
                if lvl == 0:
                    assert node.filters.complib is None
                else:
                    assert node.filters.complib == lib
            h5table.close()


def test_invalid_complib():
    df = DataFrame(np.random.rand(4, 5),
                   index=list('abcd'),
                   columns=list('ABCDE'))
    with ensure_clean_path() as path:
        with pytest.raises(ValueError):
            df.to_hdf(path, 'df', complib='foolib')
