"""
manage legacy pickle tests

How to add pickle tests:

1. Install pandas version intended to output the pickle.

2. Execute "generate_legacy_storage_files.py" to create the pickle.
$ python generate_legacy_storage_files.py <output_dir> pickle

3. Move the created pickle to "data/legacy_pickle/<version>" directory.
"""
import bz2
import glob
import gzip
import os
import pickle
import shutil
from warnings import catch_warnings, simplefilter
import zipfile

import pytest

from pandas.compat import _get_lzma_file, _import_lzma, is_platform_little_endian

import pandas as pd
from pandas import Index
import pandas._testing as tm

from pandas.tseries.offsets import Day, MonthEnd

lzma = _import_lzma()


@pytest.fixture(scope="module")
def current_pickle_data():
    # our current version pickle data
    from pandas.tests.io.generate_legacy_storage_files import create_pickle_data

    return create_pickle_data()


# ---------------------
# comparison functions
# ---------------------
def compare_element(result, expected, typ, version=None):
    if isinstance(expected, Index):
        tm.assert_index_equal(expected, result)
        return

    if typ.startswith("sp_"):
        comparator = tm.assert_equal
        comparator(result, expected)
    elif typ == "timestamp":
        if expected is pd.NaT:
            assert result is pd.NaT
        else:
            assert result == expected
            assert result.freq == expected.freq
    else:
        comparator = getattr(
            tm, "assert_{typ}_equal".format(typ=typ), tm.assert_almost_equal
        )
        comparator(result, expected)


def compare(data, vf, version):

    data = pd.read_pickle(vf)

    m = globals()
    for typ, dv in data.items():
        for dt, result in dv.items():
            expected = data[typ][dt]

            # use a specific comparator
            # if available
            comparator = "compare_{typ}_{dt}".format(typ=typ, dt=dt)

            comparator = m.get(comparator, m["compare_element"])
            comparator(result, expected, typ, version)
    return data


def compare_series_ts(result, expected, typ, version):
    # GH 7748
    tm.assert_series_equal(result, expected)
    assert result.index.freq == expected.index.freq
    assert not result.index.freq.normalize
    tm.assert_series_equal(result > 0, expected > 0)

    # GH 9291
    freq = result.index.freq
    assert freq + Day(1) == Day(2)

    res = freq + pd.Timedelta(hours=1)
    assert isinstance(res, pd.Timedelta)
    assert res == pd.Timedelta(days=1, hours=1)

    res = freq + pd.Timedelta(nanoseconds=1)
    assert isinstance(res, pd.Timedelta)
    assert res == pd.Timedelta(days=1, nanoseconds=1)


def compare_series_dt_tz(result, expected, typ, version):
    tm.assert_series_equal(result, expected)


def compare_series_cat(result, expected, typ, version):
    tm.assert_series_equal(result, expected)


def compare_frame_dt_mixed_tzs(result, expected, typ, version):
    tm.assert_frame_equal(result, expected)


def compare_frame_cat_onecol(result, expected, typ, version):
    tm.assert_frame_equal(result, expected)


def compare_frame_cat_and_float(result, expected, typ, version):
    compare_frame_cat_onecol(result, expected, typ, version)


def compare_index_period(result, expected, typ, version):
    tm.assert_index_equal(result, expected)
    assert isinstance(result.freq, MonthEnd)
    assert result.freq == MonthEnd()
    assert result.freqstr == "M"
    tm.assert_index_equal(result.shift(2), expected.shift(2))


files = glob.glob(
    os.path.join(os.path.dirname(__file__), "data", "legacy_pickle", "*", "*.pickle")
)


@pytest.fixture(params=files)
def legacy_pickle(request, datapath):
    return datapath(request.param)


# ---------------------
# tests
# ---------------------
def test_pickles(current_pickle_data, legacy_pickle):
    if not is_platform_little_endian():
        pytest.skip("known failure on non-little endian")

    version = os.path.basename(os.path.dirname(legacy_pickle))
    with catch_warnings(record=True):
        simplefilter("ignore")
        compare(current_pickle_data, legacy_pickle, version)


def test_round_trip_current(current_pickle_data):
    def python_pickler(obj, path):
        with open(path, "wb") as fh:
            pickle.dump(obj, fh, protocol=-1)

    def python_unpickler(path):
        with open(path, "rb") as fh:
            fh.seek(0)
            return pickle.load(fh)

    data = current_pickle_data
    for typ, dv in data.items():
        for dt, expected in dv.items():

            for writer in [pd.to_pickle, python_pickler]:
                if writer is None:
                    continue

                with tm.ensure_clean() as path:

                    # test writing with each pickler
                    writer(expected, path)

                    # test reading with each unpickler
                    result = pd.read_pickle(path)
                    compare_element(result, expected, typ)

                    result = python_unpickler(path)
                    compare_element(result, expected, typ)


def test_pickle_path_pathlib():
    df = tm.makeDataFrame()
    result = tm.round_trip_pathlib(df.to_pickle, pd.read_pickle)
    tm.assert_frame_equal(df, result)


def test_pickle_path_localpath():
    df = tm.makeDataFrame()
    result = tm.round_trip_localpath(df.to_pickle, pd.read_pickle)
    tm.assert_frame_equal(df, result)


def test_legacy_sparse_warning(datapath):
    """

    Generated with

    >>> df = pd.DataFrame({"A": [1, 2, 3, 4], "B": [0, 0, 1, 1]}).to_sparse()
    >>> df.to_pickle("pandas/tests/io/data/pickle/sparseframe-0.20.3.pickle.gz",
    ...              compression="gzip")

    >>> s = df['B']
    >>> s.to_pickle("pandas/tests/io/data/pickle/sparseseries-0.20.3.pickle.gz",
    ...             compression="gzip")
    """
    with tm.assert_produces_warning(FutureWarning):
        simplefilter("ignore", DeprecationWarning)  # from boto
        pd.read_pickle(
            datapath("io", "data", "pickle", "sparseseries-0.20.3.pickle.gz"),
            compression="gzip",
        )

    with tm.assert_produces_warning(FutureWarning):
        simplefilter("ignore", DeprecationWarning)  # from boto
        pd.read_pickle(
            datapath("io", "data", "pickle", "sparseframe-0.20.3.pickle.gz"),
            compression="gzip",
        )


# ---------------------
# test pickle compression
# ---------------------


@pytest.fixture
def get_random_path():
    return "__{}__.pickle".format(tm.rands(10))


class TestCompression:

    _compression_to_extension = {
        None: ".none",
        "gzip": ".gz",
        "bz2": ".bz2",
        "zip": ".zip",
        "xz": ".xz",
    }

    def compress_file(self, src_path, dest_path, compression):
        if compression is None:
            shutil.copyfile(src_path, dest_path)
            return

        if compression == "gzip":
            f = gzip.open(dest_path, "w")
        elif compression == "bz2":
            f = bz2.BZ2File(dest_path, "w")
        elif compression == "zip":
            with zipfile.ZipFile(dest_path, "w", compression=zipfile.ZIP_DEFLATED) as f:
                f.write(src_path, os.path.basename(src_path))
        elif compression == "xz":
            f = _get_lzma_file(lzma)(dest_path, "w")
        else:
            msg = "Unrecognized compression type: {}".format(compression)
            raise ValueError(msg)

        if compression != "zip":
            with open(src_path, "rb") as fh, f:
                f.write(fh.read())

    def test_write_explicit(self, compression, get_random_path):
        base = get_random_path
        path1 = base + ".compressed"
        path2 = base + ".raw"

        with tm.ensure_clean(path1) as p1, tm.ensure_clean(path2) as p2:
            df = tm.makeDataFrame()

            # write to compressed file
            df.to_pickle(p1, compression=compression)

            # decompress
            with tm.decompress_file(p1, compression=compression) as f:
                with open(p2, "wb") as fh:
                    fh.write(f.read())

            # read decompressed file
            df2 = pd.read_pickle(p2, compression=None)

            tm.assert_frame_equal(df, df2)

    @pytest.mark.parametrize("compression", ["", "None", "bad", "7z"])
    def test_write_explicit_bad(self, compression, get_random_path):
        with pytest.raises(ValueError, match="Unrecognized compression type"):
            with tm.ensure_clean(get_random_path) as path:
                df = tm.makeDataFrame()
                df.to_pickle(path, compression=compression)

    @pytest.mark.parametrize("ext", ["", ".gz", ".bz2", ".no_compress", ".xz"])
    def test_write_infer(self, ext, get_random_path):
        base = get_random_path
        path1 = base + ext
        path2 = base + ".raw"
        compression = None
        for c in self._compression_to_extension:
            if self._compression_to_extension[c] == ext:
                compression = c
                break

        with tm.ensure_clean(path1) as p1, tm.ensure_clean(path2) as p2:
            df = tm.makeDataFrame()

            # write to compressed file by inferred compression method
            df.to_pickle(p1)

            # decompress
            with tm.decompress_file(p1, compression=compression) as f:
                with open(p2, "wb") as fh:
                    fh.write(f.read())

            # read decompressed file
            df2 = pd.read_pickle(p2, compression=None)

            tm.assert_frame_equal(df, df2)

    def test_read_explicit(self, compression, get_random_path):
        base = get_random_path
        path1 = base + ".raw"
        path2 = base + ".compressed"

        with tm.ensure_clean(path1) as p1, tm.ensure_clean(path2) as p2:
            df = tm.makeDataFrame()

            # write to uncompressed file
            df.to_pickle(p1, compression=None)

            # compress
            self.compress_file(p1, p2, compression=compression)

            # read compressed file
            df2 = pd.read_pickle(p2, compression=compression)

            tm.assert_frame_equal(df, df2)

    @pytest.mark.parametrize("ext", ["", ".gz", ".bz2", ".zip", ".no_compress", ".xz"])
    def test_read_infer(self, ext, get_random_path):
        base = get_random_path
        path1 = base + ".raw"
        path2 = base + ext
        compression = None
        for c in self._compression_to_extension:
            if self._compression_to_extension[c] == ext:
                compression = c
                break

        with tm.ensure_clean(path1) as p1, tm.ensure_clean(path2) as p2:
            df = tm.makeDataFrame()

            # write to uncompressed file
            df.to_pickle(p1, compression=None)

            # compress
            self.compress_file(p1, p2, compression=compression)

            # read compressed file by inferred compression method
            df2 = pd.read_pickle(p2)

            tm.assert_frame_equal(df, df2)


# ---------------------
# test pickle compression
# ---------------------


class TestProtocol:
    @pytest.mark.parametrize("protocol", [-1, 0, 1, 2])
    def test_read(self, protocol, get_random_path):
        with tm.ensure_clean(get_random_path) as path:
            df = tm.makeDataFrame()
            df.to_pickle(path, protocol=protocol)
            df2 = pd.read_pickle(path)
            tm.assert_frame_equal(df, df2)


def test_unicode_decode_error():
    # pickle file written with py27, should be readable without raising
    #  UnicodeDecodeError, see GH#28645
    path = os.path.join(os.path.dirname(__file__), "data", "pickle", "test_py27.pkl")
    df = pd.read_pickle(path)

    # just test the columns are correct since the values are random
    excols = pd.Index(["a", "b", "c"])
    tm.assert_index_equal(df.columns, excols)
