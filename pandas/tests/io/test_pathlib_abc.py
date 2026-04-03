"""
Tests for pathlib-abc ReadablePath/WritablePath support.
Testing support via UPath.

These tests verify that pandas IO functions accept pathlib_abc compatible
path objects such as UPath from universal_pathlib.
"""

import numpy as np
import pytest

from pandas._config import using_string_dtype

from pandas.compat import HAS_PYARROW
from pandas.compat.pyarrow import pa_version_under14p0

from pandas import (
    DataFrame,
    date_range,
    read_csv,
    read_excel,
    read_feather,
    read_json,
    read_parquet,
    read_pickle,
    read_stata,
    read_table,
)
import pandas._testing as tm

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)


@pytest.fixture
def upath():
    """Import and return UPath, skipping if not available."""
    upath_module = pytest.importorskip("upath")
    return upath_module.UPath


@pytest.fixture
def cleared_fs():
    """Provide a cleared memory filesystem."""
    fsspec = pytest.importorskip("fsspec")
    memfs = fsspec.filesystem("memory")
    yield memfs
    memfs.store.clear()


@pytest.fixture
def df1():
    """Standard test DataFrame."""
    return DataFrame(
        {
            "int": [1, 3],
            "float": [2.0, np.nan],
            "str": ["t", "s"],
            "dt": date_range("2018-06-18", periods=2),
        }
    )


class TestUPathReadable:
    """Tests for reading data using UPath (ReadablePath subclass)."""

    def test_read_csv(self, cleared_fs, df1, upath):
        # Write test data using fsspec directly
        text = str(df1.to_csv(index=False)).encode()
        with cleared_fs.open("test/test.csv", "wb") as w:
            w.write(text)

        # Read using UPath object
        path = upath("memory://test/test.csv")
        df2 = read_csv(path, parse_dates=["dt"])

        tm.assert_frame_equal(df2, df1)

    def test_read_table(self, cleared_fs, df1, upath):
        # Write test data using fsspec directly
        text = str(df1.to_csv(index=False)).encode()
        with cleared_fs.open("test/test.csv", "wb") as w:
            w.write(text)

        # Read using UPath object
        path = upath("memory://test/test.csv")
        df2 = read_table(path, sep=",", parse_dates=["dt"])

        tm.assert_frame_equal(df2, df1)

    def test_read_json(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})

        # Write JSON using fsspec directly
        json_bytes = df.to_json().encode()
        with cleared_fs.open("test/test.json", "wb") as w:
            w.write(json_bytes)

        # Read using UPath object
        path = upath("memory://test/test.json")
        out = read_json(path)
        tm.assert_frame_equal(df, out)

    @pytest.mark.xfail(
        using_string_dtype() and HAS_PYARROW and not pa_version_under14p0,
        reason="TODO(infer_string) fastparquet",
    )
    def test_read_parquet_fastparquet(self, cleared_fs, df1, upath):
        pytest.importorskip("fastparquet")

        # Write parquet file using string path (known to work)
        df1.to_parquet(
            "memory://test/test.parquet",
            index=False,
            engine="fastparquet",
            compression=None,
        )

        # Read using UPath object
        path = upath("memory://test/test.parquet")
        df2 = read_parquet(path, engine="fastparquet")
        tm.assert_frame_equal(df1, df2)

    def test_read_parquet_pyarrow(self, cleared_fs, df1, upath):
        pytest.importorskip("pyarrow")

        # Write parquet file using string path
        df1.to_parquet(
            "memory://test/test.parquet",
            index=False,
            engine="pyarrow",
            compression=None,
        )

        # Read using UPath object
        path = upath("memory://test/test.parquet")
        df2 = read_parquet(path, engine="pyarrow")

        tm.assert_frame_equal(df2, df1)

    def test_read_feather(self, cleared_fs, upath):
        pytest.importorskip("pyarrow")
        df = DataFrame({"a": [0]})

        # Write feather using string path
        df.to_feather("memory://test/test.feather")

        # Read using UPath object
        path = upath("memory://test/test.feather")
        out = read_feather(path)
        tm.assert_frame_equal(df, out)

    def test_read_pickle(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})

        # Write pickle using string path
        df.to_pickle("memory://test/test.pkl")

        # Read using UPath object
        path = upath("memory://test/test.pkl")
        out = read_pickle(path)
        tm.assert_frame_equal(df, out)

    def test_read_stata(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})

        # Write stata using string path
        df.to_stata("memory://test/test.dta", write_index=False)

        # Read using UPath object
        path = upath("memory://test/test.dta")
        out = read_stata(path)
        tm.assert_frame_equal(df, out.astype("int64"))

    def test_read_excel(self, cleared_fs, df1, upath):
        pytest.importorskip("openpyxl")

        # Write excel using string path
        df1.to_excel("memory://test/test.xlsx", index=True)

        # Read using UPath object
        path = upath("memory://test/test.xlsx")
        df2 = read_excel(path, parse_dates=["dt"], index_col=0)

        tm.assert_frame_equal(df2, df1)


class TestUPathWritable:
    """Tests for writing data using UPath (WritablePath subclass)."""

    def test_to_csv(self, cleared_fs, df1, upath):
        path = upath("memory://test/test.csv")

        # Write using UPath object
        df1.to_csv(path, index=True)

        # Read back using string URL to verify
        df2 = read_csv("memory://test/test.csv", parse_dates=["dt"], index_col=0)

        tm.assert_frame_equal(df2, df1)

    def test_to_json(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.json")

        # Write using UPath object
        df.to_json(path)

        # Read back using string URL
        out = read_json("memory://test/test.json")
        tm.assert_frame_equal(df, out)

    @pytest.mark.xfail(
        using_string_dtype() and HAS_PYARROW and not pa_version_under14p0,
        reason="TODO(infer_string) fastparquet",
    )
    def test_to_parquet_fastparquet(self, cleared_fs, df1, upath):
        pytest.importorskip("fastparquet")
        path = upath("memory://test/test.parquet")

        # Write using UPath object
        df1.to_parquet(path, index=False, engine="fastparquet", compression=None)

        # Read back using string URL
        df2 = read_parquet("memory://test/test.parquet", engine="fastparquet")
        tm.assert_frame_equal(df1, df2)

    def test_to_parquet_pyarrow(self, cleared_fs, df1, upath):
        pytest.importorskip("pyarrow")
        path = upath("memory://test/test.parquet")

        # Write using UPath object
        df1.to_parquet(path, index=False, engine="pyarrow", compression=None)

        # Read back using string URL
        df2 = read_parquet("memory://test/test.parquet", engine="pyarrow")

        tm.assert_frame_equal(df2, df1)

    def test_to_feather(self, cleared_fs, upath):
        pytest.importorskip("pyarrow")
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.feather")

        # Write using UPath object
        df.to_feather(path)

        # Read back using string URL
        out = read_feather("memory://test/test.feather")
        tm.assert_frame_equal(df, out)

    def test_to_pickle(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.pkl")

        # Write using UPath object
        df.to_pickle(path)

        # Read back using string URL
        out = read_pickle("memory://test/test.pkl")
        tm.assert_frame_equal(df, out)

    def test_to_stata(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.dta")

        # Write using UPath object
        df.to_stata(path, write_index=False)

        # Read back using string URL
        out = read_stata("memory://test/test.dta")
        tm.assert_frame_equal(df, out.astype("int64"))

    def test_to_excel(self, cleared_fs, df1, upath):
        pytest.importorskip("openpyxl")
        path = upath("memory://test/test.xlsx")

        # Write using UPath object
        df1.to_excel(path, index=True)

        # Read back using string URL
        df2 = read_excel("memory://test/test.xlsx", parse_dates=["dt"], index_col=0)

        tm.assert_frame_equal(df2, df1)

    def test_to_markdown(self, cleared_fs, upath):
        pytest.importorskip("tabulate")
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.md")

        # Write using UPath object
        df.to_markdown(path)

        # Verify file exists and has content
        assert cleared_fs.cat("test/test.md")


class TestUPathRoundtrip:
    """Tests for round-trip read/write using UPath objects."""

    def test_csv_roundtrip(self, cleared_fs, df1, upath):
        path = upath("memory://test/test.csv")

        # Write using UPath
        df1.to_csv(path, index=False)

        # Read using UPath
        df2 = read_csv(path, parse_dates=["dt"])

        tm.assert_frame_equal(df2, df1)

    def test_json_roundtrip(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.json")

        # Write and read using UPath
        df.to_json(path)
        out = read_json(path)
        tm.assert_frame_equal(df, out)

    def test_parquet_roundtrip_pyarrow(self, cleared_fs, df1, upath):
        pytest.importorskip("pyarrow")
        path = upath("memory://test/test.parquet")

        # Write and read using UPath
        df1.to_parquet(path, index=False, engine="pyarrow", compression=None)
        df2 = read_parquet(path, engine="pyarrow")

        tm.assert_frame_equal(df2, df1)

    def test_pickle_roundtrip(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.pkl")

        # Write and read using UPath
        df.to_pickle(path)
        out = read_pickle(path)
        tm.assert_frame_equal(df, out)

    def test_feather_roundtrip(self, cleared_fs, upath):
        pytest.importorskip("pyarrow")
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.feather")

        # Write and read using UPath
        df.to_feather(path)
        out = read_feather(path)
        tm.assert_frame_equal(df, out)

    def test_excel_roundtrip(self, cleared_fs, df1, upath):
        pytest.importorskip("openpyxl")
        path = upath("memory://test/test.xlsx")

        # Write and read using UPath
        df1.to_excel(path, index=True)
        df2 = read_excel(path, parse_dates=["dt"], index_col=0)

        tm.assert_frame_equal(df2, df1)

    def test_stata_roundtrip(self, cleared_fs, upath):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.dta")

        # Write and read using UPath
        df.to_stata(path, write_index=False)
        out = read_stata(path)
        tm.assert_frame_equal(df, out.astype("int64"))


class TestUPathWithCompression:
    """Tests for UPath with various compression options."""

    @pytest.mark.parametrize("compression", ["gzip", "bz2", "xz"])
    def test_csv_compression_roundtrip(self, cleared_fs, df1, upath, compression):
        ext_map = {"gzip": ".gz", "bz2": ".bz2", "xz": ".xz"}
        ext = ext_map.get(compression, "")
        path = upath(f"memory://test/test.csv{ext}")

        # Write and read using UPath with compression
        df1.to_csv(path, index=False, compression=compression)
        df2 = read_csv(path, parse_dates=["dt"], compression=compression)

        tm.assert_frame_equal(df2, df1)

    @pytest.mark.parametrize("compression", [None, "gzip", "bz2"])
    def test_json_compression_roundtrip(self, cleared_fs, upath, compression):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.json")

        # Write and read using UPath with compression
        df.to_json(path, compression=compression)
        out = read_json(path, compression=compression)
        tm.assert_frame_equal(df, out)

    @pytest.mark.parametrize("compression", [None, "gzip"])
    def test_pickle_compression_roundtrip(self, cleared_fs, upath, compression):
        df = DataFrame({"a": [0]})
        path = upath("memory://test/test.pkl")

        # Write and read using UPath with compression
        df.to_pickle(path, compression=compression)
        out = read_pickle(path, compression=compression)
        tm.assert_frame_equal(df, out)


class TestUPathErrorHandling:
    """Tests for error handling with UPath objects."""

    def test_read_nonexistent_file(self, cleared_fs, upath):
        path = upath("memory://nonexistent/file.csv")
        with pytest.raises((FileNotFoundError, OSError)):
            read_csv(path)

    def test_write_creates_parent_dirs(self, cleared_fs, df1, upath):
        # UPath should be able to create parent directories
        path = upath("memory://deep/nested/path/test.csv")

        # This should work if parent creation is supported
        path.parent.mkdir(parents=True, exist_ok=True)
        df1.to_csv(path, index=False)

        # Verify the file was written
        df2 = read_csv("memory://deep/nested/path/test.csv", parse_dates=["dt"])
        tm.assert_frame_equal(df2, df1)
