import gc

import numpy as np
import pytest

from pandas import DataFrame, date_range, read_csv
import pandas._testing as tm
from pandas.util import _test_decorators as td

from pandas.io.common import is_fsspec_url

df1 = DataFrame(
    {
        "int": [1, 3],
        "float": [2.0, np.nan],
        "str": ["t", "s"],
        "dt": date_range("2018-06-18", periods=2),
    }
)
text = df1.to_csv(index=False).encode()


@pytest.fixture
@td.skip_if_installed("fsspec")
def cleared_fs():
    import fsspec

    memfs = fsspec.filesystem("memory")
    try:
        yield memfs
    finally:
        memfs.store.clear()


def test_is_fsspec_url():
    assert is_fsspec_url("gcs://pandas/somethingelse.com")
    assert is_fsspec_url("gs://pandas/somethingelse.com")
    assert not is_fsspec_url("random:pandas/somethingelse.com")


def test_read_csv(cleared_fs):
    from fsspec.implementations.memory import MemoryFile

    cleared_fs.store["test/test.csv"] = MemoryFile(data=text)
    df2 = read_csv("memory://test/test.csv", parse_dates=["dt"])

    tm.assert_frame_equal(df1, df2)


def test_to_csv(cleared_fs):
    df1.to_csv("memory://test/test.csv", index=True)
    gc.collect()  # pandas does not explicitly close file buffers
    df2 = read_csv("memory://test/test.csv", parse_dates=["dt"], index_col=0)

    tm.assert_frame_equal(df1, df2)


@td.skip_if_no("fastparquet")
def test_to_parquet_new_file(monkeypatch):
    """Regression test for writing to a not-yet-existent GCS Parquet file."""
    df1.to_parquet(
        "memory://test/test.csv", index=True, engine="fastparquet", compression=None
    )


@td.skip_if_installed("fsspec")
def test_not_present_exception():
    with pytest.raises(ImportError) as e:
        read_csv("memory://test/test.csv")
        assert "fsspec library is required" in str(e.value)
