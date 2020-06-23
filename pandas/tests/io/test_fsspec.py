import numpy as np
import pytest

from pandas import DataFrame, date_range, read_csv, read_parquet
import pandas._testing as tm
from pandas.util import _test_decorators as td

df1 = DataFrame(
    {
        "int": [1, 3],
        "float": [2.0, np.nan],
        "str": ["t", "s"],
        "dt": date_range("2018-06-18", periods=2),
    }
)
# the ignore on the following line accounts for to_csv returning Optional(str)
# in general, but always str in the case we give no filename
text = df1.to_csv(index=False).encode()  # type: ignore


@pytest.fixture
def cleared_fs():
    fsspec = pytest.importorskip("fsspec")

    memfs = fsspec.filesystem("memory")
    yield memfs
    memfs.store.clear()


def test_read_csv(cleared_fs):
    from fsspec.implementations.memory import MemoryFile

    cleared_fs.store["test/test.csv"] = MemoryFile(data=text)
    df2 = read_csv("memory://test/test.csv", parse_dates=["dt"])

    tm.assert_frame_equal(df1, df2)


def test_reasonable_error(monkeypatch, cleared_fs):
    from fsspec.registry import known_implementations
    from fsspec import registry

    registry.target.clear()
    with pytest.raises(ValueError) as e:
        read_csv("nosuchprotocol://test/test.csv")
        assert "nosuchprotocol" in str(e.value)
    err_mgs = "test error messgae"
    monkeypatch.setitem(
        known_implementations,
        "couldexist",
        {"class": "unimportable.CouldExist", "err": err_mgs},
    )
    with pytest.raises(ImportError) as e:
        read_csv("couldexist://test/test.csv")
        assert err_mgs in str(e.value)


def test_to_csv(cleared_fs):
    df1.to_csv("memory://test/test.csv", index=True)
    df2 = read_csv("memory://test/test.csv", parse_dates=["dt"], index_col=0)

    tm.assert_frame_equal(df1, df2)


@td.skip_if_no("fastparquet")
def test_to_parquet_new_file(monkeypatch, cleared_fs):
    """Regression test for writing to a not-yet-existent GCS Parquet file."""
    df1.to_parquet(
        "memory://test/test.csv", index=True, engine="fastparquet", compression=None
    )


@td.skip_if_no("s3fs")
def test_from_s3_csv(s3_resource, tips_file):
    tm.assert_equal(read_csv("s3://pandas-test/tips.csv"), read_csv(tips_file))
    # the following are decompressed by pandas, not fsspec
    tm.assert_equal(read_csv("s3://pandas-test/tips.csv.gz"), read_csv(tips_file))
    tm.assert_equal(read_csv("s3://pandas-test/tips.csv.bz2"), read_csv(tips_file))


@pytest.mark.parametrize("protocol", ["s3", "s3a", "s3n"])
@td.skip_if_no("s3fs")
def test_s3_protocols(s3_resource, tips_file, protocol):
    tm.assert_equal(
        read_csv("%s://pandas-test/tips.csv" % protocol), read_csv(tips_file)
    )


@td.skip_if_no("s3fs")
@td.skip_if_no("fastparquet")
def test_s3_parquet(s3_resource):
    fn = "s3://pandas-test/test.parquet"
    df1.to_parquet(fn, index=False, engine="fastparquet", compression=None)
    df2 = read_parquet(fn, engine="fastparquet")
    tm.assert_equal(df1, df2)


@td.skip_if_installed("fsspec")
def test_not_present_exception():
    with pytest.raises(ImportError) as e:
        read_csv("memory://test/test.csv")
        assert "fsspec library is required" in str(e.value)
