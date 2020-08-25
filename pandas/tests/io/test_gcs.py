from io import BytesIO
import os

import numpy as np
import pytest

from pandas import DataFrame, date_range, read_csv
import pandas._testing as tm
from pandas.util import _test_decorators as td


@td.skip_if_no("gcsfs")
def test_read_csv_gcs(monkeypatch):
    from fsspec import AbstractFileSystem, registry

    registry.target.clear()  # noqa  # remove state

    df1 = DataFrame(
        {
            "int": [1, 3],
            "float": [2.0, np.nan],
            "str": ["t", "s"],
            "dt": date_range("2018-06-18", periods=2),
        }
    )

    class MockGCSFileSystem(AbstractFileSystem):
        def open(*args, **kwargs):
            return BytesIO(df1.to_csv(index=False).encode())

    monkeypatch.setattr("gcsfs.GCSFileSystem", MockGCSFileSystem)
    df2 = read_csv("gs://test/test.csv", parse_dates=["dt"])

    tm.assert_frame_equal(df1, df2)


@td.skip_if_no("gcsfs")
def test_to_csv_gcs(monkeypatch):
    from fsspec import AbstractFileSystem, registry

    registry.target.clear()  # noqa  # remove state
    df1 = DataFrame(
        {
            "int": [1, 3],
            "float": [2.0, np.nan],
            "str": ["t", "s"],
            "dt": date_range("2018-06-18", periods=2),
        }
    )
    s = BytesIO()
    s.close = lambda: True

    class MockGCSFileSystem(AbstractFileSystem):
        def open(*args, **kwargs):
            s.seek(0)
            return s

    monkeypatch.setattr("gcsfs.GCSFileSystem", MockGCSFileSystem)
    df1.to_csv("gs://test/test.csv", index=True)

    def mock_get_filepath_or_buffer(*args, **kwargs):
        return BytesIO(df1.to_csv(index=True).encode()), None, None, False

    monkeypatch.setattr(
        "pandas.io.common.get_filepath_or_buffer", mock_get_filepath_or_buffer
    )

    df2 = read_csv("gs://test/test.csv", parse_dates=["dt"], index_col=0)

    tm.assert_frame_equal(df1, df2)


@td.skip_if_no("fastparquet")
@td.skip_if_no("gcsfs")
def test_to_parquet_gcs_new_file(monkeypatch, tmpdir):
    """Regression test for writing to a not-yet-existent GCS Parquet file."""
    from fsspec import AbstractFileSystem, registry

    registry.target.clear()  # noqa  # remove state
    df1 = DataFrame(
        {
            "int": [1, 3],
            "float": [2.0, np.nan],
            "str": ["t", "s"],
            "dt": date_range("2018-06-18", periods=2),
        }
    )

    class MockGCSFileSystem(AbstractFileSystem):
        def open(self, path, mode="r", *args):
            if "w" not in mode:
                raise FileNotFoundError
            return open(os.path.join(tmpdir, "test.parquet"), mode)

    monkeypatch.setattr("gcsfs.GCSFileSystem", MockGCSFileSystem)
    df1.to_parquet(
        "gs://test/test.csv", index=True, engine="fastparquet", compression=None
    )


@td.skip_if_installed("gcsfs")
def test_gcs_not_present_exception():
    with pytest.raises(ImportError) as e:
        read_csv("gs://test/test.csv")
        assert "gcsfs library is required" in str(e.value)
