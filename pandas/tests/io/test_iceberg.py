"""
Tests for the Apache Iceberg format.

Tests in this file use a simple Iceberg catalog based on SQLite, with the same
data used for Parquet tests (``pandas/tests/io/data/parquet/simple.parquet``).

The next code has been used to generate the Iceberg catalog:

```python
import pyarrow.parquet as pq
from pyiceberg.catalog import load_catalog

warehouse = "pandas/tests/io/data/iceberg"
catalog = load_catalog(
    "default",
    **{
        "type": "sql",
        "uri": f"sqlite:///{warehouse}/catalog.sqlite",
        "warehouse": f"file://{warehouse}",
    },
)
catalog.create_namespace("default")
df = pq.read_table("pandas/tests/io/data/parquet/simple.parquet")
table = catalog.create_table("default.simple", schema=df.schema)
table.append(df)
```
"""

from contextlib import contextmanager
import pathlib
import tempfile

import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.iceberg import read_iceberg

pyiceberg = pytest.importorskip("pyiceberg")
pyiceberg_catalog = pytest.importorskip("pyiceberg.catalog")
pq = pytest.importorskip("pyarrow.parquet")


@contextmanager
def create_catalog():
    # the catalog stores the full path of data files, so the catalog needs to be
    # created dynamically, and not saved in pandas/tests/io/data as other formats
    with tempfile.TemporaryDirectory("pandas-iceberg.tmp") as catalog_path:
        uri = f"sqlite:///{catalog_path}/catalog.sqlite"
        catalog = pyiceberg_catalog.load_catalog(
            "default",
            type="sql",
            uri=uri,
            warehouse=f"file://{catalog_path}",
        )
        catalog.create_namespace("default")

        df = pq.read_table(
            pathlib.Path(__file__).parent / "data" / "parquet" / "simple.parquet"
        )
        table = catalog.create_table("default.simple", schema=df.schema)
        table.append(df)

        yield uri


class TestIceberg:
    def test_read(self):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        with create_catalog() as catalog_uri:
            result = read_iceberg(
                "default.simple",
                catalog_properties={"uri": catalog_uri},
            )
        tm.assert_frame_equal(result, expected)

    def test_read_by_catalog_name(self):
        config_path = pathlib.Path.home() / ".pyiceberg.yaml"
        with create_catalog() as catalog_uri:
            with open(config_path, "w") as f:
                f.write(f"""\
    catalog:
      pandas_tests_catalog:
        uri: {catalog_uri}""")
            expected = pd.DataFrame(
                {
                    "A": [1, 2, 3],
                    "B": ["foo", "foo", "foo"],
                }
            )
            result = read_iceberg(
                "default.simple",
                catalog_name="pandas_tests_catalog",
            )
        tm.assert_frame_equal(result, expected)
        # config_path.unlink()

    def test_read_with_row_filter(self):
        expected = pd.DataFrame(
            {
                "A": [2, 3],
                "B": ["foo", "foo"],
            }
        )
        with create_catalog() as catalog_uri:
            result = read_iceberg(
                "default.simple",
                catalog_properties={"uri": catalog_uri},
                row_filter="A > 1",
            )
        tm.assert_frame_equal(result, expected)

    def test_read_with_case_sensitive(self):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
            }
        )
        with create_catalog() as catalog_uri:
            result = read_iceberg(
                "default.simple",
                catalog_properties={"uri": catalog_uri},
                selected_fields=["a"],
                case_sensitive=False,
            )
        tm.assert_frame_equal(result, expected)

        with create_catalog() as catalog_uri:
            with pytest.raises(ValueError, match="^Could not find column"):
                read_iceberg(
                    "default.simple",
                    catalog_properties={"uri": catalog_uri},
                    selected_fields=["a"],
                    case_sensitive=True,
                )

    def test_read_with_limit(self):
        expected = pd.DataFrame(
            {
                "A": [1, 2],
                "B": ["foo", "foo"],
            }
        )
        with create_catalog() as catalog_uri:
            result = read_iceberg(
                "default.simple",
                catalog_properties={"uri": catalog_uri},
                limit=2,
            )
        tm.assert_frame_equal(result, expected)
