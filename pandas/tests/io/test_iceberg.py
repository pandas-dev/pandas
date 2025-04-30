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

import pathlib

import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.iceberg import read_iceberg

pyiceberg = pytest.importorskip("pyiceberg")


@pytest.fixture
def catalog_uri(datapath):
    path = datapath("io", "data", "iceberg", "catalog.sqlite")
    return f"sqlite:///{path}"


class TestIceberg:
    def test_read(self, catalog_uri):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        result = read_iceberg(
            "default.simple",
            catalog_properties={"uri": catalog_uri},
        )
        tm.assert_frame_equal(result, expected)

    def test_read_by_catalog_name(self, catalog_uri):
        config_path = pathlib.Path.home() / ".pyiceberg.yaml"
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

    def test_read_with_row_filter(self, catalog_uri):
        expected = pd.DataFrame(
            {
                "A": [2, 3],
                "B": ["foo", "foo"],
            }
        )
        result = read_iceberg(
            "default.simple",
            catalog_properties={"uri": catalog_uri},
            row_filter="A > 1",
        )
        tm.assert_frame_equal(result, expected)

    def test_read_with_case_sensitive(self, catalog_uri):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
            }
        )
        result = read_iceberg(
            "default.simple",
            catalog_properties={"uri": catalog_uri},
            selected_fields=["a"],
            case_sensitive=False,
        )
        tm.assert_frame_equal(result, expected)

        with pytest.raises(ValueError, match="^Could not find column"):
            read_iceberg(
                "default.simple",
                catalog_properties={"uri": catalog_uri},
                selected_fields=["a"],
                case_sensitive=True,
            )

    def test_read_with_limit(self, catalog_uri):
        expected = pd.DataFrame(
            {
                "A": [1, 2],
                "B": ["foo", "foo"],
            }
        )
        result = read_iceberg(
            "default.simple",
            catalog_properties={"uri": catalog_uri},
            limit=2,
        )
        tm.assert_frame_equal(result, expected)
