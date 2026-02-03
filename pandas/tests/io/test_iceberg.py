"""
Tests for the Apache Iceberg format.

Tests in this file use a simple Iceberg catalog based on SQLite, with the same
data used for Parquet tests (``pandas/tests/io/data/parquet/simple.parquet``).
"""

import collections
import importlib
import pathlib

import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.iceberg import read_iceberg

pytestmark = pytest.mark.single_cpu

pyiceberg = pytest.importorskip("pyiceberg")
pyiceberg_catalog = pytest.importorskip("pyiceberg.catalog")
pq = pytest.importorskip("pyarrow.parquet")

Catalog = collections.namedtuple("Catalog", ["name", "uri", "warehouse"])


@pytest.fixture
def catalog(request, tmp_path):
    # the catalog stores the full path of data files, so the catalog needs to be
    # created dynamically, and not saved in pandas/tests/io/data as other formats
    uri = f"sqlite:///{tmp_path}/catalog.sqlite"
    warehouse = f"file://{tmp_path}"
    catalog_name = request.param if hasattr(request, "param") else None
    catalog = pyiceberg_catalog.load_catalog(
        catalog_name or "default",
        type="sql",
        uri=uri,
        warehouse=warehouse,
    )
    catalog.create_namespace("ns")

    df = pq.read_table(
        pathlib.Path(__file__).parent / "data" / "parquet" / "simple.parquet"
    )
    table = catalog.create_table("ns.my_table", schema=df.schema)
    table.append(df)

    if catalog_name is not None:
        config_path = pathlib.Path.home() / ".pyiceberg.yaml"
        with open(config_path, "w", encoding="utf-8") as f:
            f.write(f"""\
catalog:
  {catalog_name}:
    type: sql
    uri: {uri}
    warehouse: {warehouse}""")

        importlib.reload(pyiceberg_catalog)  # needed to reload the config file

    yield Catalog(name=catalog_name or "default", uri=uri, warehouse=warehouse)

    if catalog_name is not None:
        config_path.unlink()


class TestIceberg:
    def test_read(self, catalog):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("catalog", ["default", "pandas_tests"], indirect=True)
    def test_read_by_catalog_name(self, catalog):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_name=catalog.name,
        )
        tm.assert_frame_equal(result, expected)

    def test_read_with_row_filter(self, catalog):
        expected = pd.DataFrame(
            {
                "A": [2, 3],
                "B": ["foo", "foo"],
            }
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
            row_filter="A > 1",
        )
        tm.assert_frame_equal(result, expected)

    def test_read_with_case_sensitive(self, catalog):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
            }
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
            columns=["a"],
            case_sensitive=False,
        )
        tm.assert_frame_equal(result, expected)

        with pytest.raises(ValueError, match="^Could not find column"):
            read_iceberg(
                "ns.my_table",
                catalog_properties={"uri": catalog.uri},
                columns=["a"],
                case_sensitive=True,
            )

    def test_read_with_limit(self, catalog):
        expected = pd.DataFrame(
            {
                "A": [1, 2],
                "B": ["foo", "foo"],
            }
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
            limit=2,
        )
        tm.assert_frame_equal(result, expected)

    def test_write(self, catalog):
        df = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        df.to_iceberg(
            "ns.new_table",
            catalog_properties={"uri": catalog.uri},
            location=catalog.warehouse,
        )
        result = read_iceberg(
            "ns.new_table",
            catalog_properties={"uri": catalog.uri},
        )
        tm.assert_frame_equal(result, df)

    @pytest.mark.parametrize("catalog", ["default", "pandas_tests"], indirect=True)
    def test_write_by_catalog_name(self, catalog):
        df = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        df.to_iceberg(
            "ns.new_table",
            catalog_name=catalog.name,
        )
        result = read_iceberg(
            "ns.new_table",
            catalog_name=catalog.name,
        )
        tm.assert_frame_equal(result, df)

    def test_write_existing_table_with_append_true(self, catalog):
        original = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
        )
        new = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        expected = pd.concat([original, new], ignore_index=True)
        new.to_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
            location=catalog.warehouse,
            append=True,
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
        )
        tm.assert_frame_equal(result, expected)

    def test_write_existing_table_with_append_false(self, catalog):
        df = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        df.to_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
            location=catalog.warehouse,
            append=False,
        )
        result = read_iceberg(
            "ns.my_table",
            catalog_properties={"uri": catalog.uri},
        )
        tm.assert_frame_equal(result, df)
