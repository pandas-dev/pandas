"""
Tests for the Apache Iceberg format.

Tests in this file use a simple Iceberg catalog based on SQLite, with the same
data used for Parquet tests (``pandas/tests/io/data/parquet/simple.parquet``).
"""

from contextlib import contextmanager
import importlib
import pathlib
import tempfile

import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.iceberg import read_iceberg

pytestmark = [pytest.mark.single_cpu]

# XXX Some tests fail in the CI because an empty .pyiceberg file is found.
# Checking here before importing the pyiceberg module, to provide better
# error message
import os

def check_pyiceberg_yaml_file():
    import strictyaml

    PYICEBERG_HOME = "PYICEBERG_HOME"
    search_dirs = [os.environ.get(PYICEBERG_HOME), os.path.expanduser("~"), os.getcwd()]
    for dir_ in search_dirs:
        path = None
        exists = False
        content = None
        if dir_:
            path = pathlib.Path(dir_) / ".pyiceberg.yaml"
            exists = path.exists()
            if exists:
                with open(path, encoding="utf-8") as f:
                    yml_str = f.read()
                content = strictyaml.load(yml_str).data
                raise RuntimeError(
                    ".pyiceberg.yaml file already exists\n"
                    f"Path: {path}\nContent:\n{content}"
                )

try:
    import strictyaml
except ImportError:
    pass
else:
    check_pyiceberg_yaml_file()
# XXX End checks


pyiceberg = pytest.importorskip("pyiceberg")
pyiceberg_catalog = pytest.importorskip("pyiceberg.catalog")
pq = pytest.importorskip("pyarrow.parquet")


@contextmanager
def create_catalog(catalog_name_in_pyiceberg_config=None):
    # the catalog stores the full path of data files, so the catalog needs to be
    # created dynamically, and not saved in pandas/tests/io/data as other formats
    with tempfile.TemporaryDirectory("-pandas-iceberg.tmp") as catalog_path:
        uri = f"sqlite:///{catalog_path}/catalog.sqlite"
        warehouse = f"file://{catalog_path}"
        catalog = pyiceberg_catalog.load_catalog(
            catalog_name_in_pyiceberg_config or "default",
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

        if catalog_name_in_pyiceberg_config is not None:
            config_path = pathlib.Path.home() / ".pyiceberg.yaml"
            with open(config_path, "w", encoding="utf-8") as f:
                f.write(f"""\
catalog:
  {catalog_name_in_pyiceberg_config}:
    type: sql
    uri: {uri}
    warehouse: {warehouse}""")
        importlib.reload(pyiceberg_catalog)  # needed to reload the config file

        try:
            yield uri
        finally:
            if catalog_name_in_pyiceberg_config is not None:
                config_path.unlink()


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
                "ns.my_table",
                catalog_properties={"uri": catalog_uri},
            )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("catalog_name", ["default", "pandas_tests"])
    def test_read_by_catalog_name(self, catalog_name):
        expected = pd.DataFrame(
            {
                "A": [1, 2, 3],
                "B": ["foo", "foo", "foo"],
            }
        )
        with create_catalog(catalog_name_in_pyiceberg_config=catalog_name):
            result = read_iceberg(
                "ns.my_table",
                catalog_name=catalog_name,
            )
        tm.assert_frame_equal(result, expected)

    def test_read_with_row_filter(self):
        expected = pd.DataFrame(
            {
                "A": [2, 3],
                "B": ["foo", "foo"],
            }
        )
        with create_catalog() as catalog_uri:
            result = read_iceberg(
                "ns.my_table",
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
                "ns.my_table",
                catalog_properties={"uri": catalog_uri},
                selected_fields=["a"],
                case_sensitive=False,
            )
        tm.assert_frame_equal(result, expected)

        with create_catalog() as catalog_uri:
            with pytest.raises(ValueError, match="^Could not find column"):
                read_iceberg(
                    "ns.my_table",
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
                "ns.my_table",
                catalog_properties={"uri": catalog_uri},
                limit=2,
            )
        tm.assert_frame_equal(result, expected)
