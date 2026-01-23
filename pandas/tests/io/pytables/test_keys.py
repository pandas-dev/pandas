import numpy as np
import pytest

from pandas import (
    DataFrame,
    HDFStore,
    Index,
    Series,
    date_range,
)

tables = pytest.importorskip("tables")

pytestmark = [pytest.mark.single_cpu]


def test_keys(temp_hdfstore):
    temp_hdfstore["a"] = Series(
        np.arange(10, dtype=np.float64), index=date_range("2020-01-01", periods=10)
    )
    temp_hdfstore["b"] = Series(
        range(10), dtype="float64", index=[f"i_{i}" for i in range(10)]
    )
    temp_hdfstore["c"] = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )

    assert len(temp_hdfstore) == 3
    expected = {"/a", "/b", "/c"}
    assert set(temp_hdfstore.keys()) == expected
    assert set(temp_hdfstore) == expected


def test_non_pandas_keys(temp_h5_path):
    class Table1(tables.IsDescription):
        value1 = tables.Float32Col()

    class Table2(tables.IsDescription):
        value2 = tables.Float32Col()

    class Table3(tables.IsDescription):
        value3 = tables.Float32Col()

    with tables.open_file(temp_h5_path, mode="w") as h5file:
        group = h5file.create_group("/", "group")
        h5file.create_table(group, "table1", Table1, "Table 1")
        h5file.create_table(group, "table2", Table2, "Table 2")
        h5file.create_table(group, "table3", Table3, "Table 3")
    with HDFStore(temp_h5_path) as store:
        assert len(store.keys(include="native")) == 3
        expected = {"/group/table1", "/group/table2", "/group/table3"}
        assert set(store.keys(include="native")) == expected
        assert set(store.keys(include="pandas")) == set()
        for name in expected:
            df = store.get(name)
            assert len(df.columns) == 1


def test_keys_illegal_include_keyword_value(temp_hdfstore):
    with pytest.raises(
        ValueError,
        match="`include` should be either 'pandas' or 'native' but is 'illegal'",
    ):
        temp_hdfstore.keys(include="illegal")


def test_keys_ignore_hdf_softlink(temp_hdfstore):
    # GH 20523
    # Puts a softlink into HDF file and rereads
    df = DataFrame({"A": range(5), "B": range(5)})
    temp_hdfstore.put("df", df)

    assert temp_hdfstore.keys() == ["/df"]

    temp_hdfstore._handle.create_soft_link(temp_hdfstore._handle.root, "symlink", "df")

    # Should ignore the softlink
    assert temp_hdfstore.keys() == ["/df"]
