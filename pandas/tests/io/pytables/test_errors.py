import datetime
from io import BytesIO
import re
import uuid

import numpy as np
import pytest

import pandas as pd
from pandas import (
    CategoricalIndex,
    DataFrame,
    HDFStore,
    Index,
    MultiIndex,
    Series,
    date_range,
    read_hdf,
)
import pandas._testing as tm

from pandas.io.pytables import Term

pytestmark = [pytest.mark.single_cpu]


def test_pass_spec_to_storer(temp_hdfstore):
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )

    temp_hdfstore.put("df", df, track_times=False)
    msg = (
        "cannot pass a column specification when reading a Fixed format "
        "store. this store must be selected in its entirety"
    )
    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.select("df", columns=["A"])
    msg = (
        "cannot pass a where specification when reading from a Fixed "
        "format store. this store must be selected in its entirety"
    )
    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.select("df", where=["columns=A"])


def test_table_index_incompatible_dtypes(temp_hdfstore):
    df1 = DataFrame({"a": [1, 2, 3]})
    df2 = DataFrame(
        {"a": [4, 5, 6]}, index=date_range("1/1/2000", periods=3, unit="ns")
    )

    temp_hdfstore.put("frame", df1, format="table", track_times=False)
    msg = re.escape("incompatible kind in col [integer - datetime64[ns]]")
    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.put("frame", df2, format="table", append=True, track_times=False)


def test_unimplemented_dtypes_table_columns(temp_hdfstore):
    dtypes = [("date", datetime.date(2001, 1, 2))]

    # currently not supported dtypes ####
    for n, f in dtypes:
        df = DataFrame(
            1.1 * np.arange(120).reshape((30, 4)),
            columns=Index(list("ABCD"), dtype=object),
            index=Index([f"i-{i}" for i in range(30)], dtype=object),
        )
        df[n] = f
        msg = re.escape(f"[{n}] is not implemented as a table column")
        with pytest.raises(TypeError, match=msg):
            temp_hdfstore.append(f"df1_{n}", df)


def test_unimplemented_dtypes_table_columns2(temp_hdfstore):
    # frame
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )
    df["obj1"] = "foo"
    df["obj2"] = "bar"
    df["datetime1"] = datetime.date(2001, 1, 2)
    df = df._consolidate()

    # this fails because we have a date in the object block......
    msg = "|".join(
        [
            re.escape(
                "Cannot serialize the column [datetime1]\nbecause its data "
                "contents are not [string] but [date] object dtype"
            ),
            re.escape("[date] is not implemented as a table column"),
        ]
    )
    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.append("df_unimplemented", df)


def test_invalid_terms(temp_hdfstore):
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=date_range("2000-01-01", periods=10, freq="B", unit="ns"),
    )
    df["string"] = "foo"
    df.loc[df.index[0:4], "string"] = "bar"

    temp_hdfstore.put("df", df, format="table", track_times=False)

    # some invalid terms
    msg = re.escape("__init__() missing 1 required positional argument: 'where'")
    with pytest.raises(TypeError, match=msg):
        Term()

    # more invalid
    msg = re.escape(
        "cannot process expression [df.index[3]], "
        "[2000-01-06 00:00:00] is not a valid condition"
    )
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.select("df", "df.index[3]")

    msg = "invalid syntax"
    with pytest.raises(SyntaxError, match=msg):
        temp_hdfstore.select("df", "index>")


def test_invalid_terms_from_docs(temp_h5_path):
    # from the docs
    dfq = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=list("ABCD"),
        index=date_range("20130101", periods=10, unit="ns"),
    )
    dfq.to_hdf(temp_h5_path, key="dfq", format="table", data_columns=True)

    # check ok
    read_hdf(
        temp_h5_path, "dfq", where="index>Timestamp('20130104') & columns=['A', 'B']"
    )
    read_hdf(temp_h5_path, "dfq", where="A>0 or C>0")


def test_invalid_terms_reference(temp_h5_path):
    # catch the invalid reference
    dfq = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)),
        columns=list("ABCD"),
        index=date_range("20130101", periods=10, unit="ns"),
    )
    dfq.to_hdf(temp_h5_path, key="dfq", format="table")

    msg = (
        r"The passed where expression: A>0 or C>0\n\s*"
        r"contains an invalid variable reference\n\s*"
        r"all of the variable references must be a reference to\n\s*"
        r"an axis \(e.g. 'index' or 'columns'\), or a data_column\n\s*"
        r"The currently defined references are: index,columns\n"
    )
    with pytest.raises(ValueError, match=msg):
        read_hdf(temp_h5_path, "dfq", where="A>0 or C>0")


def test_select_too_many_conditions_raises(temp_hdfstore):
    # GH#39752 a where with too many comparisons over indexed columns hits
    # numexpr's input limit; we should raise an actionable message rather than
    # leak the opaque "too many inputs" ValueError.
    # 64 clauses x 2 comparisons each exceeds the limit (NPY_MAXARGS-1, i.e.
    # 31 or 63 depending on the numexpr build).
    n = 64
    df = DataFrame(
        {"A": range(n)},
        index=MultiIndex.from_arrays([["a"] * n, range(n)], names=("la", "lb")),
    )
    temp_hdfstore.put("df", df, format="table", track_times=False)

    where = " | ".join(f"(la == 'a' & lb == {i})" for i in range(n))
    msg = "too many comparisons"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.select("df", where=where)


def test_append_with_diff_col_name_types_raises_value_error(temp_hdfstore):
    df = DataFrame(np.random.default_rng(2).standard_normal((10, 1)))
    df2 = DataFrame({"a": np.random.default_rng(2).standard_normal(10)})
    df3 = DataFrame({(1, 2): np.random.default_rng(2).standard_normal(10)})
    df4 = DataFrame({("1", 2): np.random.default_rng(2).standard_normal(10)})
    df5 = DataFrame({("1", 2, object): np.random.default_rng(2).standard_normal(10)})

    name = "df_diff_valerror"
    temp_hdfstore.append(name, df)

    for d in (df2, df3, df4, df5):
        msg = re.escape(
            "cannot match existing table structure for [0] on appending data"
        )
        with pytest.raises(ValueError, match=msg):
            temp_hdfstore.append(name, d)


def test_invalid_complib(temp_h5_path):
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=list("abcd"),
        columns=list("ABCDE"),
    )
    msg = r"complib only supports \[.*\] compression."
    with pytest.raises(ValueError, match=msg):
        df.to_hdf(temp_h5_path, key="df", complib="foolib")


@pytest.mark.parametrize(
    "idx",
    [
        date_range("2019", freq="D", periods=3, tz="UTC", unit="ns"),
        CategoricalIndex(list("abc")),
    ],
)
def test_to_hdf_multiindex_extension_dtype(idx, temp_h5_path):
    # GH 7775
    mi = MultiIndex.from_arrays([idx, idx])
    df = DataFrame(0, index=mi, columns=["a"])
    with pytest.raises(NotImplementedError, match="Saving a MultiIndex"):
        df.to_hdf(temp_h5_path, key="df")


@pytest.mark.parametrize(
    "values, dtype_match",
    [
        # GH#42070
        (pd.arrays.SparseArray([1.0, 2.0, None, 3.0]), r"Sparse\[float64"),
        # GH#26144
        (pd.array([1, 2, None], dtype="Int32"), "Int32"),
        (pd.array([1.0, None, 3.0], dtype="Float64"), "Float64"),
        (pd.array([True, False, None], dtype="boolean"), "boolean"),
        # GH#38305
        (pd.IntervalIndex.from_arrays([0.5, 1.5], [0.9, 1.9]), r"interval\["),
    ],
)
@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_to_hdf_unsupported_extension_dtype_column(
    values, dtype_match, fmt, temp_h5_path
):
    df = DataFrame({"a": values})
    msg = rf"Cannot store a column with dtype {dtype_match}"
    with pytest.raises(NotImplementedError, match=msg):
        df.to_hdf(temp_h5_path, key="df", format=fmt)


@pytest.mark.parametrize(
    "idx, dtype_match",
    [
        # GH#38305
        (pd.IntervalIndex.from_arrays([0.5, 1.5], [0.9, 1.9]), r"interval\["),
        # GH#42070
        (Index(pd.arrays.SparseArray([1.0, 2.0])), r"Sparse\[float64"),
        # GH#26144
        (Index(pd.array([1, 2], dtype="Int32")), "Int32"),
        (Index(pd.array([1.0, 2.0], dtype="Float64")), "Float64"),
    ],
)
@pytest.mark.parametrize("fmt", ["fixed", "table"])
def test_to_hdf_unsupported_extension_dtype_index(idx, dtype_match, fmt, temp_h5_path):
    # GH#26144, GH#38305, GH#42070
    df = DataFrame({"a": [1, 2]}, index=idx)
    msg = rf"Cannot store an Index with dtype {dtype_match}"
    with pytest.raises(NotImplementedError, match=msg):
        df.to_hdf(temp_h5_path, key="df", format=fmt)


def test_to_hdf_multiindex_level_named_index_raises(temp_hdfstore):
    # GH#6208 a MultiIndex level named 'index' collides with the table
    # format's implicit row index; surface a clear error instead of the
    # confusing reshape failure that used to come from write_data
    mi = MultiIndex.from_tuples(
        [("foo", "one"), ("foo", "two"), ("bar", "one")],
        names=["index", "second"],
    )
    df = DataFrame({"A": [1, 2, 3]}, index=mi)
    msg = "cannot store a MultiIndex with a level named 'index' as a table"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put("df", df, format="table", track_times=False)
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df", df, format="table")
    series = Series([1, 2, 3], index=mi, name="vals")
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put("s", series, format="table", track_times=False)


@pytest.mark.parametrize("data_columns", [["index"], True, ["column_1", "index"]])
def test_to_hdf_data_column_named_index_raises(temp_hdfstore, data_columns):
    # GH#41437 a data column named 'index' collides with the table format's
    # implicit row index; surface a clear error instead of the confusing
    # reshape failure that used to come from write_data
    df = DataFrame({"column_1": [1, 2], "index": [3, 4]})
    df.index.name = "something_else"
    msg = "cannot use a column named 'index' as a data_column"
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.put(
            "df", df, format="table", data_columns=data_columns, track_times=False
        )
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df", df, format="table", data_columns=data_columns)


def test_to_hdf_column_named_index_without_data_columns(temp_h5_path):
    # GH#41437 a column named 'index' is fine as long as it is not a
    # data_column; it round-trips like any other column
    df = DataFrame({"column_1": [1, 2], "index": [3, 4]})
    df.index.name = "something_else"
    df.to_hdf(temp_h5_path, key="df", format="table")
    result = read_hdf(temp_h5_path, "df")
    tm.assert_frame_equal(result, df)


def test_unsupported_hdf_file_error(datapath):
    # GH 9539
    data_path = datapath("io", "data", "legacy_hdf/incompatible_dataset.h5")
    message = (
        r"Dataset\(s\) incompatible with Pandas data types, "
        "not table, or no datasets found in HDF5 file."
    )

    with pytest.raises(ValueError, match=message):
        read_hdf(data_path)


def test_read_hdf_errors(temp_h5_path):
    df = DataFrame(
        np.random.default_rng(2).random((4, 5)),
        index=list("abcd"),
        columns=list("ABCDE"),
    )

    msg = r"File [\S]* does not exist"
    with pytest.raises(OSError, match=msg):
        read_hdf(f"{uuid.uuid4()}.h5", "key")

    df.to_hdf(temp_h5_path, key="df")
    store = HDFStore(temp_h5_path, mode="r")
    store.close()

    msg = "The HDFStore must be open for reading."
    with pytest.raises(OSError, match=msg):
        read_hdf(store, "df")


def test_read_hdf_generic_buffer_errors():
    msg = "Support for generic buffers has not been implemented."
    with pytest.raises(NotImplementedError, match=msg):
        read_hdf(BytesIO(b""), "df")
