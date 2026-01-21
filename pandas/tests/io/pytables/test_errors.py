import datetime
from io import BytesIO
import re
import uuid

import numpy as np
import pytest

from pandas import (
    CategoricalIndex,
    DataFrame,
    HDFStore,
    Index,
    MultiIndex,
    date_range,
    read_hdf,
)

from pandas.io.pytables import (
    Term,
    _maybe_adjust_name,
)

pytestmark = [pytest.mark.single_cpu]


def test_pass_spec_to_storer(temp_hdfstore):
    df = DataFrame(
        1.1 * np.arange(120).reshape((30, 4)),
        columns=Index(list("ABCD"), dtype=object),
        index=Index([f"i-{i}" for i in range(30)], dtype=object),
    )

    temp_hdfstore.put("df", df)
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

    temp_hdfstore.put("frame", df1, format="table")
    msg = re.escape("incompatible kind in col [integer - datetime64[ns]]")
    with pytest.raises(TypeError, match=msg):
        temp_hdfstore.put("frame", df2, format="table", append=True)


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

    temp_hdfstore.put("df", df, format="table")

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


def test_unsuppored_hdf_file_error(datapath):
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


@pytest.mark.parametrize("bad_version", [(1, 2), (1,), [], "12", "123"])
def test_maybe_adjust_name_bad_version_raises(bad_version):
    msg = "Version is incorrect, expected sequence of 3 integers"
    with pytest.raises(ValueError, match=msg):
        _maybe_adjust_name("values_block_0", version=bad_version)
