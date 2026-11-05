import gzip
from io import (
    BytesIO,
    StringIO,
)

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestSas:
    def test_sas_buffer_format(self):
        # see gh-14947
        b = StringIO("")

        msg = (
            "If this is a buffer object rather than a string "
            "name, you must specify a format string"
        )
        with pytest.raises(ValueError, match=msg):
            pd.read_sas(b)

    def test_sas_read_no_format_or_extension(self, temp_file):
        # see gh-24548
        msg = "unable to infer format of SAS file.+"
        with pytest.raises(ValueError, match=msg):
            pd.read_sas(temp_file)


def test_sas_archive(datapath):
    fname_uncompressed = datapath("io", "sas", "data", "airline.sas7bdat")
    df_uncompressed = pd.read_sas(fname_uncompressed)
    fname_compressed = datapath("io", "sas", "data", "airline.sas7bdat.gz")
    df_compressed = pd.read_sas(fname_compressed, format="sas7bdat")
    tm.assert_frame_equal(df_uncompressed, df_compressed)


@pytest.mark.parametrize(
    "fname, usecols",
    [
        # numeric, string and date columns, in an order other than the file's
        ("productsales.sas7bdat", ["MONTH", "COUNTRY", "ACTUAL"]),
        ("productsales.sas7bdat", ["COUNTRY"]),
        ("productsales.sas7bdat", ["ACTUAL"]),
        # RLE compressed
        ("test2.sas7bdat", ["Column3", "Column1"]),
        # RDC compressed
        ("test3.sas7bdat", ["Column100", "Column2"]),
        ("DEMO_G.xpt", ["RIDAGEYR", "SEQN"]),
        ("DRXFCD_G.xpt", ["DRXFCLD", "DRXFDCD"]),
    ],
)
@pytest.mark.parametrize("encoding", [None, "latin-1", "infer"])
def test_read_sas_usecols(datapath, fname, usecols, encoding):
    # GH#37088
    path = datapath("io", "sas", "data", fname)
    full = pd.read_sas(path, encoding=encoding)

    result = pd.read_sas(path, encoding=encoding, usecols=usecols)

    expected = full[[col for col in full.columns if col in usecols]]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "fname, usecols, index",
    [
        ("productsales.sas7bdat", ["YEAR", "COUNTRY"], "COUNTRY"),
        ("DEMO_G.xpt", ["SEQN", "RIDAGEYR"], "SEQN"),
    ],
)
def test_read_sas_usecols_chunks_index(datapath, fname, usecols, index):
    # GH#37088
    path = datapath("io", "sas", "data", fname)
    expected = pd.read_sas(path, encoding="latin-1", index=index)
    expected = expected[[col for col in expected.columns if col in usecols]]

    with pd.read_sas(
        path, encoding="latin-1", index=index, usecols=usecols, chunksize=100
    ) as reader:
        result = pd.concat(reader)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "box", [list, tuple, set, np.array, pd.Index, pd.Series, dict.fromkeys, iter]
)
@pytest.mark.parametrize("fname", ["productsales.sas7bdat", "DEMO_G.xpt"])
def test_read_sas_usecols_list_like(datapath, fname, box):
    # GH#37088
    path = datapath("io", "sas", "data", fname)
    full = pd.read_sas(path, encoding="latin-1")
    expected = full.iloc[:, [0, 2]]

    usecols = box([full.columns[0], full.columns[2], full.columns[0]])
    result = pd.read_sas(path, encoding="latin-1", usecols=usecols)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("fname", ["productsales.sas7bdat", "DEMO_G.xpt"])
def test_read_sas_usecols_invalid(datapath, fname):
    # GH#37088
    path = datapath("io", "sas", "data", fname)

    with pytest.raises(TypeError, match="usecols must be list-like"):
        pd.read_sas(path, encoding="latin-1", usecols="YEAR")

    msg = r"columns expected but not found: \['nonexistent'\]"
    with pytest.raises(ValueError, match=msg):
        pd.read_sas(path, encoding="latin-1", usecols=["nonexistent"])


@pytest.mark.parametrize("fname", ["productsales.sas7bdat", "DEMO_G.xpt"])
def test_read_sas_usecols_excludes_index(datapath, fname):
    # GH#37088
    path = datapath("io", "sas", "data", fname)
    columns = pd.read_sas(path, encoding="latin-1").columns

    # iterator=True so this only passes if the check happens at open
    with pytest.raises(ValueError, match="usecols must include the index column"):
        pd.read_sas(
            path,
            encoding="latin-1",
            index=columns[0],
            usecols=[columns[1]],
            iterator=True,
        )


@pytest.mark.parametrize("fname", ["productsales.sas7bdat", "DEMO_G.xpt"])
def test_read_sas_usecols_empty(datapath, fname):
    # GH#37088
    path = datapath("io", "sas", "data", fname)
    full = pd.read_sas(path, encoding="latin-1")

    result = pd.read_sas(path, encoding="latin-1", usecols=[])
    assert result.shape == (len(full), 0)
    tm.assert_index_equal(result.columns, pd.Index([], dtype=object))

    with pd.read_sas(path, encoding="latin-1", usecols=[], chunksize=100) as reader:
        assert sum(len(chunk) for chunk in reader) == len(full)


@pytest.mark.parametrize(
    "fname, text_col, numeric_col",
    [
        ("productsales.sas7bdat", "COUNTRY", "ACTUAL"),
        ("DRXFCD_G.xpt", "DRXFCLD", "DRXFDCD"),
    ],
)
def test_read_sas_usecols_encoding_default_warning(
    datapath, fname, text_col, numeric_col
):
    # GH#37088 a text column that is not read does not trigger the warning
    path = datapath("io", "sas", "data", fname)

    with tm.assert_produces_warning(None):
        pd.read_sas(path, usecols=[numeric_col])

    with tm.assert_produces_warning(
        Pandas4Warning, match="default value of 'encoding'"
    ):
        pd.read_sas(path, usecols=[text_col])


@pytest.mark.parametrize(
    "fname, col", [("productsales.sas7bdat", "ACTUAL"), ("DEMO_G.xpt", "SEQN")]
)
def test_read_sas_usecols_only_index_chunks(datapath, fname, col):
    # GH#37088 chunks with rows but no columns must not end the iteration
    path = datapath("io", "sas", "data", fname)
    expected = pd.read_sas(path, encoding="latin-1", index=col)

    with pd.read_sas(
        path, encoding="latin-1", usecols=[col], index=col, chunksize=100
    ) as reader:
        result = pd.concat(reader)

    assert result.shape[1] == 0
    tm.assert_index_equal(result.index, expected.index)


@pytest.mark.parametrize(
    "fname, fmt, text, col",
    [
        ("productsales.sas7bdat", "sas7bdat", b"CANADA", "ACTUAL"),
        ("DRXFCD_G.xpt", "xport", b"Milk, human", "DRXFDCD"),
    ],
)
def test_read_sas_usecols_skips_other_columns(datapath, fname, fmt, text, col):
    # GH#37088 a byte the encoding rejects only raises if its column is parsed
    with open(datapath("io", "sas", "data", fname), "rb") as fh:
        data = bytearray(fh.read())
    data[data.index(text)] = 0xFF

    with pytest.raises(UnicodeDecodeError, match="ascii"):
        pd.read_sas(BytesIO(data), format=fmt, encoding="ascii")

    result = pd.read_sas(BytesIO(data), format=fmt, encoding="ascii", usecols=[col])
    assert list(result.columns) == [col]


@pytest.mark.parametrize("chunksize", [100, 2000])
def test_read_sas_usecols_empty_truncated_file(datapath, chunksize):
    # GH#37088 with no column to mismatch, running out of pages must still raise.
    #  gzip because the size of an uncompressed file already bounds the row count
    with open(datapath("io", "sas", "data", "productsales.sas7bdat"), "rb") as fh:
        data = gzip.compress(fh.read()[: 1024 + 8192 * 9])

    with pytest.raises(ValueError, match="the file is corrupt"):
        with pd.read_sas(
            BytesIO(data),
            format="sas7bdat",
            compression="gzip",
            encoding="latin-1",
            usecols=[],
            chunksize=chunksize,
        ) as reader:
            list(reader)
