import contextlib
from datetime import datetime
import io
import os
from pathlib import Path
import struct

import numpy as np
import pytest

from pandas._config import using_string_dtype

from pandas.compat import HAS_PYARROW
from pandas.errors import (
    EmptyDataError,
    Pandas4Warning,
)

import pandas as pd
import pandas._testing as tm

from pandas.io.sas.sas7bdat import (
    SAS7BDATReader,
    _utf8_translation_table,
)
import pandas.io.sas.sas_constants as const

encoding_depr_msg = "The default value of 'encoding' in read_sas is deprecated"


@pytest.fixture
def dirpath(datapath):
    return datapath("io", "sas", "data")


@pytest.fixture(params=[(1, range(1, 16)), (2, [16])])
def data_test_ix(request, dirpath):
    i, test_ix = request.param
    fname = os.path.join(dirpath, f"test_sas7bdat_{i}.csv")
    df = pd.read_csv(fname)
    epoch = datetime(1960, 1, 1)
    t1 = pd.to_timedelta(df["Column4"], unit="D")
    df["Column4"] = (epoch + t1).astype("M8[s]")
    t2 = pd.to_timedelta(df["Column12"], unit="D")
    df["Column12"] = (epoch + t2).astype("M8[s]")
    for k in range(df.shape[1]):
        col = df.iloc[:, k]
        if col.dtype == np.int64:
            df.isetitem(k, df.iloc[:, k].astype(np.float64))
    return df, test_ix


# https://github.com/cython/cython/issues/1720
class TestSAS7BDAT:
    @pytest.mark.slow
    def test_from_file(self, dirpath, data_test_ix):
        expected, test_ix = data_test_ix
        for k in test_ix:
            fname = os.path.join(dirpath, f"test{k}.sas7bdat")
            df = pd.read_sas(fname, encoding="utf-8")
            tm.assert_frame_equal(df, expected)

    @pytest.mark.slow
    def test_from_buffer(self, dirpath, data_test_ix):
        expected, test_ix = data_test_ix
        for k in test_ix:
            fname = os.path.join(dirpath, f"test{k}.sas7bdat")
            with open(fname, "rb") as f:
                byts = f.read()
            buf = io.BytesIO(byts)
            with pd.read_sas(
                buf, format="sas7bdat", iterator=True, encoding="utf-8"
            ) as rdr:
                df = rdr.read()
            tm.assert_frame_equal(df, expected)

    @pytest.mark.slow
    def test_from_iterator(self, dirpath, data_test_ix):
        expected, test_ix = data_test_ix
        for k in test_ix:
            fname = os.path.join(dirpath, f"test{k}.sas7bdat")
            with pd.read_sas(fname, iterator=True, encoding="utf-8") as rdr:
                df = rdr.read(2)
                tm.assert_frame_equal(df, expected.iloc[0:2, :])
                df = rdr.read(3)
                tm.assert_frame_equal(df, expected.iloc[2:5, :])

    @pytest.mark.slow
    def test_path_pathlib(self, dirpath, data_test_ix):
        expected, test_ix = data_test_ix
        for k in test_ix:
            fname = Path(os.path.join(dirpath, f"test{k}.sas7bdat"))
            df = pd.read_sas(fname, encoding="utf-8")
            tm.assert_frame_equal(df, expected)

    @pytest.mark.slow
    @pytest.mark.parametrize("chunksize", (3, 5, 10, 11))
    @pytest.mark.parametrize("k", range(1, 17))
    def test_iterator_loop(self, dirpath, k, chunksize):
        # github #13654
        fname = os.path.join(dirpath, f"test{k}.sas7bdat")
        with pd.read_sas(fname, chunksize=chunksize, encoding="utf-8") as rdr:
            y = 0
            for x in rdr:
                y += x.shape[0]
        assert y == rdr.row_count

    def test_iterator_read_too_much(self, dirpath):
        # github #14734
        fname = os.path.join(dirpath, "test1.sas7bdat")
        with pd.read_sas(
            fname, format="sas7bdat", iterator=True, encoding="utf-8"
        ) as rdr:
            d1 = rdr.read(rdr.row_count + 20)

        with pd.read_sas(fname, iterator=True, encoding="utf-8") as rdr:
            d2 = rdr.read(rdr.row_count + 20)
        tm.assert_frame_equal(d1, d2)


def test_encoding_options(datapath):
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    df1 = pd.read_sas(fname, encoding=None)
    df2 = pd.read_sas(fname, encoding="utf-8")
    for col in df1.columns:
        try:
            df1[col] = df1[col].str.decode("utf-8")
        except AttributeError:
            pass
    tm.assert_frame_equal(df1, df2)

    with contextlib.closing(
        SAS7BDATReader(fname, encoding=None, convert_header_text=False)
    ) as rdr:
        df3 = rdr.read()
    for x, y in zip(df1.columns, df3.columns, strict=True):
        assert x == y.decode()


def test_encoding_default_deprecated(datapath):
    # GH#66470
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    with tm.assert_produces_warning(Pandas4Warning, match=encoding_depr_msg):
        result = pd.read_sas(fname)
    assert result["Column2"].iloc[0] == b"pear"

    for encoding in [None, "infer", "cp1252"]:
        with tm.assert_produces_warning(None):
            pd.read_sas(fname, encoding=encoding)


def test_encoding_default_deprecated_sas7bdatreader(datapath):
    # GH#66470 the reader's own default moved too
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    with tm.assert_produces_warning(Pandas4Warning, match=encoding_depr_msg):
        with contextlib.closing(SAS7BDATReader(fname)) as rdr:
            rdr.read()


def test_encoding_default_not_deprecated_without_text(datapath):
    # GH#66470 nothing changes for files with no character columns and no
    #  non-ascii header text
    fname = datapath("io", "sas", "data", "test_meta2_page.sas7bdat")
    with tm.assert_produces_warning(None):
        pd.read_sas(fname)


def test_encoding_default_deprecated_header_text_only(datapath):
    # GH#66470 the default governs header text too, so a numeric-only file with
    #  non-ascii column names still changes and must warn
    fname = datapath("io", "sas", "data", "datetime.sas7bdat")
    with open(fname, "rb") as fd:
        data = bytearray(fd.read())
    assert b"s" not in SAS7BDATReader(fname, encoding=None)._column_types
    # This file declares cp1251; make a column name non-ascii so latin-1 and
    #  the declared encoding disagree
    ix = data.index(b"DateTimeHi")
    data[ix : ix + 2] = b"\xc4\xe0"

    with tm.assert_produces_warning(Pandas4Warning, match=encoding_depr_msg):
        result = pd.read_sas(io.BytesIO(data), format="sas7bdat")
    assert "ÄàteTimeHi" in result.columns

    expected = pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding="infer")
    assert b"\xc4\xe0".decode("cp1251") + "teTimeHi" in expected.columns


@pytest.mark.parametrize(
    "fname, expected",
    [
        # code 0: no encoding recorded, the most common case
        ("airline", "cp1252"),
        ("cars", "cp1252"),
        # code 28: us-ascii
        ("productsales", "ascii"),
    ],
)
def test_encoding_infer_newly_recognized_codes(datapath, fname, expected):
    # GH#66470 these codes used to be absent from encoding_names, so
    #  encoding="infer" raised LookupError instead of reading the file
    fname = datapath("io", "sas", "data", f"{fname}.sas7bdat")
    with pd.read_sas(fname, encoding="infer", iterator=True) as reader:
        assert reader.inferred_encoding == expected
        assert reader.encoding == expected


def test_encoding_names_are_valid_codecs():
    # GH#66470 a bogus name only surfaces as a LookupError when a user
    #  happens to read a file recording that encoding
    for name in const.encoding_names.values():
        b"".decode(name)


def test_encoding_infer_unrecognized_code(datapath):
    # GH#66470 encoding="infer" used to raise LookupError("unknown encoding: infer")
    #  for files whose header encoding code pandas does not map
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    with open(fname, "rb") as fd:
        data = bytearray(fd.read())
    assert data[const.encoding_offset] not in (250, 251)
    data[const.encoding_offset] = 250  # unassigned by SAS

    with pd.read_sas(
        io.BytesIO(data), format="sas7bdat", encoding="infer", iterator=True
    ) as reader:
        assert reader.inferred_encoding == "unknown (code=250)"
        result = reader.read()
    expected = pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding="latin-1")
    tm.assert_frame_equal(result, expected)


def test_header_text_encoding_contract(datapath):
    # GH#66470 "infer" decodes header text with the encoding the file declares;
    #  encoding=None opts out of decoding and keeps the latin-1 fallback
    fname = datapath("io", "sas", "data", "many_columns.sas7bdat")
    with open(fname, "rb") as fd:
        data = bytearray(fd.read())
    # This file declares utf-8; swap two ascii bytes of a column name for the
    #  two-byte utf-8 encoding of "é" so the name is non-ascii but same length
    ix = data.index(b"ecgrtxt")
    data[ix : ix + 2] = "é".encode()

    result = pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding="infer")
    assert "égrtxt" in result.columns

    # encoding=None opts out of decoding entirely, so header text keeps the
    #  latin-1 fallback, which cannot raise on arbitrary bytes
    result = pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)
    assert "Ã©grtxt" in result.columns


def test_header_text_never_raises_without_encoding(datapath):
    # GH#66470 encoding=None is the documented way to opt out of decoding, so
    #  it must not start raising on header bytes the declared encoding rejects
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    with open(fname, "rb") as fd:
        data = bytearray(fd.read())
    # test1 declares cp1252, which is undefined at 0x81
    data[data.index(b"Column1")] = 0x81

    result = pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)
    assert "\x81olumn1" in result.columns


def test_encoding_infer(datapath):
    fname = datapath("io", "sas", "data", "test1.sas7bdat")

    with pd.read_sas(fname, encoding="infer", iterator=True) as df1_reader:
        # check: is encoding inferred correctly from file
        assert df1_reader.inferred_encoding == "cp1252"
        df1 = df1_reader.read()

    with pd.read_sas(fname, encoding="cp1252", iterator=True) as df2_reader:
        df2 = df2_reader.read()

    # check: reader reads correct information
    tm.assert_frame_equal(df1, df2)


def test_productsales(datapath):
    fname = datapath("io", "sas", "data", "productsales.sas7bdat")
    df = pd.read_sas(fname, encoding="utf-8")
    fname = datapath("io", "sas", "data", "productsales.csv")
    df0 = pd.read_csv(fname, parse_dates=["MONTH"])
    vn = ["ACTUAL", "PREDICT", "QUARTER", "YEAR"]
    df0[vn] = df0[vn].astype(np.float64)

    df0["MONTH"] = df0["MONTH"].astype("M8[s]")
    tm.assert_frame_equal(df, df0)


def test_12659(datapath):
    fname = datapath("io", "sas", "data", "test_12659.sas7bdat")
    df = pd.read_sas(fname)
    fname = datapath("io", "sas", "data", "test_12659.csv")
    df0 = pd.read_csv(fname)
    df0 = df0.astype(np.float64)
    tm.assert_frame_equal(df, df0)


def test_airline(datapath):
    fname = datapath("io", "sas", "data", "airline.sas7bdat")
    df = pd.read_sas(fname)
    fname = datapath("io", "sas", "data", "airline.csv")
    df0 = pd.read_csv(fname)
    df0 = df0.astype(np.float64)
    tm.assert_frame_equal(df, df0)


def test_date_time(datapath):
    # Support of different SAS date/datetime formats (PR #15871)
    fname = datapath("io", "sas", "data", "datetime.sas7bdat")
    df = pd.read_sas(fname)
    fname = datapath("io", "sas", "data", "datetime.csv")
    df0 = pd.read_csv(
        fname, parse_dates=["Date1", "Date2", "DateTime", "DateTimeHi", "Taiw"]
    )
    # GH 19732: Timestamps imported from sas will incur floating point errors
    # See GH#56014 for discussion of the correct "expected" results
    #  We are really just testing that we are "close". This only seems to be
    #  an issue near the implementation bounds.

    df[df.columns[3]] = df.iloc[:, 3].dt.round("us")
    df0["Date1"] = df0["Date1"].astype("M8[s]")
    df0["Date2"] = df0["Date2"].astype("M8[s]")
    df0["DateTime"] = df0["DateTime"].astype("M8[ms]")
    df0["Taiw"] = df0["Taiw"].astype("M8[s]")

    res = df0["DateTimeHi"].astype("M8[us]").dt.round("ms")
    df0["DateTimeHi"] = res.astype("M8[ms]")

    tm.assert_frame_equal(df, df0)


@pytest.mark.parametrize("column", ["WGT", "CYL"])
def test_compact_numerical_values(datapath, column):
    # Regression test for #21616
    fname = datapath("io", "sas", "data", "cars.sas7bdat")
    df = pd.read_sas(fname, encoding="latin-1")
    # The two columns CYL and WGT in cars.sas7bdat have column
    # width < 8 and only contain integral values.
    # Test that pandas doesn't corrupt the numbers by adding
    # decimals.
    result = df[column]
    expected = df[column].round()
    tm.assert_series_equal(result, expected, check_exact=True)


def test_many_columns(datapath):
    # Test for looking for column information in more places (PR #22628)
    fname = datapath("io", "sas", "data", "many_columns.sas7bdat")

    df = pd.read_sas(fname, encoding="latin-1")

    fname = datapath("io", "sas", "data", "many_columns.csv")
    df0 = pd.read_csv(fname, encoding="latin-1")
    tm.assert_frame_equal(df, df0)


def test_inconsistent_number_of_rows(datapath):
    # Regression test for issue #16615. (PR #22628)
    fname = datapath("io", "sas", "data", "load_log.sas7bdat")
    df = pd.read_sas(fname, encoding="latin-1")
    assert len(df) == 2097


def test_zero_variables(datapath):
    # Check if the SAS file has zero variables (PR #18184)
    fname = datapath("io", "sas", "data", "zero_variables.sas7bdat")
    with pytest.raises(EmptyDataError, match="No columns to parse from file"):
        pd.read_sas(fname)


@pytest.mark.parametrize("encoding", [None, "utf8"])
def test_zero_rows(datapath, encoding):
    # GH 18198
    fname = datapath("io", "sas", "data", "zero_rows.sas7bdat")
    result = pd.read_sas(fname, encoding=encoding)
    str_value = b"a" if encoding is None else "a"
    expected = pd.DataFrame([{"char_field": str_value, "num_field": 1.0}]).iloc[:0]
    tm.assert_frame_equal(result, expected)


def test_corrupt_read(datapath):
    # We don't really care about the exact failure, the important thing is
    # that the resource should be cleaned up afterwards (BUG #35566)
    fname = datapath("io", "sas", "data", "corrupt.sas7bdat")
    msg = "'SAS7BDATReader' object has no attribute 'column_count'"
    with pytest.raises(AttributeError, match=msg):
        pd.read_sas(fname)


def test_max_sas_date(datapath):
    # GH 20927, GH 56014
    # NB. max datetime in SAS dataset is 31DEC9999:23:59:59.999
    # SAS uses a modified Gregorian calendar where years divisible by 4000
    # are not leap years, producing a 2-day offset for dates >= year 4000.
    fname = datapath("io", "sas", "data", "max_sas_date.sas7bdat")
    df = pd.read_sas(fname, encoding="iso-8859-1")

    expected = pd.DataFrame(
        {
            "text": ["max", "normal"],
            "dt_as_float": [253717747199.999, 1880323199.999],
            "dt_as_dt": np.array(
                [
                    datetime(9999, 12, 31, 23, 59, 59, 999000),
                    datetime(2019, 8, 1, 23, 59, 59, 999000),
                ],
                dtype="M8[ms]",
            ),
            "date_as_float": [2936547.0, 21762.0],
            "date_as_date": np.array(
                [
                    datetime(9999, 12, 31),
                    datetime(2019, 8, 1),
                ],
                dtype="M8[s]",
            ),
        },
        columns=["text", "dt_as_float", "dt_as_dt", "date_as_float", "date_as_date"],
    )
    tm.assert_frame_equal(df, expected)


def test_max_sas_date_iterator(datapath):
    # GH 20927
    # when called as an iterator, only those chunks with a date > pd.Timestamp.max
    # are returned as datetime.datetime, if this happens that whole chunk is returned
    # as datetime.datetime
    col_order = ["text", "dt_as_float", "dt_as_dt", "date_as_float", "date_as_date"]
    fname = datapath("io", "sas", "data", "max_sas_date.sas7bdat")
    results = []
    for df in pd.read_sas(fname, encoding="iso-8859-1", chunksize=1):
        # GH 19732: Timestamps imported from sas will incur floating point errors
        df.reset_index(inplace=True, drop=True)
        results.append(df)
    expected = [
        pd.DataFrame(
            {
                "text": ["max"],
                "dt_as_float": [253717747199.999],
                "dt_as_dt": np.array(
                    [datetime(9999, 12, 31, 23, 59, 59, 999000)], dtype="M8[ms]"
                ),
                "date_as_float": [2936547.0],
                "date_as_date": np.array([datetime(9999, 12, 31)], dtype="M8[s]"),
            },
            columns=col_order,
        ),
        pd.DataFrame(
            {
                "text": ["normal"],
                "dt_as_float": [1880323199.999],
                "dt_as_dt": np.array(["2019-08-01 23:59:59.999"], dtype="M8[ms]"),
                "date_as_float": [21762.0],
                "date_as_date": np.array(["2019-08-01"], dtype="M8[s]"),
            },
            columns=col_order,
        ),
    ]

    tm.assert_frame_equal(results[0], expected[0])
    tm.assert_frame_equal(results[1], expected[1])


def test_null_date(datapath):
    fname = datapath("io", "sas", "data", "dates_null.sas7bdat")
    df = pd.read_sas(fname, encoding="utf-8")

    expected = pd.DataFrame(
        {
            "datecol": np.array(
                [
                    datetime(9999, 12, 31),
                    np.datetime64("NaT", "ns"),
                ],
                dtype="M8[s]",
            ),
            "datetimecol": np.array(
                [
                    datetime(9999, 12, 31, 23, 59, 59, 999000),
                    np.datetime64("NaT", "ns"),
                ],
                dtype="M8[ms]",
            ),
        },
    )
    tm.assert_frame_equal(df, expected)


def test_meta2_page(datapath):
    # GH 35545
    fname = datapath("io", "sas", "data", "test_meta2_page.sas7bdat")
    df = pd.read_sas(fname)
    assert len(df) == 1000


@pytest.mark.parametrize(
    "test_file, override_offset, override_value, expected_msg",
    [
        ("test2.sas7bdat", 0x10000 + 55229, 0x80 | 0x0F, "Out of bounds"),
        ("test2.sas7bdat", 0x10000 + 55229, 0x10, "unknown control byte"),
        ("test3.sas7bdat", 118170, 184, "Out of bounds"),
    ],
)
def test_rle_rdc_exceptions(
    datapath, test_file, override_offset, override_value, expected_msg
):
    """Errors in RLE/RDC decompression should propagate."""
    with open(datapath("io", "sas", "data", test_file), "rb") as fd:
        data = bytearray(fd.read())
    data[override_offset] = override_value
    with pytest.raises(Exception, match=expected_msg):
        pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)


def test_rdc_overlapping_back_reference(datapath):
    # GH#47339 An RDC back-reference may point closer to the write position than
    # the run it expands is long, so the copy has to propagate byte by byte; a
    # bulk memcpy/memmove would instead read the bytes as they stood before the
    # copy began. No fixture happens to contain such a back-reference, so
    # re-point one: these two command bytes expand 8 bytes from 16 back, and are
    # rewritten to expand 8 bytes from 4 back.
    with open(datapath("io", "sas", "data", "test11.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())
    original = pd.read_sas(io.BytesIO(bytes(data)), format="sas7bdat", encoding=None)

    assert data[120922:120924] == b"\x8d\x00"
    data[120922:120924] = b"\x81\x00"
    result = pd.read_sas(io.BytesIO(bytes(data)), format="sas7bdat", encoding=None)

    # Column11 is bytes 56-63 of the row and is now copied from byte 52, so it
    # repeats the second half of Column9 (bytes 48-55) twice. The fixture is
    # big-endian.
    second_half = struct.pack(">d", original.loc[0, "Column9"])[4:]
    expected = struct.unpack(">d", second_half * 2)[0]
    assert result.loc[0, "Column11"] == expected


@pytest.mark.parametrize(
    "fname, fmt, field_offset, value",
    [
        # Fields of subheader pointers on page 0 of test9 (64-bit), which is
        # walked from Python during _parse_metadata.
        # Offset field, big enough that offset + int_length wraps in uint64.
        ("test9.sas7bdat", "<Q", 68144, (1 << 64) - 1),
        # Length field, big enough to truncate to a negative C int.
        ("test9.sas7bdat", "<Q", 68152, 1 << 31),
        # Length field of a *non-data* (ColumnName) pointer, big enough to
        # truncate to a positive C int that is wrong but in bounds -- unchecked,
        # this fabricates an extra column rather than raising.
        ("test9.sas7bdat", "<Q", 65680, (1 << 32) + 836),
        # Offset field on page 1 of test_meta2_page (32-bit, three pages), which
        # is walked from Parser.read_next_page instead.
        ("test_meta2_page.sas7bdat", "<I", 131096, 0xFFFFFFFF),
    ],
)
def test_corrupt_subheader_pointer(datapath, fname, fmt, field_offset, value):
    # GH#47339 a corrupt subheader pointer must be rejected before it is used to
    # index the page.
    with open(datapath("io", "sas", "data", fname), "rb") as fd:
        data = bytearray(fd.read())
    struct.pack_into(fmt, data, field_offset, value)
    # Rejected while walking the pointers, not later from inside the decompressor
    # once a bad offset or length has already been accepted.
    with pytest.raises(ValueError, match="cached page is too small"):
        pd.read_sas(io.BytesIO(data), format="sas7bdat")


def test_unknown_subheader_signature(datapath):
    # GH#47339 an unrecognized signature on an uncompressed file is reported with
    # the offending bytes. 121376 is where the last live subheader pointer on
    # page 0 of test1.sas7bdat points.
    with open(datapath("io", "sas", "data", "test1.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())
    data[121376:121380] = b"\xde\xad\xbe\xef"
    msg = r"Unknown subheader signature b'\\xde\\xad\\xbe\\xef'"
    with pytest.raises(ValueError, match=msg):
        pd.read_sas(io.BytesIO(data), format="sas7bdat")


@pytest.mark.parametrize(
    "test_file, offset_field, int_len, row_length",
    [
        # byte offset of the first column's data-offset field in the
        # column-attributes subheader, the file's word size, and its row length
        ("test2.sas7bdat", 126580, 4, 809),  # 32-bit, RLE compressed
        ("test1.sas7bdat", 126588, 4, 816),  # 32-bit, uncompressed
        ("test7.sas7bdat", 125532, 8, 816),  # 64-bit
    ],
)
@pytest.mark.parametrize("overshoot", [1, 8, 4096])
def test_column_offset_past_row_raises(
    datapath, test_file, offset_field, int_len, row_length, overshoot
):
    # GH#47339 a column whose data runs past the end of the row would be read
    # out of bounds; the 8-byte numeric fast path did not check.
    with open(datapath("io", "sas", "data", test_file), "rb") as fd:
        data = bytearray(fd.read())
    bad_offset = row_length + overshoot - 8
    data[offset_field : offset_field + int_len] = bad_offset.to_bytes(int_len, "little")
    with pytest.raises(ValueError, match="the file is corrupt"):
        pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)


@pytest.mark.parametrize(
    "test_file, row_length_field_offset, int_len, row_length",
    [
        # byte offset of the row-size subheader's row_length field (it lands on
        # page 1, not the header page, for these fixtures), the file's word
        # size, and its known-good row length
        ("test1.sas7bdat", 130612, 4, 816),  # 32-bit, uncompressed
        ("test2.sas7bdat", 130612, 4, 809),  # 32-bit, RLE compressed
        ("test7.sas7bdat", 130304, 8, 816),  # 64-bit
    ],
)
def test_row_length_exceeds_page_size_raises(
    datapath, test_file, row_length_field_offset, int_len, row_length
):
    # GH#66475 an inflated row_length overflowed the C `int` used in
    # process_byte_array_with_data's offset + length bounds check in sas.pyx,
    # causing a segfault instead of a clean error.
    with open(datapath("io", "sas", "data", test_file), "rb") as fd:
        data = bytearray(fd.read())
    bad_value = 0x7FFFFFF0
    data[row_length_field_offset : row_length_field_offset + int_len] = (
        bad_value.to_bytes(int_len, "little")
    )
    with pytest.raises(ValueError, match=r"row_length \(\d+\) exceeds the page size"):
        pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)


@pytest.mark.parametrize("bad_length", [9, 16, 4096])
def test_numeric_column_longer_than_8_bytes_raises(datapath, bad_length):
    # GH#47339 a numeric column is widened into 8 bytes of the output buffer, so
    # a longer one wrote past the front of it even though it fits in the row.
    with open(datapath("io", "sas", "data", "test1.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())
    # byte offset of the first column's data-length field, which is 4 bytes wide
    data[126592:126596] = bad_length.to_bytes(4, "little")
    with pytest.raises(ValueError, match="the file is corrupt"):
        pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)


@pytest.mark.parametrize("bad_offset", [2**63 - 1, 2**64 - 1])
def test_column_offset_near_int64_max_raises(datapath, bad_offset):
    # GH#47339 a 64-bit file can declare an offset large enough that adding the
    # column length wraps around in int64, which must not slip past validation.
    with open(datapath("io", "sas", "data", "test7.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())
    data[125532:125540] = bad_offset.to_bytes(8, "little")
    with pytest.raises(ValueError, match="the file is corrupt"):
        pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)


def _dates_null_with_late_metadata_page(
    datapath,
    sub_offset,
    sub_length,
    patches,
    total_rows=6,
    tail_page_type=const.page_data_type,
):
    """
    dates_null.sas7bdat, with a metadata page inserted after its data rows.

    dates_null is a 64-bit uncompressed file: a 65536-byte header followed by
    one 65536-byte mix page holding two 16-byte rows. The copy built here claims
    more rows than that page holds and continues with a metadata page carrying a
    patched copy of the subheader at ``sub_offset``, so that the patch is
    processed only once the parser has read past the mix page, and then with a
    data page for the rows it still owes.
    """
    page_length = header_length = 65536
    with open(datapath("io", "sas", "data", "dates_null.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())

    rowsize = header_length + 64728
    # row_count, and the number of rows the mix page is said to hold: enough
    # rows are owed to reach the appended pages, and the mix page's own two are
    # all it hands out.
    struct.pack_into("<q", data, rowsize + 6 * 8, total_rows)
    struct.pack_into("<q", data, rowsize + 15 * 8, 2)

    # page type, block count and subheader count are the first three fields of
    # a page header, and the pointers to its subheaders follow them
    bit_offset = const.page_bit_offset_x64
    pointers = bit_offset + const.subheader_pointers_offset

    meta_page = bytearray(page_length)
    struct.pack_into("<HHH", meta_page, bit_offset, const.page_meta_type, 0, 1)
    sub_position = 4096
    # the page's single subheader pointer: offset then length
    struct.pack_into("<qq", meta_page, pointers, sub_position, sub_length)
    subheader = bytearray(data[header_length + sub_offset :][:sub_length])
    for field_offset, value in patches:
        struct.pack_into("<q", subheader, field_offset, value)
    meta_page[sub_position : sub_position + sub_length] = subheader

    data_page = bytearray(page_length)
    struct.pack_into("<HHH", data_page, bit_offset, tail_page_type, 4, 0)

    return bytes(data) + bytes(meta_page) + bytes(data_page)


@pytest.mark.parametrize(
    "sub_offset, sub_length, patches",
    [
        # The row-size subheader, whose row_length the columns are laid out in.
        #  Shrinking it leaves them hanging off the end of the row -- and, for
        #  the last row on a page, off the page; growing it reads every column
        #  from the wrong bytes. row_count and the mix page's share of it are
        #  kept as page 0 declares them.
        (64728, 808, [(5 * 8, 8), (6 * 8, 6), (15 * 8, 2)]),
        (64728, 808, [(5 * 8, 4096), (6 * 8, 6), (15 * 8, 2)]),
        # The column-size subheader. Above the number of columns the file
        #  describes, the parser walks the offset and length arrays off their
        #  ends; below it, a column the parser did read is dropped from the
        #  frame and comes back all-NaN.
        (64704, 24, [(8, 4096)]),
        (64704, 24, [(8, 1)]),
        # row_count, which each chunk's size is taken from: too many and the
        #  reader keeps going past what the file holds, too few and it stops
        #  short of rows the file does hold. The mix page's own share of it
        #  moves the boundary between the pages the rows come from.
        (64728, 808, [(5 * 8, 16), (6 * 8, 12), (15 * 8, 2)]),
        (64728, 808, [(5 * 8, 16), (6 * 8, 2), (15 * 8, 2)]),
        (64728, 808, [(5 * 8, 16), (6 * 8, 6), (15 * 8, 9999)]),
        # Replaying the column-name and column-attributes subheaders verbatim,
        #  which append rather than overwrite: the file ends up describing each
        #  of its columns twice.
        (63960, 44, []),
        (63900, 60, []),
    ],
    ids=[
        "row_length shrunk",
        "row_length grown",
        "count grown",
        "count shrunk",
        "row_count grown",
        "row_count shrunk",
        "mix page row count",
        "column names repeated",
        "column attributes repeated",
    ],
)
def test_late_metadata_page_changing_layout_raises(
    datapath, sub_offset, sub_length, patches
):
    # GH#47339 the column ranges are validated against row_length when the file
    #  is opened, but a metadata page may follow the first data page and rewrite
    #  the layout the reader was built around.
    data = _dates_null_with_late_metadata_page(
        datapath, sub_offset, sub_length, patches
    )
    with pd.read_sas(
        io.BytesIO(data), format="sas7bdat", chunksize=2, encoding=None
    ) as reader:
        # the chunk that walks onto the metadata page is the one that reports it
        with pytest.raises(ValueError, match="changes the layout it declared"):
            reader.read()


def test_late_metadata_page_reached_mid_chunk_raises(datapath):
    # GH#47339 the parser re-read the mix page's row count from the reader for
    #  every row it took from such a page, so a metadata page reached partway
    #  through a chunk moved the page-advance boundary for the rest of that same
    #  chunk -- reading off the end of the page and failing an internal bounds
    #  assertion, rather than being reported when the chunk ended.
    data = _dates_null_with_late_metadata_page(
        datapath,
        64728,
        808,
        [(5 * 8, 16), (6 * 8, 5000), (15 * 8, 9999)],
        total_rows=5000,
        tail_page_type=const.page_mix_type,
    )
    with pytest.raises(ValueError, match="changes the layout it declared"):
        pd.read_sas(io.BytesIO(data), format="sas7bdat", encoding=None)


def test_mix_page_row_count_too_large_to_reach_raises(datapath):
    # GH#47339 the mix page's row count bounds a row index the parser holds in
    #  an int, so a count that does not fit one has to be rejected rather than
    #  leave that boundary unreachable and the rest of the page handed out as
    #  rows. What the file gets is an error, not which error.
    with open(datapath("io", "sas", "data", "dates_null.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())
    # the row-size subheader's row_count and mix-page row count
    struct.pack_into("<q", data, 65536 + 64728 + 6 * 8, 2**31)
    struct.pack_into("<q", data, 65536 + 64728 + 15 * 8, 2**31)
    with pytest.raises(OverflowError, match="too large to convert"):
        pd.read_sas(
            io.BytesIO(data), format="sas7bdat", chunksize=2, encoding=None
        ).read()


def test_late_metadata_page_repeating_layout_reads(datapath):
    # GH#47339 the check above must not fire on a metadata page that follows the
    #  data pages and restates the layout the file was opened with -- the parser
    #  advances onto the next page as soon as it takes the last row of the one
    #  it is on, so a trailing metadata page is walked routinely.
    data = _dates_null_with_late_metadata_page(
        datapath, 64728, 808, [(5 * 8, 16), (6 * 8, 6), (15 * 8, 2)]
    )
    with pd.read_sas(
        io.BytesIO(data), format="sas7bdat", chunksize=2, encoding=None
    ) as reader:
        assert len(reader.read()) == 2
        assert len(reader.read()) == 2


def _dates_null_with_overrun_data_page(datapath):
    # a data page claiming far more rows than fit on it, so that the parser
    # walks off the end of the page
    with open(datapath("io", "sas", "data", "dates_null.sas7bdat"), "rb") as fd:
        data = bytearray(fd.read())
    struct.pack_into("<q", data, 65536 + 64728 + 6 * 8, 100000)
    struct.pack_into("<q", data, 65536 + 64728 + 15 * 8, 2)
    data_page = bytearray(65536)
    struct.pack_into(
        "<HHH", data_page, const.page_bit_offset_x64, const.page_data_type, 60000, 0
    )
    return bytes(data) + bytes(data_page)


@pytest.mark.parametrize(
    "build, expected, match",
    [
        (
            lambda datapath: _dates_null_with_late_metadata_page(
                datapath, 64704, 24, [(8, 1)]
            ),
            ValueError,
            "changes the layout it declared",
        ),
        (_dates_null_with_overrun_data_page, Exception, "Out of bounds read"),
    ],
    ids=["layout redefined", "page overrun"],
)
def test_chunked_read_closes_the_file_it_rejects(
    datapath, tmp_path, build, expected, match
):
    # GH#47339 GH#56127 read_sas closes the file for the caller only when it
    #  reads the whole of it itself, so a chunked reader that gives up on a file
    #  has to close it however it gave up. Read from a path, not a buffer: a
    #  buffer the caller opened is deliberately left to the caller.
    path = tmp_path / "rejected.sas7bdat"
    path.write_bytes(build(datapath))
    reader = pd.read_sas(path, format="sas7bdat", chunksize=2, encoding=None)
    with pytest.raises(expected, match=match):
        while not reader.read().empty:
            pass
    assert reader.handles.handle.closed


def test_0x40_control_byte(datapath):
    # GH 31243
    fname = datapath("io", "sas", "data", "0x40controlbyte.sas7bdat")
    df = pd.read_sas(fname, encoding="ascii")
    fname = datapath("io", "sas", "data", "0x40controlbyte.csv")
    df0 = pd.read_csv(fname, dtype="str")
    tm.assert_frame_equal(df, df0)


def test_0x00_control_byte(datapath):
    # GH 47099
    fname = datapath("io", "sas", "data", "0x00controlbyte.sas7bdat.bz2")
    with pd.read_sas(fname, chunksize=11_000, encoding=None) as reader:
        df = next(reader)
    assert df.shape == (11_000, 20)


def _fast_string_path_available() -> bool:
    # Mirrors the environment half of SAS7BDATReader._string_mode's gate. Builds
    # without pyarrow compare the object path against itself below, so the test
    # has to know which of the two it is looking at.
    return (
        using_string_dtype()
        and HAS_PYARROW
        and pd.StringDtype(na_value=np.nan).storage == "pyarrow"
    )


# Sized per file to be the largest chunk that still lands boundaries in all four
# positions the parser distinguishes: starting mid-page, spanning a page, falling
# wholly inside a data page, and a short final chunk. Chunking these row counts
# any finer only repeats those four at a few hundred chunks apiece.
_string_path_chunksizes = {
    "test1": 7,
    "test2": 7,
    "test3": 7,
    "test16": 7,
    "load_log": 250,
    "productsales": 50,
}


# cp037 is EBCDIC: the one encoding here whose bytes below 0x80 are not their
# own utf-8, so it is what exercises the parser's non-ascii_identity path.
@pytest.mark.parametrize(
    "encoding",
    ["utf-8", "latin-1", "cp1252", "cp437", "cp1251", "latin9", "ascii", "cp037"],
)
# test2 and test3 are the RLE- and RDC-compressed copies of test1, where the
# parser reads cells out of the decompression buffer rather than the page.
@pytest.mark.parametrize("name", list(_string_path_chunksizes))
@pytest.mark.parametrize("chunked", [False, True])
def test_string_fast_path_matches_object_path(
    datapath, monkeypatch, encoding, name, chunked
):
    # GH#47339 string columns are built as pyarrow arrays directly rather than
    # via an object-dtype array of bytes; the two must agree, including on
    # which files fail to decode.
    fname = datapath("io", "sas", "data", f"{name}.sas7bdat")
    chunksize = _string_path_chunksizes[name] if chunked else None
    if not _fast_string_path_available():
        expected_mode = const.string_mode_object
    elif encoding == "utf-8":
        expected_mode = const.string_mode_utf8
    else:
        expected_mode = const.string_mode_table

    def read():
        with contextlib.closing(
            SAS7BDATReader(fname, encoding=encoding, chunksize=chunksize)
        ) as rdr:
            try:
                # Concatenated rather than first-chunk-only so that the
                # per-chunk buffers, and the short final chunk, are compared.
                frame = pd.concat(list(rdr)) if chunksize else rdr.read()
            except UnicodeDecodeError:
                # Which cell is named differs: the table mode raises from the
                # parse loop in row-major order, the object path column by
                # column. Only that both reject the file is guaranteed.
                return None, True, rdr._str_mode
            return frame, False, rdr._str_mode

    with monkeypatch.context() as patcher:
        patcher.setattr(
            SAS7BDATReader,
            "_string_mode",
            lambda self, ns: (const.string_mode_object, None),
        )
        expected, expected_raised, _ = read()
    result, result_raised, result_mode = read()
    # Without this the comparison above silently degrades to object-vs-object
    # if the gate in _string_mode ever stops opening.
    assert result_mode == expected_mode
    assert result_raised == expected_raised
    if not expected_raised:
        tm.assert_frame_equal(result, expected)


def test_string_fast_path_blank_missing(datapath, monkeypatch):
    # GH#47339 with blank_missing=False a blank cell is an empty string rather
    # than NaN, so the fast path writes a valid zero-length cell for it.
    fname = datapath("io", "sas", "data", "test1.sas7bdat")

    def read():
        with contextlib.closing(
            SAS7BDATReader(fname, encoding="latin-1", blank_missing=False)
        ) as rdr:
            return rdr.read(), rdr._str_mode

    with monkeypatch.context() as patcher:
        patcher.setattr(
            SAS7BDATReader,
            "_string_mode",
            lambda self, ns: (const.string_mode_object, None),
        )
        expected, _ = read()
    result, result_mode = read()
    if _fast_string_path_available():
        assert result_mode == const.string_mode_table
    # Column54 opens with a blank cell, so its buffer starts out empty.
    assert result["Column54"].iloc[0] == ""
    tm.assert_frame_equal(result, expected)


def test_string_fast_path_invalid_utf8(datapath, tmp_path, monkeypatch):
    # GH#47339 pyarrow rejects invalid utf-8 with ArrowInvalid; the fast path has
    # to surface the UnicodeDecodeError the object path raises instead.
    data = bytearray(
        Path(datapath("io", "sas", "data", "test16.sas7bdat")).read_bytes()
    )
    # Break the lead byte of a multi-byte character stored in a string cell.
    data[data.index(b"\xe9\xab\x98")] = 0xFF
    fname = tmp_path / "invalid_utf8.sas7bdat"
    fname.write_bytes(bytes(data))

    def read():
        with contextlib.closing(SAS7BDATReader(fname, encoding="utf-8")) as rdr:
            with pytest.raises(UnicodeDecodeError) as excinfo:
                rdr.read()
            return str(excinfo.value), rdr._str_mode

    with monkeypatch.context() as patcher:
        patcher.setattr(
            SAS7BDATReader,
            "_string_mode",
            lambda self, ns: (const.string_mode_object, None),
        )
        expected, _ = read()
    result, result_mode = read()
    if _fast_string_path_available():
        assert result_mode == const.string_mode_utf8
    assert result == expected


def test_string_fast_path_truncated_file(datapath, tmp_path, monkeypatch):
    # GH#47339 a file whose pages run out before the row count in its header
    # leaves the tail of the offsets buffer unwritten; both paths must report the
    # truncation the same way rather than pyarrow rejecting the offsets.
    fname = datapath("io", "sas", "data", "load_log.sas7bdat")
    with contextlib.closing(SAS7BDATReader(fname, encoding="utf-8")) as rdr:
        page_length = rdr._page_length
    truncated = tmp_path / "truncated.sas7bdat"
    truncated.write_bytes(Path(fname).read_bytes()[:-page_length])

    def read():
        with contextlib.closing(SAS7BDATReader(truncated, encoding="utf-8")) as rdr:
            with pytest.raises(ValueError, match="Length of values") as excinfo:
                rdr.read()
            return str(excinfo.value), rdr._str_mode

    with monkeypatch.context() as patcher:
        patcher.setattr(
            SAS7BDATReader,
            "_string_mode",
            lambda self, ns: (const.string_mode_object, None),
        )
        expected, _ = read()
    result, result_mode = read()
    if _fast_string_path_available():
        assert result_mode == const.string_mode_utf8
    assert result == expected


@pytest.mark.skipif(
    not using_string_dtype(),
    reason="string columns are object dtype, with no storage, under infer_string=0",
)
def test_string_fast_path_respects_string_storage(datapath):
    # GH#47339 the fast path builds a pyarrow-backed column, so it must not be
    # taken when the resolved storage is python.
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    with pd.option_context("mode.string_storage", "python"):
        with contextlib.closing(SAS7BDATReader(fname, encoding="cp1252")) as rdr:
            assert (
                rdr._string_mode(rdr._column_types.count(b"s"))[0]
                == const.string_mode_object
            )
            str_col = rdr.column_names[rdr._column_types.index(b"s")]
            result = rdr.read()
    assert result[str_col].dtype.storage == "python"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"encoding": None},
        {"encoding": "cp1252", "convert_text": False},
        {"encoding": "shift_jis"},
        {"encoding": "big5"},
        {"encoding": "euc_jp"},
    ],
)
def test_string_fast_path_skipped(datapath, kwargs):
    # encoding=None and convert_text=False keep raw bytes, and multi-byte
    # encodings have no per-byte translation to utf-8, so all fall back to the
    # object path.
    fname = datapath("io", "sas", "data", "test1.sas7bdat")
    with contextlib.closing(SAS7BDATReader(fname, **kwargs)) as rdr:
        assert (
            rdr._string_mode(rdr._column_types.count(b"s"))[0]
            == const.string_mode_object
        )


@pytest.mark.parametrize(
    "encoding",
    ["latin-1", "cp1252", "cp437", "cp1251", "latin9", "ascii", "cp874", "cp037"],
)
def test_utf8_translation_table(encoding):
    # GH#47339 the parser translates single-byte encodings through this table
    # instead of decoding cell by cell, so it must agree with the codec.
    table, lengths, ascii_identity = _utf8_translation_table(encoding)
    width = len(table) // 256
    assert ascii_identity == all(
        bytes([value]).decode(encoding, "replace") == chr(value)
        for value in range(0x80)
    )
    for value in range(256):
        try:
            expected = bytes([value]).decode(encoding).encode("utf-8")
        except UnicodeDecodeError:
            assert lengths[value] == const.undefined_byte
            continue
        assert lengths[value] == len(expected)
        assert bytes(table[value * width : value * width + len(expected)]) == expected


@pytest.mark.parametrize(
    "encoding",
    [
        "utf-8",
        "utf-8-sig",
        "utf-16",
        "shift_jis",
        "big5",
        "cp932",
        "euc_jp",
        "raw_unicode_escape",
    ],
)
def test_utf8_translation_table_rejects_context_dependent(encoding):
    # A byte whose meaning depends on what follows it cannot go through the
    # table. Multi-byte lead bytes are undefined in isolation but decode fine
    # once a continuation byte follows; raw_unicode_escape instead maps every
    # byte one to one yet still reads b"\\u0041" as "A".
    assert _utf8_translation_table(encoding) is None
