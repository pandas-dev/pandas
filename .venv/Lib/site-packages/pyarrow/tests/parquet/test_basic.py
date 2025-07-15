# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import os
from collections import OrderedDict
import io
import warnings
from shutil import copytree
from decimal import Decimal

import pytest

import pyarrow as pa
from pyarrow import fs
from pyarrow.tests import util
from pyarrow.tests.parquet.common import (_check_roundtrip, _roundtrip_table,
                                          _test_dataframe)

try:
    import pyarrow.parquet as pq
    from pyarrow.tests.parquet.common import _read_table, _write_table
except ImportError:
    pq = None


try:
    import pandas as pd
    import pandas.testing as tm

    from pyarrow.tests.pandas_examples import dataframe_with_lists
    from pyarrow.tests.parquet.common import alltypes_sample
except ImportError:
    pd = tm = None

try:
    import numpy as np
except ImportError:
    np = None

# Marks all of the tests in this module
# Ignore these with pytest ... -m 'not parquet'
pytestmark = pytest.mark.parquet


def test_parquet_invalid_version(tempdir):
    table = pa.table({'a': [1, 2, 3]})
    with pytest.raises(ValueError, match="Unsupported Parquet format version"):
        _write_table(table, tempdir / 'test_version.parquet', version="2.2")
    with pytest.raises(ValueError, match="Unsupported Parquet data page " +
                       "version"):
        _write_table(table, tempdir / 'test_version.parquet',
                     data_page_version="2.2")


def test_set_data_page_size():
    arr = pa.array([1, 2, 3] * 100000)
    t = pa.Table.from_arrays([arr], names=['f0'])

    # 128K, 512K
    page_sizes = [2 << 16, 2 << 18]
    for target_page_size in page_sizes:
        _check_roundtrip(t, data_page_size=target_page_size)


@pytest.mark.pandas
def test_set_write_batch_size():
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)

    _check_roundtrip(
        table, data_page_size=10, write_batch_size=1, version='2.4'
    )


@pytest.mark.pandas
def test_set_dictionary_pagesize_limit():
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)

    _check_roundtrip(table, dictionary_pagesize_limit=1,
                     data_page_size=10, version='2.4')

    with pytest.raises(TypeError):
        _check_roundtrip(table, dictionary_pagesize_limit="a",
                         data_page_size=10, version='2.4')


@pytest.mark.pandas
def test_chunked_table_write():
    # ARROW-232
    tables = []
    batch = pa.RecordBatch.from_pandas(alltypes_sample(size=10))
    tables.append(pa.Table.from_batches([batch] * 3))
    df, _ = dataframe_with_lists()
    batch = pa.RecordBatch.from_pandas(df)
    tables.append(pa.Table.from_batches([batch] * 3))

    for data_page_version in ['1.0', '2.0']:
        for use_dictionary in [True, False]:
            for table in tables:
                _check_roundtrip(
                    table, version='2.6',
                    data_page_version=data_page_version,
                    use_dictionary=use_dictionary)


@pytest.mark.pandas
def test_memory_map(tempdir):
    df = alltypes_sample(size=10)

    table = pa.Table.from_pandas(df)
    _check_roundtrip(table, read_table_kwargs={'memory_map': True},
                     version='2.6')

    filename = str(tempdir / 'tmp_file')
    with open(filename, 'wb') as f:
        _write_table(table, f, version='2.6')
    table_read = pq.read_pandas(filename, memory_map=True)
    assert table_read.equals(table)


@pytest.mark.pandas
def test_enable_buffered_stream(tempdir):
    df = alltypes_sample(size=10)

    table = pa.Table.from_pandas(df)
    _check_roundtrip(table, read_table_kwargs={'buffer_size': 1025},
                     version='2.6')

    filename = str(tempdir / 'tmp_file')
    with open(filename, 'wb') as f:
        _write_table(table, f, version='2.6')
    table_read = pq.read_pandas(filename, buffer_size=4096)
    assert table_read.equals(table)


def test_special_chars_filename(tempdir):
    table = pa.Table.from_arrays([pa.array([42])], ["ints"])
    filename = "foo # bar"
    path = tempdir / filename
    assert not path.exists()
    _write_table(table, str(path))
    assert path.exists()
    table_read = _read_table(str(path))
    assert table_read.equals(table)


def test_invalid_source():
    # Test that we provide an helpful error message pointing out
    # that None wasn't expected when trying to open a Parquet None file.
    with pytest.raises(TypeError, match="None"):
        pq.read_table(None)

    with pytest.raises(TypeError, match="None"):
        pq.ParquetFile(None)


@pytest.mark.slow
def test_file_with_over_int16_max_row_groups():
    # PARQUET-1857: Parquet encryption support introduced a INT16_MAX upper
    # limit on the number of row groups, but this limit only impacts files with
    # encrypted row group metadata because of the int16 row group ordinal used
    # in the Parquet Thrift metadata. Unencrypted files are not impacted, so
    # this test checks that it works (even if it isn't a good idea)
    t = pa.table([list(range(40000))], names=['f0'])
    _check_roundtrip(t, row_group_size=1)


@pytest.mark.pandas
def test_empty_table_roundtrip():
    df = alltypes_sample(size=10)

    # Create a non-empty table to infer the types correctly, then slice to 0
    table = pa.Table.from_pandas(df)
    table = pa.Table.from_arrays(
        [col.chunk(0)[:0] for col in table.itercolumns()],
        names=table.schema.names)

    assert table.schema.field('null').type == pa.null()
    assert table.schema.field('null_list').type == pa.list_(pa.null())
    _check_roundtrip(
        table, version='2.6')


@pytest.mark.pandas
def test_empty_table_no_columns():
    df = pd.DataFrame()
    empty = pa.Table.from_pandas(df, preserve_index=False)
    _check_roundtrip(empty)


def test_write_nested_zero_length_array_chunk_failure():
    # Bug report in ARROW-3792
    cols = OrderedDict(
        int32=pa.int32(),
        list_string=pa.list_(pa.string())
    )
    data = [[], [OrderedDict(int32=1, list_string=('G',)), ]]

    # This produces a table with a column like
    # <Column name='list_string' type=ListType(list<item: string>)>
    # [
    #   [],
    #   [
    #     [
    #       "G"
    #     ]
    #   ]
    # ]
    #
    # Each column is a ChunkedArray with 2 elements
    my_arrays = [pa.array(batch, type=pa.struct(cols)).flatten()
                 for batch in data]
    my_batches = [pa.RecordBatch.from_arrays(batch, schema=pa.schema(cols))
                  for batch in my_arrays]
    tbl = pa.Table.from_batches(my_batches, pa.schema(cols))
    _check_roundtrip(tbl)


@pytest.mark.pandas
def test_multiple_path_types(tempdir):
    # Test compatibility with PEP 519 path-like objects
    path = tempdir / 'zzz.parquet'
    df = pd.DataFrame({'x': np.arange(10, dtype=np.int64)})
    _write_table(df, path)
    table_read = _read_table(path)
    df_read = table_read.to_pandas()
    tm.assert_frame_equal(df, df_read)

    # Test compatibility with plain string paths
    path = str(tempdir) + 'zzz.parquet'
    df = pd.DataFrame({'x': np.arange(10, dtype=np.int64)})
    _write_table(df, path)
    table_read = _read_table(path)
    df_read = table_read.to_pandas()
    tm.assert_frame_equal(df, df_read)


def test_fspath(tempdir):
    # ARROW-12472 support __fspath__ objects without using str()
    path = tempdir / "test.parquet"
    table = pa.table({"a": [1, 2, 3]})
    _write_table(table, path)

    fs_protocol_obj = util.FSProtocolClass(path)

    result = _read_table(fs_protocol_obj)
    assert result.equals(table)

    # combined with non-local filesystem raises
    with pytest.raises(TypeError):
        _read_table(fs_protocol_obj, filesystem=fs.FileSystem())


@pytest.mark.parametrize("filesystem", [
    None, fs.LocalFileSystem()
])
@pytest.mark.parametrize("name", ("data.parquet", "例.parquet"))
def test_relative_paths(tempdir, filesystem, name):
    # reading and writing from relative paths
    table = pa.table({"a": [1, 2, 3]})
    path = tempdir / name

    # reading
    pq.write_table(table, str(path))
    with util.change_cwd(tempdir):
        result = pq.read_table(name, filesystem=filesystem)
    assert result.equals(table)

    path.unlink()
    assert not path.exists()

    # writing
    with util.change_cwd(tempdir):
        pq.write_table(table, name, filesystem=filesystem)
    result = pq.read_table(path)
    assert result.equals(table)


def test_read_non_existing_file():
    # ensure we have a proper error message
    with pytest.raises(FileNotFoundError):
        pq.read_table('i-am-not-existing.parquet')


def test_file_error_python_exception():
    class BogusFile(io.BytesIO):
        def read(self, *args):
            raise ZeroDivisionError("zorglub")

        def seek(self, *args):
            raise ZeroDivisionError("zorglub")

    # ensure the Python exception is restored
    with pytest.raises(ZeroDivisionError, match="zorglub"):
        pq.read_table(BogusFile(b""))


def test_parquet_read_from_buffer(tempdir):
    # reading from a buffer from python's open()
    table = pa.table({"a": [1, 2, 3]})
    pq.write_table(table, str(tempdir / "data.parquet"))

    with open(str(tempdir / "data.parquet"), "rb") as f:
        result = pq.read_table(f)
    assert result.equals(table)

    with open(str(tempdir / "data.parquet"), "rb") as f:
        result = pq.read_table(pa.PythonFile(f))
    assert result.equals(table)


def test_byte_stream_split():
    # This is only a smoke test.
    arr_float = pa.array(list(map(float, range(100))))
    arr_int = pa.array(list(map(int, range(100))))
    arr_bool = pa.array([True, False] * 50)
    data_float = [arr_float, arr_float]
    table = pa.Table.from_arrays(data_float, names=['a', 'b'])

    # Check with byte_stream_split for both columns.
    _check_roundtrip(table, expected=table, compression="gzip",
                     use_dictionary=False, use_byte_stream_split=True)

    # Check with byte_stream_split for column 'b' and dictionary
    # for column 'a'.
    _check_roundtrip(table, expected=table, compression="gzip",
                     use_dictionary=['a'],
                     use_byte_stream_split=['b'])

    # Check with a collision for both columns.
    _check_roundtrip(table, expected=table, compression="gzip",
                     use_dictionary=['a', 'b'],
                     use_byte_stream_split=['a', 'b'])

    # Check with mixed column types.
    mixed_table = pa.Table.from_arrays([arr_float, arr_float, arr_int, arr_int],
                                       names=['a', 'b', 'c', 'd'])
    _check_roundtrip(mixed_table, expected=mixed_table,
                     use_dictionary=['b', 'd'],
                     use_byte_stream_split=['a', 'c'])

    # Try to use the wrong data type with the byte_stream_split encoding.
    # This should throw an exception.
    table = pa.Table.from_arrays([arr_bool], names=['tmp'])
    with pytest.raises(IOError, match='BYTE_STREAM_SPLIT only supports'):
        _check_roundtrip(table, expected=table, use_byte_stream_split=True,
                         use_dictionary=False)


def test_store_decimal_as_integer(tempdir):
    arr_decimal_1_9 = pa.array(list(map(Decimal, range(100))),
                               type=pa.decimal128(5, 2))
    arr_decimal_10_18 = pa.array(list(map(Decimal, range(100))),
                                 type=pa.decimal128(16, 9))
    arr_decimal_gt18 = pa.array(list(map(Decimal, range(100))),
                                type=pa.decimal128(22, 2))
    arr_bool = pa.array([True, False] * 50)
    data_decimal = [arr_decimal_1_9, arr_decimal_10_18, arr_decimal_gt18]
    table = pa.Table.from_arrays(data_decimal, names=['a', 'b', 'c'])

    # Check with store_decimal_as_integer.
    _check_roundtrip(table,
                     expected=table,
                     compression="gzip",
                     use_dictionary=False,
                     store_decimal_as_integer=True)

    # Check physical type in parquet schema
    pqtestfile_path = os.path.join(tempdir, 'test.parquet')
    pq.write_table(table, pqtestfile_path,
                   compression="gzip",
                   use_dictionary=False,
                   store_decimal_as_integer=True)

    pqtestfile = pq.ParquetFile(pqtestfile_path)
    pqcol_decimal_1_9 = pqtestfile.schema.column(0)
    pqcol_decimal_10_18 = pqtestfile.schema.column(1)

    assert pqcol_decimal_1_9.physical_type == 'INT32'
    assert pqcol_decimal_10_18.physical_type == 'INT64'

    # Check with store_decimal_as_integer and delta-int encoding.
    # DELTA_BINARY_PACKED requires parquet physical type to be INT64 or INT32
    _check_roundtrip(table,
                     expected=table,
                     compression="gzip",
                     use_dictionary=False,
                     store_decimal_as_integer=True,
                     column_encoding={
                         'a': 'DELTA_BINARY_PACKED',
                         'b': 'DELTA_BINARY_PACKED'
                     })

    # Check with mixed column types.
    mixed_table = pa.Table.from_arrays(
        [arr_decimal_1_9, arr_decimal_10_18, arr_decimal_gt18, arr_bool],
        names=['a', 'b', 'c', 'd'])
    _check_roundtrip(mixed_table,
                     expected=mixed_table,
                     use_dictionary=False,
                     store_decimal_as_integer=True)


def test_column_encoding():
    arr_float = pa.array(list(map(float, range(100))))
    arr_int = pa.array(list(map(int, range(100))))
    arr_bin = pa.array([str(x) for x in range(100)], type=pa.binary())
    arr_flba = pa.array(
        [str(x).zfill(10) for x in range(100)], type=pa.binary(10))
    arr_bool = pa.array([False, True, False, False] * 25)
    mixed_table = pa.Table.from_arrays(
        [arr_float, arr_int, arr_bin, arr_flba, arr_bool],
        names=['a', 'b', 'c', 'd', 'e'])

    # Check "BYTE_STREAM_SPLIT" for columns 'a', 'b', 'd'
    # and "PLAIN" column_encoding for column 'c'.
    _check_roundtrip(mixed_table, expected=mixed_table, use_dictionary=False,
                     column_encoding={'a': "BYTE_STREAM_SPLIT",
                                      'b': "BYTE_STREAM_SPLIT",
                                      'c': "PLAIN",
                                      'd': "BYTE_STREAM_SPLIT"})

    # Check "PLAIN" for all columns.
    _check_roundtrip(mixed_table, expected=mixed_table,
                     use_dictionary=False,
                     column_encoding="PLAIN")

    # Check "DELTA_BINARY_PACKED" for integer columns.
    _check_roundtrip(mixed_table, expected=mixed_table,
                     use_dictionary=False,
                     column_encoding={'a': "PLAIN",
                                      'b': "DELTA_BINARY_PACKED",
                                      'c': "PLAIN"})

    # Check "DELTA_LENGTH_BYTE_ARRAY" for byte columns.
    _check_roundtrip(mixed_table, expected=mixed_table,
                     use_dictionary=False,
                     column_encoding={'a': "PLAIN",
                                      'b': "DELTA_BINARY_PACKED",
                                      'c': "DELTA_LENGTH_BYTE_ARRAY"})

    # Check "DELTA_BYTE_ARRAY" for byte columns.
    _check_roundtrip(mixed_table, expected=mixed_table,
                     use_dictionary=False,
                     column_encoding={'a': "PLAIN",
                                      'b': "DELTA_BINARY_PACKED",
                                      'c': "DELTA_BYTE_ARRAY",
                                      'd': "DELTA_BYTE_ARRAY"})

    # Check "RLE" for boolean columns.
    _check_roundtrip(mixed_table, expected=mixed_table,
                     use_dictionary=False,
                     column_encoding={'e': "RLE"})

    # Try to pass "BYTE_STREAM_SPLIT" column encoding for boolean column 'e'.
    # This should throw an error as it is does not support BOOLEAN.
    with pytest.raises(IOError,
                       match="BYTE_STREAM_SPLIT only supports"):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         column_encoding={'a': "PLAIN",
                                          'c': "PLAIN",
                                          'e': "BYTE_STREAM_SPLIT"})

    # Try to pass use "DELTA_BINARY_PACKED" encoding on float column.
    # This should throw an error as only integers are supported.
    with pytest.raises(OSError,
                       match="DELTA_BINARY_PACKED encoder only supports"):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         column_encoding={'a': "DELTA_BINARY_PACKED",
                                          'b': "PLAIN",
                                          'c': "PLAIN"})

    # Try to pass "RLE_DICTIONARY".
    # This should throw an error as dictionary encoding is already used by
    # default and not supported to be specified as "fallback" encoding
    with pytest.raises(ValueError,
                       match="'RLE_DICTIONARY' is already used by default"):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         column_encoding="RLE_DICTIONARY")

    # Try to pass unsupported encoding.
    with pytest.raises(ValueError,
                       match="Unsupported column encoding: 'MADE_UP_ENCODING'"):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         column_encoding={'a': "MADE_UP_ENCODING"})

    # Try to pass column_encoding and use_dictionary.
    # This should throw an error.
    with pytest.raises(ValueError):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=['b'],
                         column_encoding={'b': "PLAIN"})

    # Try to pass column_encoding and use_dictionary=True (default value).
    # This should throw an error.
    with pytest.raises(ValueError):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         column_encoding={'b': "PLAIN"})

    # Try to pass column_encoding and use_byte_stream_split on same column.
    # This should throw an error.
    with pytest.raises(ValueError):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         use_byte_stream_split=['a'],
                         column_encoding={'a': "RLE",
                                          'b': "BYTE_STREAM_SPLIT",
                                          'c': "PLAIN"})

    # Try to pass column_encoding and use_byte_stream_split=True.
    # This should throw an error.
    with pytest.raises(ValueError):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         use_byte_stream_split=True,
                         column_encoding={'a': "RLE",
                                          'b': "BYTE_STREAM_SPLIT",
                                          'c': "PLAIN"})

    # Try to pass column_encoding=True.
    # This should throw an error.
    with pytest.raises(TypeError):
        _check_roundtrip(mixed_table, expected=mixed_table,
                         use_dictionary=False,
                         column_encoding=True)


def test_compression_level():
    arr = pa.array(list(map(int, range(1000))))
    data = [arr, arr]
    table = pa.Table.from_arrays(data, names=['a', 'b'])

    # Check one compression level.
    _check_roundtrip(table, expected=table, compression="gzip",
                     compression_level=1)

    # Check another one to make sure that compression_level=1 does not
    # coincide with the default one in Arrow.
    _check_roundtrip(table, expected=table, compression="gzip",
                     compression_level=5)

    # Check that the user can provide a compression per column
    _check_roundtrip(table, expected=table,
                     compression={'a': "gzip", 'b': "snappy"})

    # Check that the user can provide a compression level per column
    _check_roundtrip(table, expected=table, compression="gzip",
                     compression_level={'a': 2, 'b': 3})

    # Check if both LZ4 compressors are working
    # (level < 3 -> fast, level >= 3 -> HC)
    _check_roundtrip(table, expected=table, compression="lz4",
                     compression_level=1)

    _check_roundtrip(table, expected=table, compression="lz4",
                     compression_level=9)

    # Check that specifying a compression level for a codec which does allow
    # specifying one, results into an error.
    # Uncompressed, snappy and lzo do not support specifying a compression
    # level.
    # GZIP (zlib) allows for specifying a compression level but as of up
    # to version 1.2.11 the valid range is [-1, 9].
    invalid_combinations = [("snappy", 4), ("gzip", -1337),
                            ("None", 444), ("lzo", 14)]
    buf = io.BytesIO()
    for (codec, level) in invalid_combinations:
        with pytest.raises((ValueError, OSError)):
            _write_table(table, buf, compression=codec,
                         compression_level=level)


def test_sanitized_spark_field_names():
    a0 = pa.array([0, 1, 2, 3, 4])
    name = 'prohib; ,\t{}'
    table = pa.Table.from_arrays([a0], [name])

    result = _roundtrip_table(table, write_table_kwargs={'flavor': 'spark'})

    expected_name = 'prohib______'
    assert result.schema[0].name == expected_name


@pytest.mark.pandas
def test_multithreaded_read():
    df = alltypes_sample(size=10000)

    table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(table, buf, compression='SNAPPY', version='2.6')

    buf.seek(0)
    table1 = _read_table(buf, use_threads=True)

    buf.seek(0)
    table2 = _read_table(buf, use_threads=False)

    assert table1.equals(table2)


@pytest.mark.pandas
def test_min_chunksize():
    data = pd.DataFrame([np.arange(4)], columns=['A', 'B', 'C', 'D'])
    table = pa.Table.from_pandas(data.reset_index())

    buf = io.BytesIO()
    _write_table(table, buf, chunk_size=-1)

    buf.seek(0)
    result = _read_table(buf)

    assert result.equals(table)

    with pytest.raises(ValueError):
        _write_table(table, buf, chunk_size=0)


@pytest.mark.pandas
def test_write_error_deletes_incomplete_file(tempdir):
    # ARROW-1285
    df = pd.DataFrame({'a': list('abc'),
                       'b': list(range(1, 4)),
                       'c': np.arange(3, 6).astype('u1'),
                       'd': np.arange(4.0, 7.0, dtype='float64'),
                       'e': [True, False, True],
                       'f': pd.Categorical(list('abc')),
                       'g': pd.date_range('20130101', periods=3),
                       'h': pd.date_range('20130101', periods=3,
                                          tz='US/Eastern'),
                       'i': pd.date_range('20130101', periods=3, freq='ns')})

    pdf = pa.Table.from_pandas(df)

    filename = tempdir / 'tmp_file'
    try:
        # Test relies on writing nanoseconds to raise an error
        # true for Parquet 2.4
        _write_table(pdf, filename, version="2.4")
    except pa.ArrowException:
        pass

    assert not filename.exists()


def test_read_non_existent_file(tempdir):
    path = 'nonexistent-file.parquet'
    try:
        pq.read_table(path)
    except Exception as e:
        assert path in e.args[0]


def test_read_table_doesnt_warn(datadir):
    with warnings.catch_warnings():
        warnings.simplefilter(action="error")
        pq.read_table(datadir / 'v0.7.1.parquet')


@pytest.mark.pandas
def test_zlib_compression_bug():
    # ARROW-3514: "zlib deflate failed, output buffer too small"
    table = pa.Table.from_arrays([pa.array(['abc', 'def'])], ['some_col'])
    f = io.BytesIO()
    pq.write_table(table, f, compression='gzip')

    f.seek(0)
    roundtrip = pq.read_table(f)
    tm.assert_frame_equal(roundtrip.to_pandas(), table.to_pandas())


def test_parquet_file_too_small(tempdir):
    path = str(tempdir / "test.parquet")
    # TODO(dataset) with datasets API it raises OSError instead
    with pytest.raises((pa.ArrowInvalid, OSError),
                       match='size is 0 bytes'):
        with open(path, 'wb') as f:
            pass
        pq.read_table(path)

    with pytest.raises((pa.ArrowInvalid, OSError),
                       match='size is 4 bytes'):
        with open(path, 'wb') as f:
            f.write(b'ffff')
        pq.read_table(path)


@pytest.mark.pandas
@pytest.mark.fastparquet
@pytest.mark.filterwarnings("ignore:RangeIndex:FutureWarning")
@pytest.mark.filterwarnings("ignore:tostring:DeprecationWarning:fastparquet")
def test_fastparquet_cross_compatibility(tempdir):
    fp = pytest.importorskip('fastparquet')

    df = pd.DataFrame(
        {
            "a": list("abc"),
            "b": list(range(1, 4)),
            "c": np.arange(4.0, 7.0, dtype="float64"),
            "d": [True, False, True],
            "e": pd.date_range("20130101", periods=3),
            "f": pd.Categorical(["a", "b", "a"]),
            # fastparquet writes list as BYTE_ARRAY JSON, so no roundtrip
            # "g": [[1, 2], None, [1, 2, 3]],
        }
    )
    table = pa.table(df)

    # Arrow -> fastparquet
    file_arrow = str(tempdir / "cross_compat_arrow.parquet")
    pq.write_table(table, file_arrow, compression=None)

    fp_file = fp.ParquetFile(file_arrow)
    df_fp = fp_file.to_pandas()
    tm.assert_frame_equal(df, df_fp)

    # Fastparquet -> arrow
    file_fastparquet = str(tempdir / "cross_compat_fastparquet.parquet")
    fp.write(file_fastparquet, df)

    table_fp = pq.read_pandas(file_fastparquet)
    # for fastparquet written file, categoricals comes back as strings
    # (no arrow schema in parquet metadata)
    df['f'] = df['f'].astype(object)
    tm.assert_frame_equal(table_fp.to_pandas(), df)


@pytest.mark.parametrize('array_factory', [
    lambda: pa.array([0, None] * 10),
    lambda: pa.array([0, None] * 10).dictionary_encode(),
    lambda: pa.array(["", None] * 10),
    lambda: pa.array(["", None] * 10).dictionary_encode(),
])
@pytest.mark.parametrize('read_dictionary', [False, True])
def test_buffer_contents(
        array_factory, read_dictionary
):
    # Test that null values are deterministically initialized to zero
    # after a roundtrip through Parquet.
    # See ARROW-8006 and ARROW-8011.
    orig_table = pa.Table.from_pydict({"col": array_factory()})
    bio = io.BytesIO()
    pq.write_table(orig_table, bio, use_dictionary=True)
    bio.seek(0)
    read_dictionary = ['col'] if read_dictionary else None
    table = pq.read_table(bio, use_threads=False,
                          read_dictionary=read_dictionary)

    for col in table.columns:
        [chunk] = col.chunks
        buf = chunk.buffers()[1]
        assert buf.to_pybytes() == buf.size * b"\0"


def test_parquet_compression_roundtrip(tempdir):
    # ARROW-10480: ensure even with nonstandard Parquet file naming
    # conventions, writing and then reading a file works. In
    # particular, ensure that we don't automatically double-compress
    # the stream due to auto-detecting the extension in the filename
    table = pa.table([pa.array(range(4))], names=["ints"])
    path = tempdir / "arrow-10480.pyarrow.gz"
    pq.write_table(table, path, compression="GZIP")
    result = pq.read_table(path)
    assert result.equals(table)


def test_empty_row_groups(tempdir):
    # ARROW-3020
    table = pa.Table.from_arrays([pa.array([], type='int32')], ['f0'])

    path = tempdir / 'empty_row_groups.parquet'

    num_groups = 3
    with pq.ParquetWriter(path, table.schema) as writer:
        for i in range(num_groups):
            writer.write_table(table)

    reader = pq.ParquetFile(path)
    assert reader.metadata.num_row_groups == num_groups

    for i in range(num_groups):
        assert reader.read_row_group(i).equals(table)


def test_reads_over_batch(tempdir):
    data = [None] * (1 << 20)
    data.append([1])
    # Large list<int64> with mostly nones and one final
    # value.  This should force batched reads when
    # reading back.
    table = pa.Table.from_arrays([data], ['column'])

    path = tempdir / 'arrow-11607.parquet'
    pq.write_table(table, path)
    table2 = pq.read_table(path)
    assert table == table2


def test_permutation_of_column_order(tempdir):
    # ARROW-2366
    case = tempdir / "dataset_column_order_permutation"
    case.mkdir(exist_ok=True)

    data1 = pa.table([[1, 2, 3], [.1, .2, .3]], names=['a', 'b'])
    pq.write_table(data1, case / "data1.parquet")

    data2 = pa.table([[.4, .5, .6], [4, 5, 6]], names=['b', 'a'])
    pq.write_table(data2, case / "data2.parquet")

    table = pq.read_table(str(case))
    table2 = pa.table([[1, 2, 3, 4, 5, 6],
                       [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]],
                      names=['a', 'b'])

    assert table == table2


def test_thrift_size_limits(tempdir):
    path = tempdir / 'largethrift.parquet'

    array = pa.array(list(range(10)))
    num_cols = 1000
    table = pa.table(
        [array] * num_cols,
        names=[f'some_long_column_name_{i}' for i in range(num_cols)])
    pq.write_table(table, path)

    with pytest.raises(
            OSError,
            match="Couldn't deserialize thrift:.*Exceeded size limit"):
        pq.read_table(path, thrift_string_size_limit=50 * num_cols)
    with pytest.raises(
            OSError,
            match="Couldn't deserialize thrift:.*Exceeded size limit"):
        pq.read_table(path, thrift_container_size_limit=num_cols)

    got = pq.read_table(path, thrift_string_size_limit=100 * num_cols)
    assert got == table
    got = pq.read_table(path, thrift_container_size_limit=2 * num_cols)
    assert got == table
    got = pq.read_table(path)
    assert got == table


def test_page_checksum_verification_write_table(tempdir):
    """Check that checksum verification works for datasets created with
    pq.write_table()"""

    # Write some sample data into a parquet file with page checksum enabled
    original_path = tempdir / 'correct.parquet'
    table_orig = pa.table({'a': [1, 2, 3, 4]})
    pq.write_table(table_orig, original_path, write_page_checksum=True)

    # Read file and verify that the data is correct
    table_check = pq.read_table(original_path, page_checksum_verification=True)
    assert table_orig == table_check

    # Read the original file as binary and swap the 31-th and 36-th bytes. This
    # should be equivalent to storing the following data:
    #    pa.table({'a': [1, 3, 2, 4]})
    bin_data = bytearray(original_path.read_bytes())

    # Swap two bytes to emulate corruption. Also, check that the two bytes are
    # different, otherwise no corruption occurs
    assert bin_data[31] != bin_data[36]
    bin_data[31], bin_data[36] = bin_data[36], bin_data[31]

    # Write the corrupted data to another parquet file
    corrupted_path = tempdir / 'corrupted.parquet'
    corrupted_path.write_bytes(bin_data)

    # Case 1: Reading the corrupted file with read_table() and without page
    # checksum verification succeeds but yields corrupted data
    table_corrupt = pq.read_table(corrupted_path,
                                  page_checksum_verification=False)
    # The read should complete without error, but the table has different
    # content than the original file!
    assert table_corrupt != table_orig
    assert table_corrupt == pa.table({'a': [1, 3, 2, 4]})

    # Case 2: Reading the corrupted file with read_table() and with page
    # checksum verification enabled raises an exception
    with pytest.raises(OSError, match="CRC checksum verification"):
        _ = pq.read_table(corrupted_path, page_checksum_verification=True)

    # Case 3: Reading the corrupted file with ParquetFile.read() and without
    # page checksum verification succeeds but yields corrupted data
    corrupted_pq_file = pq.ParquetFile(corrupted_path,
                                       page_checksum_verification=False)
    table_corrupt2 = corrupted_pq_file.read()
    assert table_corrupt2 != table_orig
    assert table_corrupt2 == pa.table({'a': [1, 3, 2, 4]})

    # Case 4: Reading the corrupted file with ParquetFile.read() and with page
    # checksum verification enabled raises an exception
    corrupted_pq_file = pq.ParquetFile(corrupted_path,
                                       page_checksum_verification=True)
    # Accessing the data should result in an error
    with pytest.raises(OSError, match="CRC checksum verification"):
        _ = corrupted_pq_file.read()


@pytest.mark.dataset
def test_checksum_write_to_dataset(tempdir):
    """Check that checksum verification works for datasets created with
    pq.write_to_dataset"""

    table_orig = pa.table({'a': [1, 2, 3, 4]})

    # Write a sample dataset with page checksum enabled
    original_dir_path = tempdir / 'correct_dir'
    pq.write_to_dataset(table_orig,
                        original_dir_path,
                        write_page_checksum=True)

    # Read file and verify that the data is correct
    original_file_path_list = list(original_dir_path.iterdir())
    assert len(original_file_path_list) == 1
    original_path = original_file_path_list[0]
    table_check = pq.read_table(original_path, page_checksum_verification=True)
    assert table_orig == table_check

    # Read the original file as binary and swap the 31-th and 36-th bytes. This
    # should be equivalent to storing the following data:
    #    pa.table({'a': [1, 3, 2, 4]})
    bin_data = bytearray(original_path.read_bytes())

    # Swap two bytes to emulate corruption. Also, check that the two bytes are
    # different, otherwise no corruption occurs
    assert bin_data[31] != bin_data[36]
    bin_data[31], bin_data[36] = bin_data[36], bin_data[31]

    # Write the corrupted data to another parquet dataset
    # Copy dataset dir (which should be just one file)
    corrupted_dir_path = tempdir / 'corrupted_dir'
    copytree(original_dir_path, corrupted_dir_path)
    # Corrupt just the one file with the dataset
    corrupted_file_path = corrupted_dir_path / original_path.name
    corrupted_file_path.write_bytes(bin_data)

    # Case 1: Reading the corrupted file with read_table() and without page
    # checksum verification succeeds but yields corrupted data
    table_corrupt = pq.read_table(corrupted_file_path,
                                  page_checksum_verification=False)
    # The read should complete without error, but the table has different
    # content than the original file!
    assert table_corrupt != table_orig
    assert table_corrupt == pa.table({'a': [1, 3, 2, 4]})

    # Case 2: Reading the corrupted file with read_table() and with page
    # checksum verification enabled raises an exception
    with pytest.raises(OSError, match="CRC checksum verification"):
        _ = pq.read_table(corrupted_file_path, page_checksum_verification=True)
