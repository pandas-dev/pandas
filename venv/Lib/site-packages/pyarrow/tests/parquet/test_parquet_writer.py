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

import pytest

import pyarrow as pa
from pyarrow import fs

try:
    import pyarrow.parquet as pq
    from pyarrow.tests.parquet.common import (_read_table, _test_dataframe,
                                              _test_table, _range_integers)
except ImportError:
    pq = None


try:
    import pandas as pd
    import pandas.testing as tm

except ImportError:
    pd = tm = None


# Marks all of the tests in this module
# Ignore these with pytest ... -m 'not parquet'
pytestmark = pytest.mark.parquet


@pytest.mark.pandas
def test_parquet_incremental_file_build(tempdir):
    df = _test_dataframe(100)
    df['unique_id'] = 0

    arrow_table = pa.Table.from_pandas(df, preserve_index=False)
    out = pa.BufferOutputStream()

    writer = pq.ParquetWriter(out, arrow_table.schema, version='2.6')

    frames = []
    for i in range(10):
        df['unique_id'] = i
        arrow_table = pa.Table.from_pandas(df, preserve_index=False)
        writer.write_table(arrow_table)

        frames.append(df.copy())

    writer.close()

    buf = out.getvalue()
    result = _read_table(pa.BufferReader(buf))

    expected = pd.concat(frames, ignore_index=True)
    tm.assert_frame_equal(result.to_pandas(), expected)


def test_validate_schema_write_table(tempdir):
    # ARROW-2926
    simple_fields = [
        pa.field('POS', pa.uint32()),
        pa.field('desc', pa.string())
    ]

    simple_schema = pa.schema(simple_fields)

    # simple_table schema does not match simple_schema
    simple_from_array = [pa.array([1]), pa.array(['bla'])]
    simple_table = pa.Table.from_arrays(simple_from_array, ['POS', 'desc'])

    path = tempdir / 'simple_validate_schema.parquet'

    with pq.ParquetWriter(path, simple_schema,
                          version='2.6',
                          compression='snappy', flavor='spark') as w:
        with pytest.raises(ValueError):
            w.write_table(simple_table)


def test_parquet_invalid_writer(tempdir):
    # avoid segfaults with invalid construction
    with pytest.raises(TypeError):
        some_schema = pa.schema([pa.field("x", pa.int32())])
        pq.ParquetWriter(None, some_schema)

    with pytest.raises(TypeError):
        pq.ParquetWriter(tempdir / "some_path", None)


@pytest.mark.pandas
def test_parquet_writer_context_obj(tempdir):
    df = _test_dataframe(100)
    df['unique_id'] = 0

    arrow_table = pa.Table.from_pandas(df, preserve_index=False)
    out = pa.BufferOutputStream()

    with pq.ParquetWriter(out, arrow_table.schema, version='2.6') as writer:

        frames = []
        for i in range(10):
            df['unique_id'] = i
            arrow_table = pa.Table.from_pandas(df, preserve_index=False)
            writer.write_table(arrow_table)

            frames.append(df.copy())

    buf = out.getvalue()
    result = _read_table(pa.BufferReader(buf))

    expected = pd.concat(frames, ignore_index=True)
    tm.assert_frame_equal(result.to_pandas(), expected)


@pytest.mark.pandas
def test_parquet_writer_context_obj_with_exception(tempdir):
    df = _test_dataframe(100)
    df['unique_id'] = 0

    arrow_table = pa.Table.from_pandas(df, preserve_index=False)
    out = pa.BufferOutputStream()
    error_text = 'Artificial Error'

    try:
        with pq.ParquetWriter(out,
                              arrow_table.schema,
                              version='2.6') as writer:

            frames = []
            for i in range(10):
                df['unique_id'] = i
                arrow_table = pa.Table.from_pandas(df, preserve_index=False)
                writer.write_table(arrow_table)
                frames.append(df.copy())
                if i == 5:
                    raise ValueError(error_text)
    except Exception as e:
        assert str(e) == error_text

    buf = out.getvalue()
    result = _read_table(pa.BufferReader(buf))

    expected = pd.concat(frames, ignore_index=True)
    tm.assert_frame_equal(result.to_pandas(), expected)


@pytest.mark.pandas
@pytest.mark.parametrize("filesystem", [
    None,
    fs.LocalFileSystem(),
])
def test_parquet_writer_write_wrappers(tempdir, filesystem):
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)
    batch = pa.RecordBatch.from_pandas(df, preserve_index=False)
    path_table = str(tempdir / 'data_table.parquet')
    path_batch = str(tempdir / 'data_batch.parquet')

    with pq.ParquetWriter(
        path_table, table.schema, filesystem=filesystem, version='2.6'
    ) as writer:
        writer.write_table(table)

    result = _read_table(path_table).to_pandas()
    tm.assert_frame_equal(result, df)

    with pq.ParquetWriter(
        path_batch, table.schema, filesystem=filesystem, version='2.6'
    ) as writer:
        writer.write_batch(batch)

    result = _read_table(path_batch).to_pandas()
    tm.assert_frame_equal(result, df)

    with pq.ParquetWriter(
        path_table, table.schema, filesystem=filesystem, version='2.6'
    ) as writer:
        writer.write(table)

    result = _read_table(path_table).to_pandas()
    tm.assert_frame_equal(result, df)

    with pq.ParquetWriter(
        path_batch, table.schema, filesystem=filesystem, version='2.6'
    ) as writer:
        writer.write(batch)

    result = _read_table(path_batch).to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.large_memory
@pytest.mark.pandas
def test_parquet_writer_chunk_size(tempdir):
    default_chunk_size = 1024 * 1024
    abs_max_chunk_size = 64 * 1024 * 1024

    def check_chunk_size(data_size, chunk_size, expect_num_chunks):
        table = pa.Table.from_arrays([
            _range_integers(data_size, 'b')
        ], names=['x'])
        if chunk_size is None:
            pq.write_table(table, tempdir / 'test.parquet')
        else:
            pq.write_table(table, tempdir / 'test.parquet', row_group_size=chunk_size)
        metadata = pq.read_metadata(tempdir / 'test.parquet')
        expected_chunk_size = default_chunk_size if chunk_size is None else chunk_size
        assert metadata.num_row_groups == expect_num_chunks
        latched_chunk_size = min(expected_chunk_size, abs_max_chunk_size)
        # First chunks should be full size
        for chunk_idx in range(expect_num_chunks - 1):
            assert metadata.row_group(chunk_idx).num_rows == latched_chunk_size
        # Last chunk may be smaller
        remainder = data_size - (expected_chunk_size * (expect_num_chunks - 1))
        if remainder == 0:
            assert metadata.row_group(
                expect_num_chunks - 1).num_rows == latched_chunk_size
        else:
            assert metadata.row_group(expect_num_chunks - 1).num_rows == remainder

    check_chunk_size(default_chunk_size * 2, default_chunk_size - 100, 3)
    check_chunk_size(default_chunk_size * 2, default_chunk_size, 2)
    check_chunk_size(default_chunk_size * 2, default_chunk_size + 100, 2)
    check_chunk_size(default_chunk_size + 100, default_chunk_size + 100, 1)
    # Even though the chunk size requested is large enough it will be capped
    # by the absolute max chunk size
    check_chunk_size(abs_max_chunk_size * 2, abs_max_chunk_size * 2, 2)

    # These tests don't pass a chunk_size to write_table and so the chunk size
    # should be default_chunk_size
    check_chunk_size(default_chunk_size, None, 1)
    check_chunk_size(default_chunk_size + 1, None, 2)


@pytest.mark.pandas
@pytest.mark.parametrize("filesystem", [
    None,
    fs.LocalFileSystem(),
])
def test_parquet_writer_filesystem_local(tempdir, filesystem):
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)
    path = str(tempdir / 'data.parquet')

    with pq.ParquetWriter(
        path, table.schema, filesystem=filesystem, version='2.6'
    ) as writer:
        writer.write_table(table)

    result = _read_table(path).to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.pandas
@pytest.mark.s3
def test_parquet_writer_filesystem_s3(s3_example_fs):
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)

    fs, uri, path = s3_example_fs

    with pq.ParquetWriter(
        path, table.schema, filesystem=fs, version='2.6'
    ) as writer:
        writer.write_table(table)

    result = _read_table(uri).to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.pandas
@pytest.mark.s3
def test_parquet_writer_filesystem_s3_uri(s3_example_fs):
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)

    fs, uri, path = s3_example_fs

    with pq.ParquetWriter(uri, table.schema, version='2.6') as writer:
        writer.write_table(table)

    result = _read_table(path, filesystem=fs).to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.pandas
@pytest.mark.s3
def test_parquet_writer_filesystem_s3fs(s3_example_s3fs):
    df = _test_dataframe(100)
    table = pa.Table.from_pandas(df, preserve_index=False)

    fs, directory = s3_example_s3fs
    path = directory + "/test.parquet"

    with pq.ParquetWriter(
        path, table.schema, filesystem=fs, version='2.6'
    ) as writer:
        writer.write_table(table)

    result = _read_table(path, filesystem=fs).to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.numpy
def test_parquet_writer_filesystem_buffer_raises():
    table = _test_table(100)
    filesystem = fs.LocalFileSystem()

    # Should raise ValueError when filesystem is passed with file-like object
    with pytest.raises(ValueError, match="specified path is file-like"):
        pq.ParquetWriter(
            pa.BufferOutputStream(), table.schema, filesystem=filesystem
        )


def test_parquet_writer_store_schema(tempdir):
    table = pa.table({'a': [1, 2, 3]})

    # default -> write schema information
    path1 = tempdir / 'test_with_schema.parquet'
    with pq.ParquetWriter(path1, table.schema) as writer:
        writer.write_table(table)

    meta = pq.read_metadata(path1)
    assert b'ARROW:schema' in meta.metadata
    assert meta.metadata[b'ARROW:schema']

    # disable adding schema information
    path2 = tempdir / 'test_without_schema.parquet'
    with pq.ParquetWriter(path2, table.schema, store_schema=False) as writer:
        writer.write_table(table)

    meta = pq.read_metadata(path2)
    assert meta.metadata is None


def test_parquet_writer_append_key_value_metadata(tempdir):
    table = pa.Table.from_arrays([pa.array([], type='int32')], ['f0'])
    path = tempdir / 'metadata.parquet'

    with pq.ParquetWriter(path, table.schema) as writer:
        writer.write_table(table)
        writer.add_key_value_metadata({'key1': '1', 'key2': 'x'})
        writer.add_key_value_metadata({'key2': '2', 'key3': '3'})
    reader = pq.ParquetFile(path)
    metadata = reader.metadata.metadata
    assert metadata[b'key1'] == b'1'
    assert metadata[b'key2'] == b'2'
    assert metadata[b'key3'] == b'3'


def test_parquet_content_defined_chunking(tempdir):
    table = pa.table({'a': range(100_000)})

    # use PLAIN encoding because we compare the overall size of the row groups
    # which would vary depending on the encoding making the assertions wrong
    pq.write_table(table, tempdir / 'unchunked.parquet',
                   use_dictionary=False,
                   column_encoding="PLAIN")
    pq.write_table(table, tempdir / 'chunked-default.parquet',
                   use_dictionary=False,
                   column_encoding="PLAIN",
                   use_content_defined_chunking=True)
    pq.write_table(table, tempdir / 'chunked-custom.parquet',
                   use_dictionary=False,
                   column_encoding="PLAIN",
                   use_content_defined_chunking={"min_chunk_size": 32_768,
                                                 "max_chunk_size": 65_536})

    # the data must be the same
    unchunked = pq.read_table(tempdir / 'unchunked.parquet')
    chunked_default = pq.read_table(tempdir / 'chunked-default.parquet')
    chunked_custom = pq.read_table(tempdir / 'chunked-custom.parquet')
    assert unchunked.equals(chunked_default)
    assert unchunked.equals(chunked_custom)

    # number of row groups and their sizes are not affected by content defined chunking
    unchunked_metadata = pq.read_metadata(tempdir / 'unchunked.parquet')
    chunked_default_metadata = pq.read_metadata(tempdir / 'chunked-default.parquet')
    chunked_custom_metadata = pq.read_metadata(tempdir / 'chunked-custom.parquet')

    assert unchunked_metadata.num_row_groups == chunked_default_metadata.num_row_groups
    assert unchunked_metadata.num_row_groups == chunked_custom_metadata.num_row_groups

    for i in range(unchunked_metadata.num_row_groups):
        rg_unchunked = unchunked_metadata.row_group(i)
        rg_chunked_default = chunked_default_metadata.row_group(i)
        rg_chunked_custom = chunked_custom_metadata.row_group(i)
        assert rg_unchunked.num_rows == rg_chunked_default.num_rows
        assert rg_unchunked.num_rows == rg_chunked_custom.num_rows
        # since PageReader is not exposed we cannot inspect the page sizes
        # so just check that the total byte size is different
        assert rg_unchunked.total_byte_size < rg_chunked_default.total_byte_size
        assert rg_unchunked.total_byte_size < rg_chunked_custom.total_byte_size
        assert rg_chunked_default.total_byte_size < rg_chunked_custom.total_byte_size


def test_parquet_content_defined_chunking_parameters(tempdir):
    table = pa.table({'a': range(100)})
    path = tempdir / 'chunked-invalid.parquet'

    # it raises OSError, not ideal but this is how parquet exceptions are handled
    # currently
    msg = "max_chunk_size must be greater than min_chunk_size"
    with pytest.raises(Exception, match=msg):
        cdc_options = {"min_chunk_size": 65_536, "max_chunk_size": 32_768}
        pq.write_table(table, path, use_content_defined_chunking=cdc_options)

    cases = [
        (
            {"min_chunk_size": 64 * 1024, "unknown_option": True},
            "Unknown options in 'use_content_defined_chunking': {'unknown_option'}"
        ),
        (
            {"min_chunk_size": 64 * 1024},
            "Missing options in 'use_content_defined_chunking': {'max_chunk_size'}"
        ),
        (
            {"max_chunk_size": 64 * 1024},
            "Missing options in 'use_content_defined_chunking': {'min_chunk_size'}"
        )
    ]
    for cdc_options, msg in cases:
        with pytest.raises(ValueError, match=msg):
            pq.write_table(table, path, use_content_defined_chunking=cdc_options)

    # using the default parametrization
    pq.write_table(table, path, use_content_defined_chunking=True)

    # using min_chunk_size and max_chunk_size
    cdc_options = {"min_chunk_size": 32_768, "max_chunk_size": 65_536}
    pq.write_table(table, path, use_content_defined_chunking=cdc_options)

    # using min_chunk_size, max_chunk_size and norm_level
    cdc_options = {"min_chunk_size": 32_768, "max_chunk_size": 65_536, "norm_level": 1}
    pq.write_table(table, path, use_content_defined_chunking=cdc_options)


@pytest.mark.parametrize("time_type, time_unit", [
    (pa.time32, "s"),
    (pa.time32, "ms"),
    (pa.time64, "us"),
    (pa.time64, "ns"),
])
@pytest.mark.parametrize("utc_flag_val", [False, True])
def test_arrow_writer_props_time_adjusted_to_utc(
    tempdir,
    utc_flag_val,
    time_type,
    time_unit,
):
    # GH-47441
    filename = tempdir / "time_adjusted_to_utc.parquet"

    time_values = [0, 123, 10_000, 86_399]

    table = pa.table({
        "time_col": pa.array(time_values, type=time_type(time_unit)),
    })

    schema = pa.schema([
        ("time_col", time_type(time_unit)),
    ])

    with pq.ParquetWriter(
        where=filename,
        schema=schema,
        write_time_adjusted_to_utc=utc_flag_val,
    ) as writer:
        writer.write_table(table)

    result = pq.read_table(filename, schema=schema)

    result.validate(full=True)

    assert result.equals(table)


@pytest.mark.parametrize(
    "max_rows_per_page",
    [1, 10, 100, 1_000, None],
)
def test_writer_props_max_rows_per_page(tempdir, max_rows_per_page):
    # GH-48096
    filename = tempdir / "max_rows_per_page.parquet"

    table = pa.table({
        "x": pa.array([1, 2, 3, 4, 5, 6, 7], type=pa.int8()),
        "y": pa.array([11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0], type=pa.float16()),
    })

    schema = pa.schema([
        ("x", pa.int8()),
        ("y", pa.float16()),
    ])

    with pq.ParquetWriter(
        where=filename,
        schema=schema,
        max_rows_per_page=max_rows_per_page,
    ) as writer:
        writer.write_table(table)

    result = pq.read_table(filename, schema=schema)

    result.validate(full=True)

    assert result.equals(table)


def test_writer_props_max_rows_per_page_file_size(tempdir):
    # GH-48096
    table = pa.table({
        "x": pa.array(range(1_000_000))
    })

    local = fs.LocalFileSystem()
    file_infos = []

    for max_rows in (1_000, 10_000):
        path = f"{tempdir}/max_rows_per_page_{max_rows}.parquet"

        with pq.ParquetWriter(
            where=path,
            schema=table.schema,
            max_rows_per_page=max_rows,
        ) as writer:
            writer.write_table(table)

        file_infos.append(local.get_file_info(path))

    # A smaller maximum rows parameter should produce a larger file
    assert file_infos[0].size > file_infos[1].size
