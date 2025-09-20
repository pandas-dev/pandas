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

import io
import os
import re
import sys
import types

import pytest
from unittest import mock

import pyarrow as pa

try:
    import pyarrow.parquet as pq
    from pyarrow.tests.parquet.common import _write_table
except ImportError:
    pq = None

try:
    import pandas as pd
    import pandas.testing as tm

    from pyarrow.tests.parquet.common import alltypes_sample
except ImportError:
    pd = tm = None


# Marks all of the tests in this module
# Ignore these with pytest ... -m 'not parquet'
pytestmark = pytest.mark.parquet


@pytest.mark.pandas
def test_pass_separate_metadata():
    # ARROW-471
    df = alltypes_sample(size=10000)

    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, compression='snappy', version='2.6')

    buf.seek(0)
    metadata = pq.read_metadata(buf)

    buf.seek(0)

    fileh = pq.ParquetFile(buf, metadata=metadata)

    tm.assert_frame_equal(df, fileh.read().to_pandas())


@pytest.mark.pandas
def test_read_single_row_group():
    # ARROW-471
    N, K = 10000, 4
    df = alltypes_sample(size=N)

    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, row_group_size=N / K,
                 compression='snappy', version='2.6')

    buf.seek(0)

    pf = pq.ParquetFile(buf)

    assert pf.num_row_groups == K

    row_groups = [pf.read_row_group(i) for i in range(K)]
    result = pa.concat_tables(row_groups)
    tm.assert_frame_equal(df, result.to_pandas())


@pytest.mark.pandas
def test_read_single_row_group_with_column_subset():
    N, K = 10000, 4
    df = alltypes_sample(size=N)
    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, row_group_size=N / K,
                 compression='snappy', version='2.6')

    buf.seek(0)
    pf = pq.ParquetFile(buf)

    cols = list(df.columns[:2])
    row_groups = [pf.read_row_group(i, columns=cols) for i in range(K)]
    result = pa.concat_tables(row_groups)
    tm.assert_frame_equal(df[cols], result.to_pandas())

    # ARROW-4267: Selection of duplicate columns still leads to these columns
    # being read uniquely.
    row_groups = [pf.read_row_group(i, columns=cols + cols) for i in range(K)]
    result = pa.concat_tables(row_groups)
    tm.assert_frame_equal(df[cols], result.to_pandas())


@pytest.mark.pandas
def test_read_multiple_row_groups():
    N, K = 10000, 4
    df = alltypes_sample(size=N)

    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, row_group_size=N / K,
                 compression='snappy', version='2.6')

    buf.seek(0)

    pf = pq.ParquetFile(buf)

    assert pf.num_row_groups == K

    result = pf.read_row_groups(range(K))
    tm.assert_frame_equal(df, result.to_pandas())


@pytest.mark.pandas
def test_read_multiple_row_groups_with_column_subset():
    N, K = 10000, 4
    df = alltypes_sample(size=N)
    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, row_group_size=N / K,
                 compression='snappy', version='2.6')

    buf.seek(0)
    pf = pq.ParquetFile(buf)

    cols = list(df.columns[:2])
    result = pf.read_row_groups(range(K), columns=cols)
    tm.assert_frame_equal(df[cols], result.to_pandas())

    # ARROW-4267: Selection of duplicate columns still leads to these columns
    # being read uniquely.
    result = pf.read_row_groups(range(K), columns=cols + cols)
    tm.assert_frame_equal(df[cols], result.to_pandas())


@pytest.mark.pandas
def test_scan_contents():
    N, K = 10000, 4
    df = alltypes_sample(size=N)
    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, row_group_size=N / K,
                 compression='snappy', version='2.6')

    buf.seek(0)
    pf = pq.ParquetFile(buf)

    assert pf.scan_contents() == 10000
    assert pf.scan_contents(df.columns[:4]) == 10000


def test_parquet_file_pass_directory_instead_of_file(tempdir):
    # ARROW-7208
    path = tempdir / 'directory'
    os.mkdir(str(path))

    msg = f"Cannot open for reading: path '{str(path)}' is a directory"
    with pytest.raises(IOError) as exc:
        pq.ParquetFile(path)
    if exc.errisinstance(PermissionError) and sys.platform == 'win32':
        return  # Windows CI can get a PermissionError here.
    exc.match(msg)


def test_read_column_invalid_index():
    table = pa.table([pa.array([4, 5]), pa.array(["foo", "bar"])],
                     names=['ints', 'strs'])
    bio = pa.BufferOutputStream()
    pq.write_table(table, bio)
    f = pq.ParquetFile(bio.getvalue())
    assert f.reader.read_column(0).to_pylist() == [4, 5]
    assert f.reader.read_column(1).to_pylist() == ["foo", "bar"]
    for index in (-1, 2):
        with pytest.raises((ValueError, IndexError)):
            f.reader.read_column(index)


@pytest.mark.pandas
@pytest.mark.parametrize('batch_size', [300, 1000, 1300])
def test_iter_batches_columns_reader(tempdir, batch_size):
    total_size = 3000
    chunk_size = 1000
    # TODO: Add categorical support
    df = alltypes_sample(size=total_size)

    filename = tempdir / 'pandas_roundtrip.parquet'
    arrow_table = pa.Table.from_pandas(df)
    _write_table(arrow_table, filename, version='2.6',
                 chunk_size=chunk_size)

    file_ = pq.ParquetFile(filename)
    for columns in [df.columns[:10], df.columns[10:]]:
        batches = file_.iter_batches(batch_size=batch_size, columns=columns)
        batch_starts = range(0, total_size+batch_size, batch_size)
        for batch, start in zip(batches, batch_starts):
            end = min(total_size, start + batch_size)
            tm.assert_frame_equal(
                batch.to_pandas(),
                df.iloc[start:end, :].loc[:, columns].reset_index(drop=True)
            )


@pytest.mark.pandas
@pytest.mark.parametrize('chunk_size', [1000])
def test_iter_batches_reader(tempdir, chunk_size):
    df = alltypes_sample(size=10000, categorical=True)

    filename = tempdir / 'pandas_roundtrip.parquet'
    arrow_table = pa.Table.from_pandas(df)
    assert arrow_table.schema.pandas_metadata is not None

    _write_table(arrow_table, filename, version='2.6',
                 chunk_size=chunk_size)

    file_ = pq.ParquetFile(filename)

    def get_all_batches(f):
        for row_group in range(f.num_row_groups):
            batches = f.iter_batches(
                batch_size=900,
                row_groups=[row_group],
            )

            for batch in batches:
                yield batch

    batches = list(get_all_batches(file_))
    batch_no = 0

    for i in range(file_.num_row_groups):
        tm.assert_frame_equal(
            batches[batch_no].to_pandas(),
            file_.read_row_groups([i]).to_pandas().head(900)
        )

        batch_no += 1

        tm.assert_frame_equal(
            batches[batch_no].to_pandas().reset_index(drop=True),
            file_.read_row_groups([i]).to_pandas().iloc[900:].reset_index(
                drop=True
            )
        )

        batch_no += 1


@pytest.mark.pandas
@pytest.mark.parametrize('pre_buffer', [False, True])
def test_pre_buffer(pre_buffer):
    N, K = 10000, 4
    df = alltypes_sample(size=N)
    a_table = pa.Table.from_pandas(df)

    buf = io.BytesIO()
    _write_table(a_table, buf, row_group_size=N / K,
                 compression='snappy', version='2.6')

    buf.seek(0)
    pf = pq.ParquetFile(buf, pre_buffer=pre_buffer)
    assert pf.read().num_rows == N


def test_parquet_file_explicitly_closed(tempdir):
    """
    Unopened files should be closed explicitly after use,
    and previously opened files should be left open.
    Applies to read_table, ParquetDataset, and ParquetFile
    """
    # create test parquet file
    fn = tempdir.joinpath('file.parquet')
    table = pa.table({'col1': [0, 1], 'col2': [0, 1]})
    pq.write_table(table, fn)

    # ParquetFile with opened file (will leave open)
    with open(fn, 'rb') as f:
        with pq.ParquetFile(f) as p:
            p.read()
            assert not f.closed
            assert not p.closed
        assert not f.closed  # opened input file was not closed
        assert not p.closed  # parquet file obj reports as not closed
    assert f.closed
    assert p.closed  # parquet file being closed reflects underlying file

    # ParquetFile with unopened file (will close)
    with pq.ParquetFile(fn) as p:
        p.read()
        assert not p.closed
    assert p.closed  # parquet file obj reports as closed


@pytest.mark.s3
@pytest.mark.parametrize("use_uri", (True, False))
def test_parquet_file_with_filesystem(s3_example_fs, use_uri):
    s3_fs, s3_uri, s3_path = s3_example_fs

    args = (s3_uri if use_uri else s3_path,)
    kwargs = {} if use_uri else dict(filesystem=s3_fs)

    table = pa.table({"a": range(10)})
    pq.write_table(table, s3_path, filesystem=s3_fs)

    parquet_file = pq.ParquetFile(*args, **kwargs)
    assert parquet_file.read() == table
    assert not parquet_file.closed
    parquet_file.close()
    assert parquet_file.closed

    with pq.ParquetFile(*args, **kwargs) as f:
        assert f.read() == table
        assert not f.closed
    assert f.closed


def test_read_statistics():
    table = pa.table({"value": pa.array([-1, None, 3])})
    buf = io.BytesIO()
    _write_table(table, buf)
    buf.seek(0)

    statistics = pq.ParquetFile(buf).read().columns[0].chunks[0].statistics
    assert statistics.null_count == 1
    assert statistics.distinct_count is None
    assert statistics.min == -1
    assert statistics.is_min_exact
    assert statistics.max == 3
    assert statistics.is_max_exact
    assert repr(statistics) == ("arrow.ArrayStatistics<"
                                "null_count=1, distinct_count=None, "
                                "min=-1, is_min_exact=True, "
                                "max=3, is_max_exact=True>")


def test_read_undefined_logical_type(parquet_test_datadir):
    test_file = f"{parquet_test_datadir}/unknown-logical-type.parquet"

    table = pq.ParquetFile(test_file).read()
    assert table.column_names == ["column with known type", "column with unknown type"]
    assert table["column with unknown type"].to_pylist() == [
        b"unknown string 1",
        b"unknown string 2",
        b"unknown string 3"
    ]


def test_parquet_file_fsspec_support():
    pytest.importorskip("fsspec")

    table = pa.table({"a": range(10)})
    pq.write_table(table, "fsspec+memory://example.parquet")
    table2 = pq.read_table("fsspec+memory://example.parquet")
    assert table.equals(table2)

    msg = "Unrecognized filesystem type in URI"
    with pytest.raises(pa.ArrowInvalid, match=msg):
        pq.read_table("non-existing://example.parquet")


def test_parquet_file_fsspec_support_through_filesystem_argument():
    try:
        from fsspec.implementations.memory import MemoryFileSystem
    except ImportError:
        pytest.skip("fsspec is not installed, skipping test")

    table = pa.table({"b": range(10)})

    fs = MemoryFileSystem()
    fs.mkdir("/path/to/prefix", create_parents=True)
    assert fs.exists("/path/to/prefix")

    fs_str = "fsspec+memory://path/to/prefix"
    pq.write_table(table, "b.parquet", filesystem=fs_str)
    table2 = pq.read_table("fsspec+memory://path/to/prefix/b.parquet")
    assert table.equals(table2)


def test_parquet_file_hugginface_support():
    try:
        from fsspec.implementations.memory import MemoryFileSystem
    except ImportError:
        pytest.skip("fsspec is not installed, skipping Hugging Face test")

    fake_hf_module = types.ModuleType("huggingface_hub")
    fake_hf_module.HfFileSystem = MemoryFileSystem
    with mock.patch.dict("sys.modules", {"huggingface_hub": fake_hf_module}):
        uri = "hf://datasets/apache/arrow/test.parquet"
        table = pa.table({"a": range(10)})
        pq.write_table(table, uri)
        table2 = pq.read_table(uri)
        assert table.equals(table2)


def test_fsspec_uri_raises_if_fsspec_is_not_available():
    # sadly cannot patch sys.modules because cython will still be able to import fsspec
    try:
        import fsspec  # noqa: F401
    except ImportError:
        pass
    else:
        pytest.skip("fsspec is available, skipping test")

    msg = re.escape(
        "`fsspec` is required to handle `fsspec+<filesystem>://` and `hf://` URIs.")
    with pytest.raises(ImportError, match=msg):
        pq.read_table("fsspec+memory://example.parquet")


def test_iter_batches_raises_batch_size_zero(tempdir):
    # See https://github.com/apache/arrow/issues/46811
    schema = pa.schema([])
    empty_table = pa.Table.from_batches([], schema=schema)
    parquet_file_path = tempdir / "empty_file.parquet"
    pq.write_table(empty_table, parquet_file_path)
    parquet_file = pq.ParquetFile(parquet_file_path)
    with pytest.raises(ValueError):
        parquet_file.iter_batches(batch_size=0)
