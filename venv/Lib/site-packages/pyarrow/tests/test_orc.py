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

import decimal
import datetime
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

import pyarrow as pa
from pyarrow import fs
from pyarrow.tests import util


# Marks all of the tests in this module
# Ignore these with pytest ... -m 'not orc'
pytestmark = pytest.mark.orc


try:
    from pandas.testing import assert_frame_equal
    import pandas as pd
except ImportError:
    pass


@pytest.fixture(scope="module")
def datadir(base_datadir):
    return base_datadir / "orc"


def fix_example_values(actual_cols, expected_cols):
    """
    Fix type of expected values (as read from JSON) according to
    actual ORC datatype.
    """
    for name in expected_cols:
        expected = expected_cols[name]
        actual = actual_cols[name]
        if (name == "map" and
                [d.keys() == {'key', 'value'} for m in expected for d in m]):
            # convert [{'key': k, 'value': v}, ...] to [(k, v), ...]
            col = expected_cols[name].copy()
            for i, m in enumerate(expected):
                col[i] = [(d['key'], d['value']) for d in m]
            expected_cols[name] = col
            continue

        typ = actual[0].__class__
        if issubclass(typ, datetime.datetime):
            # timestamp fields are represented as strings in JSON files
            expected = pd.to_datetime(expected)
        elif issubclass(typ, datetime.date):
            # date fields are represented as strings in JSON files
            expected = expected.dt.date
        elif typ is decimal.Decimal:
            converted_decimals = [None] * len(expected)
            # decimal fields are represented as reals in JSON files
            for i, (d, v) in enumerate(zip(actual, expected)):
                if not pd.isnull(v):
                    exp = d.as_tuple().exponent
                    factor = 10 ** -exp
                    converted_decimals[i] = (
                        decimal.Decimal(round(v * factor)).scaleb(exp))
            expected = pd.Series(converted_decimals)

        expected_cols[name] = expected


def check_example_values(orc_df, expected_df, start=None, stop=None):
    if start is not None or stop is not None:
        expected_df = expected_df[start:stop].reset_index(drop=True)
    assert_frame_equal(orc_df, expected_df, check_dtype=False)


def check_example_file(orc_path, expected_df, need_fix=False):
    """
    Check a ORC file against the expected columns dictionary.
    """
    from pyarrow import orc

    orc_file = orc.ORCFile(orc_path)
    # Exercise ORCFile.read()
    table = orc_file.read()
    assert isinstance(table, pa.Table)
    table.validate()

    # This workaround needed because of ARROW-3080
    orc_df = pd.DataFrame(table.to_pydict())

    assert set(expected_df.columns) == set(orc_df.columns)

    # reorder columns if necessary
    if not orc_df.columns.equals(expected_df.columns):
        expected_df = expected_df.reindex(columns=orc_df.columns)

    if need_fix:
        fix_example_values(orc_df, expected_df)

    check_example_values(orc_df, expected_df)
    # Exercise ORCFile.read_stripe()
    json_pos = 0
    for i in range(orc_file.nstripes):
        batch = orc_file.read_stripe(i)
        check_example_values(pd.DataFrame(batch.to_pydict()),
                             expected_df,
                             start=json_pos,
                             stop=json_pos + len(batch))
        json_pos += len(batch)
    assert json_pos == orc_file.nrows


@pytest.mark.pandas
@pytest.mark.parametrize('filename', [
    'TestOrcFile.test1.orc',
    'TestOrcFile.testDate1900.orc',
    'decimal.orc'
])
def test_example_using_json(filename, datadir):
    """
    Check a ORC file example against the equivalent JSON file, as given
    in the Apache ORC repository (the JSON file has one JSON object per
    line, corresponding to one row in the ORC file).
    """
    # Read JSON file
    path = datadir / filename
    table = pd.read_json(str(path.with_suffix('.jsn.gz')), lines=True)
    check_example_file(path, table, need_fix=True)


def test_timezone_database_absent(datadir):
    # Example file relies on the timezone "US/Pacific". It should gracefully
    # fail, not crash, if the timezone database is not found.
    path = datadir / 'TestOrcFile.testDate1900.orc'
    code = f"""if 1:
        import os
        os.environ['TZDIR'] = '/tmp/non_existent'

        from pyarrow import orc
        try:
            orc_file = orc.ORCFile({str(path)!r})
            orc_file.read()
        except Exception as e:
            assert "time zone database" in str(e).lower(), e
        else:
            assert False, "Should have raised exception"
    """
    subprocess.run([sys.executable, "-c", code], check=True)


def test_timezone_absent(datadir, tmpdir):
    # Example file relies on the timezone "US/Pacific". It should gracefully
    # fail, not crash, if the timezone database is present but the timezone
    # is not found (GH-40633).
    source_tzdir = Path('/usr/share/zoneinfo')
    if not source_tzdir.exists():
        pytest.skip(f"Test needs timezone database in {source_tzdir}")
    tzdir = Path(tmpdir / 'zoneinfo')
    try:
        shutil.copytree(source_tzdir, tzdir, symlinks=True)
    except OSError as e:
        pytest.skip(f"Failed to copy timezone database: {e}")
    # ORC 2.1.1 Creates an alias between some legacy Timezones
    # https://github.com/apache/orc/pull/2422
    # Example US/Pacific -> America/Los_Angeles
    # Remove both to simulate missing timezone and avoid alias resolution
    (tzdir / 'US' / 'Pacific').unlink(missing_ok=True)
    (tzdir / 'America' / 'Los_Angeles').unlink(missing_ok=True)

    path = datadir / 'TestOrcFile.testDate1900.orc'
    code = f"""if 1:
        import os
        os.environ['TZDIR'] = {str(tzdir)!r}

        from pyarrow import orc
        orc_file = orc.ORCFile({str(path)!r})
        try:
            orc_file.read()
        except Exception as e:
            timezones = ["zoneinfo/US/Pacific", "zoneinfo/America/Los_Angeles"]
            assert any(tz in str(e) for tz in timezones), e
        else:
            assert False, "Should have raised exception"
    """
    subprocess.run([sys.executable, "-c", code], check=True)


def test_orcfile_empty(datadir):
    from pyarrow import orc

    table = orc.ORCFile(datadir / "TestOrcFile.emptyFile.orc").read()
    assert table.num_rows == 0

    expected_schema = pa.schema([
        ("boolean1", pa.bool_()),
        ("byte1", pa.int8()),
        ("short1", pa.int16()),
        ("int1", pa.int32()),
        ("long1", pa.int64()),
        ("float1", pa.float32()),
        ("double1", pa.float64()),
        ("bytes1", pa.binary()),
        ("string1", pa.string()),
        ("middle", pa.struct(
            [("list", pa.list_(
                pa.struct([("int1", pa.int32()),
                           ("string1", pa.string())])))
             ])),
        ("list", pa.list_(
            pa.struct([("int1", pa.int32()),
                       ("string1", pa.string())])
        )),
        ("map", pa.map_(pa.string(),
                        pa.struct([("int1", pa.int32()),
                                   ("string1", pa.string())])
                        )),
    ])
    assert table.schema == expected_schema


def test_filesystem_uri(tmpdir):
    from pyarrow import orc
    table = pa.table({"a": [1, 2, 3]})

    directory = tmpdir / "data_dir"
    directory.mkdir()
    path = directory / "data.orc"
    orc.write_table(table, str(path))

    # filesystem object
    result = orc.read_table(path, filesystem=fs.LocalFileSystem())
    assert result.equals(table)

    # filesystem URI
    result = orc.read_table(
        "data_dir/data.orc", filesystem=util._filesystem_uri(tmpdir))
    assert result.equals(table)

    # use the path only
    result = orc.read_table(
        util._filesystem_uri(path))
    assert result.equals(table)


def test_orcfile_readwrite(tmpdir):
    from pyarrow import orc
    a = pa.array([1, None, 3, None])
    b = pa.array([None, "Arrow", None, "ORC"])
    table = pa.table({"int64": a, "utf8": b})
    file = tmpdir.join("test.orc")
    orc.write_table(table, file)
    output_table = orc.read_table(file)
    assert table.equals(output_table)

    output_table = orc.read_table(file, [])
    assert 4 == output_table.num_rows
    assert 0 == output_table.num_columns

    output_table = orc.read_table(file, columns=["int64"])
    assert 4 == output_table.num_rows
    assert 1 == output_table.num_columns


def test_bytesio_readwrite():
    from pyarrow import orc
    from io import BytesIO

    buf = BytesIO()
    a = pa.array([1, None, 3, None])
    b = pa.array([None, "Arrow", None, "ORC"])
    table = pa.table({"int64": a, "utf8": b})
    orc.write_table(table, buf)
    buf.seek(0)
    orc_file = orc.ORCFile(buf)
    output_table = orc_file.read()
    assert table.equals(output_table)


def test_buffer_readwrite():
    from pyarrow import orc

    buffer_output_stream = pa.BufferOutputStream()
    a = pa.array([1, None, 3, None])
    b = pa.array([None, "Arrow", None, "ORC"])
    table = pa.table({"int64": a, "utf8": b})
    orc.write_table(table, buffer_output_stream)
    buffer_reader = pa.BufferReader(buffer_output_stream.getvalue())
    orc_file = orc.ORCFile(buffer_reader)
    output_table = orc_file.read()
    assert table.equals(output_table)
    # Check for default WriteOptions
    assert orc_file.compression == 'UNCOMPRESSED'
    assert orc_file.file_version == '0.12'
    assert orc_file.row_index_stride == 10000
    assert orc_file.compression_size == 65536

    # deprecated keyword order
    buffer_output_stream = pa.BufferOutputStream()
    with pytest.warns(FutureWarning):
        orc.write_table(buffer_output_stream, table)
    buffer_reader = pa.BufferReader(buffer_output_stream.getvalue())
    orc_file = orc.ORCFile(buffer_reader)
    output_table = orc_file.read()
    assert table.equals(output_table)
    # Check for default WriteOptions
    assert orc_file.compression == 'UNCOMPRESSED'
    assert orc_file.file_version == '0.12'
    assert orc_file.row_index_stride == 10000
    assert orc_file.compression_size == 65536


@pytest.mark.snappy
def test_buffer_readwrite_with_writeoptions():
    from pyarrow import orc

    buffer_output_stream = pa.BufferOutputStream()
    a = pa.array([1, None, 3, None])
    b = pa.array([None, "Arrow", None, "ORC"])
    table = pa.table({"int64": a, "utf8": b})
    orc.write_table(
        table,
        buffer_output_stream,
        compression='snappy',
        file_version='0.11',
        row_index_stride=5000,
        compression_block_size=65536,
    )
    buffer_reader = pa.BufferReader(buffer_output_stream.getvalue())
    orc_file = orc.ORCFile(buffer_reader)
    output_table = orc_file.read()
    assert table.equals(output_table)
    # Check for modified WriteOptions
    assert orc_file.compression == 'SNAPPY'
    assert orc_file.file_version == '0.11'
    assert orc_file.row_index_stride == 5000
    assert orc_file.compression_size == 65536

    # deprecated keyword order
    buffer_output_stream = pa.BufferOutputStream()
    with pytest.warns(FutureWarning):
        orc.write_table(
            buffer_output_stream,
            table,
            compression='uncompressed',
            file_version='0.11',
            row_index_stride=20000,
            compression_block_size=65536,
        )
    buffer_reader = pa.BufferReader(buffer_output_stream.getvalue())
    orc_file = orc.ORCFile(buffer_reader)
    output_table = orc_file.read()
    assert table.equals(output_table)
    # Check for default WriteOptions
    assert orc_file.compression == 'UNCOMPRESSED'
    assert orc_file.file_version == '0.11'
    assert orc_file.row_index_stride == 20000
    assert orc_file.compression_size == 65536


def test_buffer_readwrite_with_bad_writeoptions():
    from pyarrow import orc
    buffer_output_stream = pa.BufferOutputStream()
    a = pa.array([1, None, 3, None])
    table = pa.table({"int64": a})

    # batch_size must be a positive integer
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            batch_size=0,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            batch_size=-100,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            batch_size=1024.23,
        )

    # file_version must be 0.11 or 0.12
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            file_version=0.13,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            file_version='1.1',
        )

    # stripe_size must be a positive integer
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            stripe_size=0,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            stripe_size=-400,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            stripe_size=4096.73,
        )

    # compression must be among the given options
    with pytest.raises(TypeError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression=0,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression='none',
        )
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression='zlid',
        )

    # compression_block_size must be a positive integer
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression_block_size=0,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression_block_size=-200,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression_block_size=1096.73,
        )

    # compression_strategy must be among the given options
    with pytest.raises(TypeError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression_strategy=0,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression_strategy='no',
        )
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            compression_strategy='large',
        )

    # row_index_stride must be a positive integer
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            row_index_stride=0,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            row_index_stride=-800,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            row_index_stride=3096.29,
        )

    # padding_tolerance must be possible to cast to float
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            padding_tolerance='cat',
        )

    # dictionary_key_size_threshold must be possible to cast to
    # float between 0.0 and 1.0
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            dictionary_key_size_threshold='arrow',
        )
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            dictionary_key_size_threshold=1.2,
        )
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            dictionary_key_size_threshold=-3.2,
        )

    # bloom_filter_columns must be convertible to a list containing
    # nonnegative integers
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            bloom_filter_columns="string",
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            bloom_filter_columns=[0, 1.4],
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            bloom_filter_columns={0, 2, -1},
        )

    # bloom_filter_fpp must be convertible to a float between 0.0 and 1.0
    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            bloom_filter_fpp='arrow',
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            bloom_filter_fpp=1.1,
        )

    with pytest.raises(ValueError):
        orc.write_table(
            table,
            buffer_output_stream,
            bloom_filter_fpp=-0.1,
        )


def test_column_selection(tempdir):
    from pyarrow import orc

    # create a table with nested types
    inner = pa.field('inner', pa.int64())
    middle = pa.field('middle', pa.struct([inner]))
    fields = [
        pa.field('basic', pa.int32()),
        pa.field(
            'list', pa.list_(pa.field('item', pa.int32()))
        ),
        pa.field(
            'struct', pa.struct([middle, pa.field('inner2', pa.int64())])
        ),
        pa.field(
            'list-struct', pa.list_(pa.field(
                'item', pa.struct([
                    pa.field('inner1', pa.int64()),
                    pa.field('inner2', pa.int64())
                ])
            ))
        ),
        pa.field('basic2', pa.int64()),
    ]
    arrs = [
        [0], [[1, 2]], [{"middle": {"inner": 3}, "inner2": 4}],
        [[{"inner1": 5, "inner2": 6}, {"inner1": 7, "inner2": 8}]], [9]]
    table = pa.table(arrs, schema=pa.schema(fields))

    path = str(tempdir / 'test.orc')
    orc.write_table(table, path)
    orc_file = orc.ORCFile(path)

    # default selecting all columns
    result1 = orc_file.read()
    assert result1.equals(table)

    # selecting with columns names
    result2 = orc_file.read(columns=["basic", "basic2"])
    assert result2.equals(table.select(["basic", "basic2"]))

    result3 = orc_file.read(columns=["list", "struct", "basic2"])
    assert result3.equals(table.select(["list", "struct", "basic2"]))

    # using dotted paths
    result4 = orc_file.read(columns=["struct.middle.inner"])
    expected4 = pa.table({"struct": [{"middle": {"inner": 3}}]})
    assert result4.equals(expected4)

    result5 = orc_file.read(columns=["struct.inner2"])
    expected5 = pa.table({"struct": [{"inner2": 4}]})
    assert result5.equals(expected5)

    result6 = orc_file.read(
        columns=["list", "struct.middle.inner", "struct.inner2"]
    )
    assert result6.equals(table.select(["list", "struct"]))

    result7 = orc_file.read(columns=["list-struct.inner1"])
    expected7 = pa.table({"list-struct": [[{"inner1": 5}, {"inner1": 7}]]})
    assert result7.equals(expected7)

    # selecting with (Arrow-based) field indices
    result2 = orc_file.read(columns=[0, 4])
    assert result2.equals(table.select(["basic", "basic2"]))

    result3 = orc_file.read(columns=[1, 2, 3])
    assert result3.equals(table.select(["list", "struct", "list-struct"]))

    # error on non-existing name or index
    with pytest.raises(IOError):
        # liborc returns ParseError, which gets translated into IOError
        # instead of ValueError
        orc_file.read(columns=["wrong"])

    with pytest.raises(ValueError):
        orc_file.read(columns=[5])


def test_wrong_usage_orc_writer(tempdir):
    from pyarrow import orc

    path = str(tempdir / 'test.orc')
    with orc.ORCWriter(path) as writer:
        with pytest.raises(AttributeError):
            writer.test()


def test_orc_writer_with_null_arrays(tempdir):
    from pyarrow import orc

    path = str(tempdir / 'test.orc')
    a = pa.array([1, None, 3, None])
    b = pa.array([None, None, None, None])
    table = pa.table({"int64": a, "utf8": b})
    with pytest.raises(pa.ArrowNotImplementedError):
        orc.write_table(table, path)
