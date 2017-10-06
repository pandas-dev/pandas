import pytest
import moto

import pandas as pd
from pandas import compat
import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal, assert_raises_regex


COMPRESSION_TYPES = [None, 'bz2', 'gzip', 'xz']


def decompress_file(path, compression):
    if compression is None:
        f = open(path, 'rb')
    elif compression == 'gzip':
        import gzip
        f = gzip.GzipFile(path, 'rb')
    elif compression == 'bz2':
        import bz2
        f = bz2.BZ2File(path, 'rb')
    elif compression == 'xz':
        lzma = compat.import_lzma()
        f = lzma.open(path, 'rb')
    else:
        msg = 'Unrecognized compression type: {}'.format(compression)
        raise ValueError(msg)

    result = f.read().decode('utf8')
    f.close()
    return result


@pytest.mark.parametrize('compression', COMPRESSION_TYPES)
def test_compression_roundtrip(compression):
    if compression == 'xz':
        tm._skip_if_no_lzma()

    df = pd.DataFrame([[0.123456, 0.234567, 0.567567],
                       [12.32112, 123123.2, 321321.2]],
                      index=['A', 'B'], columns=['X', 'Y', 'Z'])

    with tm.ensure_clean() as path:
        df.to_json(path, compression=compression)
        assert_frame_equal(df, pd.read_json(path, compression=compression))

        # explicitly ensure file was compressed.
        uncompressed_content = decompress_file(path, compression)
        assert_frame_equal(df, pd.read_json(uncompressed_content))


def test_compress_zip_value_error():
    df = pd.DataFrame([[0.123456, 0.234567, 0.567567],
                       [12.32112, 123123.2, 321321.2]],
                      index=['A', 'B'], columns=['X', 'Y', 'Z'])

    with tm.ensure_clean() as path:
        import zipfile
        pytest.raises(zipfile.BadZipfile, df.to_json, path, compression="zip")


def test_read_zipped_json():
    uncompressed_path = tm.get_data_path("tsframe_v012.json")
    uncompressed_df = pd.read_json(uncompressed_path)

    compressed_path = tm.get_data_path("tsframe_v012.json.zip")
    compressed_df = pd.read_json(compressed_path, compression='zip')

    assert_frame_equal(uncompressed_df, compressed_df)


@pytest.mark.parametrize('compression', COMPRESSION_TYPES)
def test_with_s3_url(compression):
    boto3 = pytest.importorskip('boto3')
    pytest.importorskip('s3fs')
    if compression == 'xz':
        tm._skip_if_no_lzma()

    df = pd.read_json('{"a": [1, 2, 3], "b": [4, 5, 6]}')
    with moto.mock_s3():
        conn = boto3.resource("s3", region_name="us-east-1")
        bucket = conn.create_bucket(Bucket="pandas-test")

        with tm.ensure_clean() as path:
            df.to_json(path, compression=compression)
            with open(path, 'rb') as f:
                bucket.put_object(Key='test-1', Body=f)

        roundtripped_df = pd.read_json('s3://pandas-test/test-1',
                                       compression=compression)
        assert_frame_equal(df, roundtripped_df)


@pytest.mark.parametrize('compression', COMPRESSION_TYPES)
def test_lines_with_compression(compression):
    if compression == 'xz':
        tm._skip_if_no_lzma()

    with tm.ensure_clean() as path:
        df = pd.read_json('{"a": [1, 2, 3], "b": [4, 5, 6]}')
        df.to_json(path, orient='records', lines=True, compression=compression)
        roundtripped_df = pd.read_json(path, lines=True,
                                       compression=compression)
        assert_frame_equal(df, roundtripped_df)


@pytest.mark.parametrize('compression', COMPRESSION_TYPES)
def test_chunksize_with_compression(compression):
    if compression == 'xz':
        tm._skip_if_no_lzma()

    with tm.ensure_clean() as path:
        df = pd.read_json('{"a": ["foo", "bar", "baz"], "b": [4, 5, 6]}')
        df.to_json(path, orient='records', lines=True, compression=compression)

        roundtripped_df = pd.concat(pd.read_json(path, lines=True, chunksize=1,
                                                 compression=compression))
        assert_frame_equal(df, roundtripped_df)


def test_write_unsupported_compression_type():
    df = pd.read_json('{"a": [1, 2, 3], "b": [4, 5, 6]}')
    with tm.ensure_clean() as path:
        msg = "Unrecognized compression type: unsupported"
        assert_raises_regex(ValueError, msg, df.to_json,
                            path, compression="unsupported")


def test_read_unsupported_compression_type():
    with tm.ensure_clean() as path:
        msg = "Unrecognized compression type: unsupported"
        assert_raises_regex(ValueError, msg, pd.read_json,
                            path, compression="unsupported")
