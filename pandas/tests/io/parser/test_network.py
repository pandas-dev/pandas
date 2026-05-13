"""
Tests parsers ability to read and parse non-local files
and hence require a network connection to be read.
"""

from io import BytesIO
import logging
import re

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.feather_format import read_feather
from pandas.io.parsers import read_csv

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)


@pytest.mark.network
@pytest.mark.single_cpu
@pytest.mark.parametrize("mode", ["explicit", "infer"])
@pytest.mark.parametrize("engine", ["python", "c"])
def test_compressed_urls(
    httpserver,
    datapath,
    salaries_table,
    mode,
    engine,
    compression_only,
    compression_to_extension,
):
    # test reading compressed urls with various engines and
    # extension inference
    if compression_only == "tar":
        pytest.skip("TODO: Add tar salaries.csv to pandas/io/parsers/data")

    extension = compression_to_extension[compression_only]
    with open(datapath("io", "parser", "data", "salaries.csv" + extension), "rb") as f:
        httpserver.serve_content(content=f.read())

    url = httpserver.url + "/salaries.csv" + extension

    if mode != "explicit":
        compression_only = mode

    url_table = read_csv(url, sep="\t", compression=compression_only, engine=engine)
    tm.assert_frame_equal(url_table, salaries_table)


@pytest.mark.network
@pytest.mark.single_cpu
def test_url_encoding_csv(httpserver, datapath):
    """
    read_csv should honor the requested encoding for URLs.

    GH 10424
    """
    with open(datapath("io", "parser", "data", "unicode_series.csv"), "rb") as f:
        httpserver.serve_content(content=f.read())
        df = read_csv(httpserver.url, encoding="latin-1", header=None)
    assert df.loc[15, 1] == "Á köldum klaka (Cold Fever) (1994)"


@pytest.fixture
def tips_df(datapath):
    """DataFrame with the tips dataset."""
    return read_csv(datapath("io", "data", "csv", "tips.csv"))


@pytest.mark.single_cpu
@pytest.mark.network
@td.skip_if_not_us_locale()
class TestS3:
    @pytest.mark.parametrize(
        "suffix, compression",
        [
            ("", None),
            (".gz", "gzip"),
            (".bz2", "bz2"),
        ],
    )
    @pytest.mark.parametrize("nrows", [None, 10])
    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_parse_public_s3_bucket(
        self,
        s3_bucket_public_with_data,
        s3so,
        tips_df,
        suffix,
        compression,
        nrows,
        engine,
    ):
        # more of an integration test due to the not-public contents portion
        # can probably mock this though.
        pytest.importorskip("s3fs")
        df = read_csv(
            f"s3://{s3_bucket_public_with_data.name}/tips.csv{suffix}",
            nrows=nrows,
            compression=compression,
            storage_options=s3so,
            engine=engine,
        )
        tm.assert_frame_equal(df, tips_df.iloc[:nrows])

    def test_parse_private_s3_bucket(self, s3_bucket_private_with_data, s3so, tips_df):
        # Read public file from bucket with not-public contents
        pytest.importorskip("s3fs")
        df = read_csv(
            f"s3://{s3_bucket_private_with_data.name}/tips.csv", storage_options=s3so
        )
        tm.assert_frame_equal(df, tips_df)

    @pytest.mark.parametrize("scheme", ["s3n", "s3a"])
    def test_parse_public_bucket_s3n_s3a(
        self, s3_bucket_public_with_data, s3so, tips_df, scheme
    ):
        nrows = 10
        df = read_csv(
            f"{scheme}://{s3_bucket_public_with_data.name}/tips.csv",
            nrows=nrows,
            storage_options=s3so,
        )
        tm.assert_frame_equal(df, tips_df.iloc[:nrows])

    @pytest.mark.parametrize(
        "suffix, compression",
        [
            ("", None),
            (".gz", "gzip"),
            (".bz2", "bz2"),
        ],
    )
    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_parse_public_s3_bucket_chunked(
        self, s3_bucket_public_with_data, s3so, tips_df, suffix, compression, engine
    ):
        # Read with a chunksize
        chunksize = 5
        with read_csv(
            f"s3://{s3_bucket_public_with_data.name}/tips.csv{suffix}",
            chunksize=chunksize,
            compression=compression,
            storage_options=s3so,
            engine=engine,
        ) as df_reader:
            assert df_reader.chunksize == chunksize
            for i_chunk in [0, 1, 2]:
                # Read a couple of chunks and make sure we see them
                # properly.
                df = df_reader.get_chunk()
                assert isinstance(df, DataFrame)
                assert not df.empty
                true_df = tips_df.iloc[chunksize * i_chunk : chunksize * (i_chunk + 1)]
                tm.assert_frame_equal(true_df, df)

    @pytest.mark.parametrize("suffix", ["", ".gz", ".bz2"])
    def test_infer_s3_compression(
        self, s3_bucket_public_with_data, s3so, tips_df, suffix
    ):
        df = read_csv(
            f"s3://{s3_bucket_public_with_data.name}/tips.csv{suffix}",
            engine="python",
            compression="infer",
            storage_options=s3so,
        )
        tm.assert_frame_equal(df, tips_df)

    def test_read_s3_fails(self, s3_bucket_public_with_data, s3so):
        msg = "The specified bucket does not exist"
        with pytest.raises(OSError, match=msg):
            read_csv("s3://nyqpug/asdf.csv", storage_options=s3so)

    def test_read_s3_fails_private(self, s3_bucket_private_with_data, s3so):
        s3_url = f"{s3_bucket_private_with_data.name}/file.csv"
        msg = rf"{s3_url}"
        # Receive a permission error when trying to read a private bucket.
        # It's irrelevant here that this isn't actually a table.
        with pytest.raises(FileNotFoundError, match=msg):
            read_csv(
                f"s3://{s3_url}",
                storage_options=s3so,
            )

    @pytest.mark.single_cpu
    def test_read_csv_handles_boto_s3_object(
        self, s3_bucket_public_with_data, tips_file
    ):
        # see gh-16135

        s3_object = s3_bucket_public_with_data.Object("tips.csv")

        with BytesIO(s3_object.get()["Body"].read()) as buffer:
            result = read_csv(buffer, encoding="utf8")
        assert isinstance(result, DataFrame)
        assert not result.empty

        expected = read_csv(tips_file)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.single_cpu
    def test_read_csv_chunked_download(self, s3_bucket_public, s3so, caplog):
        # 8 MB, S3FS uses 5MB chunks
        df = DataFrame(np.zeros((100000, 4)), columns=list("abcd"))
        with BytesIO(df.to_csv().encode("utf-8")) as buf:
            s3_bucket_public.put_object(Key="large-file.csv", Body=buf)
            uri = f"{s3_bucket_public.name}/large-file.csv"
            match_re = re.compile(rf"^Fetch: {uri}, 0-(?P<stop>\d+)$")
            with caplog.at_level(logging.DEBUG, logger="s3fs"):
                read_csv(
                    f"s3://{uri}",
                    nrows=5,
                    storage_options=s3so,
                )
                for log in caplog.messages:
                    if match := re.match(match_re, log):
                        # Less than 8 MB
                        assert int(match.group("stop")) < 8000000

    def test_read_s3_with_hash_in_key(self, s3_bucket_public_with_data, s3so, tips_df):
        # GH 25945
        result = read_csv(
            f"s3://{s3_bucket_public_with_data.name}/tips#1.csv", storage_options=s3so
        )
        tm.assert_frame_equal(tips_df, result)

    def test_read_feather_s3_file_path(
        self, s3_bucket_public_with_data, s3so, feather_file
    ):
        # GH 29055
        pytest.importorskip("pyarrow")
        expected = read_feather(feather_file)
        res = read_feather(
            f"s3://{s3_bucket_public_with_data.name}/simple_dataset.feather",
            storage_options=s3so,
        )
        tm.assert_frame_equal(expected, res)
