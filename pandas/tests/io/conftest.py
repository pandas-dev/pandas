import uuid

import pytest

from pandas.compat import (
    is_ci_environment,
    is_platform_arm,
    is_platform_mac,
    is_platform_windows,
)
import pandas.util._test_decorators as td

import pandas.io.common as icom
from pandas.io.parsers import read_csv


@pytest.fixture
def compression_to_extension():
    return {value: key for key, value in icom.extension_to_compression.items()}


@pytest.fixture
def tips_file(datapath):
    """Path to the tips dataset"""
    return datapath("io", "data", "csv", "tips.csv")


@pytest.fixture
def jsonl_file(datapath):
    """Path to a JSONL dataset"""
    return datapath("io", "parser", "data", "items.jsonl")


@pytest.fixture
def salaries_table(datapath):
    """DataFrame with the salaries dataset"""
    return read_csv(datapath("io", "parser", "data", "salaries.csv"), sep="\t")


@pytest.fixture
def feather_file(datapath):
    return datapath("io", "data", "feather", "feather-0_3_1.feather")


@pytest.fixture
def xml_file(datapath):
    return datapath("io", "data", "xml", "books.xml")


@pytest.fixture(scope="session")
def aws_credentials(monkeysession):
    """Mocked AWS Credentials for moto."""
    monkeysession.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeysession.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeysession.setenv("AWS_SECURITY_TOKEN", "testing")
    monkeysession.setenv("AWS_SESSION_AWS_SESSION_TOKEN", "testing")
    monkeysession.setenv("AWS_DEFAULT_REGION", "us-east-1")


@pytest.fixture(scope="session")
def moto_server(aws_credentials):
    # use service container for Linux on GitHub Actions
    if is_ci_environment() and not (
        is_platform_mac() or is_platform_arm() or is_platform_windows()
    ):
        yield "http://localhost:5000"
    else:
        moto_server = pytest.importorskip("moto.server")
        server = moto_server.ThreadedMotoServer(port=0)
        server.start()
        host, port = server.get_host_and_port()
        yield f"http://{host}:{port}"
        server.stop()


@pytest.fixture
def moto_s3_resource(moto_server):
    boto3 = pytest.importorskip("boto3")
    s3 = boto3.resource("s3", endpoint_url=moto_server)
    return s3


@pytest.fixture(scope="session")
def s3so(moto_server):
    return {
        "client_kwargs": {
            "endpoint_url": moto_server,
        }
    }


@pytest.fixture
def s3_bucket_public(moto_s3_resource):
    """
    Create a public S3 bucket using moto.
    """
    bucket_name = f"pandas-test-{uuid.uuid4()}"
    bucket = moto_s3_resource.Bucket(bucket_name)
    bucket.create(ACL="public-read")
    yield bucket
    bucket.objects.delete()
    bucket.delete()


@pytest.fixture
def s3_bucket_private(moto_s3_resource):
    """
    Create a private S3 bucket using moto.
    """
    bucket_name = f"cant_get_it-{uuid.uuid4()}"
    bucket = moto_s3_resource.Bucket(bucket_name)
    bucket.create(ACL="private")
    yield bucket
    bucket.objects.delete()
    bucket.delete()


@pytest.fixture
def s3_bucket_public_with_data(
    s3_bucket_public, tips_file, jsonl_file, feather_file, xml_file
):
    """
    The following datasets
    are loaded.

    - tips.csv
    - tips.csv.gz
    - tips.csv.bz2
    - items.jsonl
    """
    test_s3_files = [
        ("tips#1.csv", tips_file),
        ("tips.csv", tips_file),
        ("tips.csv.gz", tips_file + ".gz"),
        ("tips.csv.bz2", tips_file + ".bz2"),
        ("items.jsonl", jsonl_file),
        ("simple_dataset.feather", feather_file),
        ("books.xml", xml_file),
    ]
    for s3_key, file_name in test_s3_files:
        with open(file_name, "rb") as f:
            s3_bucket_public.put_object(Key=s3_key, Body=f)
    return s3_bucket_public


@pytest.fixture
def s3_bucket_private_with_data(
    s3_bucket_private, tips_file, jsonl_file, feather_file, xml_file
):
    """
    The following datasets
    are loaded.

    - tips.csv
    - tips.csv.gz
    - tips.csv.bz2
    - items.jsonl
    """
    test_s3_files = [
        ("tips#1.csv", tips_file),
        ("tips.csv", tips_file),
        ("tips.csv.gz", tips_file + ".gz"),
        ("tips.csv.bz2", tips_file + ".bz2"),
        ("items.jsonl", jsonl_file),
        ("simple_dataset.feather", feather_file),
        ("books.xml", xml_file),
    ]
    for s3_key, file_name in test_s3_files:
        with open(file_name, "rb") as f:
            s3_bucket_private.put_object(Key=s3_key, Body=f)
    return s3_bucket_private


_compression_formats_params = [
    (".no_compress", None),
    ("", None),
    (".gz", "gzip"),
    (".GZ", "gzip"),
    (".bz2", "bz2"),
    (".BZ2", "bz2"),
    (".zip", "zip"),
    (".ZIP", "zip"),
    (".xz", "xz"),
    (".XZ", "xz"),
    pytest.param((".zst", "zstd"), marks=td.skip_if_no("zstandard")),
    pytest.param((".ZST", "zstd"), marks=td.skip_if_no("zstandard")),
]


@pytest.fixture(params=_compression_formats_params[1:])
def compression_format(request):
    return request.param


@pytest.fixture(params=_compression_formats_params)
def compression_ext(request):
    return request.param[0]
