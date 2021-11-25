import time

import pytest

from pandas.io.parsers import read_csv


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
def s3so():
    return {"client_kwargs": {"endpoint_url": "http://localhost:5000/"}}


@pytest.fixture()
def s3_resource(tips_file, jsonl_file, feather_file):
    """
    Sets up S3 bucket with contents

    The primary bucket name is "pandas-test". The following datasets
    are loaded.

    - tips.csv
    - tips.csv.gz
    - tips.csv.bz2
    - items.jsonl

    A private bucket "cant_get_it" is also created. The boto3 s3 resource
    is yielded by the fixture.
    """
    boto3 = pytest.importorskip("boto3")
    s3fs = pytest.importorskip("s3fs")

    # motoserver/moto CI service
    s3_base = "http://localhost:5000"

    test_s3_files = [
        ("tips#1.csv", tips_file),
        ("tips.csv", tips_file),
        ("tips.csv.gz", tips_file + ".gz"),
        ("tips.csv.bz2", tips_file + ".bz2"),
        ("items.jsonl", jsonl_file),
        ("simple_dataset.feather", feather_file),
    ]

    def add_tips_files(bucket_name):
        for s3_key, file_name in test_s3_files:
            with open(file_name, "rb") as f:
                cli.put_object(Bucket=bucket_name, Key=s3_key, Body=f)

    bucket = "pandas-test"
    conn = boto3.resource("s3", endpoint_url=s3_base)
    cli = boto3.client("s3", endpoint_url=s3_base)

    try:
        cli.create_bucket(Bucket=bucket)
    except Exception:
        # OK is bucket already exists
        pass
    try:
        cli.create_bucket(Bucket="cant_get_it", ACL="private")
    except Exception:
        # OK is bucket already exists
        pass
    timeout = 2
    while not cli.list_buckets()["Buckets"] and timeout > 0:
        time.sleep(0.1)
        timeout -= 0.1

    add_tips_files(bucket)
    add_tips_files("cant_get_it")
    s3fs.S3FileSystem.clear_instance_cache()
    yield conn

    s3 = s3fs.S3FileSystem(client_kwargs={"endpoint_url": s3_base})

    try:
        s3.rm(bucket, recursive=True)
    except Exception:
        pass
    try:
        s3.rm("cant_get_it", recursive=True)
    except Exception:
        pass
    timeout = 2
    while cli.list_buckets()["Buckets"] and timeout > 0:
        time.sleep(0.1)
        timeout -= 0.1
