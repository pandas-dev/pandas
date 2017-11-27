import os

import moto
import pytest
from pandas.io.parsers import read_table

HERE = os.path.dirname(__file__)


@pytest.fixture(scope='module')
def tips_file():
    """Path to the tips dataset"""
    return os.path.join(HERE, 'parser', 'data', 'tips.csv')


@pytest.fixture(scope='module')
def jsonl_file():
    """Path a JSONL dataset"""
    return os.path.join(HERE, 'parser', 'data', 'items.jsonl')


@pytest.fixture(scope='module')
def salaries_table():
    """DataFrame with the salaries dataset"""
    path = os.path.join(HERE, 'parser', 'data', 'salaries.csv')
    return read_table(path)


@pytest.fixture(scope='module')
def s3_resource(tips_file, jsonl_file):
    """Fixture for mocking S3 interaction.

    The primary bucket name is "pandas-test". The following datasets
    are loaded.

    - tips.csv
    - tips.csv.gz
    - tips.csv.bz2
    - items.jsonl

    A private bucket "cant_get_it" is also created. The boto3 s3 resource
    is yielded by the fixture.
    """
    pytest.importorskip('s3fs')
    moto.mock_s3().start()

    test_s3_files = [
        ('tips.csv', tips_file),
        ('tips.csv.gz', tips_file + '.gz'),
        ('tips.csv.bz2', tips_file + '.bz2'),
        ('items.jsonl', jsonl_file),
    ]

    def add_tips_files(bucket_name):
        for s3_key, file_name in test_s3_files:
            with open(file_name, 'rb') as f:
                conn.Bucket(bucket_name).put_object(
                    Key=s3_key,
                    Body=f)

    boto3 = pytest.importorskip('boto3')
    # see gh-16135
    bucket = 'pandas-test'

    conn = boto3.resource("s3", region_name="us-east-1")
    conn.create_bucket(Bucket=bucket)
    add_tips_files(bucket)

    conn.create_bucket(Bucket='cant_get_it', ACL='private')
    add_tips_files('cant_get_it')

    yield conn

    moto.mock_s3().stop()
