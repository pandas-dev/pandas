import pytest
import pandas as pd

from pandas.io.common import _is_s3_url
from pandas.util import testing as tm

class TestS3URL(object):

    def test_is_s3_url(self):
        assert _is_s3_url("s3://pandas/somethingelse.com")
        assert not _is_s3_url("s4://pandas/somethingelse.com")

class TestIntegration(object):
    def test_s3_roundtrip(self):
        expected = pd.DataFrame({'A': [1, 2, 3], 'B': 'foo'})

        boto3 = pytest.importorskip('boto3')
        moto = pytest.importorskip('moto')

        with moto.mock_s3():
            conn = boto3.resource("s3", region_name="us-east-1")
            conn.create_bucket(Bucket="pandas-test")

            expected.to_parquet('s3://pandas-test/test.parquet')
            actual = pd.read_parquet('s3://pandas-test/test.parquet')

            tm.assert_frame_equal(actual, expected)
