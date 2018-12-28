import pytest

from pandas import read_csv
from pandas.compat import BytesIO
from pandas.io.common import is_s3_url


class TestS3URL(object):

    def test_is_s3_url(self):
        assert is_s3_url("s3://pandas/somethingelse.com")
        assert not is_s3_url("s4://pandas/somethingelse.com")


def test_streaming_s3_objects():
    # GH17135
    # botocore gained iteration support in 1.10.47, can now be used in read_*
    pytest.importorskip('botocore', minversion='1.10.47')
    from botocore.response import StreamingBody

    data = [
        b'foo,bar,baz\n1,2,3\n4,5,6\n',
        b'just,the,header\n',
    ]
    for el in data:
        body = StreamingBody(BytesIO(el), content_length=len(el))
        read_csv(body)
