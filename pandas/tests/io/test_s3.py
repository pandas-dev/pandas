from io import BytesIO
import os

import pytest

from pandas import read_csv
import pandas._testing as tm

from pandas.io.common import is_s3_url


class TestS3URL:
    def test_is_s3_url(self):
        assert is_s3_url("s3://pandas/somethingelse.com")
        assert not is_s3_url("s4://pandas/somethingelse.com")


def test_streaming_s3_objects():
    # GH17135
    # botocore gained iteration support in 1.10.47, can now be used in read_*
    pytest.importorskip("botocore", minversion="1.10.47")
    from botocore.response import StreamingBody

    data = [b"foo,bar,baz\n1,2,3\n4,5,6\n", b"just,the,header\n"]
    for el in data:
        body = StreamingBody(BytesIO(el), content_length=len(el))
        read_csv(body)


@tm.network
@pytest.mark.slow
def test_read_s3_public():
    # ensure we can read from a public bucket with credentials
    pytest.importorskip("s3fs")

    with tm.ensure_safe_environment_variables():
        # temporary workaround as moto fails for botocore >= 1.11 otherwise,
        # see https://github.com/spulec/moto/issues/1924 & 1952
        os.environ.setdefault("AWS_ACCESS_KEY_ID", "foobar_key")
        os.environ.setdefault("AWS_SECRET_ACCESS_KEY", "foobar_secret")
        df = read_csv(
            "s3://gdelt-open-data/events/20130420.export.csv",
            nrows=5,
            sep="\t",
            header=None,
        )
        assert len(df) == 5
