from io import BytesIO

import pytest

from pandas import read_csv


@pytest.mark.parametrize("data", [b"foo,bar,baz\n1,2,3\n4,5,6\n", b"just,the,header\n"])
def test_streaming_s3_objects(data):
    # GH 17135
    # botocore gained iteration support in 1.10.47, can now be used in read_*
    pytest.importorskip("botocore", minversion="1.10.47")
    from botocore.response import StreamingBody

    body = StreamingBody(BytesIO(data), content_length=len(data))
    read_csv(body)


@pytest.mark.single_cpu
@pytest.mark.parametrize("header", ["infer", None])
def test_read_with_and_without_creds_from_pub_bucket(
    s3_bucket_public_with_data, s3so, header
):
    # GH 34626
    pytest.importorskip("s3fs")
    nrows = 5
    df = read_csv(
        f"s3://{s3_bucket_public_with_data.name}/tips.csv",
        nrows=nrows,
        header=header,
        storage_options=s3so,
    )
    assert len(df) == nrows
