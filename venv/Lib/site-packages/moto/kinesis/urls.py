from .responses import KinesisResponse

url_bases = [
    # Need to avoid conflicting with kinesisvideo
    r"https?://kinesis\.(.+)\.amazonaws\.com",
    # Somewhere around boto3-1.26.31 botocore-1.29.31, AWS started using a new endpoint:
    # 	111122223333.control-kinesis.us-east-1.amazonaws.com
    r"https?://(.+)\.control-kinesis\.(.+)\.amazonaws\.com",
    # When passing in the StreamARN to get_shard_iterator/get_records, this endpoint is called:
    r"https?://(.+)\.data-kinesis\.(.+)\.amazonaws\.com",
]

url_paths = {"{0}/$": KinesisResponse.dispatch}
