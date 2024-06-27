from .responses import TimestreamWriteResponse

url_bases = [
    r"https?://ingest\.timestream\.(.+)\.amazonaws\.com",
    r"https?://ingest\.timestream\.(.+)\.amazonaws\.com/",
]

response = TimestreamWriteResponse()

# Boto3 sends a request to 'https://ingest.timestream.amazonaws.com'
# Which means we need two url_paths - one without slash (for boto3), and one with (for other SDK's/API's)
url_paths = {"{0}$": response.dispatch, "{0}/$": response.dispatch}
