"""s3control base URL and path."""

from .responses import S3ControlResponse

url_bases = [
    r"https?://([0-9]+)\.s3-control\.(.+)\.amazonaws\.com",
]


url_paths = {
    r"{0}/v20180820/configuration/publicAccessBlock$": S3ControlResponse.dispatch,
    r"{0}/v20180820/accesspoint/(?P<name>[\w_:%-]+)$": S3ControlResponse.dispatch,
    r"{0}/v20180820/accesspoint/(?P<name>[\w_:%-]+)/policy$": S3ControlResponse.dispatch,
    r"{0}/v20180820/accesspoint/(?P<name>[\w_:%-]+)/policyStatus$": S3ControlResponse.dispatch,
}
