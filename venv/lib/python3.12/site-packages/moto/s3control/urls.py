"""s3control base URL and path."""

from .responses import S3ControlResponse

url_bases = [
    r"https?://([0-9]+)\.s3-control\.(.+)\.amazonaws\.com",
]


url_paths = {
    r"{0}/v20180820/accesspoint$": S3ControlResponse.dispatch,
    r"{0}/v20180820/configuration/publicAccessBlock$": S3ControlResponse.dispatch,
    r"{0}/v20180820/accesspoint/(?P<name>[\w_:%-]+)$": S3ControlResponse.dispatch,
    r"{0}/v20180820/accesspoint/(?P<name>[\w_:%-]+)/policy$": S3ControlResponse.dispatch,
    r"{0}/v20180820/accesspoint/(?P<name>[\w_:%-]+)/policyStatus$": S3ControlResponse.dispatch,
    r"{0}/v20180820/async-requests/mrap.*$": S3ControlResponse.dispatch,
    r"{0}/v20180820/mrap/instances.*$": S3ControlResponse.dispatch,
    "{0}/v20180820/storagelens/(?P<storagelensid>[^/]+)$": S3ControlResponse.dispatch,
    "{0}/v20180820/storagelens/(?P<storagelensid>[^/]+)/tagging$": S3ControlResponse.dispatch,
    "{0}/v20180820/storagelens$": S3ControlResponse.dispatch,
    r"{0}/v20180820/tags/(?P<arn>[\w_:%-]+)$": S3ControlResponse.dispatch,
    "/v20180820/accesspoint": S3ControlResponse.dispatch,
    "/v20180820/configuration/publicAccessBlock": S3ControlResponse.dispatch,
    "/v20180820/accesspoint/<name>": S3ControlResponse.dispatch,
    "/v20180820/accesspoint/<name>/policy": S3ControlResponse.dispatch,
    "/v20180820/accesspoint/<name>/policyStatus": S3ControlResponse.dispatch,
    "/v20180820/async-requests/mrap/create": S3ControlResponse.dispatch,
    "/v20180820/async-requests/mrap/delete": S3ControlResponse.dispatch,
    "/v20180820/async-requests/mrap/put-policy": S3ControlResponse.dispatch,
    "/v20180820/async-requests/mrap/<path:request_token_arn>": S3ControlResponse.dispatch,
    "/v20180820/mrap/instances": S3ControlResponse.dispatch,
    "/v20180820/mrap/instances/<mrap_name>": S3ControlResponse.dispatch,
    "/v20180820/mrap/instances/<mrap_name>/policy": S3ControlResponse.dispatch,
    "/v20180820/mrap/instances/<mrap_name>/policystatus": S3ControlResponse.dispatch,
    "/v20180820/storagelens/<storagelensid>": S3ControlResponse.dispatch,
    "/v20180820/storagelens/<storagelensid>/tagging": S3ControlResponse.dispatch,
    "/v20180820/storagelens": S3ControlResponse.dispatch,
    "/v20180820/tags/<arn>": S3ControlResponse.dispatch,
}
