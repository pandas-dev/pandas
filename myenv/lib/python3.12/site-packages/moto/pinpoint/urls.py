"""pinpoint base URL and path."""

from .responses import PinpointResponse

url_bases = [
    r"https?://pinpoint\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/v1/apps$": PinpointResponse.dispatch,
    "{0}/v1/apps/(?P<app_id>[^/]+)$": PinpointResponse.dispatch,
    "{0}/v1/apps/(?P<app_id>[^/]+)/eventstream": PinpointResponse.dispatch,
    "{0}/v1/apps/(?P<app_id>[^/]+)/settings$": PinpointResponse.dispatch,
    "{0}/v1/tags/(?P<app_arn>[^/]+)$": PinpointResponse.method_dispatch(
        PinpointResponse.tags
    ),
    "{0}/v1/tags/(?P<app_arn_pt_1>[^/]+)/(?P<app_arn_pt_2>[^/]+)$": PinpointResponse.method_dispatch(
        PinpointResponse.tags
    ),
}
