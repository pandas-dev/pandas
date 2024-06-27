"""appsync base URL and path."""

from .responses import AppSyncResponse

url_bases = [
    r"https?://appsync\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/v1/apis$": AppSyncResponse.dispatch,
    "{0}/v1/apis/(?P<api_id>[^/]+)$": AppSyncResponse.dispatch,
    "{0}/v1/apis/(?P<api_id>[^/]+)/apikeys$": AppSyncResponse.dispatch,
    "{0}/v1/apis/(?P<api_id>[^/]+)/apikeys/(?P<api_key_id>[^/]+)$": AppSyncResponse.dispatch,
    "{0}/v1/apis/(?P<api_id>[^/]+)/schemacreation$": AppSyncResponse.dispatch,
    "{0}/v1/apis/(?P<api_id>[^/]+)/schema$": AppSyncResponse.dispatch,
    "{0}/v1/tags/(?P<resource_arn>.+)$": AppSyncResponse.dispatch,
    "{0}/v1/tags/(?P<resource_arn_pt1>.+)/(?P<resource_arn_pt2>.+)$": AppSyncResponse.dispatch,
    "{0}/v1/apis/(?P<api_id>[^/]+)/types/(?P<type_name>.+)$": AppSyncResponse.dispatch,
}
