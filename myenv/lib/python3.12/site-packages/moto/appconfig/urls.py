"""appconfig base URL and path."""

from .responses import AppConfigResponse

url_bases = [
    r"https?://appconfig\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/applications$": AppConfigResponse.dispatch,
    "{0}/applications/(?P<app_id>[^/]+)$": AppConfigResponse.dispatch,
    "{0}/applications/(?P<app_id>[^/]+)/configurationprofiles$": AppConfigResponse.dispatch,
    "{0}/applications/(?P<app_id>[^/]+)/configurationprofiles/(?P<config_profile_id>[^/]+)$": AppConfigResponse.dispatch,
    "{0}/applications/(?P<app_id>[^/]+)/configurationprofiles/(?P<config_profile_id>[^/]+)/hostedconfigurationversions$": AppConfigResponse.dispatch,
    "{0}/applications/(?P<app_id>[^/]+)/configurationprofiles/(?P<config_profile_id>[^/]+)/hostedconfigurationversions/(?P<version>[^/]+)$": AppConfigResponse.dispatch,
    "{0}/tags/(?P<app_id>.+)$": AppConfigResponse.dispatch,
    "{0}/tags/(?P<arn_part_1>[^/]+)/(?P<app_id>[^/]+)$": AppConfigResponse.method_dispatch(
        AppConfigResponse.tags  # type: ignore
    ),
    "{0}/tags/(?P<arn_part_1>[^/]+)/(?P<app_id>[^/]+)/configurationprofile/(?P<cp_id>[^/]+)$": AppConfigResponse.method_dispatch(
        AppConfigResponse.tags  # type: ignore
    ),
}
