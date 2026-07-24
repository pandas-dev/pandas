"""amp base URL and path."""

from .responses import PrometheusServiceResponse

url_bases = [
    r"https?://aps\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/workspaces$": PrometheusServiceResponse.dispatch,
    "{0}/workspaces/(?P<workspace_id>[^/]+)$": PrometheusServiceResponse.dispatch,
    "{0}/workspaces/(?P<workspace_id>[^/]+)/alias$": PrometheusServiceResponse.dispatch,
    "{0}/workspaces/(?P<workspace_id>[^/]+)/logging$": PrometheusServiceResponse.dispatch,
    "{0}/workspaces/(?P<workspace_id>[^/]+)/rulegroupsnamespaces$": PrometheusServiceResponse.dispatch,
    "{0}/workspaces/(?P<workspace_id>[^/]+)/rulegroupsnamespaces/(?P<name>[^/]+)$": PrometheusServiceResponse.dispatch,
    "{0}/tags/(?P<resource_arn>[^/]+)$": PrometheusServiceResponse.dispatch,
    "{0}/tags/(?P<arn_prefix>[^/]+)/(?P<workspace_id>[^/]+)$": PrometheusServiceResponse.dispatch,
    "{0}/tags/(?P<arn_prefix>[^/]+)/(?P<workspace_id>[^/]+)/(?P<ns_name>[^/]+)$": PrometheusServiceResponse.dispatch,
}
