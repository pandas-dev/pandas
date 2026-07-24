"""emrcontainers base URL and path."""

from .responses import EMRContainersResponse

url_bases = [
    r"https?://emr-containers\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/virtualclusters$": EMRContainersResponse.dispatch,
    "{0}/virtualclusters/(?P<virtualClusterId>[^/]+)$": EMRContainersResponse.dispatch,
    "{0}/virtualclusters/(?P<virtualClusterId>[^/]+)/jobruns$": EMRContainersResponse.dispatch,
    "{0}/virtualclusters/(?P<virtualClusterId>[^/]+)/jobruns/(?P<jobRunId>[^/]+)$": EMRContainersResponse.dispatch,
}
