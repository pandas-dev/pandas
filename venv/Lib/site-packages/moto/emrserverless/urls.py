"""emrserverless base URL and path."""

from .responses import EMRServerlessResponse

url_bases = [
    r"https?://emr-serverless\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/applications$": EMRServerlessResponse.dispatch,
    "{0}/applications/(?P<applicationId>[^/]+)$": EMRServerlessResponse.dispatch,
    "{0}/applications/(?P<applicationId>[^/]+)/start$": EMRServerlessResponse.dispatch,
    "{0}/applications/(?P<applicationId>[^/]+)/stop$": EMRServerlessResponse.dispatch,
    "{0}/applications/(?P<applicationId>[^/]+)/jobruns$": EMRServerlessResponse.dispatch,
    "{0}/applications/(?P<applicationId>[^/]+)/jobruns/(?P<jobRunId>[^/]+)$": EMRServerlessResponse.dispatch,
}
