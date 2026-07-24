"""synthetics base URL and path."""

from .responses import SyntheticsResponse

url_bases = [
    r"https?://synthetics\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": SyntheticsResponse.dispatch,
    "{0}/canary$": SyntheticsResponse.dispatch,
    "{0}/canary/(?P<name>[^/]+)$": SyntheticsResponse.dispatch,
    "{0}/canaries$": SyntheticsResponse.dispatch,
    "{0}/tags/(?P<resourceArn>[^/]+)$": SyntheticsResponse.dispatch,
}
