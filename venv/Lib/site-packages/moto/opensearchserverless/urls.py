"""opensearchserverless base URL and path."""

from .responses import OpenSearchServiceServerlessResponse

url_bases = [
    r"https?://aoss\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": OpenSearchServiceServerlessResponse.dispatch,
}
