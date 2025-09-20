"""comprehend base URL and path."""

from .responses import ComprehendResponse

url_bases = [
    r"https?://comprehend\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": ComprehendResponse.dispatch,
}
