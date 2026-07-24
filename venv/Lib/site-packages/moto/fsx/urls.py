"""fsx base URL and path."""

from .responses import FSxResponse

url_bases = [
    r"https?://fsx\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": FSxResponse.dispatch,
}
