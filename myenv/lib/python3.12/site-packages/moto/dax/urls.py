"""dax base URL and path."""

from .responses import DAXResponse

url_bases = [
    r"https?://dax\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": DAXResponse.dispatch,
}
