"""shield base URL and path."""

from .responses import ShieldResponse

url_bases = [
    r"https?://shield\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": ShieldResponse.dispatch,
}
