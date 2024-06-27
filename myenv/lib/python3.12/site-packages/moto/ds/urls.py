"""ds base URL and path."""

from .responses import DirectoryServiceResponse

url_bases = [
    r"https?://ds\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": DirectoryServiceResponse.dispatch,
}
