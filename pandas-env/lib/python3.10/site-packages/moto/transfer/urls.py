"""transfer base URL and path."""

from .responses import TransferResponse

url_bases = [
    r"https?://transfer\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": TransferResponse.dispatch,
}
