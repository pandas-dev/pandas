"""acmpca base URL and path."""

from .responses import ACMPCAResponse

url_bases = [
    r"https?://acm-pca\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": ACMPCAResponse.dispatch,
}
