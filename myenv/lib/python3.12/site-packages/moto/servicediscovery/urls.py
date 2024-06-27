"""servicediscovery base URL and path."""

from .responses import ServiceDiscoveryResponse

url_bases = [
    r"https?://servicediscovery\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": ServiceDiscoveryResponse.dispatch,
}
