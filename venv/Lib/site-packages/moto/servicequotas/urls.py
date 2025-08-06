"""servicequotas base URL and path."""

from .responses import ServiceQuotasResponse

url_bases = [
    r"https?://servicequotas\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": ServiceQuotasResponse.dispatch,
}
