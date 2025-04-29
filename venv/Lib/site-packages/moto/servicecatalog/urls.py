"""servicecatalog base URL and path."""

from .responses import ServiceCatalogResponse

url_bases = [
    r"https?://servicecatalog\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": ServiceCatalogResponse.dispatch,
}
