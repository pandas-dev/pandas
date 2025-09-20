"""servicecatalog base URL and path."""

from moto.servicecatalogappregistry.responses import AppRegistryResponse
from moto.servicecatalogappregistry.urls import url_paths as catalog_paths

from .responses import ServiceCatalogResponse

url_bases = [
    r"https?://servicecatalog\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": ServiceCatalogResponse.dispatch,
    **{path: AppRegistryResponse.dispatch for path in catalog_paths},
}
