"""servicecatalogappregistry base URL and path."""

from .responses import AppRegistryResponse

url_bases = [
    r"https?://servicecatalog-appregistry\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/applications$": AppRegistryResponse.dispatch,
    "{0}/applications/(?P<application>[^/]+)/resources/(?P<resource_type>[^/]+)/(?P<resource>.+)$": AppRegistryResponse.dispatch,
    "{0}/applications/(?P<application>.+)/resources$": AppRegistryResponse.dispatch,
}
