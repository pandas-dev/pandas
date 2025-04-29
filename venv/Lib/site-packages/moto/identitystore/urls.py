"""identitystore base URL and path."""

from .responses import IdentityStoreResponse

url_bases = [
    r"https?://identitystore\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": IdentityStoreResponse.dispatch,
}
