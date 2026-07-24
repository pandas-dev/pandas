"""ssoadmin base URL and path."""

from .responses import SSOAdminResponse

url_bases = [
    r"https?://sso\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": SSOAdminResponse.dispatch,
}
