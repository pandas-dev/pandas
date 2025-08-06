"""workspaces base URL and path."""

from .responses import WorkSpacesResponse

url_bases = [
    r"https?://workspaces\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": WorkSpacesResponse.dispatch,
}
