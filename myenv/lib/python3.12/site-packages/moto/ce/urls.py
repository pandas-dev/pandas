"""ce base URL and path."""

from .responses import CostExplorerResponse

url_bases = [
    r"https?://ce\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": CostExplorerResponse.dispatch,
}
