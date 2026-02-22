"""pipes base URL and path."""

from .responses import EventBridgePipesResponse

url_bases = [
    r"https?://pipes\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/v1/pipes/(?P<Name>[^/]+)$": EventBridgePipesResponse.dispatch,
    "{0}/tags/(?P<resourceArn>.+)$": EventBridgePipesResponse.dispatch,
    "{0}/v1/pipes$": EventBridgePipesResponse.dispatch,
}
