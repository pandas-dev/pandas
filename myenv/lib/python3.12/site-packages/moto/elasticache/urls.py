"""elasticache base URL and path."""

from .responses import ElastiCacheResponse

url_bases = [
    r"https?://elasticache\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": ElastiCacheResponse.dispatch,
}
