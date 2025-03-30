"""memorydb base URL and path."""

from .responses import MemoryDBResponse

url_bases = [
    r"https?://memory-db\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": MemoryDBResponse.dispatch,
}
