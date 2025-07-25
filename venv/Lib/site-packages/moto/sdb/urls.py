from .responses import SimpleDBResponse

url_bases = [
    r"https?://sdb\.(.+)\.amazonaws\.com",
]

url_paths = {"{0}/$": SimpleDBResponse.dispatch}
