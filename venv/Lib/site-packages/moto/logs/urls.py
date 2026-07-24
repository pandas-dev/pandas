from .responses import LogsResponse

url_bases = [
    r"https?://logs\.(.+)\.amazonaws\.com",
    r"https?://stream-logs\.(.+)\.amazonaws\.com",
    r"https?://streaming-logs\.(.+)\.amazonaws\.com",
]

url_paths = {"{0}/$": LogsResponse.dispatch}
