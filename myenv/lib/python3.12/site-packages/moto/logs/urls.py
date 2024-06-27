from .responses import LogsResponse

url_bases = [r"https?://logs\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": LogsResponse.dispatch}
