from .responses import CloudWatchResponse

url_bases = [r"https?://monitoring\.(.+)\.amazonaws.com"]

url_paths = {"{0}/$": CloudWatchResponse.dispatch}
