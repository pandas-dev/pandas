from .responses import ConfigResponse

url_bases = [r"https?://config\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": ConfigResponse.dispatch}
