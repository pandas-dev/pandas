from .responses import GlueResponse

url_bases = [r"https?://glue\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": GlueResponse.dispatch}
