from .responses import MediaStoreResponse

url_bases = [r"https?://mediastore\.(.+)\.amazonaws\.com"]

response = MediaStoreResponse()

url_paths = {"{0}/$": response.dispatch, "{0}/(?P<Path>[^/.]+)$": response.dispatch}
