from .responses import DataSyncResponse

url_bases = [r"https?://(.*\.)?(datasync)\.(.+)\.amazonaws.com"]

url_paths = {"{0}/$": DataSyncResponse.dispatch}
