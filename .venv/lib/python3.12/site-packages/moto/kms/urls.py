from .responses import KmsResponse

url_bases = [r"https?://kms\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": KmsResponse.dispatch}
