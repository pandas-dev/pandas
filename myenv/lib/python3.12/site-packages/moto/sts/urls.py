from .responses import TokenResponse

url_bases = [r"https?://sts\.(.*\.)?amazonaws\.com"]

url_paths = {"{0}/$": TokenResponse.dispatch}
