from .responses import IamResponse

url_bases = [r"https?://iam\.(.*\.)?amazonaws\.com"]

url_paths = {"{0}/$": IamResponse.dispatch}
