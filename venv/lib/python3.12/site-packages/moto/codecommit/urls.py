from .responses import CodeCommitResponse

url_bases = [r"https?://codecommit\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": CodeCommitResponse.dispatch}
