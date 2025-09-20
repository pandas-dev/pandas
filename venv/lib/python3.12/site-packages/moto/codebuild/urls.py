from .responses import CodeBuildResponse

url_bases = [r"https?://codebuild\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": CodeBuildResponse.dispatch}
