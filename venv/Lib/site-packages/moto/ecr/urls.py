from .responses import ECRResponse

url_bases = [
    r"https?://ecr\.(.+)\.amazonaws\.com",
    r"https?://api\.ecr\.(.+)\.amazonaws\.com",
]

url_paths = {"{0}/$": ECRResponse.dispatch}
