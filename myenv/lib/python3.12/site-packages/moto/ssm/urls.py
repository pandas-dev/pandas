from .responses import SimpleSystemManagerResponse

url_bases = [
    r"https?://ssm\.(.+)\.amazonaws\.com",
    r"https?://ssm\.(.+)\.amazonaws\.com\.cn",
]

url_paths = {"{0}/$": SimpleSystemManagerResponse.dispatch}
