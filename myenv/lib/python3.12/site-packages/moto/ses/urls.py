from .responses import EmailResponse

url_bases = [
    r"https?://email\.(.+)\.amazonaws\.com",
    r"https?://ses\.(.+)\.amazonaws\.com",
]

url_paths = {"{0}/$": EmailResponse.dispatch}
