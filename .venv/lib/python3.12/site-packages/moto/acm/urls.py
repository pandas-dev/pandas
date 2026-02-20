from .responses import AWSCertificateManagerResponse

url_bases = [r"https?://acm\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": AWSCertificateManagerResponse.dispatch}
