from .responses import EC2Response

url_bases = [r"https?://ec2\.(.+)\.amazonaws\.com(|\.cn)"]

url_paths = {"{0}/": EC2Response.dispatch}
