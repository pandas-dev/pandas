from .responses import OpsWorksResponse

# AWS OpsWorks has a single endpoint: opsworks.us-east-1.amazonaws.com
# and only supports HTTPS requests.
url_bases = [r"https?://opsworks\.us-east-1\.amazonaws.com"]

url_paths = {"{0}/$": OpsWorksResponse.dispatch}
