"""cloudtrail base URL and path."""

from .responses import CloudTrailResponse

response = CloudTrailResponse()

url_bases = [
    r"https?://cloudtrail\.(.+)\.amazonaws\.com",
]


url_paths = {"{0}/$": response.dispatch}
