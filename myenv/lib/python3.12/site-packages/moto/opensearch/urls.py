"""opensearch base URL and path."""

from .responses import OpenSearchServiceResponse

url_bases = [r"https?://es\.(.+)\.amazonaws\.com"]


response = OpenSearchServiceResponse()


url_paths = {"{0}/.*$": response.dispatch}
