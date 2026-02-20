"""Firehose base URL and path."""

from .responses import FirehoseResponse

url_bases = [r"https?://firehose\.(.+)\.amazonaws\.com"]
url_paths = {"{0}/$": FirehoseResponse.dispatch}
