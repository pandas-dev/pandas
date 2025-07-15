"""textract base URL and path."""

from .responses import TextractResponse

url_bases = [
    r"https?://textract\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": TextractResponse.dispatch,
}
