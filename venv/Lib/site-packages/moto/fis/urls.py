"""fis base URL and path."""

from .responses import FISResponse

url_bases = [
    r"https?://fis\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/experimentTemplates$": FISResponse.dispatch,
    "{0}/experimentTemplates/(?P<id>.+)$": FISResponse.dispatch,
    "{0}/tags/(?P<resourceArn>.+)$": FISResponse.dispatch,
}
