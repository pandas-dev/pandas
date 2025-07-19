"""personalize base URL and path."""

from .responses import PersonalizeResponse

url_bases = [
    r"https?://personalize\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": PersonalizeResponse.dispatch,
}
