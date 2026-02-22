"""account base URL and path."""

from .responses import AccountResponse

url_bases = [
    r"https?://account\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/putAlternateContact$": AccountResponse.dispatch,
    "{0}/getAlternateContact$": AccountResponse.dispatch,
    "{0}/deleteAlternateContact$": AccountResponse.dispatch,
}
