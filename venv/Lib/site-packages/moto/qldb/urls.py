"""qldb base URL and path."""

from .responses import QLDBResponse

url_bases = [
    r"https?://qldb\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/ledgers$": QLDBResponse.dispatch,
    "{0}/ledgers/(?P<ledger_name>[^/]+)$": QLDBResponse.dispatch,
    "{0}/tags/(?P<resource_arn>[^/]+)$": QLDBResponse.dispatch,
    "{0}/tags/(?P<resource_arn>[^/]+)/mock_ledger$": QLDBResponse.dispatch,
}
