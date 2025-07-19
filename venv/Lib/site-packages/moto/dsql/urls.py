"""dsql base URL and path."""

from .responses import AuroraDSQLResponse

url_bases = [
    r"https?://dsql\.(.+)\.api\.aws",
]

url_paths = {
    "{0}/cluster$": AuroraDSQLResponse.dispatch,
    "{0}/cluster/(?P<identifier>[^/]+)$": AuroraDSQLResponse.dispatch,
}
