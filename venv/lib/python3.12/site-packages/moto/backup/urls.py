"""backup base URL and path."""

from .responses import BackupResponse

url_bases = [
    r"https?://backup\.(.+)\.amazonaws\.com",
]


response = BackupResponse()


url_paths = {
    "{0}/backup/plans/?$": response.dispatch,
    "{0}/backup/plans/(?P<name>.+)/?$": response.dispatch,
    "{0}/backup-vaults/$": response.dispatch,
    "{0}/backup-vaults/(?P<name>[^/]+)$": response.dispatch,
    "{0}/tags/(?P<resource_arn>.+)$": response.dispatch,
    "{0}/untag/(?P<resource_arn>.+)$": response.dispatch,
}
