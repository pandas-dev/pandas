from .responses import GlacierResponse

url_bases = [r"https?://glacier\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/(?P<account_number>.+)/vaults$": GlacierResponse.dispatch,
    "{0}/(?P<account_number>.+)/vaults/(?P<vault_name>[^/]+)$": GlacierResponse.dispatch,
    "{0}/(?P<account_number>.+)/vaults/(?P<vault_name>.+)/archives$": GlacierResponse.dispatch,
    "{0}/(?P<account_number>.+)/vaults/(?P<vault_name>.+)/archives/(?P<archive_id>.+)$": GlacierResponse.dispatch,
    "{0}/(?P<account_number>.+)/vaults/(?P<vault_name>.+)/jobs$": GlacierResponse.dispatch,
    "{0}/(?P<account_number>.+)/vaults/(?P<vault_name>.+)/jobs/(?P<job_id>[^/.]+)$": GlacierResponse.dispatch,
    "{0}/(?P<account_number>.+)/vaults/(?P<vault_name>.+)/jobs/(?P<job_id>.+)/output$": GlacierResponse.dispatch,
}
