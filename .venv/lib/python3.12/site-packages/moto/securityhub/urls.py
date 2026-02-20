"""securityhub base URL and path."""

from .responses import SecurityHubResponse

url_bases = [
    r"https?://securityhub\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": SecurityHubResponse.dispatch,
    "{0}/accounts$": SecurityHubResponse.dispatch,
    "{0}/accounts/describe$": SecurityHubResponse.dispatch,
    "{0}/findings$": SecurityHubResponse.dispatch,
    "{0}/findings/import$": SecurityHubResponse.dispatch,
    "{0}/organization/admin/enable$": SecurityHubResponse.dispatch,
    "{0}/organization/configuration$": SecurityHubResponse.dispatch,
    "{0}/administrator$": SecurityHubResponse.dispatch,
    "{0}/hub$": SecurityHubResponse.dispatch,
    "{0}/hub/enable$": SecurityHubResponse.dispatch,
    "{0}/hub/disable$": SecurityHubResponse.dispatch,
    "{0}/members$": SecurityHubResponse.dispatch,
    "{0}/members/get$": SecurityHubResponse.dispatch,
    "{0}/master$": SecurityHubResponse.dispatch,
}
