"""connectcampaigns base URL and path."""

from .responses import ConnectCampaignServiceResponse

url_bases = [
    r"https?://connect-campaigns\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/tags/(?P<arn>.+)$": ConnectCampaignServiceResponse.dispatch,
    "{0}/tags/(?P<arn>[^/]+)?tagKeys=(?P<tagKeys>[^/]+)$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns-summary$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns/(?P<id>[^/]+)$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns/(?P<id>[^/]+)/start$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns/(?P<id>[^/]+)/stop$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns/(?P<id>[^/]+)/pause$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns/(?P<id>[^/]+)/resume$": ConnectCampaignServiceResponse.dispatch,
    "{0}/campaigns/(?P<id>[^/]+)/state$": ConnectCampaignServiceResponse.dispatch,
    "{0}/connect-instance/(?P<connectInstanceId>[^/]+)/config$": ConnectCampaignServiceResponse.dispatch,
    "{0}/connect-instance/(?P<connectInstanceId>[^/]+)/onboarding$": ConnectCampaignServiceResponse.dispatch,
}
