"""networkmanager base URL and path."""

from .responses import NetworkManagerResponse

url_bases = [
    r"https?://networkmanager\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": NetworkManagerResponse.dispatch,
    "{0}/global-networks$": NetworkManagerResponse.dispatch,
    "{0}/core-networks$": NetworkManagerResponse.dispatch,
    "{0}/core-networks/(?P<networkid>[^/.]+)$": NetworkManagerResponse.dispatch,
    "{0}/global-networks/(?P<networkid>[^/.]+)$": NetworkManagerResponse.dispatch,
    "{0}/tags$": NetworkManagerResponse.dispatch,
    "{0}/tags/(?P<resourcearn>[^/.]+)$": NetworkManagerResponse.dispatch,
    "{0}/tags/(?P<arn_prefix>[^/]+)/(?P<resource_id>[^/]+)$": NetworkManagerResponse.dispatch,
}
