from .responses import MediaPackageResponse

url_bases = [
    r"https?://mediapackage\.(.+)\.amazonaws.com",
]


response = MediaPackageResponse()


url_paths = {
    "{0}/channels$": response.dispatch,
    "{0}/channels/(?P<channelid>[^/.]+)$": response.dispatch,
    "{0}/origin_endpoints$": response.dispatch,
    "{0}/origin_endpoints/(?P<id>[^/.]+)$": response.dispatch,
}
