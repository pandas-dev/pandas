from .responses import MediaLiveResponse

url_bases = [
    r"https?://medialive\.(.+)\.amazonaws.com",
]


response = MediaLiveResponse()


url_paths = {
    "{0}/prod/channels$": response.dispatch,
    "{0}/prod/channels/(?P<channelid>[^/.]+)$": response.dispatch,
    "{0}/prod/channels/(?P<channelid>[^/.]+)/start$": response.dispatch,
    "{0}/prod/channels/(?P<channelid>[^/.]+)/stop$": response.dispatch,
    "{0}/prod/inputs$": response.dispatch,
    "{0}/prod/inputs/(?P<inputid>[^/.]+)$": response.dispatch,
}
