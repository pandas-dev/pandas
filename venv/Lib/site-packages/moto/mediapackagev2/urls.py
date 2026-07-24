"""mediapackagev2 base URL and path."""

from .responses import mediapackagev2Response

url_bases = [
    r"https?://mediapackagev2\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/channelGroup$": mediapackagev2Response.dispatch,
    "{0}/channelGroup/(?P<ChannelGroupName>[^/]+)/channel$": mediapackagev2Response.dispatch,
    "{0}/channelGroup/(?P<ChannelGroupName>[^/]+)$": mediapackagev2Response.dispatch,
    "{0}/channelGroup/(?P<ChannelGroupName>[^/]+)/channel/(?P<ChannelName>[^/]+)/$": mediapackagev2Response.dispatch,
}
