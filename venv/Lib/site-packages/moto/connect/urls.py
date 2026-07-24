from .responses import ConnectResponse

url_bases = [
    r"https?://connect\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/analytics-data/instance/(?P<InstanceId>[^/]+)/association$": ConnectResponse.dispatch,
    "{0}/instance/(?P<InstanceId>[^/]+)$": ConnectResponse.dispatch,
    "{0}/instance$": ConnectResponse.dispatch,
    "{0}/tags/(?P<resourceArn>.+)$": ConnectResponse.dispatch,
}
