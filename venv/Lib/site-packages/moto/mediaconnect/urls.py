from .responses import MediaConnectResponse

url_bases = [
    r"https?://mediaconnect\.(.+)\.amazonaws.com",
]


response = MediaConnectResponse()


url_paths = {
    "{0}/v1/flows$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/vpcInterfaces$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/vpcInterfaces/(?P<vpcinterfacename>[^/.]+)$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/source$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/source/(?P<sourcearn>[^/.]+)$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/outputs$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/outputs/(?P<outputarn>[^/.]+)$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/entitlements$": response.dispatch,
    "{0}/v1/flows/(?P<flowarn>[^/.]+)/entitlements/(?P<entitlementarn>[^/.]+)$": response.dispatch,
    "{0}/v1/flows/start/(?P<flowarn>[^/.]+)$": response.dispatch,
    "{0}/v1/flows/stop/(?P<flowarn>[^/.]+)$": response.dispatch,
    "{0}/tags/(?P<resourcearn>[^/.]+)$": response.dispatch,
}
