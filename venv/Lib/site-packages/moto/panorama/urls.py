from .responses import PanoramaResponse

url_bases = [
    r"https?://panorama\.(.+)\.amazonaws.com",
]

url_paths = {
    "{0}/$": PanoramaResponse.dispatch,
    "{0}/devices$": PanoramaResponse.dispatch,
    "{0}/devices/(?P<DeviceId>[^/]+)$": PanoramaResponse.dispatch,
    "{0}/packages/template-job$": PanoramaResponse.dispatch,
    "{0}/packages/template-job/(?P<JobId>[^/]+)$": PanoramaResponse.dispatch,
    "{0}/nodes$": PanoramaResponse.dispatch,
    "{0}/application-instances$": PanoramaResponse.dispatch,
    "{0}/application-instances/(?P<ApplicationInstanceId>[^/]+)$": PanoramaResponse.dispatch,
    "{0}/application-instances/(?P<ApplicationInstanceId>[^/]+)/details$": PanoramaResponse.dispatch,
}
