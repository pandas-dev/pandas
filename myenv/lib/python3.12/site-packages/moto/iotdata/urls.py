from .responses import IoTDataPlaneResponse

url_bases = [
    r"https?://data\.iot\.(.+)\.amazonaws.com",
    r"https?://data-ats\.iot\.(.+)\.amazonaws.com",
]


response = IoTDataPlaneResponse()


url_paths = {
    "{0}/.*$": response.dispatch,
    #
    # (Flask) Paths for :class:`moto.core.models.ServerModeMockAWS`
    #
    "{0}/<path:route>$": response.dispatch,
}
