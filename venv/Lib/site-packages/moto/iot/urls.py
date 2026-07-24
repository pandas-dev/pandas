from .responses import IoTResponse

url_bases = [r"https?://iot\.(.+)\.amazonaws\.com"]


url_paths = {
    "{0}/.*$": IoTResponse.dispatch,
}
