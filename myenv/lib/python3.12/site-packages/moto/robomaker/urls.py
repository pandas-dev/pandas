from .responses import RoboMakerResponse

url_bases = [
    r"https?://robomaker\.(.+)\.amazonaws\.com",
]


response = RoboMakerResponse()


url_paths = {
    "{0}/createRobotApplication$": response.dispatch,
    "{0}/deleteRobotApplication$": response.dispatch,
    "{0}/describeRobotApplication$": response.dispatch,
    "{0}/listRobotApplications$": response.dispatch,
}
