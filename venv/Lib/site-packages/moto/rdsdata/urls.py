"""rdsdata base URL and path."""

from .responses import RDSDataServiceResponse

url_bases = [
    r"https?://rds-data\.(.+)\.amazonaws\.com",
]


response = RDSDataServiceResponse()


url_paths = {
    "{0}/Execute$": response.dispatch,
}
