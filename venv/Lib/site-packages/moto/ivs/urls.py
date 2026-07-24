"""ivs base URL and path."""

from .responses import IVSResponse

url_bases = [
    r"https?://ivs\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/CreateChannel": IVSResponse.dispatch,
    "{0}/ListChannels": IVSResponse.dispatch,
    "{0}/GetChannel": IVSResponse.dispatch,
    "{0}/BatchGetChannel": IVSResponse.dispatch,
    "{0}/UpdateChannel": IVSResponse.dispatch,
    "{0}/DeleteChannel": IVSResponse.dispatch,
}
