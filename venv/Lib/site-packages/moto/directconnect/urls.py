"""directconnect base URL and path."""

from .responses import DirectConnectResponse

url_bases = [
    r"https?://directconnect\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": DirectConnectResponse.dispatch,
}
