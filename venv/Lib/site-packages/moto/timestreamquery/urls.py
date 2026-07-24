"""timestreamquery base URL and path."""

from moto.timestreamwrite.responses import TimestreamWriteResponse

url_bases = [
    r"https?://query\.timestream\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/?$": TimestreamWriteResponse.dispatch,
}
