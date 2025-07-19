from .responses import KinesisVideoResponse

url_bases = [
    r"https?://kinesisvideo\.(.+)\.amazonaws.com",
]


response = KinesisVideoResponse()


url_paths = {
    "{0}/createStream$": response.dispatch,
    "{0}/describeStream$": response.dispatch,
    "{0}/deleteStream$": response.dispatch,
    "{0}/listStreams$": response.dispatch,
    "{0}/getDataEndpoint$": response.dispatch,
}
