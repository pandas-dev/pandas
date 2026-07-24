from .responses import EBResponse

url_bases = [
    r"https?://elasticbeanstalk\.(?P<region>[a-zA-Z0-9\-_]+)\.amazonaws.com",
]

url_paths = {
    "{0}/$": EBResponse.dispatch,
}
