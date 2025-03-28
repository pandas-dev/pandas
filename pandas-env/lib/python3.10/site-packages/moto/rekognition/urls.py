"""rekognition base URL and path."""

from .responses import RekognitionResponse

url_bases = [
    r"https?://rekognition\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": RekognitionResponse.dispatch,
}
