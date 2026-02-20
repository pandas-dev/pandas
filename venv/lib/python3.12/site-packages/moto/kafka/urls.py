"""kafka base URL and path."""

from .responses import KafkaResponse

url_bases = [
    r"https?://kafka\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/api/v2/clusters$": KafkaResponse.dispatch,
    "{0}/api/v2/clusters/(?P<clusterArn>.+)$": KafkaResponse.dispatch,
    "{0}/v1/tags/(?P<resourceArn>.+)$": KafkaResponse.dispatch,
    "{0}/v1/clusters$": KafkaResponse.dispatch,
    "{0}/v1/clusters/(?P<clusterArn>.+)$": KafkaResponse.dispatch,
}
