"""mq base URL and path."""

from .responses import MQResponse

url_bases = [
    r"https?://mq\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/v1/brokers/(?P<broker_id>[^/]+)$": MQResponse.dispatch,
    "{0}/v1/brokers/(?P<broker_id>[^/]+)/reboot$": MQResponse.dispatch,
    "{0}/v1/brokers/(?P<broker_id>[^/]+)/users$": MQResponse.dispatch,
    "{0}/v1/brokers/(?P<broker_id>[^/]+)/users/(?P<user_name>[^/]+)$": MQResponse.dispatch,
    "{0}/v1/brokers$": MQResponse.dispatch,
    "{0}/v1/configurations$": MQResponse.dispatch,
    "{0}/v1/configurations/(?P<config_id>[^/]+)$": MQResponse.dispatch,
    "{0}/v1/configurations/(?P<config_id>[^/]+)/revisions/(?P<revision_id>[^/]+)$": MQResponse.dispatch,
    "{0}/v1/tags/(?P<resource_arn>[^/]+)$": MQResponse.dispatch,
}
