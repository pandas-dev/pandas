"""scheduler base URL and path."""

from .responses import EventBridgeSchedulerResponse

url_bases = [
    r"https?://scheduler\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/schedules$": EventBridgeSchedulerResponse.dispatch,
    "{0}/schedules/(?P<name>[^/]+)$": EventBridgeSchedulerResponse.dispatch,
    "{0}/schedule-groups$": EventBridgeSchedulerResponse.dispatch,
    "{0}/schedule-groups/(?P<name>[^/]+)$": EventBridgeSchedulerResponse.dispatch,
    "{0}/tags/(?P<ResourceArn>.+)$": EventBridgeSchedulerResponse.dispatch,
    "{0}/tags/arn:aws:scheduler:(?P<region_name>[^/]+):(?P<account_id>[^/]+):schedule/(?P<group_name>[^/]+)/(?P<schedule_name>[^/]+)/?$": EventBridgeSchedulerResponse.dispatch,
}
