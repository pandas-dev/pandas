from .responses import GreengrassResponse

url_bases = [
    r"https?://greengrass\.(.+)\.amazonaws.com",
]

url_paths = {
    "{0}/greengrass/definition/cores$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/cores/(?P<definition_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/cores/(?P<definition_id>[^/]+)/versions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/cores/(?P<definition_id>[^/]+)/versions/(?P<definition_version_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/devices$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/devices/(?P<definition_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/devices/(?P<definition_id>[^/]+)/versions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/devices/(?P<definition_id>[^/]+)/versions/(?P<definition_version_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/functions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/functions/(?P<definition_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/functions/(?P<definition_id>[^/]+)/versions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/functions/(?P<definition_id>[^/]+)/versions/(?P<definition_version_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/resources$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/resources/(?P<definition_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/resources/(?P<definition_id>[^/]+)/versions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/subscriptions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/subscriptions/(?P<definition_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/subscriptions/(?P<definition_id>[^/]+)/versions$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/subscriptions/(?P<definition_id>[^/]+)/versions/(?P<definition_version_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/definition/resources/(?P<definition_id>[^/]+)/versions/(?P<definition_version_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<group_id>[^/]+)/?$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<group_id>[^/]+)/role$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<group_id>[^/]+)/versions$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<group_id>[^/]+)/deployments$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<GroupId>[^/]+)/deployments/\\$reset$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<group_id>[^/]+)/deployments/(?P<group_version_id>[^/]+)/status$": GreengrassResponse.dispatch,
    "{0}/greengrass/groups/(?P<group_id>[^/]+)/versions/(?P<group_version_id>[^/]+)/?$": GreengrassResponse.dispatch,
}
