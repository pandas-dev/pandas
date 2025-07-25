from .responses import ResourceGroupsResponse

url_bases = [r"https?://resource-groups(-fips)?\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/delete-group$": ResourceGroupsResponse.dispatch,
    "{0}/get-group$": ResourceGroupsResponse.dispatch,
    "{0}/get-group-configuration$": ResourceGroupsResponse.dispatch,
    "{0}/put-group-configuration$": ResourceGroupsResponse.dispatch,
    "{0}/get-group-query$": ResourceGroupsResponse.dispatch,
    "{0}/groups$": ResourceGroupsResponse.dispatch,
    "{0}/groups/(?P<resource_group_name>[^/]+)$": ResourceGroupsResponse.dispatch,
    "{0}/groups/(?P<resource_group_name>[^/]+)/query$": ResourceGroupsResponse.dispatch,
    "{0}/groups-list$": ResourceGroupsResponse.dispatch,
    "{0}/resources/(?P<resource_arn>[^/]+)/tags$": ResourceGroupsResponse.dispatch,
    "{0}/update-group$": ResourceGroupsResponse.dispatch,
    "{0}/update-group-query$": ResourceGroupsResponse.dispatch,
}
