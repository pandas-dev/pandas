"""quicksight base URL and path."""

from .responses import QuickSightResponse

url_bases = [
    r"https?://quicksight\.(.+)\.amazonaws\.com",
]


url_paths = {
    r"{0}/accounts/(?P<account_id>[\d]+)/data-sets$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/data-sets/(?P<datasetid>[^/.]+)/ingestions/(?P<ingestionid>[^/.]+)$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/namespaces/(?P<namespace>[a-zA-Z0-9._-]+)/groups$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/namespaces/(?P<namespace>[a-zA-Z0-9._-]+)/groups/(?P<groupname>[^/]+)$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/namespaces/(?P<namespace>[a-zA-Z0-9._-]+)/groups/(?P<groupname>[^/]+)/members$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/namespaces/(?P<namespace>[a-zA-Z0-9._-]+)/groups/(?P<groupname>[^/]+)/members/(?P<username>[^/]+)$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/namespaces/(?P<namespace>[a-zA-Z0-9._-]+)/users$": QuickSightResponse.dispatch,
    r"{0}/accounts/(?P<account_id>[\d]+)/namespaces/(?P<namespace>[a-zA-Z0-9._-]+)/users/(?P<username>[^/]+)$": QuickSightResponse.dispatch,
}
