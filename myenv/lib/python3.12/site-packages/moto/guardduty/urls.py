from .responses import GuardDutyResponse

url_bases = [
    "https?://guardduty\\.(.+)\\.amazonaws\\.com",
]


url_paths = {
    "{0}/detector$": GuardDutyResponse.dispatch,
    "{0}/detector/(?P<detector_id>[^/]+)$": GuardDutyResponse.dispatch,
    "{0}/detector/(?P<detector_id>[^/]+)/filter$": GuardDutyResponse.dispatch,
    "{0}/detector/(?P<detector_id>[^/]+)/filter/(?P<filter_name>[^/]+)$": GuardDutyResponse.dispatch,
    "{0}/admin/enable$": GuardDutyResponse.dispatch,
    "{0}/admin$": GuardDutyResponse.dispatch,
}
