from .responses import MacieResponse

url_bases = [
    r"https?://macie2\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/invitations$": MacieResponse.dispatch,
    "{0}/invitations/decline$": MacieResponse.dispatch,
    "{0}/invitations/accept$": MacieResponse.dispatch,
    "{0}/members$": MacieResponse.dispatch,
    "{0}/members/[^/]+$": MacieResponse.dispatch,
    "{0}/members/(?P<id>[^/]+)$": MacieResponse.dispatch,
    "{0}/administrator$": MacieResponse.dispatch,
    "{0}/macie$": MacieResponse.dispatch,
    "{0}/session$": MacieResponse.dispatch,
    "{0}/members/disassociate/(?P<id>[^/]+)$": MacieResponse.dispatch,
    "{0}/admin$": MacieResponse.dispatch,
}
