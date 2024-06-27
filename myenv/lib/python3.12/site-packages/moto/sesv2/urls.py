"""sesv2 base URL and path."""

from .responses import SESV2Response

url_bases = [
    r"https?://email\.(.+)\.amazonaws\.com",
]


response = SESV2Response()


url_paths = {
    "{0}/v2/email/outbound-emails$": response.dispatch,
    "{0}/v2/email/contact-lists/(?P<name>[^/]+)$": response.dispatch,
    "{0}/v2/email/contact-lists/(?P<name>[^/]+)/contacts$": response.dispatch,
    "{0}/v2/email/contact-lists/(?P<name>[^/]+)/contacts/(?P<email>[^/]+)$": response.dispatch,
    "{0}/v2/email/contact-lists$": response.dispatch,
    "{0}/v2/.*$": response.dispatch,
}
