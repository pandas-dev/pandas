"""route53domains base URL and path."""

from .responses import Route53DomainsResponse

url_bases = [
    r"https?://route53domains\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": Route53DomainsResponse.dispatch,
}
