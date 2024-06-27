"""route53resolver base URL and path."""

from .responses import Route53ResolverResponse

url_bases = [
    r"https?://route53resolver\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": Route53ResolverResponse.dispatch,
}
