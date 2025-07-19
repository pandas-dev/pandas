from .responses import WAFV2Response

url_bases = [
    r"https?://wafv2\.(.+)\.amazonaws.com",
]

url_paths = {
    "{0}/": WAFV2Response.dispatch,
    "{0}/$": WAFV2Response.dispatch,
}
