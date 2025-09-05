"""cloudhsmv2 base URL and path."""

from .responses import CloudHSMV2Response

url_bases = [
    r"https?://cloudhsm\.(.+)\.amazonaws\.com",
    r"https?://cloudhsmv2\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": CloudHSMV2Response.dispatch,
}
