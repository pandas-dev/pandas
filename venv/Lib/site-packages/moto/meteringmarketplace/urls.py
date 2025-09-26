from .responses import MarketplaceMeteringResponse

url_bases = [
    r"https?://metering\.marketplace\.(.+)\.amazonaws\.com",
    r"https?://aws-marketplace\.(.+)\.amazonaws\.com",
]

url_paths = {"{0}/$": MarketplaceMeteringResponse.dispatch}
