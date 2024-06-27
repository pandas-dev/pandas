from .responses import MarketplaceMeteringResponse

url_bases = [
    "https?://metering.marketplace.(.+).amazonaws.com",
    "https?://aws-marketplace.(.+).amazonaws.com",
]

url_paths = {"{0}/$": MarketplaceMeteringResponse.dispatch}
