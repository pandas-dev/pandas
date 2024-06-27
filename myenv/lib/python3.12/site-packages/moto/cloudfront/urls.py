"""cloudfront base URL and path."""

from .responses import CloudFrontResponse

url_bases = [
    r"https?://cloudfront\.amazonaws\.com",
    r"https?://cloudfront\.(.+)\.amazonaws\.com",
]
url_paths = {
    "{0}/2020-05-31/distribution$": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/distribution/(?P<distribution_id>[^/]+)$": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/distribution/(?P<distribution_id>[^/]+)/config$": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/distribution/(?P<distribution_id>[^/]+)/invalidation": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/tagging$": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/origin-access-control$": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/origin-access-control/(?P<oac_id>[^/]+)$": CloudFrontResponse.dispatch,
    "{0}/2020-05-31/origin-access-control/(?P<oac_id>[^/]+)/config$": CloudFrontResponse.dispatch,
}
