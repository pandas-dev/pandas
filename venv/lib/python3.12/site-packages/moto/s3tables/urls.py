"""s3tables base URL and path."""

from .responses import S3TablesResponse

url_bases = [
    r"https?://s3tables\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/buckets$": S3TablesResponse.dispatch,
    "{0}/buckets/(?P<tableBucketARN>.+)$": S3TablesResponse.dispatch,
    "{0}/buckets/(?P<tableBucketARN_pt_1>[^/]+)/(?P<tableBucketARN_pt_2>[^/]+)$": S3TablesResponse.dispatch,
    "{0}/get-table$": S3TablesResponse.dispatch,
    "{0}/namespaces/(?P<tableBucketARN>.+)$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN>[^/]+)/(?P<namespace>[^/]+)$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN>[^/]+)/(?P<namespace>[^/]+)/(?P<name>[^/]+)$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN>[^/]+)$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN_pt_1>[^/]+)/(?P<tableBucketARN_pt_2>[^/]+)/(?P<namespace>[^/]+)/(?P<name>[^/]+)$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN>[^/]+)/(?P<namespace>[^/]+)/(?P<name>[^/]+)/metadata-location$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN_pt_1>[^/]+)/(?P<tableBucketARN_pt_2>[^/]+)/(?P<namespace>[^/]+)/(?P<name>[^/]+)/metadata-location$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN>[^/]+)/(?P<namespace>[^/]+)/(?P<name>[^/]+)/rename$": S3TablesResponse.dispatch,
    "{0}/tables/(?P<tableBucketARN_pt_1>[^/]+)/(?P<tableBucketARN_pt_2>[^/]+)/(?P<namespace>[^/]+)/(?P<name>[^/]+)/rename$": S3TablesResponse.dispatch,
}
