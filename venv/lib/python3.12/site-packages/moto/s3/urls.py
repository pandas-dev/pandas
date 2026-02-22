from moto import settings

from .responses import S3Response

# Catch s3.amazonaws.com, but not s3-control.amazonaws.com
url_bases = [
    r"https?://s3(?!(-control|tables|vectors))(.*)\.amazonaws.com",
    r"https?://(?P<bucket_name>[a-zA-Z0-9\-_.]*)\.?s3(?!(-control|tables|vectors))(.*)\.amazonaws.com",
]

url_bases.extend(settings.get_s3_custom_endpoints())

url_paths = {
    # subdomain bucket
    "{0}/$": S3Response.bucket_dispatch,
    # subdomain key of path-based bucket
    "{0}/(?P<key_or_bucket_name>[^/]+)$": S3Response.ambiguous_dispatch,
    "{0}/(?P<key_or_bucket_name>[^/]+)/$": S3Response.ambiguous_dispatch,
    # path-based bucket + key
    "{0}/(?P<bucket_name_path>[^/]+)/(?P<key_name>.+)": S3Response.key_dispatch,
    # subdomain bucket + key with empty first part of path
    "{0}/(?P<key_name>/.*)$": S3Response.key_dispatch,
}
