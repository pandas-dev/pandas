"""s3vectors base URL and path."""

from .responses import S3VectorsResponse

url_bases = [
    r"https?://s3vectors\.(.+)\.api\.aws",
]

url_paths = {
    "{0}/CreateVectorBucket$": S3VectorsResponse.dispatch,
    "{0}/DeleteVectorBucket$": S3VectorsResponse.dispatch,
    "{0}/GetVectorBucket$": S3VectorsResponse.dispatch,
    "{0}/ListVectorBuckets$": S3VectorsResponse.dispatch,
    "{0}/CreateIndex$": S3VectorsResponse.dispatch,
    "{0}/DeleteIndex$": S3VectorsResponse.dispatch,
    "{0}/GetIndex$": S3VectorsResponse.dispatch,
    "{0}/ListIndexes$": S3VectorsResponse.dispatch,
    "{0}/DeleteVectorBucketPolicy$": S3VectorsResponse.dispatch,
    "{0}/GetVectorBucketPolicy$": S3VectorsResponse.dispatch,
    "{0}/PutVectorBucketPolicy$": S3VectorsResponse.dispatch,
    "{0}/ListVectors$": S3VectorsResponse.dispatch,
    "{0}/DeleteVectors$": S3VectorsResponse.dispatch,
    "{0}/GetVectors$": S3VectorsResponse.dispatch,
    "{0}/PutVectors$": S3VectorsResponse.dispatch,
}
