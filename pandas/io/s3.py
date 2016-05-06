""" s3 support for remote file interactivity """
from pandas import compat
try:
    import s3fs
except:
    raise ImportError("The s3fs library is required to handle s3 files")

if compat.PY3:
    from urllib.parse import urlparse as parse_url
else:
    from urlparse import urlparse as parse_url


def _strip_schema(url):
    """Returns the url without the s3:// part"""
    result = parse_url(url)
    return result.netloc + result.path


def get_filepath_or_buffer(filepath_or_buffer, encoding=None,
                           compression=None):

    # Assuming AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY and AWS_S3_HOST
    # are environment variables
    fs = s3fs.S3FileSystem(anon=False)
    try:
        filepath_or_buffer = fs.open(_strip_schema(filepath_or_buffer))
    except OSError:
        # boto3 has troubles when trying to access a public file
        # when credentialed...
        fs = s3fs.S3FileSystem(anon=True)
        filepath_or_buffer = fs.open(_strip_schema(filepath_or_buffer))
    return filepath_or_buffer, None, compression
