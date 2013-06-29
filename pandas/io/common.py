"""Common IO api utilities"""

import sys
import urlparse
import urllib2
import zipfile
from contextlib import contextmanager, closing
from StringIO import StringIO

from pandas.util import py3compat

_VALID_URLS = set(urlparse.uses_relative + urlparse.uses_netloc +
                  urlparse.uses_params)
_VALID_URLS.discard('')


class PerformanceWarning(Warning):
    pass


def _is_url(url):
    """Check to see if a URL has a valid protocol.

    Parameters
    ----------
    url : str or unicode

    Returns
    -------
    isurl : bool
        If `url` has a valid protocol return True otherwise False.
    """
    try:
        return urlparse.urlparse(url).scheme in _VALID_URLS
    except:
        return False


def _is_s3_url(url):
    """Check for an s3 url"""
    try:
        return urlparse.urlparse(url).scheme == 's3'
    except:
        return False


def get_filepath_or_buffer(filepath_or_buffer, encoding=None):
    """
    If the filepath_or_buffer is a url, translate and return the buffer
    passthru otherwise.

    Parameters
    ----------
    filepath_or_buffer : a url, filepath, or buffer
    encoding : the encoding to use to decode py3 bytes, default is 'utf-8'

    Returns
    -------
    a filepath_or_buffer, the encoding
    """

    if _is_url(filepath_or_buffer):
        from urllib2 import urlopen
        filepath_or_buffer = urlopen(filepath_or_buffer)
        if py3compat.PY3:  # pragma: no cover
            if encoding:
                errors = 'strict'
            else:
                errors = 'replace'
                encoding = 'utf-8'
            bytes = filepath_or_buffer.read().decode(encoding, errors)
            filepath_or_buffer = StringIO(bytes)
            return filepath_or_buffer, encoding
        return filepath_or_buffer, None

    if _is_s3_url(filepath_or_buffer):
        try:
            import boto
        except:
            raise ImportError("boto is required to handle s3 files")
        # Assuming AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY
        # are environment variables
        parsed_url = urlparse.urlparse(filepath_or_buffer)
        conn = boto.connect_s3()
        b = conn.get_bucket(parsed_url.netloc)
        k = boto.s3.key.Key(b)
        k.key = parsed_url.path
        filepath_or_buffer = StringIO(k.get_contents_as_string())
        return filepath_or_buffer, None

    return filepath_or_buffer, None


# ----------------------
# Prevent double closing
if py3compat.PY3:
    urlopen = urllib2.urlopen
else:
    @contextmanager
    def urlopen(*args, **kwargs):
        with closing(urllib2.urlopen(*args, **kwargs)) as f:
            yield f

# ZipFile is not a context manager for <= 2.6
# must be tuple index here since 2.6 doesn't use namedtuple for version_info
if sys.version_info[1] <= 6:
    @contextmanager
    def ZipFile(*args, **kwargs):
        with closing(zipfile.ZipFile(*args, **kwargs)) as zf:
            yield zf
else:
    ZipFile = zipfile.ZipFile
