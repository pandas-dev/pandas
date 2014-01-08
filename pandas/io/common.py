"""Common IO api utilities"""

import sys
import zipfile
from contextlib import contextmanager, closing

from pandas.compat import StringIO
from pandas import compat


if compat.PY3:
    from urllib.request import urlopen, pathname2url
    _urlopen = urlopen
    from urllib.parse import urlparse as parse_url
    import urllib.parse as compat_parse
    from urllib.parse import (uses_relative, uses_netloc, uses_params,
                              urlencode, urljoin)
    from urllib.error import URLError
    from http.client import HTTPException
else:
    from urllib2 import urlopen as _urlopen
    from urllib import urlencode, pathname2url
    from urlparse import urlparse as parse_url
    from urlparse import uses_relative, uses_netloc, uses_params, urljoin
    from urllib2 import URLError
    from httplib import HTTPException
    from contextlib import contextmanager, closing
    from functools import wraps

    # @wraps(_urlopen)
    @contextmanager
    def urlopen(*args, **kwargs):
        with closing(_urlopen(*args, **kwargs)) as f:
            yield f


_VALID_URLS = set(uses_relative + uses_netloc + uses_params)
_VALID_URLS.discard('')


class PerformanceWarning(Warning):
    pass


class DtypeWarning(Warning):
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
        return parse_url(url).scheme in _VALID_URLS
    except:
        return False


def _is_s3_url(url):
    """Check for an s3 url"""
    try:
        return parse_url(url).scheme == 's3'
    except:
        return False


def maybe_read_encoded_stream(reader, encoding=None):
    """read an encoded stream from the reader and transform the bytes to
    unicode if required based on the encoding

        Parameters
        ----------
        reader : a streamable file-like object
        encoding : optional, the encoding to attempt to read

        Returns
        -------
        a tuple of (a stream of decoded bytes, the encoding which was used)

    """

    if compat.PY3 or encoding is not None:  # pragma: no cover
        if encoding:
            errors = 'strict'
        else:
            errors = 'replace'
            encoding = 'utf-8'
        reader = StringIO(reader.read().decode(encoding, errors))
    else:
        encoding = None
    return reader, encoding


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
        req = _urlopen(str(filepath_or_buffer))
        return maybe_read_encoded_stream(req, encoding)

    if _is_s3_url(filepath_or_buffer):
        try:
            import boto
        except:
            raise ImportError("boto is required to handle s3 files")
        # Assuming AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY
        # are environment variables
        parsed_url = parse_url(filepath_or_buffer)
        conn = boto.connect_s3()
        b = conn.get_bucket(parsed_url.netloc)
        k = boto.s3.key.Key(b)
        k.key = parsed_url.path
        filepath_or_buffer = StringIO(k.get_contents_as_string())
        return filepath_or_buffer, None

    return filepath_or_buffer, None


def file_path_to_url(path):
    """
    converts an absolute native path to a FILE URL.

    Parameters
    ----------
    path : a path in native format

    Returns
    -------
    a valid FILE URL
    """
    return urljoin('file:', pathname2url(path))


# ZipFile is not a context manager for <= 2.6
# must be tuple index here since 2.6 doesn't use namedtuple for version_info
if sys.version_info[1] <= 6:
    @contextmanager
    def ZipFile(*args, **kwargs):
        with closing(zipfile.ZipFile(*args, **kwargs)) as zf:
            yield zf
else:
    ZipFile = zipfile.ZipFile
