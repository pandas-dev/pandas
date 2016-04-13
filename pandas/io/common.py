"""Common IO api utilities"""

import sys
import os
import csv
import codecs
import zipfile
from contextlib import contextmanager, closing

from pandas.compat import StringIO, BytesIO, string_types, text_type
from pandas import compat
from pandas.formats.printing import pprint_thing
from pandas.core.common import is_number, AbstractMethodError


try:
    import pathlib
    _PATHLIB_INSTALLED = True
except ImportError:
    _PATHLIB_INSTALLED = False


try:
    from py.path import local as LocalPath
    _PY_PATH_INSTALLED = True
except:
    _PY_PATH_INSTALLED = False


if compat.PY3:
    from urllib.request import urlopen, pathname2url
    _urlopen = urlopen
    from urllib.parse import urlparse as parse_url
    from urllib.parse import (uses_relative, uses_netloc, uses_params,
                              urlencode, urljoin)
    from urllib.error import URLError
    from http.client import HTTPException  # noqa
else:
    from urllib2 import urlopen as _urlopen
    from urllib import urlencode, pathname2url  # noqa
    from urlparse import urlparse as parse_url
    from urlparse import uses_relative, uses_netloc, uses_params, urljoin
    from urllib2 import URLError  # noqa
    from httplib import HTTPException  # noqa
    from contextlib import contextmanager, closing  # noqa
    from functools import wraps  # noqa

    # @wraps(_urlopen)
    @contextmanager
    def urlopen(*args, **kwargs):
        with closing(_urlopen(*args, **kwargs)) as f:
            yield f


_VALID_URLS = set(uses_relative + uses_netloc + uses_params)
_VALID_URLS.discard('')


class CParserError(ValueError):
    """
    Exception that is thrown by the C engine when it encounters
    a parsing error in `pd.read_csv`
    """
    pass


class DtypeWarning(Warning):
    """
    Warning that is raised whenever `pd.read_csv` encounters non-
    uniform dtypes in a column(s) of a given CSV file
    """
    pass


class EmptyDataError(ValueError):
    """
    Exception that is thrown in `pd.read_csv` (by both the C and
    Python engines) when empty data or header is encountered
    """
    pass


class ParserWarning(Warning):
    """
    Warning that is raised in `pd.read_csv` whenever it is necessary
    to change parsers (generally from 'c' to 'python') contrary to the
    one specified by the user due to lack of support or functionality for
    parsing particular attributes of a CSV file with the requsted engine
    """
    pass


class BaseIterator(object):
    """Subclass this and provide a "__next__()" method to obtain an iterator.
    Useful only when the object being iterated is non-reusable (e.g. OK for a
    parser, not for an in-memory table, yes for its iterator)."""
    def __iter__(self):
        return self

    def __next__(self):
        raise AbstractMethodError(self)

if not compat.PY3:
    BaseIterator.next = lambda self: self.__next__()


try:
    from boto.s3 import key

    class BotoFileLikeReader(key.Key):
        """boto Key modified to be more file-like

        This modification of the boto Key will read through a supplied
        S3 key once, then stop. The unmodified boto Key object will repeatedly
        cycle through a file in S3: after reaching the end of the file,
        boto will close the file. Then the next call to `read` or `next` will
        re-open the file and start reading from the beginning.

        Also adds a `readline` function which will split the returned
        values by the `\n` character.
        """

        def __init__(self, *args, **kwargs):
            encoding = kwargs.pop("encoding", None)  # Python 2 compat
            super(BotoFileLikeReader, self).__init__(*args, **kwargs)
            # Add a flag to mark the end of the read.
            self.finished_read = False
            self.buffer = ""
            self.lines = []
            if encoding is None and compat.PY3:
                encoding = "utf-8"
            self.encoding = encoding
            self.lines = []

        def next(self):
            return self.readline()

        __next__ = next

        def read(self, *args, **kwargs):
            if self.finished_read:
                return b'' if compat.PY3 else ''
            return super(BotoFileLikeReader, self).read(*args, **kwargs)

        def close(self, *args, **kwargs):
            self.finished_read = True
            return super(BotoFileLikeReader, self).close(*args, **kwargs)

        def seekable(self):
            """Needed for reading by bz2"""
            return False

        def readline(self):
            """Split the contents of the Key by '\n' characters."""
            if self.lines:
                retval = self.lines[0]
                self.lines = self.lines[1:]
                return retval
            if self.finished_read:
                if self.buffer:
                    retval, self.buffer = self.buffer, ""
                    return retval
                else:
                    raise StopIteration

            if self.encoding:
                self.buffer = "{}{}".format(
                    self.buffer, self.read(8192).decode(self.encoding))
            else:
                self.buffer = "{}{}".format(self.buffer, self.read(8192))

            split_buffer = self.buffer.split("\n")
            self.lines.extend(split_buffer[:-1])
            self.buffer = split_buffer[-1]

            return self.readline()
except ImportError:
    # boto is only needed for reading from S3.
    pass
except TypeError:
    # boto/boto3 issues
    # GH11915
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
    """Check for an s3, s3n, or s3a url"""
    try:
        return parse_url(url).scheme in ['s3', 's3n', 's3a']
    except:
        return False


def maybe_read_encoded_stream(reader, encoding=None, compression=None):
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

        if compression == 'gzip':
            reader = BytesIO(reader.read())
        else:
            reader = StringIO(reader.read().decode(encoding, errors))
    else:
        if compression == 'gzip':
            reader = BytesIO(reader.read())
        encoding = None
    return reader, encoding


def _expand_user(filepath_or_buffer):
    """Return the argument with an initial component of ~ or ~user
       replaced by that user's home directory.

    Parameters
    ----------
    filepath_or_buffer : object to be converted if possible

    Returns
    -------
    expanded_filepath_or_buffer : an expanded filepath or the
                                  input if not expandable
    """
    if isinstance(filepath_or_buffer, string_types):
        return os.path.expanduser(filepath_or_buffer)
    return filepath_or_buffer


def _validate_header_arg(header):
    if isinstance(header, bool):
        raise TypeError("Passing a bool to header is invalid. "
                        "Use header=None for no header or "
                        "header=int or list-like of ints to specify "
                        "the row(s) making up the column names")


def _stringify_path(filepath_or_buffer):
    """Return the argument coerced to a string if it was a pathlib.Path
       or a py.path.local

    Parameters
    ----------
    filepath_or_buffer : object to be converted

    Returns
    -------
    str_filepath_or_buffer : a the string version of the input path
    """
    if _PATHLIB_INSTALLED and isinstance(filepath_or_buffer, pathlib.Path):
        return text_type(filepath_or_buffer)
    if _PY_PATH_INSTALLED and isinstance(filepath_or_buffer, LocalPath):
        return filepath_or_buffer.strpath
    return filepath_or_buffer


def get_filepath_or_buffer(filepath_or_buffer, encoding=None,
                           compression=None):
    """
    If the filepath_or_buffer is a url, translate and return the buffer
    passthru otherwise.

    Parameters
    ----------
    filepath_or_buffer : a url, filepath (str, py.path.local or pathlib.Path),
                         or buffer
    encoding : the encoding to use to decode py3 bytes, default is 'utf-8'

    Returns
    -------
    a filepath_or_buffer, the encoding, the compression
    """

    if _is_url(filepath_or_buffer):
        req = _urlopen(str(filepath_or_buffer))
        if compression == 'infer':
            content_encoding = req.headers.get('Content-Encoding', None)
            if content_encoding == 'gzip':
                compression = 'gzip'
            else:
                compression = None
        # cat on the compression to the tuple returned by the function
        to_return = (list(maybe_read_encoded_stream(req, encoding,
                                                    compression)) +
                     [compression])
        return tuple(to_return)

    if _is_s3_url(filepath_or_buffer):
        try:
            import boto
        except:
            raise ImportError("boto is required to handle s3 files")
        # Assuming AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY and AWS_S3_HOST
        # are environment variables
        parsed_url = parse_url(filepath_or_buffer)
        s3_host = os.environ.get('AWS_S3_HOST', 's3.amazonaws.com')

        try:
            conn = boto.connect_s3(host=s3_host)
        except boto.exception.NoAuthHandlerFound:
            conn = boto.connect_s3(host=s3_host, anon=True)

        b = conn.get_bucket(parsed_url.netloc, validate=False)
        if compat.PY2 and (compression == 'gzip' or
                           (compression == 'infer' and
                            filepath_or_buffer.endswith(".gz"))):
            k = boto.s3.key.Key(b, parsed_url.path)
            filepath_or_buffer = BytesIO(k.get_contents_as_string(
                encoding=encoding))
        else:
            k = BotoFileLikeReader(b, parsed_url.path, encoding=encoding)
            k.open('r')  # Expose read errors immediately
            filepath_or_buffer = k
        return filepath_or_buffer, None, compression

    # It is a pathlib.Path/py.path.local or string
    filepath_or_buffer = _stringify_path(filepath_or_buffer)
    return _expand_user(filepath_or_buffer), None, compression


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


def _get_handle(path, mode, encoding=None, compression=None):
    """Gets file handle for given path and mode.
    """
    if compression is not None:
        if encoding is not None and not compat.PY3:
            msg = 'encoding + compression not yet supported in Python 2'
            raise ValueError(msg)

        if compression == 'gzip':
            import gzip
            f = gzip.GzipFile(path, mode)
        elif compression == 'bz2':
            import bz2
            f = bz2.BZ2File(path, mode)
        elif compression == 'zip':
            import zipfile
            zip_file = zipfile.ZipFile(path)
            zip_names = zip_file.namelist()

            if len(zip_names) == 1:
                file_name = zip_names.pop()
                f = zip_file.open(file_name)
            elif len(zip_names) == 0:
                raise ValueError('Zero files found in ZIP file {}'
                                 .format(path))
            else:
                raise ValueError('Multiple files found in ZIP file.'
                                 ' Only one file per ZIP :{}'
                                 .format(zip_names))
        elif compression == 'xz':
            lzma = compat.import_lzma()
            f = lzma.LZMAFile(path, mode)
        else:
            raise ValueError('Unrecognized compression type: %s' %
                             compression)
        if compat.PY3:
            from io import TextIOWrapper
            f = TextIOWrapper(f, encoding=encoding)
        return f
    else:
        if compat.PY3:
            if encoding:
                f = open(path, mode, encoding=encoding)
            else:
                f = open(path, mode, errors='replace')
        else:
            f = open(path, mode)

    return f


class UTF8Recoder(BaseIterator):

    """
    Iterator that reads an encoded stream and reencodes the input to UTF-8
    """

    def __init__(self, f, encoding):
        self.reader = codecs.getreader(encoding)(f)

    def read(self, bytes=-1):
        return self.reader.read(bytes).encode("utf-8")

    def readline(self):
        return self.reader.readline().encode("utf-8")

    def next(self):
        return next(self.reader).encode("utf-8")


if compat.PY3:  # pragma: no cover
    def UnicodeReader(f, dialect=csv.excel, encoding="utf-8", **kwds):
        # ignore encoding
        return csv.reader(f, dialect=dialect, **kwds)

    def UnicodeWriter(f, dialect=csv.excel, encoding="utf-8", **kwds):
        return csv.writer(f, dialect=dialect, **kwds)
else:
    class UnicodeReader(BaseIterator):

        """
        A CSV reader which will iterate over lines in the CSV file "f",
        which is encoded in the given encoding.

        On Python 3, this is replaced (below) by csv.reader, which handles
        unicode.
        """

        def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
            f = UTF8Recoder(f, encoding)
            self.reader = csv.reader(f, dialect=dialect, **kwds)

        def __next__(self):
            row = next(self.reader)
            return [compat.text_type(s, "utf-8") for s in row]

    class UnicodeWriter:

        """
        A CSV writer which will write rows to CSV file "f",
        which is encoded in the given encoding.
        """

        def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
            # Redirect output to a queue
            self.queue = StringIO()
            self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
            self.stream = f
            self.encoder = codecs.getincrementalencoder(encoding)()
            self.quoting = kwds.get("quoting", None)

        def writerow(self, row):
            def _check_as_is(x):
                return (self.quoting == csv.QUOTE_NONNUMERIC and
                        is_number(x)) or isinstance(x, str)

            row = [x if _check_as_is(x)
                   else pprint_thing(x).encode("utf-8") for x in row]

            self.writer.writerow([s for s in row])
            # Fetch UTF-8 output from the queue ...
            data = self.queue.getvalue()
            data = data.decode("utf-8")
            # ... and reencode it into the target encoding
            data = self.encoder.encode(data)
            # write to the target stream
            self.stream.write(data)
            # empty queue
            self.queue.truncate(0)

        def writerows(self, rows):
            def _check_as_is(x):
                return (self.quoting == csv.QUOTE_NONNUMERIC and
                        is_number(x)) or isinstance(x, str)

            for i, row in enumerate(rows):
                rows[i] = [x if _check_as_is(x)
                           else pprint_thing(x).encode("utf-8") for x in row]

            self.writer.writerows([[s for s in row] for row in rows])
            # Fetch UTF-8 output from the queue ...
            data = self.queue.getvalue()
            data = data.decode("utf-8")
            # ... and reencode it into the target encoding
            data = self.encoder.encode(data)
            # write to the target stream
            self.stream.write(data)
            # empty queue
            self.queue.truncate(0)
