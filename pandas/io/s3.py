""" s3 support for remote file interactivity """

import os
from pandas import compat
from pandas.compat import BytesIO

try:
    import boto
    from boto.s3 import key
except:
    raise ImportError("boto is required to handle s3 files")

if compat.PY3:
    from urllib.parse import urlparse as parse_url
else:
    from urlparse import urlparse as parse_url


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


def get_filepath_or_buffer(filepath_or_buffer, encoding=None,
                           compression=None):

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
