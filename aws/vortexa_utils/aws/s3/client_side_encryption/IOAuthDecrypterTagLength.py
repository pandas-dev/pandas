# @Author: richard
# @Date:   2018-11-28T17:01:36+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T16:00:31+00:00
import logging
from .IODecrypter import IODecrypter
from io import BytesIO, IOBase
logger = logging.getLogger(__name__)


class StreamChunker(IOBase):
    """StreamChunker a class to keep the last tag bites of a file

    keeps hold of the last `tag_length` bytes in `self.tag`
    when reading from a `BytesIO` object.
    """

    def __init__(self, io: BytesIO, tag_length: int):
        self.io = io
        self.tag_length = tag_length
        # get the first chunk if this is the first read
        self.tag = self.io.read(self.tag_length)

    def read(self, chunk=None):
        bytes = self.tag + self.io.read(chunk)
        bytes, self.tag = bytes[:-self.tag_length], bytes[-self.tag_length:]
        return bytes

    def close(self):
        """Close the underlying http response stream."""
        self.io.close()

    def readable(self):
        return True

    def seekable(self):
        return False

    def writable(self):
        return False


class IOAuthDecrypterTagLength(IODecrypter):
    def __init__(self, cipher, io, tag_length, chunk_size=16*1024):
        super().__init__(cipher, StreamChunker(io, tag_length))

    def verify(self):
        # the remaining bytes should be the auth tag
        tag = self.io.tag
        logger.debug("Verifing Tag %s", tag)
        self.cipher.verify(tag)

    def iter_chunks(self, chunk_size=None):
        """Return an iterator to yield chunks of chunk_size bytes from the raw
        stream.
        """
        if chunk_size is None:
            chunk_size = self._DEFAULT_CHUNK_SIZE

        while True:
            bytes = self.read(chunk_size)
            if bytes == b'':
                break
            yield bytes
        self.verify()
