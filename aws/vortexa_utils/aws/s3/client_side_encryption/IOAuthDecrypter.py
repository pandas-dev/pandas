# @Author: richard
# @Date:   2018-11-28T17:01:36+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T16:00:31+00:00
import logging
from .IODecrypter import IODecrypter

logger = logging.getLogger(__name__)


class IOAuthDecrypter(IODecrypter):
    def __init__(self, cipher, io, content_length, chunk_size=16*1024):
        super().__init__(cipher, io)
        self.bytes_read = 0
        self.content_length = content_length

    def read(self, chunk=None):
        chunk = min(chunk, self.content_length - self.bytes_read)
        bytes = super().read(chunk)
        logger.debug("Bytes Read %s/%s", self.bytes_read, self.content_length)
        self.bytes_read += len(bytes)
        return bytes

    def verify(self):
        # the remaining bytes should be the auth tag
        tag = self.io.read()
        logger.debug("Verifing Tag %s", tag)
        self.cipher.verify(tag)

    def iter_chunks(self, chunk_size=None):
        """Return an iterator to yield chunks of chunk_size bytes from the raw
        stream.
        """
        if chunk_size is None:
            chunk_size = self._DEFAULT_CHUNK_SIZE

        while self.bytes_read < self.content_length:
            bytes = self.read(chunk_size)
            yield bytes
        self.verify()
