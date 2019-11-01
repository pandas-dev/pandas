# @Author: richard
# @Date:   2018-11-28T17:01:20+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T16:57:10+00:00
# from typing import Iterable

from io import IOBase
from botocore.response import StreamingBody

import logging

logger = logging.getLogger(__name__)


class IODecrypter(IOBase):
    _DEFAULT_CHUNK_SIZE = 1024

    def __init__(self, cipher, io: StreamingBody):
        self.cipher: object = cipher
        self.io: StreamingBody = io

    def read(self, chunk=None):
        bytes = self.io.read(chunk)
        return self.cipher.decrypt(bytes)

    def __iter__(self):
        """Return an iterator to yield 1k chunks from the raw stream."""
        return self.iter_chunks(self._DEFAULT_CHUNK_SIZE)

    def iter_chunks(self, chunk_size: int = _DEFAULT_CHUNK_SIZE):
        # type: (...) -> Iterable[bytes]
        """Return an iterator to yield chunks bytes from the raw `io` stream.

        Parameters
        ----------
        chunk_size : int
            iterates over no more than Chunk size bytes. If `None` use
            `self._DEFAULT_CHUNK_SIZE`.

        Returns
        -------
        Iterator[bytes]

        """
        decrypt = self.cipher.decrypt
        chunks = self.io.iter_chunks(chunk_size)

        return (decrypt(bytes) for bytes in chunks)

    def close(self):
        """Close the underlying http response stream."""
        self.io.close()

    def readable(self):
        return True

    def seekable(self):
        return False

    def writable(self):
        return False
