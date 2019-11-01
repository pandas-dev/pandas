# @Author: richard
# @Date:   2018-11-28T17:01:20+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T16:57:10+00:00
from typing import Iterable
from botocore.response import StreamingBody
from .IODecrypter import IODecrypter

import logging

logger = logging.getLogger(__name__)


class IONocrypter(IODecrypter):

    def __init__(self, io):
        self.io: StreamingBody = io

    def read(self, chunk=None):
        return self.io.read(chunk)

    def iter_chunks(self, chunk_size: int = None) -> Iterable[bytes]:
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
        if chunk_size is None:
            chunk_size = self._DEFAULT_CHUNK_SIZE
        return self.io.iter_chunks(chunk_size)
