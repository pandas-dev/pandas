# @Author: richard
# @Date:   2018-11-28T18:11:28+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T13:07:14+00:00
from io import IOBase

from vortexa_utils.aws.s3.client_side_encryption.IODecrypter import IODecrypter
import unittest
from nose2.tools import params


class DummyCipher(object):
    def __init__(self, valid: bool = True):
        self.valid = valid

    def decrypt(self, bytes):
        return bytes

    def verify(self, tag):
        if not self.valid:
            raise ValueError("MAC check failed")
        pass


class DummyChunksIO(IOBase):
    _DEFAULT_CHUNK_SIZE = 1024

    def __init__(self, size):
        self.bytes_read = 0
        self.size = size

    def read(self, chunk=-1):
        if chunk < 0:
            chunk = self.size - self.bytes_read
        else:
            chunk = min(chunk, abs(self.size - self.bytes_read))
        self.bytes_read += chunk
        return b' ' * chunk

    def __iter__(self):
        """Return an iterator to yield 1k chunks from the raw stream.
        """
        return self.iter_chunks(self._DEFAULT_CHUNK_SIZE)

    def iter_chunks(self, chunk_size=_DEFAULT_CHUNK_SIZE):
        """Return an iterator to yield chunks of chunk_size bytes from the raw
        stream.
        """
        while True:
            bytes = self.read(chunk_size)
            if bytes == b'':
                break
            yield bytes

    def close(self):
        pass

    def readable(self):
        return True

    def seekable(self):
        return False

    def writable(self):
        return False


class IODecrypterTestCase(unittest.TestCase):
    io_decrypter_class = IODecrypter

    def get_decrypter(self, cypher, io, content_length):
        return self.io_decrypter_class(cypher, io)

    def get_io(self, content_length):
        return DummyChunksIO(content_length)

    def make_decrypter(self, content_length, valid=True):
        io = DummyChunksIO(content_length)
        cypher = DummyCipher(valid=valid)
        return self.get_decrypter(cypher, io, content_length)

    @params(123, 1024, 1024*3, 1024*3+123, 1, 0)
    def test_read(self, content_length):
        with self.make_decrypter(content_length) as decrypter:
            bytes = list(decrypter)
        self.assertEqual(b''.join(bytes),  b' ' * content_length)

    @params(123, 1024, 1024*3, 1024*3+123, 1, 0)
    def test_invalid(self, content_length):
        self.invalid_decryption(content_length)

    def invalid_decryption(self, content_length):
        with self.make_decrypter(content_length, valid=False) as decrypter:
            list(decrypter)
