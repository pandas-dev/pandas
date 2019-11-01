import unittest
import io
from vortexa_utils.aws.s3.client_side_encryption.IOAuthDecrypterTagLength \
        import StreamChunker
from nose2.tools import params


class StreamChunkerTestCase(unittest.TestCase):

    def get_chunker(self, io, tag_length):
        return StreamChunker(io, tag_length)

    def test_tagged(self):
        fixture = io.BytesIO(b'1234567890')
        chunker = StreamChunker(fixture, 3)
        bytes = chunker.read()
        self.assertEqual(chunker.tag, b'890')
        self.assertEqual(bytes, b'1234567')

    @params(*range(1, 11))
    def test_read_in_chunks(self, chunk):
        bytes = b'1234567890'
        fixture = io.BytesIO(bytes)
        tag_length = 3
        chunker = StreamChunker(fixture, tag_length)
        result = []
        index = 0
        while True:
            byte = chunker.read(chunk)
            if byte == b'':
                break
            result.append(byte)
            self.assertEqual(bytes[index:index + len(byte)], byte)
            index += len(byte)
        print(result)
        self.assertEqual(bytes[-tag_length:], chunker.tag)
        self.assertEqual(b''.join(result), bytes[:-tag_length])
        # check that subsuquent reads return nothing and tag is correct
        for i in range(10):
            byte = chunker.read(chunk)
            self.assertEqual(b'', byte)
            self.assertEqual(bytes[-tag_length:], chunker.tag)
