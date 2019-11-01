# @Author: richard
# @Date:   2018-11-28T18:11:28+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T13:06:33+00:00
from vortexa_utils.aws.s3.client_side_encryption import IOAuthDecrypter
from nose2.tools import params
from .test_IODecrypter import DummyChunksIO, IODecrypterTestCase


class IOAuthDecrypter(IODecrypterTestCase):
    io_decrypter_class = IOAuthDecrypter.IOAuthDecrypter

    def get_decrypter(self, cypher, io, content_length):
        return self.io_decrypter_class(cypher, io, content_length)

    def get_io(self, content_length):
        tag_length = 128
        return DummyChunksIO(content_length + tag_length)

    def invalid_decryption(self, content_length):
        with self.assertRaises(ValueError):
            super().invalid_decryption(content_length)
