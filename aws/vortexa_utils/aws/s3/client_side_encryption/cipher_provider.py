# @Author: richard
# @Date:   2018-11-27T18:22:34+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-05T16:25:56+00:00
from Cryptodome.Cipher import AES  # pycryptodomex


class CipherProvider(object):
    def __init__(self, key):
        self.key = key

    def decryptor(self, envelope):
        pass

    def encryptor(self):
        cipher = AES.new(self.key, AES.MODE_GCM)
        return cipher
