# @Author: richard
# @Date:   2018-12-05T16:23:13+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-05T19:43:28+00:00
import unittest
from vortexa_utils.aws.s3.client_side_encryption import kms_cipher_provider
import logging


logger = logging.getLogger(__name__)


def log_bytes(*bytes):
    logger.info(f' bytes: {bytes}')


class KMSCipherProviderTest(unittest.TestCase):
    test_key_id = 'alias/python_utils_test_key'

    def get_cipher(self):
        return kms_cipher_provider.KMSCipherProvider(self.test_key_id)

    def test_encrypt(self):
        envelope, cipher = self.get_cipher().encryptor()
        plaintext = b"The quick brown fox jumped over the lazy dog"
        self.plaintext = plaintext
        ciphertext, tag = cipher.encrypt_and_digest(plaintext)
        log_bytes(ciphertext, tag)
        self.assertNotEqual(ciphertext, plaintext)
        package = (envelope, ciphertext, tag)
        return package

    def test_decrypt(self):
        envelope, ciphertext, tag = self.test_encrypt()
        cipher = kms_cipher_provider.KMSCipherProvider().decryptor(envelope)
        plaintext = cipher.decrypt(ciphertext)
        log_bytes(ciphertext, tag, plaintext)
        self.assertEqual(plaintext, self.plaintext)
        cipher.verify(tag)
