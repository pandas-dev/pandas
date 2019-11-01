# @Author: richard
# @Date:   2018-11-27T18:20:28+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-05T17:20:59+00:00
import base64
import boto3
import json

from Cryptodome.Cipher import AES  # pycryptodomex
from .cipher_provider import CipherProvider


class KMSCipherProvider(CipherProvider):
    aes_mode_map = {
        'AES/GCM/NoPadding': AES.MODE_GCM,
        'AES/CBC/PKCS5Padding': AES.MODE_CBC,
        'AES/CBC/PKCS7Padding': AES.MODE_CBC
    }

    def __init__(self, key_id=None, **kwargs):
        self.kms = boto3.client('kms', **kwargs)
        self.key_id = key_id

    def decryptor(self, envelope):
        key_alg = envelope['x-amz-cek-alg']
        aes_mode = self.aes_mode_map.get(key_alg)
        if aes_mode is None:
            raise Exception(f'unknown encryption algorythem {key_alg}')

        envelope_key = base64.b64decode(envelope['x-amz-key-v2'])
        iv = base64.b64decode(envelope['x-amz-iv'])
        encryption_context = json.loads(envelope['x-amz-matdesc'])

        decrypted_envelope = self.kms.decrypt(
            CiphertextBlob=envelope_key,
            EncryptionContext=encryption_context
        )
        key = decrypted_envelope['Plaintext']
        cipher = AES.new(key, aes_mode, iv)
        return cipher

    def encryptor(self):
        encryption_context = {"kms_cmk_id": self.key_id}

        key_data = self.kms.generate_data_key(
            KeyId=self.key_id,
            EncryptionContext=encryption_context,
            KeySpec='AES_256'
        )

        key = key_data['Plaintext']
        cipher = AES.new(key, AES.MODE_GCM)

        envelope = {
            'x-amz-key-v2': base64.encodebytes(key_data['CiphertextBlob']),
            'x-amz-iv': base64.encodebytes(cipher.nonce),
            'x-amz-cek-alg': 'AES/GCM/NoPadding',
            'x-amz-wrap-alg': 'kms',
            'x-amz-matdesc': json.dumps(encryption_context)
        }
        return envelope, cipher
