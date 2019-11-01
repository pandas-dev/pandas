# @Author: richard
# @Date:   2018-11-27T17:24:50+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T16:38:19+00:00
import boto3
import logging

from .IODecrypter import IODecrypter
from .IONocrypter import IONocrypter
from .IOAuthDecrypter import IOAuthDecrypter
from .IOAuthDecrypterTagLength import IOAuthDecrypterTagLength

logger = logging.getLogger(__name__)
kms = boto3.client('kms')


class DecryptionException(Exception):
    pass


class DecryptHandeler(object):

    V1_ENVELOPE_KEYS = {
        'x-amz-key',
        'x-amz-iv',
        'x-amz-matdesc'
    }

    V2_ENVELOPE_KEYS = {
      'x-amz-key-v2',
      'x-amz-iv',
      'x-amz-cek-alg',
      'x-amz-wrap-alg',
      'x-amz-matdesc'
    }

    POSSIBLE_ENVELOPE_KEYS = V1_ENVELOPE_KEYS | V2_ENVELOPE_KEYS

    POSSIBLE_ENCRYPTION_FORMATS = {
      'AES/GCM/NoPadding',
      'AES/CBC/PKCS5Padding',
      'AES/CBC/PKCS7Padding'
    }

    def __init__(self, obj, context, strict=False):
        self.obj = obj
        self.context = context
        self.metadata = obj['Metadata']
        self.body = obj['Body']
        self.strict = strict

    def decrypt(self):
        cipher = self.decryption_cipher()
        logger.debug(self.metadata)
        if cipher:
            logger.debug(cipher)
            if self.auth_tag():
                return self.decrypt_auth(cipher)
            return IODecrypter(cipher=cipher, io=self.body)
        # Object not encrypted with an envelope
        mesg = f"Unencrypted Object at {self.obj['ETag']}"
        if self.strict:
            logger.error(mesg)
            raise ValueError(mesg)
        else:
            logger.warning(mesg)
            return IONocrypter(io=self.body)

    def auth_tag(self):
        return 'x-amz-tag-len' in self.metadata

    def decryption_cipher(self):
        envelope = self.extract_envelop(self.metadata)
        if envelope:
            return self.context.cipher_provider.decryptor(envelope)

    def extract_envelop(self, meta):
        if 'x-amz-key' in meta:
            return self.envelope_v1(meta)
        elif 'x-amz-key-v2' in meta:
            return self.envelope_v2(meta)

        key_prefix = 'x-amz-key'
        key = next((k for k in meta.keys() if k.startswith(key_prefix)), None)
        if key is not None:
            key_version = key[len(key_prefix):]
            mesg = f'Unknown envelope encryption version {key_version}'
            raise DecryptionException(mesg)
        # no envelope found
        return None

    def envelope_v2(self, meta):
        if meta['x-amz-cek-alg'] not in self.POSSIBLE_ENCRYPTION_FORMATS:
            alg = meta['x-amz-cek-alg']
            msg = f'unsuported content encrypting key format: {alg}'
            raise DecryptionException(msg)
        if meta['x-amz-wrap-alg'] != 'kms':
            alg = meta['x-amz-wrap-alg']
            msg = f'unsupported key wrapping algorithm: {alg}'
            raise DecryptionException(msg)
        if not self.V2_ENVELOPE_KEYS <= set(meta.keys()):
            msg = "incomplete v2 encryption envelope:\n"
            msg += f"  expected: #{', '.join(self.V2_ENVELOPE_KEYS)}\n"
            msg += f"  got: #{', '.join(meta.keys)}"
        return meta

    def envelope_v1(self, meta):
        return meta

    def decrypt_auth(self, cipher):
        meta = self.metadata

        content_length_string = meta.get(
            'x-amz-unencrypted-content-length',
            None
        )
        if content_length_string is not None:
            content_length = int(content_length_string)
            return IOAuthDecrypter(cipher, self.body, content_length)
        tag_length = int(meta['x-amz-tag-len'])//8
        return IOAuthDecrypterTagLength(cipher, self.body, tag_length)
