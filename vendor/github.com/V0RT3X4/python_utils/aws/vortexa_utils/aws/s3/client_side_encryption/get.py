# @Author: richard
# @Date:   2018-11-27T14:58:39+00:00
# @Last modified by:   richard
# @Last modified time: 2018-11-30T12:09:27+00:00

# see https://medium.com/@samnco/reading-aws-ses-encrypted-emails-with-boto3-9c177f8ba130
# and https://github.com/boto/boto3/issues/38

import base64
import json
from Cryptodome.Cipher import AES  # pycryptodomex
import boto3


s3 = boto3.client('s3')
kms = boto3.client('kms')


def chunker(length, chunk_size):
    index = 0
    while index < length:
        chunk = min(chunk_size, length - index)
        index += chunk
        yield chunk

list(chunker(2, 3))


def content_streamer(bytes_io, content_length, chunk_size=16*1024):
    for chunk in chunker(content_length, chunk_size):
        yield bytes_io.read(chunk)





def decrypt_object(obj):
    metadata = obj['Metadata']
    key_alg = metadata['x-amz-cek-alg']

    envelope_key = base64.b64decode(metadata['x-amz-key-v2'])
    envelope_iv = base64.b64decode(metadata['x-amz-iv'])
    encrypt_ctx = json.loads(metadata['x-amz-matdesc'])

    # x-amz-tag-len in is in bits so /8 to get bytes
    tag_len = int(metadata['x-amz-tag-len'])/8
    original_size = int(metadata['x-amz-unencrypted-content-length'])

    decrypted_envelope_key = kms.decrypt(
        CiphertextBlob=envelope_key,
        EncryptionContext=encrypt_ctx
    )
    key = decrypted_envelope_key['Plaintext']

    if key_alg == 'AES/GCM/NoPadding':
        # x-amz-tag-len in is in bits so /8 to get bytes
        cipher = AES.new(key, AES.MODE_GCM, envelope_iv)
    elif key_alg == 'AES/CBC/PKCS5Padding':
        cipher = AES.new(key, AES.MODE_CBC, envelope_iv)
    else:
        raise Exception('unknown encryption algorythem')

    body = obj['Body']

    body = body.read()
    body, tag = body[:original_size], body[original_size:]
    email = cipher.decrypt(body)
    cipher.verify(tag)
    return email


def get_object(bucket, key):
    obj = s3.get_object(Bucket=bucket_name, Key=key)
    location_info = s3.get_bucket_location(Bucket=bucket_name)
    bucket_region = location_info['LocationConstraint']
