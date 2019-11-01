# @Author: richard
# @Date:   2018-12-06T13:27:47+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T17:24:43+00:00
import logging
import sys
# cd aws/vortexa_utils/
# import aws.s3.client_side_encryption.client as client
import vortexa_utils.aws.s3.client_side_encryption.client as client
import email.parser
from nose2.tools.such import helper


logger = logging.getLogger(__name__)

Bucket = 'ops-data.incoming-emails'
Key = 'incoming_email/4pnlhtml86pobumjn9d59mbkcq3to1i43sjbd201'


def test_get_obj():
    self = client.Client()
    location_info = self.s3.get_bucket_location(Bucket=Bucket)
    logger.info('location %s', location_info)

    obj = self.s3.get_object(Bucket=Bucket, Key=Key)
    handeler = client.DecryptHandeler(obj, self)
    envelop = handeler.envelope_v2(handeler.metadata)
    cipher = self.cipher_provider.decryptor(envelop)
    assert handeler.auth_tag()
    io = handeler.decrypt_auth(cipher)

    bytes = []
    while True:
        byte = io.read(1024)
        if byte == b'':
            break
        logger.info("Bytes Read %s/%s", io.bytes_read, io.content_length)
        logger.debug("Bytes    %s", byte)
        bytes.append(byte)
    io.verify()
    io.close()
    # logger.info('bytes %s', str(bytes))


def test_get_obj_io():
    cl = client.Client()
    with cl.get_object(Bucket, Key) as io:
        list(io)


def test_get_obj_mime():
    cl = client.Client()
    parser = email.parser.BytesParser()
    with cl.get_object(Bucket, Key) as io:
        parsed = parser.parse(io)

    target = parsed['to']
    source = parsed['from']
    helper.assertEqual(target, 'test@opsdata.vortexa.com')
    helper.assertEqual(source, 'Richard Mathie <richard.mathie@cantab.net>')

    logger.info('\twalking message')
    for part in parsed.walk():
        if part.get_content_type().startswith('text'):
            logger.info('\t%s', part.get_payload())
