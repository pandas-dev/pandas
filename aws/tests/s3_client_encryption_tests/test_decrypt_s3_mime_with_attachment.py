# @Author: richard
# @Date:   2018-12-06T17:26:08+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-06T19:36:16+00:00
# cd aws/vortexa_utils/
# import aws.s3.client_side_encryption.client as client
import logging
import vortexa_utils.aws.s3.client_side_encryption.client as client
import io
import email.parser
from email import policy
from email.iterators import _structure
import base64
from nose2.tools.such import helper

import pandas as pd

logger = logging.getLogger(__name__)

Bucket = 'ops-data.incoming-emails'
Key = 'incoming_email/akrk0l8sq4lm7qkgj8hpurfshpnj8frgqpqe9mg1'
Key = 'incoming_email/8ej2ldqnsmako2tgsbdpqg8tdi6tdnduoscojdo1'


def test_get_attachment():
    cl = client.Client()
    parser = email.parser.BytesParser(policy=policy.default)
    with cl.get_object(Bucket, Key) as io:
        parsed = parser.parse(io)
    _structure(parsed)

    # with open("/home/richard/an_email", 'wb') as f:
    #     for b in io:
    #         f.write(b)
    #
    # atts = list(parsed.iter_attachments())
    # [a.get_filename() for a in atts]
    # [a.get_content_type() for a in atts]
    # att = atts[2]
    # att
    # att.get_content_type()
    # pd.read_excel(io.BytesIO(att.get_content()))

    target = parsed['to']
    source = parsed['from']
    helper.assertEqual(target, 'test@opsdata.vortexa.com')
    helper.assertEqual(source, 'Richard Mathie <richard.mathie@vortexa.com>')

    parsed['subject']

    for part in parsed.walk():
        print(part.get_content_type())
    att = parsed.get_payload()
    att[0].get_content_type()
    att[0].get_payload()[1].get_payload()

    logger.debug('\nwalking message')
    for part in parsed.walk():
        content_type = part.get_content_type()
        if content_type.startswith('text'):
            logger.debug(content_type)
            payload = part.get_payload()
            if content_type == 'text/csv':
                csv = base64.decodebytes(payload.encode('utf-8'))
                for line in csv.splitlines():
                    logger.debug(line)
            else:
                logger.debug('\n%s', payload)
