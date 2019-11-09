import base64
import boto3
import json
import mimetypes
import sendgrid
from sendgrid.helpers.mail import *
from typing import List


secretsmanager = boto3.client('secretsmanager')


def create_sendgrid_client():
    secret = secretsmanager.get_secret_value(SecretId='prod/sendgrid')
    api_key = json.loads(secret['SecretString'])['SENDGRID_API_KEY']

    return sendgrid.SendGridAPIClient(apikey=api_key)


def build_attachment(buf: bytes, filename: str, disposition: str = "attachment", content_id: str = None):
    encoded = base64.b64encode(buf).decode()

    mime_type, encoding = mimetypes.guess_type(filename)

    attachment = Attachment()
    attachment.content = encoded
    attachment.type = mime_type
    attachment.filename = filename
    attachment.disposition = disposition
    attachment.content_id = content_id

    return attachment


def add_recipients(recipients: List[str], mail: Mail):
    personalization = Personalization()

    for rec in recipients:
        personalization.add_to(Email(rec))

    mail.add_personalization(personalization)

    return mail
