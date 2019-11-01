from vortexa_utils.aws.ses.notification.types import Record
from json import loads
from .test_mail import mail_json
from .test_receipt import receipt_json


ses = dict(
    receipt=receipt_json,
    mail=mail_json
)


record_json = loads("""
{
    "eventSource": "aws:ses",
    "eventVersion": "1.0",
    "ses": {
       "receipt": {
        },
       "mail": {
       }
     }
}
""")

record_json.update(ses=ses)


def test_record():
    record = Record(**record_json)
    record.ses
    assert record.eventSource == "aws:ses"
