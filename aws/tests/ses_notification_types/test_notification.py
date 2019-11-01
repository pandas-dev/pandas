from vortexa_utils.aws.ses.notification.types import Notification
from json import loads
from .test_mail import mail_json
from .test_action import action_json_sns
from .test_receipt import receipt_json


nodification_json = loads("""
{
"notificationType": "Received",
"content": "blarblarblar"
}
"""
)

nodification_json.update(
    mail=mail_json,
    receipt=receipt_json
)


def test_init():
    Notification(**nodification_json)
