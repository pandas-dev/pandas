from vortexa_utils.aws.utils.dataclasses import nested_dataclass
from . import Mail, Receipt


@nested_dataclass
class Notification:
    """Notification Object.

    Attributes
    ----------
    notificationType: str
        The notification type. For this type of notification,
        the value is always Received.
    receipt : Recipt
        Object that contains information about the email delivery.
    mail : Mail
        Object that contains information about the email
        associated with the notification.
    content : str
        String that contains the raw, unmodified email, which is typically
        in Multipurpose Internet Mail Extensions (MIME) format.
        *Only if the notification was triggered by an SNS action.*

    """

    notificationType: str
    receipt: Receipt
    mail: Mail
    content: str
