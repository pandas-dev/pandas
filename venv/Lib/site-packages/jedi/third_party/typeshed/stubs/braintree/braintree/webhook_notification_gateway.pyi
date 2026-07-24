from _typeshed import Incomplete

from braintree.webhook_notification import WebhookNotification

text_type = str

class WebhookNotificationGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def parse(self, signature: str, payload: str) -> WebhookNotification: ...
    def verify(self, challenge: str) -> str: ...
