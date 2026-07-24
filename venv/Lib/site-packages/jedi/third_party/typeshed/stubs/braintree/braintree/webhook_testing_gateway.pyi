from _typeshed import Incomplete

class WebhookTestingGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def sample_notification(self, kind: str, id: str, source_merchant_id: str | None = None) -> dict[str, str | bytes]: ...
