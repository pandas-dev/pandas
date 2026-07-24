from _typeshed import Incomplete

from braintree.braintree_gateway import BraintreeGateway

class ClientTokenGateway:
    gateway: BraintreeGateway
    config: Incomplete
    def __init__(self, gateway: BraintreeGateway) -> None: ...
    def generate(self, params: dict[str, Incomplete] | None = None) -> str: ...
