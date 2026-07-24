from _typeshed import Incomplete

from braintree.braintree_gateway import BraintreeGateway

class ClientToken:
    @staticmethod
    def generate(params: dict[str, Incomplete] | None = None, gateway: BraintreeGateway | None = None) -> str: ...
    @staticmethod
    def generate_signature() -> list[str | dict[str, list[str]]]: ...
