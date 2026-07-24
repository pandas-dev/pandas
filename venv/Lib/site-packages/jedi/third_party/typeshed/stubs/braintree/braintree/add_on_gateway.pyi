from _typeshed import Incomplete

from braintree.add_on import AddOn

class AddOnGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def all(self) -> list[AddOn]: ...
