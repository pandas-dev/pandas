from _typeshed import Incomplete

from braintree.discount import Discount

class DiscountGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def all(self) -> list[Discount]: ...
