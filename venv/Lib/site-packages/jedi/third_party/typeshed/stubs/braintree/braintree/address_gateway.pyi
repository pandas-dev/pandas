from _typeshed import Incomplete

from braintree.address import Address
from braintree.braintree_gateway import BraintreeGateway
from braintree.error_result import ErrorResult
from braintree.successful_result import SuccessfulResult

class AddressGateway:
    gateway: BraintreeGateway
    config: Incomplete
    def __init__(self, gateway: BraintreeGateway) -> None: ...
    def create(self, params: dict[str, Incomplete] | None = None) -> SuccessfulResult | ErrorResult | None: ...
    def delete(self, customer_id: str, address_id: str) -> SuccessfulResult: ...
    def find(self, customer_id: str, address_id: str) -> Address: ...
    def update(
        self, customer_id: str, address_id: str, params: dict[str, Incomplete] | None = None
    ) -> SuccessfulResult | ErrorResult | None: ...
