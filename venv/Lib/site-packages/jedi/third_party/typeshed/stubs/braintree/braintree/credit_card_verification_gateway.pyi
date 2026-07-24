from _typeshed import Incomplete

from braintree.credit_card_verification import CreditCardVerification
from braintree.error_result import ErrorResult
from braintree.resource_collection import ResourceCollection
from braintree.successful_result import SuccessfulResult

class CreditCardVerificationGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def find(self, verification_id: str) -> CreditCardVerification: ...
    def search(self, *query) -> ResourceCollection: ...
    def create(self, params: dict[str, Incomplete] | None) -> SuccessfulResult | ErrorResult | None: ...
