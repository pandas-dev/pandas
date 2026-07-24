from _typeshed import Incomplete

from braintree.error_result import ErrorResult
from braintree.payment_method_nonce import PaymentMethodNonce
from braintree.successful_result import SuccessfulResult

class PaymentMethodNonceGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def create(self, payment_method_token: str, params=...) -> SuccessfulResult | ErrorResult: ...
    def find(self, payment_method_nonce: str) -> PaymentMethodNonce: ...
