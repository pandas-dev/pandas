from braintree.braintree_gateway import BraintreeGateway
from braintree.configuration import Configuration
from braintree.error_result import ErrorResult
from braintree.successful_result import SuccessfulResult

class PayPalPaymentResourceGateway:
    config: Configuration
    gateway: BraintreeGateway
    def __init__(self, gateway: BraintreeGateway) -> None: ...
    def update(self, params) -> SuccessfulResult | ErrorResult: ...
