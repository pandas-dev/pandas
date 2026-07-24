from _typeshed import Incomplete

from braintree.error_result import ErrorResult
from braintree.resource import Resource
from braintree.successful_result import SuccessfulResult

class PayPalPaymentResource(Resource):
    def __init__(self, gateway, attributes) -> None: ...
    @staticmethod
    def update(request) -> SuccessfulResult | ErrorResult: ...
    @staticmethod
    def update_signature() -> list[Incomplete]: ...
