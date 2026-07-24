from typing_extensions import deprecated

from braintree.address import Address
from braintree.credit_card_verification import CreditCardVerification
from braintree.resource import Resource
from braintree.subscription import Subscription

@deprecated("Visa Checkout is no longer supported for creating new transactions.")
class VisaCheckoutCard(Resource):
    billing_address: Address | None
    subscriptions: list[Subscription]
    verification: CreditCardVerification
    def __init__(self, gateway, attributes): ...
    @property
    def expiration_date(self) -> str: ...
    @property
    def masked_number(self) -> str: ...
