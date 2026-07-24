from braintree.address import Address
from braintree.resource import Resource
from braintree.subscription import Subscription

class MasterpassCard(Resource):
    billing_address: Address | None
    subscriptions: list[Subscription]
    def __init__(self, gateway, attributes) -> None: ...
    @property
    def expiration_date(self) -> str: ...
    @property
    def masked_number(self) -> str: ...
