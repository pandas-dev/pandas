from braintree.resource import Resource
from braintree.subscription import Subscription

class AmexExpressCheckoutCard(Resource):
    subscriptions: list[Subscription]
    def __init__(self, gateway, attributes) -> None: ...
    @property
    def expiration_date(self): ...
