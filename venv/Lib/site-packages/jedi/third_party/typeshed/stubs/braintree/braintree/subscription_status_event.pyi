from decimal import Decimal

from braintree.resource import Resource

class SubscriptionStatusEvent(Resource):
    balance: Decimal
    price: Decimal
    def __init__(self, gateway, attributes) -> None: ...
