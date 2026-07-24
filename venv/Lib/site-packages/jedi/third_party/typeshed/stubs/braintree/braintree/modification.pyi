from decimal import Decimal

from braintree.resource import Resource

class Modification(Resource):
    amount: Decimal
    def __init__(self, gateway, attributes) -> None: ...
