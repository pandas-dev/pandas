from _typeshed import Incomplete

from braintree.resource import Resource

class LocalPaymentFunded(Resource):
    transaction: Incomplete
    def __init__(self, gateway, attributes) -> None: ...
