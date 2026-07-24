from braintree.resource import Resource

class LocalPaymentExpired(Resource):
    def __init__(self, gateway, attributes) -> None: ...
