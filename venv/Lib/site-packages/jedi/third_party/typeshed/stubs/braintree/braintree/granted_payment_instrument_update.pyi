from _typeshed import Incomplete

from braintree.resource import Resource

class GrantedPaymentInstrumentUpdate(Resource):
    payment_method_nonce: Incomplete
    def __init__(self, gateway, attributes) -> None: ...
