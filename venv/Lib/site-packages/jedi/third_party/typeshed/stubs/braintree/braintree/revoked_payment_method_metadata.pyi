from _typeshed import Incomplete

from braintree.resource import Resource

class RevokedPaymentMethodMetadata(Resource):
    revoked_payment_method: Incomplete
    customer_id: Incomplete
    token: Incomplete
    def __init__(self, gateway, attributes) -> None: ...
