from _typeshed import Incomplete

from braintree.resource import Resource

class PaymentMethodCustomerDataUpdatedMetadata(Resource):
    payment_method: Incomplete
    enriched_customer_data: Incomplete
    def __init__(self, gateway, attributes) -> None: ...
