from _typeshed import Incomplete

from braintree.resource import Resource

class EnrichedCustomerData(Resource):
    profile_data: Incomplete
    def __init__(self, gateway, attributes) -> None: ...
