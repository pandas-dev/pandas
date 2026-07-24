from braintree.resource import Resource
from braintree.transaction import Transaction

class LocalPaymentCompleted(Resource):
    transaction: Transaction
    def __init__(self, gateway, attributes) -> None: ...
