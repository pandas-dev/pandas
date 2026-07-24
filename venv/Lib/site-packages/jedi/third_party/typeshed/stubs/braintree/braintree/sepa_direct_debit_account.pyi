from braintree.resource import Resource
from braintree.subscription import Subscription

class SepaDirectDebitAccount(Resource):
    @staticmethod
    def find(sepa_direct_debit_account_token): ...
    @staticmethod
    def delete(sepa_direct_debit_account_token): ...
    subscriptions: list[Subscription]
    def __init__(self, gateway, attributes) -> None: ...
