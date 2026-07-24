from braintree.merchant_account import MerchantAccount
from braintree.resource import Resource

class Merchant(Resource):
    merchant_accounts: list[MerchantAccount]
    def __init__(self, gateway, attributes) -> None: ...
