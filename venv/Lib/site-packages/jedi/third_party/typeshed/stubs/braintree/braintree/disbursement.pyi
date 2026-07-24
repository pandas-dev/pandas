from decimal import Decimal
from typing import Final

from braintree.merchant_account import MerchantAccount
from braintree.resource import Resource

class Disbursement(Resource):
    class Type:
        Credit: Final = "credit"
        Debit: Final = "debit"

    amount: Decimal
    merchant_account: MerchantAccount
    def __init__(self, gateway, attributes) -> None: ...
    def transactions(self): ...
    def is_credit(self) -> bool: ...
    def is_debit(self) -> bool: ...
