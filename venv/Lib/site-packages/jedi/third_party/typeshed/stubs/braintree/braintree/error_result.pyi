from _typeshed import Incomplete
from typing import Literal

from braintree.credit_card_verification import CreditCardVerification
from braintree.errors import Errors
from braintree.plan import Plan
from braintree.subscription import Subscription
from braintree.transaction import Transaction

class ErrorResult:
    params: Incomplete
    errors: Errors
    message: Incomplete
    credit_card_verification: CreditCardVerification | None
    transaction: Transaction
    subscription: Subscription
    merchant_account: Plan
    def __init__(self, gateway, attributes: dict[str, Incomplete]) -> None: ...
    @property
    def is_success(self) -> Literal[False]: ...
