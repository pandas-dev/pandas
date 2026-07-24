from _typeshed import Incomplete

from braintree.us_bank_account import UsBankAccount

class UsBankAccountGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def find(self, us_bank_account_token: str) -> UsBankAccount | None: ...
