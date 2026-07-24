from _typeshed import Incomplete

from braintree.transaction_line_item import TransactionLineItem

class TransactionLineItemGateway:
    gateway: Incomplete
    config: Incomplete
    def __init__(self, gateway) -> None: ...
    def find_all(self, transaction_id: str) -> list[TransactionLineItem]: ...
