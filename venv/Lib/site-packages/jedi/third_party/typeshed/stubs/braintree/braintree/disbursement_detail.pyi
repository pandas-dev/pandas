from decimal import Decimal

from braintree.attribute_getter import AttributeGetter

class DisbursementDetail(AttributeGetter):
    settlement_amount: Decimal | None
    settlement_currency_exchange_rate: Decimal | None
    def __init__(self, attributes) -> None: ...
    @property
    def is_valid(self) -> bool: ...
