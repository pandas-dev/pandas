from decimal import Decimal

from braintree.attribute_getter import AttributeGetter

class MonetaryAmount(AttributeGetter):
    value: Decimal
    def __init__(self, attributes) -> None: ...
