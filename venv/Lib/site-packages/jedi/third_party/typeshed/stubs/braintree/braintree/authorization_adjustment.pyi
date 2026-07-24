from decimal import Decimal

from braintree.attribute_getter import AttributeGetter

class AuthorizationAdjustment(AttributeGetter):
    amount: Decimal | None
    def __init__(self, attributes) -> None: ...
