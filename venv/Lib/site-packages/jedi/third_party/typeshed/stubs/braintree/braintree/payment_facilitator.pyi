from braintree.attribute_getter import AttributeGetter
from braintree.sub_merchant import SubMerchant

class PaymentFacilitator(AttributeGetter):
    sub_merchant: SubMerchant
    def __init__(self, attributes) -> None: ...
