from typing import ClassVar

from braintree.attribute_getter import AttributeGetter

class AddressDetails(AttributeGetter):
    detail_list: ClassVar[list[str]]
    def __init__(self, attributes) -> None: ...
