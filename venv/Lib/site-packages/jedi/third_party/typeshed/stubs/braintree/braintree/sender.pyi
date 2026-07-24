from typing import Any

from braintree.attribute_getter import AttributeGetter

class Sender(AttributeGetter):
    def __init__(self, attributes: dict[str, Any] | None) -> None: ...
