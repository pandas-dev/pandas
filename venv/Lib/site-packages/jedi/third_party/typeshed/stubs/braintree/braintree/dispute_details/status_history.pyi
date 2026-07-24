from typing import Any

from braintree.attribute_getter import AttributeGetter

class DisputeStatusHistory(AttributeGetter):
    def __init__(self, attributes: dict[str, Any] | None) -> None: ...
