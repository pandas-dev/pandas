from typing import Any

class AttributeGetter:
    def __init__(self, attributes: dict[str, Any] | None = None) -> None: ...
    # This doesn't exist at runtime, but subclasses should define their own fields populated by attributes in __init__
    # Until that's done, keep __getattribute__ to fill in the gaps
    def __getattribute__(self, name: str, /) -> Any: ...
