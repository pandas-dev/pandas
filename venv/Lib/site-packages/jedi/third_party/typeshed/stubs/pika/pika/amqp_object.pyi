from typing import ClassVar

class AMQPObject:
    NAME: ClassVar[str]
    INDEX: ClassVar[int | None]
    def __eq__(self, other: AMQPObject | None) -> bool: ...  # type: ignore[override]

class Class(AMQPObject): ...

class Method(AMQPObject):
    # This is a class attribute in the implementation, but subclasses use @property,
    # so it's more convenient to use that here as well.
    @property
    def synchronous(self) -> bool: ...
    def get_properties(self) -> Properties: ...
    def get_body(self) -> str: ...

class Properties(AMQPObject): ...
