from typing import Any

class GeneratedServiceType(type):
    def __init__(cls, name: str, bases: tuple[type, ...], dictionary: dict[str, Any]) -> None: ...

class GeneratedServiceStubType(GeneratedServiceType):
    def __init__(cls, name: str, bases: tuple[type, ...], dictionary: dict[str, Any]) -> None: ...
