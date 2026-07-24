from typing import Final

PY3: Final[bool]

class DefusedXmlException(ValueError): ...

class DTDForbidden(DefusedXmlException):
    name: str
    sysid: str | None
    pubid: str | None
    def __init__(self, name: str, sysid: str | None, pubid: str | None) -> None: ...

class EntitiesForbidden(DefusedXmlException):
    name: str
    value: str | None
    base: str | None
    sysid: str | None
    pubid: str | None
    notation_name: str | None
    def __init__(
        self, name: str, value: str | None, base: str | None, sysid: str | None, pubid: str | None, notation_name: str | None
    ) -> None: ...

class ExternalReferenceForbidden(DefusedXmlException):
    context: str
    base: str | None
    sysid: str | None
    pubid: str | None
    def __init__(self, context: str, base: str | None, sysid: str | None, pubid: str | None) -> None: ...

class NotSupportedError(DefusedXmlException): ...
