import datetime
import decimal
from collections.abc import Iterable, Mapping
from typing_extensions import TypeAlias

integer_types = int
text_type = str
binary_type = bytes

_XMLValue: TypeAlias = (
    str
    | bytes
    | int
    | bool
    | decimal.Decimal
    | Iterable[_XMLValue]
    | Mapping[str, _XMLValue]
    | datetime.datetime
    | datetime.date
    | None
)
_XML: TypeAlias = Mapping[str, _XMLValue]

class Generator:
    dict: _XML
    def __init__(self, dict: _XML) -> None: ...
    def generate(self) -> str: ...
