import re
from _typeshed import SupportsRichComparison
from collections.abc import Callable, Iterator
from typing import Any, Literal, NoReturn

ESCAPE: re.Pattern[str]
ESCAPE_ASCII: re.Pattern[str]
HAS_UTF8: re.Pattern[str]
ESCAPE_DCT: dict[str, str]
FLOAT_REPR: Callable[[object], str]

class JSONEncoder:
    item_separator: str
    key_separator: str
    skipkeys: bool
    ensure_ascii: bool
    check_circular: bool
    allow_nan: bool
    sort_keys: bool
    indent: str
    encoding: str
    use_decimal: bool
    namedtuple_as_object: bool
    tuple_as_array: bool
    bigint_as_string: bool
    item_sort_key: Callable[[Any], SupportsRichComparison] | None
    for_json: bool
    ignore_nan: bool
    int_as_string_bitcount: int | None
    iterable_as_array: bool

    def __init__(
        self,
        skipkeys: bool = False,
        ensure_ascii: bool = True,
        check_circular: bool = True,
        allow_nan: bool = False,
        sort_keys: bool = False,
        indent: str | int | None = None,
        separators: tuple[str, str] | None = None,
        encoding: str = "utf-8",
        default: Callable[[Any], Any] | None = None,
        use_decimal: bool = True,
        namedtuple_as_object: bool = True,
        tuple_as_array: bool = True,
        bigint_as_string: bool = False,
        item_sort_key: Callable[[Any], SupportsRichComparison] | None = None,
        for_json: bool = False,
        ignore_nan: bool = False,
        int_as_string_bitcount: int | None = None,
        iterable_as_array: bool = False,
    ) -> None: ...
    def encode(self, o: Any) -> str: ...
    def default(self, o: Any) -> NoReturn: ...
    def iterencode(self, o: Any) -> Iterator[str]: ...

class JSONEncoderForHTML(JSONEncoder): ...

def encode_basestring(s: str | bytes, _PY3: Literal[True] = True, _q: str = '"') -> str: ...
def encode_basestring_ascii(s: str | bytes, /) -> str: ...
