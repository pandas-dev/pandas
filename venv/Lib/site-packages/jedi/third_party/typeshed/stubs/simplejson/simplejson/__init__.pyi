from _typeshed import SupportsRichComparison
from collections import OrderedDict
from collections.abc import Callable
from typing import IO, Any, TypeVar, overload
from typing_extensions import TypeAlias

from simplejson.decoder import JSONDecoder as JSONDecoder
from simplejson.encoder import JSONEncoder as JSONEncoder, JSONEncoderForHTML as JSONEncoderForHTML
from simplejson.raw_json import RawJSON as RawJSON
from simplejson.scanner import JSONDecodeError as JSONDecodeError

_LoadsString: TypeAlias = str | bytes | bytearray
_T = TypeVar("_T")

@overload
def dumps(
    obj: Any,
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = False,
    *,
    cls: type[JSONEncoder],
    indent: str | int | None = None,
    separators: tuple[str, str] | None = None,
    encoding: str | None = "utf-8",
    default: Callable[[Any], Any] | None = None,
    use_decimal: bool = True,
    namedtuple_as_object: bool = True,
    tuple_as_array: bool = True,
    bigint_as_string: bool = False,
    sort_keys: bool = False,
    item_sort_key: Callable[[Any], SupportsRichComparison] | None = None,
    for_json: bool = False,
    ignore_nan: bool = False,
    int_as_string_bitcount: int | None = None,
    iterable_as_array: bool = False,
    **kw: Any,
) -> str: ...
@overload
def dumps(
    obj: Any,
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = False,
    cls: type[JSONEncoder] | None = None,
    indent: str | int | None = None,
    separators: tuple[str, str] | None = None,
    encoding: str | None = "utf-8",
    default: Callable[[Any], Any] | None = None,
    use_decimal: bool = True,
    namedtuple_as_object: bool = True,
    tuple_as_array: bool = True,
    bigint_as_string: bool = False,
    sort_keys: bool = False,
    item_sort_key: Callable[[Any], SupportsRichComparison] | None = None,
    for_json: bool = False,
    ignore_nan: bool = False,
    int_as_string_bitcount: int | None = None,
    iterable_as_array: bool = False,
) -> str: ...
@overload
def dump(
    obj: Any,
    fp: IO[str],
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = False,
    *,
    cls: type[JSONEncoder],
    indent: str | int | None = None,
    separators: tuple[str, str] | None = None,
    encoding: str | None = "utf-8",
    default: Callable[[Any], Any] | None = None,
    use_decimal: bool = True,
    namedtuple_as_object: bool = True,
    tuple_as_array: bool = True,
    bigint_as_string: bool = False,
    sort_keys: bool = False,
    item_sort_key: Callable[[Any], SupportsRichComparison] | None = None,
    for_json: bool = False,
    ignore_nan: bool = False,
    int_as_string_bitcount: int | None = None,
    iterable_as_array: bool = False,
    **kw: Any,
) -> None: ...
@overload
def dump(
    obj: Any,
    fp: IO[str],
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = False,
    cls: type[JSONEncoder] | None = None,
    indent: str | int | None = None,
    separators: tuple[str, str] | None = None,
    encoding: str | None = "utf-8",
    default: Callable[[Any], Any] | None = None,
    use_decimal: bool = True,
    namedtuple_as_object: bool = True,
    tuple_as_array: bool = True,
    bigint_as_string: bool = False,
    sort_keys: bool = False,
    item_sort_key: Callable[[Any], SupportsRichComparison] | None = None,
    for_json: bool = False,
    ignore_nan: bool = False,
    int_as_string_bitcount: int | None = None,
    iterable_as_array: bool = False,
) -> None: ...
@overload
def loads(
    s: _LoadsString,
    encoding: str | None = None,
    *,
    cls: type[JSONDecoder],
    object_hook: Callable[[dict[Any, Any]], Any] | None = None,
    parse_float: Callable[[str], Any] | None = None,
    parse_int: Callable[[str], Any] | None = None,
    parse_constant: Callable[[str], Any] | None = None,
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    use_decimal: bool = False,
    allow_nan: bool = False,
    **kw: Any,
) -> Any: ...
@overload
def loads(
    s: _LoadsString,
    encoding: str | None = None,
    cls: type[JSONDecoder] | None = None,
    object_hook: Callable[[dict[Any, Any]], Any] | None = None,
    parse_float: Callable[[str], Any] | None = None,
    parse_int: Callable[[str], Any] | None = None,
    parse_constant: Callable[[str], Any] | None = None,
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    use_decimal: bool = False,
    allow_nan: bool = False,
) -> Any: ...
@overload
def load(
    fp: IO[str],
    encoding: str | None = None,
    *,
    cls: type[JSONDecoder],
    object_hook: Callable[[dict[Any, Any]], Any] | None = None,
    parse_float: Callable[[str], Any] | None = None,
    parse_int: Callable[[str], Any] | None = None,
    parse_constant: Callable[[str], Any] | None = None,
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    use_decimal: bool = False,
    allow_nan: bool = False,
    **kw: Any,
) -> Any: ...
@overload
def load(
    fp: IO[str],
    encoding: str | None = None,
    cls: type[JSONDecoder] | None = None,
    object_hook: Callable[[dict[Any, Any]], Any] | None = None,
    parse_float: Callable[[str], Any] | None = None,
    parse_int: Callable[[str], Any] | None = None,
    parse_constant: Callable[[str], Any] | None = None,
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    use_decimal: bool = False,
    allow_nan: bool = False,
) -> Any: ...
def simple_first(kv: tuple[_T, object]) -> tuple[bool, _T]: ...

__all__ = [
    "dump",
    "dumps",
    "load",
    "loads",
    "JSONDecoder",
    "JSONDecodeError",
    "JSONEncoder",
    "OrderedDict",
    "simple_first",
    "RawJSON",
]
