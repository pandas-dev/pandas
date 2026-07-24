from collections.abc import Callable, Collection, Iterable, Mapping
from typing import Any, TypeVar, overload
from typing_extensions import TypeAlias

from ._utils import NO_DEFAULT, ExtractorError

_Traversable: TypeAlias = Mapping[str, Any] | Iterable[Any]
_PathArg: TypeAlias = str | int

def traverse_obj(
    obj: _Traversable,
    *paths: _PathArg,
    default: Any = ...,  # Anything or type[NO_DEFAULT]
    expected_type: type[Any] | None = None,
    get_all: bool = True,
    casesense: bool = True,
    is_user_input: bool | type[NO_DEFAULT] = ...,
    traverse_string: bool = False,
) -> Any: ...  # Unknown return type

_T = TypeVar("_T")

def value(value: _T, /) -> _T: ...
def require(name: str, /, *, expected: bool = False) -> Callable[[_T], _T]: ...

class _RequiredError(ExtractorError): ...

@overload
def subs_list_to_dict(
    *, lang: str | None = "und", ext: str | None = None
) -> Callable[[list[dict[str, Any]]], dict[str, list[dict[str, Any]]]]: ...
@overload
def subs_list_to_dict(
    subs: list[dict[str, Any]] | None, /, *, lang: str | None = "und", ext: str | None = None
) -> dict[str, list[dict[str, Any]]]: ...
@overload
def find_element(*, attr: str, value: str, tag: str | None = None, html: bool = False, regex: bool = False) -> str: ...
@overload
def find_element(*, cls: str, html: bool = False) -> str: ...
@overload
def find_element(*, id: str, tag: str | None = None, html: bool = False, regex: bool = False) -> str: ...
@overload
def find_element(*, tag: str, html: bool = False, regex: bool = False) -> str: ...
@overload
def find_element(
    *,
    tag: str | None = None,
    id: str | None = None,
    cls: str | None = None,
    attr: str | None = None,
    value: str | None = None,
    html: bool = False,
    regex: bool = False,
) -> str: ...
@overload
def find_elements(*, cls: str, html: bool = False) -> list[str]: ...
@overload
def find_elements(*, attr: str, value: str, tag: str | None = None, html: bool = False, regex: bool = False) -> list[str]: ...
@overload
def find_elements(
    *,
    tag: str | None = None,
    cls: str | None = None,
    attr: str | None = None,
    value: str | None = None,
    html: bool = False,
    regex: bool = False,
) -> list[str]: ...
def trim_str(*, start: str | None = None, end: str | None = None) -> Callable[[str], str]: ...

# Returns a callable f(items) which calls func(*items, **kwargs).
def unpack(func: Callable[..., Any], **kwargs: Any) -> Callable[..., Any]: ...
def get_first(
    obj: _Traversable,
    *paths: _PathArg,
    default: Any = ...,  # Anything or type[NO_DEFAULT]
    expected_type: type[Any] | None = None,
    get_all: bool = True,
    casesense: bool = True,
    is_user_input: bool | type[NO_DEFAULT] = ...,
    traverse_string: bool = False,
) -> Any: ...
@overload
def dict_get(d: str, key_or_keys: str | Collection[str]) -> Any | None: ...
@overload
def dict_get(
    d: str, key_or_keys: str | Collection[str], default: Any | None = None, skip_false_values: bool = True
) -> Any | None: ...
