from collections.abc import Callable
from typing import (
    Any,
    Literal,
    overload,
)

@overload
def ujson_dumps(
    obj: Any,
    ensure_ascii: bool = ...,
    double_precision: int | None = ...,
    encode_html_chars: bool = ...,
    orient: str = ...,
    date_unit: str = ...,
    iso_dates: bool = ...,
    default_handler: Callable[[Any], str | float | bool | list | dict | None]
    | None = ...,
    indent: int = ...,
    detect_float_format_change: Literal[False] = ...,
) -> str: ...
@overload
def ujson_dumps(
    obj: Any,
    ensure_ascii: bool = ...,
    double_precision: int | None = ...,
    encode_html_chars: bool = ...,
    orient: str = ...,
    date_unit: str = ...,
    iso_dates: bool = ...,
    default_handler: Callable[[Any], str | float | bool | list | dict | None]
    | None = ...,
    indent: int = ...,
    *,
    detect_float_format_change: Literal[True],
) -> tuple[str, bool]: ...
@overload
def ujson_loads(
    s: str,
    precise_float: bool = ...,
    detect_float_parse_change: Literal[False] = ...,
) -> Any: ...
@overload
def ujson_loads(
    s: str,
    precise_float: bool = ...,
    *,
    detect_float_parse_change: Literal[True],
) -> tuple[Any, bool]: ...
