from collections.abc import Callable
from typing import (
    Any,
)

def ujson_dumps(
    obj: Any,
    ensure_ascii: bool = ...,
    double_precision: int = ...,
    indent: int = ...,
    orient: str = ...,
    date_unit: str = ...,
    iso_dates: bool = ...,
    default_handler: None
    | Callable[[Any], str | float | bool | list | dict | None] = ...,
) -> str: ...
