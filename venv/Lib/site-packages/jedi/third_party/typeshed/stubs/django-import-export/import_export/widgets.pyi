from collections.abc import Mapping
from datetime import datetime
from typing import Any, ClassVar, Generic, TypeVar, overload
from typing_extensions import deprecated

from django.db.models import Model, QuerySet

def format_datetime(value: datetime, datetime_format: str) -> str: ...

class Widget:
    coerce_to_string: bool
    def __init__(self, coerce_to_string: bool = True) -> None: ...
    def clean(self, value: Any, row: Mapping[str, Any] | None = None, **kwargs: Any) -> Any: ...
    @overload
    @deprecated("The 'obj' parameter is deprecated and will be removed in a future release.")
    def render(self, value: Any, obj: Model, **kwargs: Any) -> Any: ...
    @overload
    def render(self, value: Any, obj: None = None, **kwargs: Any) -> Any: ...

class NumberWidget(Widget):
    def is_empty(self, value: Any) -> bool: ...

class FloatWidget(NumberWidget): ...
class IntegerWidget(NumberWidget): ...
class DecimalWidget(NumberWidget): ...

class CharWidget(Widget):
    allow_blank: bool
    def __init__(self, coerce_to_string: bool = True, allow_blank: bool = True) -> None: ...

class BooleanWidget(Widget):
    TRUE_VALUES: ClassVar[list[str | int | bool]]
    FALSE_VALUES: ClassVar[list[str | int | bool]]
    NULL_VALUES: ClassVar[list[str | None]]
    def __init__(self, coerce_to_string: bool = True) -> None: ...

class DateWidget(Widget):
    formats: tuple[str, ...]
    def __init__(self, format: str | None = None, coerce_to_string: bool = True) -> None: ...

class DateTimeWidget(Widget):
    formats: tuple[str, ...]
    def __init__(self, format: str | None = None, coerce_to_string: bool = True) -> None: ...

class TimeWidget(Widget):
    formats: tuple[str, ...]
    def __init__(self, format: str | None = None, coerce_to_string: bool = True) -> None: ...

class DurationWidget(Widget): ...

class SimpleArrayWidget(Widget):
    separator: str
    def __init__(self, separator: str | None = None, coerce_to_string: bool = True) -> None: ...

class JSONWidget(Widget): ...

_ModelT = TypeVar("_ModelT", bound=Model)

class ForeignKeyWidget(Widget, Generic[_ModelT]):
    model: type[_ModelT]
    field: str
    key_is_id: bool
    use_natural_foreign_keys: bool
    def __init__(
        self,
        model: type[_ModelT],
        field: str = "pk",
        use_natural_foreign_keys: bool = False,
        key_is_id: bool = False,
        **kwargs: Any,
    ) -> None: ...
    def get_queryset(self, value: Any, row: Mapping[str, Any], *args: Any, **kwargs: Any) -> QuerySet[_ModelT]: ...
    def get_instance_by_natural_key(self, value: str | bytes | bytearray) -> _ModelT: ...
    def get_instance_by_lookup_fields(self, value: Any, row: Mapping[str, Any], **kwargs: Any) -> _ModelT: ...
    def get_lookup_kwargs(self, value: Any, row: Mapping[str, Any] | None = None, **kwargs: Any) -> dict[str, Any]: ...

class _CachedQuerySetWrapper(Generic[_ModelT]):
    queryset: QuerySet[_ModelT]
    model: type[_ModelT]
    def __init__(self, queryset: QuerySet[_ModelT]) -> None: ...
    def get(self, **lookup_fields: Any) -> _ModelT: ...  # instance can have different fields

class CachedForeignKeyWidget(ForeignKeyWidget[_ModelT]):
    def get_instance_by_lookup_fields(self, value: Any, row: Mapping[str, Any], **kwargs: Any) -> _ModelT: ...

class ManyToManyWidget(Widget, Generic[_ModelT]):
    model: _ModelT
    separator: str
    field: str
    def __init__(self, model: _ModelT, separator: str | None = None, field: str | None = None, **kwargs: Any) -> None: ...
