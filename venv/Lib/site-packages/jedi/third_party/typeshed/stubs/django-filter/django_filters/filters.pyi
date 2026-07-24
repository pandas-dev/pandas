from collections.abc import Callable, Iterable
from typing import Any

from django import forms
from django.db.models import QuerySet
from django.forms import Field
from django_stubs_ext import StrOrPromise

from .fields import (
    BaseCSVField,
    BaseRangeField,
    DateRangeField,
    DateTimeRangeField,
    IsoDateTimeField,
    IsoDateTimeRangeField,
    Lookup,
    LookupChoiceField,
    ModelChoiceField,
    ModelMultipleChoiceField,
    RangeField,
    TimeRangeField,
)

__all__ = [
    "AllValuesFilter",
    "AllValuesMultipleFilter",
    "BaseCSVFilter",
    "BaseInFilter",
    "BaseRangeFilter",
    "BooleanFilter",
    "CharFilter",
    "ChoiceFilter",
    "DateFilter",
    "DateFromToRangeFilter",
    "DateRangeFilter",
    "DateTimeFilter",
    "DateTimeFromToRangeFilter",
    "DurationFilter",
    "Filter",
    "IsoDateTimeFilter",
    "IsoDateTimeFromToRangeFilter",
    "LookupChoiceFilter",
    "ModelChoiceFilter",
    "ModelMultipleChoiceFilter",
    "MultipleChoiceFilter",
    "NumberFilter",
    "NumericRangeFilter",
    "OrderingFilter",
    "RangeFilter",
    "TimeFilter",
    "TimeRangeFilter",
    "TypedChoiceFilter",
    "TypedMultipleChoiceFilter",
    "UUIDFilter",
]

class Filter:
    creation_counter: int
    field_class: type[Any]  # Subclasses specify more specific field types
    field_name: str | None
    lookup_expr: str
    distinct: bool
    exclude: bool
    extra: dict[str, Any]  # Field kwargs can include various types of parameters
    def __init__(
        self,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        distinct: bool = False,
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...
    def get_method(self, qs: QuerySet[Any]) -> Callable[..., QuerySet[Any]]: ...  # Returns QuerySet filtering methods
    method: Callable[..., Any] | str | None  # Custom filter methods return various types
    label: StrOrPromise | None  # Filter label for display
    @property
    def field(self) -> Field: ...
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...  # Filter value can be any user input type

class CharFilter(Filter):
    field_class: type[forms.CharField]

class BooleanFilter(Filter):
    field_class: type[forms.NullBooleanField]

class ChoiceFilter(Filter):
    field_class: type[Any]  # Base class for choice-based filters
    null_value: Any  # Null value can be any type (None, empty string, etc.)
    def __init__(
        self,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        null_value: Any = ...,  # Null value can be any type (None, empty string, etc.)
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        distinct: bool = False,
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...

class TypedChoiceFilter(Filter):
    field_class: type[forms.TypedChoiceField]

class UUIDFilter(Filter):
    field_class: type[forms.UUIDField]

class MultipleChoiceFilter(Filter):
    field_class: type[Any]  # Base class for multiple choice filters
    always_filter: bool
    conjoined: bool
    null_value: Any  # Multiple choice null values vary by implementation
    def __init__(
        self,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        distinct: bool = True,  # Overrides distinct default
        conjoined: bool = False,
        null_value: Any = ...,  # Multiple choice null values vary by implementation
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...
    def is_noop(self, qs: QuerySet[Any], value: Any) -> bool: ...  # Value can be any filter input
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...
    def get_filter_predicate(self, v: Any) -> dict[str, Any]: ...  # Predicate value can be any filter input type

class TypedMultipleChoiceFilter(MultipleChoiceFilter):
    field_class: type[forms.TypedMultipleChoiceField]  # More specific than parent MultipleChoiceField

class DateFilter(Filter):
    field_class: type[forms.DateField]

class DateTimeFilter(Filter):
    field_class: type[forms.DateTimeField]

class IsoDateTimeFilter(DateTimeFilter):
    field_class: type[IsoDateTimeField]

class TimeFilter(Filter):
    field_class: type[forms.TimeField]

class DurationFilter(Filter):
    field_class: type[forms.DurationField]

class QuerySetRequestMixin:
    queryset: QuerySet[Any] | None
    def __init__(self, *, queryset: QuerySet[Any] | None) -> None: ...
    def get_request(self) -> Any: ...  # Request can be HttpRequest or other request types
    def get_queryset(self, request: Any) -> QuerySet[Any]: ...  # Request parameter accepts various request types
    @property
    def field(self) -> Field: ...

class ModelChoiceFilter(QuerySetRequestMixin, ChoiceFilter):
    field_class: type[ModelChoiceField]  # More specific than parent ChoiceField
    def __init__(
        self,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        # Inherited from QuerySetRequestMixin
        queryset: QuerySet[Any] | None = None,
        # Inherited from ChoiceFilter
        null_value: Any = ...,  # Null value can be any type (None, empty string, etc.)
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        distinct: bool = False,
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...

class ModelMultipleChoiceFilter(QuerySetRequestMixin, MultipleChoiceFilter):
    field_class: type[ModelMultipleChoiceField]  # More specific than parent MultipleChoiceField
    def __init__(
        self,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        # Inherited from QuerySetRequestMixin
        queryset: QuerySet[Any] | None = None,
        # Inherited from MultipleChoiceFilter
        distinct: bool = True,  # Overrides distinct default
        conjoined: bool = False,
        null_value: Any = ...,  # Multiple choice null values vary by implementation
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...

class NumberFilter(Filter):
    field_class: type[forms.DecimalField]
    def get_max_validator(self) -> Any: ...  # Validator can be various Django validator types
    @property
    def field(self) -> Field: ...

class NumericRangeFilter(Filter):
    field_class: type[RangeField]
    lookup_expr: str
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...

class RangeFilter(Filter):
    field_class: type[RangeField]
    lookup_expr: str
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...

class DateRangeFilter(ChoiceFilter):
    choices: list[tuple[str, str]] | None
    filters: dict[str, Filter] | None
    def __init__(
        self,
        choices: list[tuple[str, str]] | None = None,
        filters: dict[str, Filter] | None = None,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        # Inherited from ChoiceFilter
        null_value: Any = ...,  # Null value can be any type (None, empty string, etc.)
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        distinct: bool = False,
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...  # Uses args/kwargs for choice and filter configuration
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...

class DateFromToRangeFilter(RangeFilter):
    field_class: type[DateRangeField]

class DateTimeFromToRangeFilter(RangeFilter):
    field_class: type[DateTimeRangeField]

class IsoDateTimeFromToRangeFilter(RangeFilter):
    field_class: type[IsoDateTimeRangeField]

class TimeRangeFilter(RangeFilter):
    field_class: type[TimeRangeField]

class AllValuesFilter(ChoiceFilter):
    @property
    def field(self) -> Field: ...

class AllValuesMultipleFilter(MultipleChoiceFilter):
    @property
    def field(self) -> Field: ...

class BaseCSVFilter(Filter):
    base_field_class: type[BaseCSVField] = ...
    field_class: type[Any]  # Base class for CSV-based filters

class BaseInFilter(BaseCSVFilter): ...

class BaseRangeFilter(BaseCSVFilter):
    base_field_class: type[BaseRangeField] = ...

class LookupChoiceFilter(Filter):
    field_class: type[forms.CharField]
    outer_class: type[LookupChoiceField] = ...
    empty_label: StrOrPromise | None
    lookup_choices: list[tuple[str, StrOrPromise]] | None
    def __init__(
        self,
        field_name: str | None = None,
        lookup_choices: list[tuple[str, StrOrPromise]] | None = None,
        field_class: type[Field] | None = None,
        *,
        empty_label: StrOrPromise = ...,
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        distinct: bool = False,
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...
    @classmethod
    def normalize_lookup(cls, lookup: str | tuple[str, StrOrPromise]) -> tuple[str, StrOrPromise]: ...
    def get_lookup_choices(self) -> list[tuple[str, StrOrPromise]]: ...
    @property
    def field(self) -> Field: ...
    lookup_expr: str
    def filter(self, qs: QuerySet[Any], lookup: Lookup) -> QuerySet[Any]: ...

class OrderingFilter(BaseCSVFilter, ChoiceFilter):
    # Inherits CSV field behavior for comma-separated ordering.
    # BaseCSVFilter constructs a custom ConcreteCSVField class that derives
    # from BaseCSVField.
    field_class: type[BaseCSVField]
    descending_fmt: str
    param_map: dict[str, str] | None
    def __init__(
        self,
        field_name: str | None = None,
        lookup_expr: str | None = None,
        *,
        fields: dict[str, str] | Iterable[str] | Iterable[tuple[str, str]] = ...,
        field_labels: dict[str, StrOrPromise] = ...,
        # Inherited from ChoiceFilter
        null_value: Any = ...,  # Null value can be any type (None, empty string, etc.)
        # Inherited from Filter
        label: StrOrPromise | None = None,
        method: Callable[..., Any] | str | None = None,  # Filter methods can return various types
        distinct: bool = False,
        exclude: bool = False,
        **kwargs: Any,  # Field kwargs stored as extra (required, help_text, etc.)
    ) -> None: ...
    def get_ordering_value(self, param: str) -> str: ...
    def filter(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...
    @classmethod
    def normalize_fields(cls, fields: Any) -> list[str]: ...
    def build_choices(self, fields: Any, labels: dict[str, StrOrPromise] | None) -> list[tuple[str, str]]: ...

class FilterMethod:
    f: Filter
    def __init__(self, filter_instance: Filter) -> None: ...
    def __call__(self, qs: QuerySet[Any], value: Any) -> QuerySet[Any]: ...
    @property
    def method(self) -> Callable[..., Any]: ...
