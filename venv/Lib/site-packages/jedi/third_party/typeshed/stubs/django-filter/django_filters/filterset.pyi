from collections import OrderedDict
from collections.abc import Sequence
from enum import Enum
from typing import Any, ClassVar

from django.db import models
from django.db.models import Model, QuerySet
from django.forms import Form
from django.http import HttpRequest, QueryDict

from .filters import Filter

def remote_queryset(field: models.Field[Any, Any]) -> QuerySet[Any]: ...  # Field type params vary by model definition

class UnknownFieldBehavior(Enum):
    RAISE = "raise"
    WARN = "warn"
    IGNORE = "ignore"

class FilterSetOptions:
    model: type[Model] | None
    fields: Sequence[str] | dict[str, Sequence[str]] | str | None
    exclude: Sequence[str] | None
    filter_overrides: dict[type[models.Field[Any, Any]], dict[str, Any]]  # Field override mapping
    form: type[Form]
    unknown_field_behavior: UnknownFieldBehavior
    def __init__(self, options: Any | None = None) -> None: ...  # Meta options can be various configuration types

class FilterSetMetaclass(type):
    # Class attrs vary by definition
    def __new__(cls, name: str, bases: tuple[type, ...], attrs: dict[str, Any]) -> FilterSetMetaclass: ...

    # Class attrs vary by definition
    @classmethod
    def get_declared_filters(cls, bases: tuple[type, ...], attrs: dict[str, Any]) -> OrderedDict[str, Filter]: ...

# Django field types vary widely - Any allows mapping all field types to their filters
FILTER_FOR_DBFIELD_DEFAULTS: dict[type[models.Field[Any, Any]], dict[str, Any]]

class BaseFilterSet:
    FILTER_DEFAULTS: ClassVar[dict[type[models.Field[Any, Any]], dict[str, Any]]] = ...  # Field type mapping
    is_bound: bool
    base_filters: OrderedDict[str, Filter]
    declared_filters: OrderedDict[str, Filter]
    data: QueryDict | dict[str, Any] | None  # Filter input data values vary
    queryset: QuerySet[Any] | None  # Base queryset for any model type
    request: HttpRequest | None
    form_prefix: str | None
    filters: OrderedDict[str, Filter]
    def __init__(
        self,
        data: QueryDict | dict[str, Any] | None = None,  # Filter data values vary
        queryset: QuerySet[Any] | None = None,  # Base queryset for any model
        *,
        request: HttpRequest | None = None,
        prefix: str | None = None,
    ) -> None: ...
    def is_valid(self) -> bool: ...
    @property
    def errors(self) -> dict[str, list[str]]: ...
    def filter_queryset(self, queryset: QuerySet[Any]) -> QuerySet[Any]: ...  # Works with any model type
    @property
    def qs(self) -> QuerySet[Any]: ...  # Filtered queryset of any model
    def get_form_class(self) -> type[Form]: ...
    @property
    def form(self) -> Form: ...
    @classmethod
    def get_fields(cls) -> dict[str, models.Field[Any, Any]]: ...  # Model fields have varying type params
    @classmethod
    def get_filter_name(cls, field_name: str, lookup_expr: str) -> str: ...
    @classmethod
    def get_filters(cls) -> OrderedDict[str, Filter]: ...
    @classmethod
    def handle_unrecognized_field(cls, field_name: str, message: str) -> None: ...
    @classmethod
    def filter_for_field(
        cls, field: models.Field[Any, Any], field_name: str, lookup_expr: str | None = None
    ) -> Filter: ...  # Accepts any Django field type
    @classmethod
    def filter_for_lookup(
        cls, field: models.Field[Any, Any], lookup_type: str  # Field type varies by model
    ) -> tuple[type[Filter], dict[str, Any]]: ...

class FilterSet(BaseFilterSet, metaclass=FilterSetMetaclass): ...

def filterset_factory(
    model: type[Model], filterset: FilterSetMetaclass = ..., fields: Sequence[str] | dict[str, Sequence[str]] | str | None = None
) -> type[FilterSet]: ...
