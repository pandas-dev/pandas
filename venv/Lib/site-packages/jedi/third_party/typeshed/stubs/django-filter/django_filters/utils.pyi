from collections.abc import Callable
from datetime import datetime
from typing import Any

from django.db import models
from django.db.models import Model

def deprecate(msg: str, level_modifier: int = 0) -> None: ...

class MigrationNotice(DeprecationWarning):
    url: str
    def __init__(self, message: str) -> None: ...

class RenameAttributesBase(type):
    renamed_attributes: tuple[tuple[str, str, DeprecationWarning], ...] = ()

    # Class attrs vary by definition
    def __new__(metacls, name: str, bases: tuple[type, ...], attrs: dict[str, Any]) -> RenameAttributesBase: ...
    def get_name(metacls, name: str) -> str: ...
    def __getattr__(metacls, name: str) -> Any: ...  # Attribute values vary by name and class
    def __setattr__(metacls, name: str, value: Any) -> None: ...  # Attribute values can be any type

def try_dbfield(
    fn: Callable[[models.Field[Any, Any]], Any], field_class: type[models.Field[Any, Any]]
) -> Any: ...  # Generic field operation
def get_all_model_fields(model: type[Model]) -> dict[str, models.Field[Any, Any]]: ...  # Fields vary by model definition
def get_model_field(model: type[Model], field_name: str) -> models.Field[Any, Any]: ...  # Field type unknown at static time
def get_field_parts(model: type[Model], field_name: str) -> list[models.Field[Any, Any]]: ...  # Relationship fields vary
def resolve_field(
    model_field: models.Field[Any, Any], lookup_expr: str
) -> tuple[models.Field[Any, Any], str]: ...  # Generic field resolution
def handle_timezone(value: datetime, is_dst: bool | None = None) -> datetime: ...
def verbose_field_name(model: type[Model], field_name: str) -> str: ...
def verbose_lookup_expr(lookup_expr: str) -> str: ...
def label_for_filter(model: type[Model], field_name: str, lookup_expr: str, exclude: bool = False) -> str: ...
def translate_validation(error_dict: dict[str, list[str]]) -> dict[str, list[str]]: ...
