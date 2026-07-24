from typing import Any

from django.core.exceptions import FieldError
from django.db import models

class FieldLookupError(FieldError):
    # Field type params are runtime-determined
    def __init__(self, model_field: models.Field[Any, Any], lookup_expr: str) -> None: ...
