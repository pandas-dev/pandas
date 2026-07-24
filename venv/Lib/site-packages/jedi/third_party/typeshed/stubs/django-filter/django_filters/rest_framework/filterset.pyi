from collections import OrderedDict
from typing import Any, ClassVar

from django.db import models
from django.forms import Form
from django_filters import filterset
from django_filters.filters import Filter

# REST framework field mappings support all Django field types
FILTER_FOR_DBFIELD_DEFAULTS: dict[type[models.Field[Any, Any]], dict[str, Any]]

class FilterSet(filterset.FilterSet):
    FILTER_DEFAULTS: ClassVar[dict[type[models.Field[Any, Any]], dict[str, Any]]] = ...  # DRF field mappings
    base_filters: OrderedDict[str, Filter]
    declared_filters: OrderedDict[str, Filter]
    @property
    def form(self) -> Form: ...
