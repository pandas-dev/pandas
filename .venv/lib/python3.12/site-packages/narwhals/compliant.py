"""Public protocols for plugin authors."""

from __future__ import annotations

from narwhals._compliant import (
    CompliantDataFrame,
    CompliantExpr,
    CompliantGroupBy,
    CompliantLazyFrame,
    CompliantNamespace,
    CompliantSeries,
)
from narwhals._compliant.any_namespace import (
    CatNamespace,
    DateTimeNamespace,
    ListNamespace,
    StringNamespace,
    StructNamespace,
)
from narwhals._compliant.expr import CompliantExprNameNamespace
from narwhals._compliant.selectors import CompliantSelector, CompliantSelectorNamespace

__all__ = [
    "CatNamespace",
    "CompliantDataFrame",
    "CompliantExpr",
    "CompliantExprNameNamespace",
    "CompliantGroupBy",
    "CompliantLazyFrame",
    "CompliantNamespace",
    "CompliantSelector",
    "CompliantSelectorNamespace",
    "CompliantSeries",
    "DateTimeNamespace",
    "ListNamespace",
    "StringNamespace",
    "StructNamespace",
]
