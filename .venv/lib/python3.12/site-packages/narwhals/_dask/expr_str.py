from __future__ import annotations

from typing import TYPE_CHECKING

import dask.dataframe as dd

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import StringNamespace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    import dask.dataframe.dask_expr as dx

    from narwhals._dask.expr import DaskExpr


class DaskExprStringNamespace(LazyExprNamespace["DaskExpr"], StringNamespace["DaskExpr"]):
    def len_chars(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.len())

    def replace(
        self, value: DaskExpr, pattern: str, *, literal: bool, n: int
    ) -> DaskExpr:
        if not value._metadata.is_literal:
            msg = "dask backed `Expr.str.replace` only supports str replacement values"
            raise TypeError(msg)

        def _replace(expr: dx.Series, value: dx.Series) -> dx.Series:
            # OK to call `compute` here as `value` is just a literal expression.
            return expr.str.replace(  # pyright: ignore[reportAttributeAccessIssue]
                pattern, value.compute(), regex=not literal, n=n
            )

        return self.compliant._with_callable(_replace, value=value)

    def replace_all(self, value: DaskExpr, pattern: str, *, literal: bool) -> DaskExpr:
        return self.replace(value, pattern, literal=literal, n=-1)

    def strip_chars(self, characters: str | None) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.strip(characters))

    def starts_with(self, prefix: str) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.startswith(prefix))

    def ends_with(self, suffix: str) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.endswith(suffix))

    def contains(self, pattern: str, *, literal: bool) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.contains(pat=pattern, regex=not literal)
        )

    def slice(self, offset: int, length: int | None) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.slice(
                start=offset, stop=offset + length if length else None
            )
        )

    def split(self, by: str) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.split(pat=by))

    def to_datetime(self, format: str | None) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: dd.to_datetime(expr, format=format)
        )

    def to_uppercase(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.upper())

    def to_lowercase(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.lower())

    def to_titlecase(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.title())

    def zfill(self, width: int) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.zfill(width))

    def pad_start(self, length: int, fill_char: str) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.rjust(width=length, fillchar=fill_char)
        )

    def pad_end(self, length: int, fill_char: str) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.ljust(width=length, fillchar=fill_char)
        )

    to_date = not_implemented()
