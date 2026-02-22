from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace
from narwhals._duckdb.utils import F, col, lambda_expr, lit, when
from narwhals._utils import requires

if TYPE_CHECKING:
    from duckdb import Expression

    from narwhals._duckdb.expr import DuckDBExpr
    from narwhals.typing import NonNestedLiteral


class DuckDBExprListNamespace(
    LazyExprNamespace["DuckDBExpr"], ListNamespace["DuckDBExpr"]
):
    def len(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("len", expr))

    @requires.backend_version((1, 3))  # bugged before 1.3
    def unique(self) -> DuckDBExpr:
        def func(expr: Expression) -> Expression:
            expr_distinct = F("list_distinct", expr)
            return when(
                F("array_position", expr, lit(None)).isnotnull(),
                F("list_append", expr_distinct, lit(None)),
            ).otherwise(expr_distinct)

        return self.compliant._with_callable(func)

    def contains(self, item: NonNestedLiteral) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("list_contains", expr, lit(item))
        )

    def get(self, index: int) -> DuckDBExpr:
        return self.compliant._with_elementwise(
            lambda expr: F("list_extract", expr, lit(index + 1))
        )

    def min(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("list_min", expr))

    def max(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("list_max", expr))

    def mean(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("list_avg", expr))

    def median(self) -> DuckDBExpr:
        return self.compliant._with_elementwise(lambda expr: F("list_median", expr))

    @requires.backend_version((1, 2))
    def sum(self) -> DuckDBExpr:
        def func(expr: Expression) -> Expression:
            elem = col("_")
            expr_no_nulls = F("list_filter", expr, lambda_expr(elem, elem.isnotnull()))
            expr_sum = F("list_sum", expr_no_nulls)
            return when(F("array_length", expr_no_nulls) == lit(0), lit(0)).otherwise(
                expr_sum
            )

        return self.compliant._with_callable(func)

    def sort(self, *, descending: bool, nulls_last: bool) -> DuckDBExpr:
        sort_direction = "DESC" if descending else "ASC"
        nulls_position = "NULLS LAST" if nulls_last else "NULLS FIRST"
        return self.compliant._with_elementwise(
            lambda expr: F("list_sort", expr, lit(sort_direction), lit(nulls_position))
        )
