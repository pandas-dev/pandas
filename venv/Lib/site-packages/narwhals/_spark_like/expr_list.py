from __future__ import annotations

import operator
from typing import TYPE_CHECKING

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace

if TYPE_CHECKING:
    from sqlframe.base.column import Column

    from narwhals._spark_like.expr import SparkLikeExpr
    from narwhals.typing import NonNestedLiteral


class SparkLikeExprListNamespace(
    LazyExprNamespace["SparkLikeExpr"], ListNamespace["SparkLikeExpr"]
):
    def len(self) -> SparkLikeExpr:
        return self.compliant._with_elementwise(self.compliant._F.array_size)

    def unique(self) -> SparkLikeExpr:
        return self.compliant._with_elementwise(self.compliant._F.array_distinct)

    def contains(self, item: NonNestedLiteral) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            return F.array_contains(expr, F.lit(item))

        return self.compliant._with_elementwise(func)

    def get(self, index: int) -> SparkLikeExpr:
        def _get(expr: Column) -> Column:
            return expr.getItem(index)

        return self.compliant._with_elementwise(_get)

    def min(self) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            return F.array_min(expr)

        return self.compliant._with_elementwise(func)

    def max(self) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            return F.array_max(F.array_compact(expr))

        return self.compliant._with_elementwise(func)

    def sum(self) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            drop_nulls = F.array_compact(expr)
            len = F.array_size(drop_nulls)
            sum = F.aggregate(drop_nulls, F.lit(0.0), operator.add)
            return F.when((len.isNotNull()) & (len == 0), F.lit(0)).otherwise(sum)

        return self.compliant._with_elementwise(func)

    def mean(self) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            return F.try_divide(
                F.aggregate(F.array_compact(expr), F.lit(0.0), operator.add),
                F.array_size(F.array_compact(expr)),
            )

        return self.compliant._with_elementwise(func)

    def median(self) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            sorted_expr = F.array_compact(F.sort_array(expr))
            size = F.array_size(sorted_expr)
            mid_index = (size / 2).cast("int")
            impl = self.compliant._implementation
            if impl.is_sqlframe() and impl._backend_version() < (3, 44, 1):
                # indexing is different in sqlframe
                # https://github.com/eakmanrq/sqlframe/issues/568
                mid_index = mid_index + 1  # pragma: no cover
            odd_case = sorted_expr[mid_index]
            even_case = (sorted_expr[mid_index - 1] + sorted_expr[mid_index]) / 2
            return (
                F.when((size.isNull()) | (size == 0), F.lit(None))
                .when(size % 2 == 1, odd_case)
                .otherwise(even_case)
            )

        return self.compliant._with_elementwise(func)

    def sort(self, *, descending: bool, nulls_last: bool) -> SparkLikeExpr:
        def func(expr: Column) -> Column:
            F = self.compliant._F
            if not descending and nulls_last:
                return F.array_sort(expr)
            if descending and not nulls_last:
                impl = self.compliant._implementation
                rev = F.array_reverse if impl.is_sqlframe() else F.reverse
                return rev(F.array_sort(expr))
            return F.sort_array(expr, asc=not descending)

        return self.compliant._with_elementwise(func)
