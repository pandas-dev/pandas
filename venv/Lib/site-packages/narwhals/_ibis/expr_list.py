from __future__ import annotations

from typing import TYPE_CHECKING

from ibis import cases, literal

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import ListNamespace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    import ibis.expr.types as ir

    from narwhals._ibis.expr import IbisExpr
    from narwhals.typing import NonNestedLiteral


class IbisExprListNamespace(LazyExprNamespace["IbisExpr"], ListNamespace["IbisExpr"]):
    def len(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.length())

    def unique(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.unique())

    def contains(self, item: NonNestedLiteral) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.contains(item))

    def get(self, index: int) -> IbisExpr:
        def _get(expr: ir.ArrayColumn) -> ir.Column:
            return expr[index]

        return self.compliant._with_callable(_get)

    def min(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.mins())

    def max(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.maxs())

    def mean(self) -> IbisExpr:
        return self.compliant._with_callable(lambda expr: expr.means())

    def sum(self) -> IbisExpr:
        def func(expr: ir.ArrayColumn) -> ir.Value:
            expr_no_nulls = expr.filter(lambda x: x.notnull())
            len = expr_no_nulls.length()
            return cases(
                (len.isnull(), literal(None)),
                (len == literal(0), literal(0)),
                else_=expr.sums(),
            )

        return self.compliant._with_callable(func)

    def sort(self, *, descending: bool, nulls_last: bool) -> IbisExpr:
        if descending:
            msg = "Descending sort is not currently supported for Ibis."
            raise NotImplementedError(msg)

        def func(expr: ir.ArrayColumn) -> ir.ArrayValue:
            if nulls_last:
                return expr.sort()
            expr_no_nulls = expr.filter(lambda x: x.notnull())
            expr_nulls = expr.filter(lambda x: x.isnull())
            return expr_nulls.concat(expr_no_nulls.sort())

        return self.compliant._with_callable(func)

    median = not_implemented()
