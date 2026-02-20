from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, Any

import ibis
import ibis.expr.types as ir

from narwhals._compliant.namespace import AlignDiagonal
from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
)
from narwhals._ibis.dataframe import IbisLazyFrame
from narwhals._ibis.expr import IbisExpr
from narwhals._ibis.selectors import IbisSelectorNamespace
from narwhals._ibis.utils import function, lit, narwhals_to_native_dtype
from narwhals._sql.namespace import SQLNamespace
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from narwhals._utils import Version
    from narwhals.typing import ConcatMethod, IntoDType, PythonLiteral


class IbisNamespace(
    SQLNamespace[IbisLazyFrame, IbisExpr, "ir.Table", "ir.Value"],
    AlignDiagonal[IbisLazyFrame, IbisExpr],
):
    _implementation: Implementation = Implementation.IBIS

    def __init__(self, *, version: Version) -> None:
        self._version = version

    @property
    def selectors(self) -> IbisSelectorNamespace:
        return IbisSelectorNamespace.from_namespace(self)

    @property
    def _expr(self) -> type[IbisExpr]:
        return IbisExpr

    @property
    def _lazyframe(self) -> type[IbisLazyFrame]:
        return IbisLazyFrame

    def _function(self, name: str, *args: ir.Value | PythonLiteral) -> ir.Value:
        return function(name, *args)

    def _lit(self, value: Any) -> ir.Value:
        return lit(value)

    def _when(
        self, condition: ir.Value, value: ir.Value, otherwise: ir.Expr | None = None
    ) -> ir.Value:
        if otherwise is None:
            return ibis.cases((condition, value))
        return ibis.cases((condition, value), else_=otherwise)  # pragma: no cover

    def _coalesce(self, *exprs: ir.Value) -> ir.Value:
        return ibis.coalesce(*exprs)

    def concat(
        self, items: Iterable[IbisLazyFrame], *, how: ConcatMethod
    ) -> IbisLazyFrame:
        frames: Sequence[IbisLazyFrame] = tuple(items)
        if how == "diagonal":
            frames = self.align_diagonal(frames)
        try:
            result = ibis.union(*(lf.native for lf in frames))
        except ibis.IbisError:
            first = frames[0].schema
            if not all(x.schema == first for x in frames[1:]):
                msg = "inputs should all have the same schema"
                raise TypeError(msg) from None
            raise
        return frames[0]._with_native(result)

    def concat_str(
        self, *exprs: IbisExpr, separator: str, ignore_nulls: bool
    ) -> IbisExpr:
        def func(df: IbisLazyFrame) -> list[ir.Value]:
            cols = chain.from_iterable(expr(df) for expr in exprs)
            cols_casted = [s.cast("string") for s in cols]

            if ignore_nulls:
                result = lit(separator).join(cols_casted)
            else:
                result = reduce(
                    lambda acc, col: acc.concat(separator, col),
                    cols_casted[1:],
                    cols_casted[0],
                )

            return [result]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def mean_horizontal(self, *exprs: IbisExpr) -> IbisExpr:
        def func(cols: Iterable[ir.Value]) -> ir.Value:
            cols = list(cols)
            return reduce(operator.add, (col.fill_null(lit(0)) for col in cols)) / reduce(
                operator.add, (col.isnull().ifelse(lit(0), lit(1)) for col in cols)
            )

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def lit(self, value: PythonLiteral, dtype: IntoDType | None) -> IbisExpr:
        def func(_df: IbisLazyFrame) -> Sequence[ir.Value]:
            ibis_dtype = narwhals_to_native_dtype(dtype, self._version) if dtype else None
            if not isinstance(value, dict):
                return [lit(value, ibis_dtype)]
            if value:
                return [ibis.struct(value, type=ibis_dtype)]
            msg = "Cannot create an empty struct type for Ibis backend"
            raise NotImplementedError(msg)

        return self._expr(
            func,
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
        )

    def len(self) -> IbisExpr:
        def func(_df: IbisLazyFrame) -> list[ir.Value]:
            return [_df.native.count()]

        return self._expr(
            call=func,
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            version=self._version,
        )
