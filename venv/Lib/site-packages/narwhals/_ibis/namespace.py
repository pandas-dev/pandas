from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, Any, cast

import ibis
import ibis.expr.types as ir

from narwhals._compliant.namespace import AlignDiagonal
from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
    evaluate_output_names_and_aliases,
)
from narwhals._ibis.dataframe import IbisLazyFrame
from narwhals._ibis.expr import IbisExpr
from narwhals._ibis.selectors import IbisSelectorNamespace
from narwhals._ibis.utils import function, lit, narwhals_to_native_dtype
from narwhals._sql.namespace import SQLNamespace
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

    from narwhals._utils import Version
    from narwhals.typing import ConcatMethod, CorrelationMethod, IntoDType, PythonLiteral


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
            numerator = reduce(
                operator.add,
                (cast("ir.NumericValue", col.fill_null(lit(0))) for col in cols),
            )
            denominator = reduce(
                operator.add,
                (
                    cast("ir.NumericValue", col.isnull().ifelse(lit(0), lit(1)))
                    for col in cols
                ),
            )
            return numerator / denominator

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

    def corr(self, a: IbisExpr, b: IbisExpr, *, method: CorrelationMethod) -> IbisExpr:
        if method != "pearson":
            msg = "Only 'pearson' correlation is supported for Ibis."
            raise NotImplementedError(msg)

        def func(_df: IbisLazyFrame) -> list[ir.Value]:
            a_ = _df._evaluate_single_output_expr(a)
            b_ = _df._evaluate_single_output_expr(b)
            return [a_.corr(b_, how="pop")]  # pyright: ignore[reportAttributeAccessIssue]

        return self._expr(
            func,
            evaluate_output_names=combine_evaluate_output_names(a, b),
            alias_output_names=combine_alias_output_names(a, b),
            version=self._version,
        )

    def cov(self, a: IbisExpr, b: IbisExpr, *, ddof: int) -> IbisExpr:
        def func(_df: IbisLazyFrame) -> list[ir.Value]:
            a_ = cast("ir.NumericColumn", _df._evaluate_single_output_expr(a))
            b_ = cast("ir.NumericColumn", _df._evaluate_single_output_expr(b))
            if ddof == 0:
                return [a_.cov(b_, how="pop")]
            if ddof == 1:
                return [a_.cov(b_, how="sample")]
            n_samples = a_.count(where=a_.notnull() & b_.notnull())
            denominator = n_samples - ddof
            rescaled = a_.cov(b_, how="sample") * ((n_samples - 1) / denominator)
            return [self._when(denominator <= lit(0), lit(None), rescaled)]

        return self._expr(
            func,
            evaluate_output_names=combine_evaluate_output_names(a, b),
            alias_output_names=combine_alias_output_names(a, b),
            version=self._version,
        )

    def struct(self, *exprs: IbisExpr) -> IbisExpr:
        version = self._version

        def func(df: IbisLazyFrame) -> list[ir.Value]:
            names_to_cols: Mapping[str, ir.Value] = {
                alias: native_expr
                for expr in exprs
                for native_expr, _, alias in zip(
                    expr(df),
                    *evaluate_output_names_and_aliases(expr, df, []),
                    strict=True,
                )
            }
            return [ibis.struct(names_to_cols)]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=version,
        )

    def list(self, *exprs: IbisExpr) -> IbisExpr:
        version = self._version

        def func(df: IbisLazyFrame) -> list[ir.Value]:
            cols = [native_expr for expr in exprs for native_expr in expr(df)]
            return [ibis.array(cols)]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=version,
        )
