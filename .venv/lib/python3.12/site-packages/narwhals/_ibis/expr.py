from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, Callable, TypeVar, cast

import ibis

from narwhals._ibis.expr_dt import IbisExprDateTimeNamespace
from narwhals._ibis.expr_list import IbisExprListNamespace
from narwhals._ibis.expr_str import IbisExprStringNamespace
from narwhals._ibis.expr_struct import IbisExprStructNamespace
from narwhals._ibis.utils import (
    IntoColumn,
    asc_nulls_first,
    asc_nulls_last,
    desc_nulls_first,
    desc_nulls_last,
    is_floating,
    lit,
    narwhals_to_native_dtype,
)
from narwhals._sql.expr import SQLExpr
from narwhals._utils import (
    Implementation,
    Version,
    extend_bool,
    no_default,
    not_implemented,
    zip_strict,
)

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    import ibis.expr.types as ir
    from typing_extensions import Self

    from narwhals._compliant import WindowInputs
    from narwhals._compliant.typing import (
        AliasNames,
        EvalNames,
        EvalSeries,
        WindowFunction,
    )
    from narwhals._ibis.dataframe import IbisLazyFrame
    from narwhals._ibis.namespace import IbisNamespace
    from narwhals._typing import NoDefault
    from narwhals._utils import _LimitedContext
    from narwhals.typing import IntoDType, RankMethod, RollingInterpolationMethod

    ExprT = TypeVar("ExprT", bound=ir.Value)
    IbisWindowFunction = WindowFunction[IbisLazyFrame, ir.Value]
    IbisWindowInputs = WindowInputs[ir.Value]


class IbisExpr(SQLExpr["IbisLazyFrame", "ir.Value"]):
    _implementation = Implementation.IBIS

    def __init__(
        self,
        call: EvalSeries[IbisLazyFrame, ir.Value],
        window_function: IbisWindowFunction | None = None,
        *,
        evaluate_output_names: EvalNames[IbisLazyFrame],
        alias_output_names: AliasNames | None,
        version: Version,
        implementation: Implementation = Implementation.IBIS,
    ) -> None:
        self._call = call
        self._evaluate_output_names = evaluate_output_names
        self._alias_output_names = alias_output_names
        self._version = version
        self._window_function: IbisWindowFunction | None = window_function

    @property
    def window_function(self) -> IbisWindowFunction:
        def default_window_func(
            df: IbisLazyFrame, window_inputs: IbisWindowInputs
        ) -> Sequence[ir.Value]:
            return [
                expr.over(
                    ibis.window(
                        group_by=window_inputs.partition_by,
                        order_by=self._sort(*window_inputs.order_by),
                    )
                )
                for expr in self(df)
            ]

        return self._window_function or default_window_func

    def _window_expression(
        self,
        expr: ir.Value,
        partition_by: Sequence[str | ir.Value] = (),
        order_by: Sequence[IntoColumn] = (),
        rows_start: int | None = None,
        rows_end: int | None = None,
        *,
        descending: Sequence[bool] | None = None,
        nulls_last: Sequence[bool] | None = None,
    ) -> ir.Value:
        if rows_start is not None and rows_end is not None:
            rows_between = {"preceding": -rows_start, "following": rows_end}
        elif rows_end is not None:
            rows_between = {"following": rows_end}
        elif rows_start is not None:  # pragma: no cover
            rows_between = {"preceding": -rows_start}
        else:
            rows_between = {}
        desc = descending or False
        last = nulls_last or False
        window = ibis.window(
            group_by=partition_by,
            order_by=self._sort(*order_by, descending=desc, nulls_last=last),
            **rows_between,
        )
        return expr.over(window)

    def _first(self, expr: ir.Value, *order_by: str) -> ir.Value:
        return cast("ir.Column", expr).first(
            order_by=self._sort(*order_by), include_null=True
        )

    def _last(self, expr: ir.Value, *order_by: str) -> ir.Value:
        return cast("ir.Column", expr).last(
            order_by=self._sort(*order_by), include_null=True
        )

    def _any_value(self, expr: ir.Value, *, ignore_nulls: bool) -> ir.Value:
        # !NOTE: ibis arbitrary returns a random non-null value
        # See: https://ibis-project.org/reference/expression-generic.html#ibis.expr.types.generic.Column.arbitrary
        expr = cast("ir.Column", expr)
        return expr.arbitrary() if ignore_nulls else expr.first(include_null=True)

    def __narwhals_namespace__(self) -> IbisNamespace:  # pragma: no cover
        from narwhals._ibis.namespace import IbisNamespace

        return IbisNamespace(version=self._version)

    def broadcast(self) -> Self:
        # Ibis does its own broadcasting.
        return self

    @staticmethod
    def _sort(
        *cols: IntoColumn,
        descending: Sequence[bool] | bool = False,
        nulls_last: Sequence[bool] | bool = False,
    ) -> Iterator[ir.Column]:
        n = len(cols)
        descending = extend_bool(descending, n)
        nulls_last = extend_bool(nulls_last, n)
        mapping = {
            (False, False): asc_nulls_first,
            (False, True): asc_nulls_last,
            (True, False): desc_nulls_first,
            (True, True): desc_nulls_last,
        }
        for col, _desc, _nulls_last in zip_strict(cols, descending, nulls_last):
            yield mapping[(_desc, _nulls_last)](col)

    @classmethod
    def from_column_names(
        cls: type[Self],
        evaluate_column_names: EvalNames[IbisLazyFrame],
        /,
        *,
        context: _LimitedContext,
    ) -> Self:
        def func(df: IbisLazyFrame) -> Sequence[ir.Column]:
            return [df.native[name] for name in evaluate_column_names(df)]

        return cls(
            func,
            evaluate_output_names=evaluate_column_names,
            alias_output_names=None,
            version=context._version,
        )

    @classmethod
    def from_column_indices(cls, *column_indices: int, context: _LimitedContext) -> Self:
        def func(df: IbisLazyFrame) -> Sequence[ir.Column]:
            return [df.native[i] for i in column_indices]

        return cls(
            func,
            evaluate_output_names=cls._eval_names_indices(column_indices),
            alias_output_names=None,
            version=context._version,
        )

    def _with_binary(self, op: Callable[..., ir.Value], other: Self) -> Self:
        return self._with_callable(op, other=other)

    def _with_elementwise(
        self, op: Callable[..., ir.Value], /, **expressifiable_args: Self
    ) -> Self:
        return self._with_callable(op, **expressifiable_args)

    @classmethod
    def _alias_native(cls, expr: ExprT, name: str, /) -> ExprT:
        return cast("ExprT", expr.name(name))

    def __invert__(self) -> Self:
        invert = cast("Callable[..., ir.Value]", operator.invert)
        return self._with_callable(invert)

    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> Self:
        if interpolation != "linear":
            msg = "Only linear interpolation methods are supported for Ibis quantile."
            raise NotImplementedError(msg)
        return self._with_callable(lambda expr: expr.quantile(quantile))

    def n_unique(self) -> Self:
        return self._with_callable(
            lambda expr: expr.nunique() + expr.isnull().any().cast("int8")
        )

    def len(self) -> Self:
        def func(df: IbisLazyFrame) -> Sequence[ir.IntegerScalar]:
            return [df.native.count() for _ in self._evaluate_output_names(df)]

        return self.__class__(
            func,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )

    def null_count(self) -> Self:
        return self._with_callable(lambda expr: expr.isnull().sum())

    def is_nan(self) -> Self:
        def func(expr: ir.FloatingValue) -> ir.Value:
            otherwise = expr.isnan() if is_floating(expr.type()) else False
            return ibis.ifelse(expr.isnull(), None, otherwise)

        return self._with_callable(func)

    def is_finite(self) -> Self:
        def func(expr: ir.IntegerValue | ir.FloatingValue) -> ir.Value:
            if is_floating(expr.type()):
                expr = cast("ir.FloatingValue", expr)
                return ~(expr.isinf() | expr.isnan())
            return ibis.ifelse(expr.isnull(), None, lit(True))

        return self._with_callable(func)

    def is_in(self, other: Sequence[Any]) -> Self:
        return self._with_callable(lambda expr: expr.isin(other))

    def fill_null(self, value: Self | None, strategy: Any, limit: int | None) -> Self:
        # Ibis doesn't yet allow ignoring nulls in first/last with window functions, which makes forward/backward
        # strategies inconsistent when there are nulls present: https://github.com/ibis-project/ibis/issues/9539
        if strategy is not None:
            msg = "`strategy` is not supported for the Ibis backend"
            raise NotImplementedError(msg)
        if limit is not None:
            msg = "`limit` is not supported for the Ibis backend"  # pragma: no cover
            raise NotImplementedError(msg)

        def _fill_null(expr: ir.Value, value: ir.Scalar) -> ir.Value:
            return expr.fill_null(value)

        assert value is not None  # noqa: S101
        return self._with_callable(_fill_null, value=value)

    def cast(self, dtype: IntoDType) -> Self:
        def _func(expr: ir.Column) -> ir.Value:
            native_dtype = narwhals_to_native_dtype(dtype, self._version)
            # ibis `cast` overloads do not include DataType, only literals
            return expr.cast(native_dtype)  # type: ignore[unused-ignore]

        return self._with_callable(_func)

    def is_unique(self) -> Self:
        return self._with_callable(
            lambda expr: expr.isnull().count().over(ibis.window(group_by=(expr))) == 1
        )

    def rank(self, method: RankMethod, *, descending: bool) -> Self:
        def _rank(expr: ir.Column) -> ir.Value:
            order_by = next(self._sort(expr, descending=descending, nulls_last=True))
            window = ibis.window(order_by=order_by)

            if method == "dense":
                rank_ = order_by.dense_rank()
            elif method == "ordinal":
                rank_ = ibis.row_number().over(window)
            else:
                rank_ = order_by.rank()

            # Ibis uses 0-based ranking. Add 1 to match polars 1-based rank.
            rank_ = rank_ + lit(1)

            # For "max" and "average", adjust using the count of rows in the partition.
            if method == "max":
                # Define a window partitioned by expr (i.e. each distinct value)
                partition = ibis.window(group_by=[expr])
                cnt = expr.count().over(partition)
                rank_ = rank_ + cnt - lit(1)
            elif method == "average":
                partition = ibis.window(group_by=[expr])
                cnt = expr.count().over(partition)
                avg = cast("ir.NumericValue", (cnt - lit(1)) / lit(2.0))
                rank_ = rank_ + avg

            return ibis.cases((expr.notnull(), rank_))

        def window_f(df: IbisLazyFrame, inputs: WindowInputs[ir.Value]) -> list[ir.Value]:
            if inputs.order_by:
                msg = "`rank` followed by `over` with `order_by` specified is not supported for Ibis backend."
                raise NotImplementedError(msg)
            return [
                _rank(cast("ir.Column", expr)).over(
                    ibis.window(group_by=inputs.partition_by)
                )
                for expr in self(df)
            ]

        return self._with_callable(_rank, window_f)

    def replace_strict(
        self,
        default: IbisExpr | NoDefault,
        old: Sequence[Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> Self:
        if default is no_default:
            msg = "`replace_strict` requires an explicit value for `default` for ibis backend."
            raise ValueError(msg)
        ns = self.__narwhals_namespace__()

        keys = list(old)
        values = list(new)
        mapping_expr = ibis.map(keys, values)

        def func(df: IbisLazyFrame) -> list[ir.Value]:
            default_col = df._evaluate_single_output_expr(default)

            results = [
                ns._when(expr.isin(keys), mapping_expr[expr], default_col)
                for expr in self(df)
            ]

            if return_dtype:
                native_dtype = narwhals_to_native_dtype(return_dtype, self._version)
                return [res.cast(native_dtype) for res in results]  # pyright: ignore[reportArgumentType, reportCallIssue]
            return results

        return self.__class__(
            func,
            None,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )

    @property
    def str(self) -> IbisExprStringNamespace:
        return IbisExprStringNamespace(self)

    @property
    def dt(self) -> IbisExprDateTimeNamespace:
        return IbisExprDateTimeNamespace(self)

    @property
    def list(self) -> IbisExprListNamespace:
        return IbisExprListNamespace(self)

    @property
    def struct(self) -> IbisExprStructNamespace:
        return IbisExprStructNamespace(self)

    # NOTE: https://github.com/ibis-project/ibis/issues/10542
    cum_prod = not_implemented()

    # NOTE: https://github.com/ibis-project/ibis/issues/11176
    skew = not_implemented()
    kurtosis = not_implemented()

    _count_star = not_implemented()

    # Intentionally not implemented, as Ibis does its own expression rewriting.
    _push_down_window_function = not_implemented()
