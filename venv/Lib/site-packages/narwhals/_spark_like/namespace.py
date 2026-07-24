from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, Any

from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
    evaluate_output_names_and_aliases,
)
from narwhals._spark_like.dataframe import SparkLikeLazyFrame
from narwhals._spark_like.expr import SparkLikeExpr
from narwhals._spark_like.selectors import SparkLikeSelectorNamespace
from narwhals._spark_like.utils import (
    import_functions,
    import_native_dtypes,
    narwhals_to_native_dtype,
    true_divide,
)
from narwhals._sql.namespace import SQLNamespace

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Mapping

    from sqlframe.base.column import Column

    from narwhals._compliant.window import WindowInputs
    from narwhals._spark_like.dataframe import SQLFrameDataFrame  # noqa: F401
    from narwhals._utils import Implementation, Version
    from narwhals.typing import ConcatMethod, CorrelationMethod, IntoDType, PythonLiteral

# Adjust slight SQL vs PySpark differences
FUNCTION_REMAPPINGS = {
    "starts_with": "startswith",
    "ends_with": "endswith",
    "trim": "btrim",
    "str_split": "split",
    "regexp_matches": "regexp",
}

_BINARY_OPS = {
    "add": operator.add,
    "subtract": operator.sub,
    "multiply": operator.mul,
    "divide": operator.truediv,
    "and": operator.and_,
}


class SparkLikeNamespace(
    SQLNamespace[SparkLikeLazyFrame, SparkLikeExpr, "SQLFrameDataFrame", "Column"]
):
    def __init__(self, *, version: Version, implementation: Implementation) -> None:
        self._version = version
        self._implementation = implementation

    @property
    def selectors(self) -> SparkLikeSelectorNamespace:
        return SparkLikeSelectorNamespace.from_namespace(self)

    @property
    def _expr(self) -> type[SparkLikeExpr]:
        return SparkLikeExpr

    @property
    def _lazyframe(self) -> type[SparkLikeLazyFrame]:
        return SparkLikeLazyFrame

    @property
    def _F(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            from sqlframe.base import functions

            return functions
        return import_functions(self._implementation)

    @property
    def _native_dtypes(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            from sqlframe.base import types

            return types
        return import_native_dtypes(self._implementation)

    def _function(self, name: str, *args: Column | PythonLiteral) -> Column:
        if name in _BINARY_OPS:
            return _BINARY_OPS[name](*args)
        if name == "isnotnull":
            return args[0].isNotNull()  # type: ignore[union-attr]
        return getattr(self._F, FUNCTION_REMAPPINGS.get(name, name))(*args)

    def _lit(self, value: Any) -> Column:
        return self._F.lit(value)

    def _when(
        self, condition: Column, value: Column, otherwise: Column | None = None
    ) -> Column:
        if otherwise is None:
            return self._F.when(condition, value)
        return self._F.when(condition, value).otherwise(otherwise)

    def _coalesce(self, *exprs: Column) -> Column:
        return self._F.coalesce(*exprs)

    def lit(self, value: PythonLiteral, dtype: IntoDType | None) -> SparkLikeExpr:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            F = df._F

            if isinstance(value, (list, tuple)):
                lit_values = [F.lit(v) for v in value]
                column = F.lit(F.array(lit_values))
            elif isinstance(value, dict):
                if not value:
                    impl = self._implementation
                    if impl.is_sqlframe():
                        msg = f"Cannot create an empty struct type for {self._implementation} backend"
                        raise NotImplementedError(msg)

                    column = (  # pragma: no cover
                        F.from_json(F.lit("{}"), df._native_dtypes.StructType())
                        if impl._backend_version() < (4, 0)
                        else F.struct([])
                    )
                else:
                    lit_values = [F.lit(v).alias(k) for k, v in value.items()]
                    column = F.struct(*lit_values)
            else:
                column = F.lit(value)

            if dtype:
                native_dtype = narwhals_to_native_dtype(
                    dtype, self._version, df._native_dtypes, df.native.sparkSession
                )
                column = column.cast(native_dtype)

            return [column]

        def window_func(
            df: SparkLikeLazyFrame, _window_inputs: WindowInputs[Column]
        ) -> list[Column]:
            return func(df)

        return self._expr(
            func,
            window_func,
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
            implementation=self._implementation,
        )

    def len(self) -> SparkLikeExpr:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            return [df._F.count("*")]

        return self._expr(
            func,
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            version=self._version,
            implementation=self._implementation,
        )

    def mean_horizontal(self, *exprs: SparkLikeExpr) -> SparkLikeExpr:
        def func(cols: Iterable[Column]) -> Column:
            cols = tuple(cols)
            F = exprs[0]._F
            numerator = reduce(
                operator.add, (self._F.coalesce(col, self._F.lit(0)) for col in cols)
            )
            denominator = reduce(
                operator.add,
                (col.isNotNull().cast(self._native_dtypes.IntegerType()) for col in cols),
            )
            return true_divide(F, numerator, denominator)

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def concat(
        self, items: Iterable[SparkLikeLazyFrame], *, how: ConcatMethod
    ) -> SparkLikeLazyFrame:
        dfs = [item._native_frame for item in items]
        if how == "vertical":
            cols_0 = dfs[0].columns
            for i, df in enumerate(dfs[1:], start=1):
                cols_current = df.columns
                if not ((len(cols_current) == len(cols_0)) and (cols_current == cols_0)):
                    msg = (
                        "unable to vstack, column names don't match:\n"
                        f"   - dataframe 0: {cols_0}\n"
                        f"   - dataframe {i}: {cols_current}\n"
                    )
                    raise TypeError(msg)

            return SparkLikeLazyFrame(
                native_dataframe=reduce(lambda x, y: x.union(y), dfs),
                version=self._version,
                implementation=self._implementation,
            )

        if how == "diagonal":
            return SparkLikeLazyFrame(
                native_dataframe=reduce(
                    lambda x, y: x.unionByName(y, allowMissingColumns=True), dfs
                ),
                version=self._version,
                implementation=self._implementation,
            )
        raise NotImplementedError

    def concat_str(
        self, *exprs: SparkLikeExpr, separator: str, ignore_nulls: bool
    ) -> SparkLikeExpr:
        def func(df: SparkLikeLazyFrame) -> list[Column]:
            F = self._F
            cols = tuple(chain.from_iterable(e(df) for e in exprs))
            result = F.concat_ws(separator, *cols)

            if not ignore_nulls:
                null_mask = reduce(operator.or_, (F.isnull(s) for s in cols))
                result = F.when(~null_mask, result).otherwise(F.lit(None))

            return [result]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
            implementation=self._implementation,
        )

    def corr(
        self, a: SparkLikeExpr, b: SparkLikeExpr, *, method: CorrelationMethod
    ) -> SparkLikeExpr:
        if method != "pearson":
            msg = "Only 'pearson' correlation is supported for Spark."
            raise NotImplementedError(msg)

        def func(df: SparkLikeLazyFrame) -> list[Column]:
            F = self._F
            a_ = df._evaluate_single_output_expr(a)
            b_ = df._evaluate_single_output_expr(b)
            return [F.corr(a_, b_)]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(a, b),
            alias_output_names=combine_alias_output_names(a, b),
            version=self._version,
            implementation=self._implementation,
        )

    def cov(self, a: SparkLikeExpr, b: SparkLikeExpr, *, ddof: int) -> SparkLikeExpr:
        F = self._F

        def _cov(
            a_: Column, b_: Column, wrap: Callable[[Column], Column]
        ) -> list[Column]:
            # `wrap` adds the window frame to each aggregate in a window context and
            # is the identity otherwise. A compound expression can't be wrapped as a
            # single window function, so each aggregate is wrapped individually.
            if ddof == 0:
                return [wrap(F.covar_pop(a_, b_))]
            if ddof == 1:
                return [wrap(F.covar_samp(a_, b_))]
            is_valid = a_.isNotNull() & b_.isNotNull()
            n_samples = wrap(F.sum(F.when(is_valid, F.lit(1)).otherwise(F.lit(0))))
            denominator = n_samples - F.lit(ddof)
            rescaled = wrap(F.covar_samp(a_, b_)) * ((n_samples - F.lit(1)) / denominator)
            return [F.when(denominator <= F.lit(0), F.lit(None)).otherwise(rescaled)]

        def func(df: SparkLikeLazyFrame) -> list[Column]:
            a_ = df._evaluate_single_output_expr(a)
            b_ = df._evaluate_single_output_expr(b)
            return _cov(a_, b_, lambda e: e)

        def window_f(
            df: SparkLikeLazyFrame, inputs: WindowInputs[Column]
        ) -> list[Column]:
            assert not inputs.order_by  # noqa: S101
            a_ = df._evaluate_single_output_expr(a)
            b_ = df._evaluate_single_output_expr(b)
            window = df._Window.partitionBy(*(inputs.partition_by or [F.lit(1)]))
            return _cov(a_, b_, lambda e: e.over(window))

        return self._expr(
            call=func,
            window_function=window_f,
            evaluate_output_names=combine_evaluate_output_names(a, b),
            alias_output_names=combine_alias_output_names(a, b),
            version=self._version,
            implementation=self._implementation,
        )

    def struct(self, *exprs: SparkLikeExpr) -> SparkLikeExpr:
        version = self._version

        def func(df: SparkLikeLazyFrame) -> list[Column]:
            F = self._F
            names_to_cols: Mapping[str, Column] = {
                alias: native_expr
                for expr in exprs
                for native_expr, _, alias in zip(
                    expr(df),
                    *evaluate_output_names_and_aliases(expr, df, []),
                    strict=True,
                )
            }
            aliased = (col.alias(name) for name, col in names_to_cols.items())
            return [F.struct(*aliased)]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=version,
            implementation=self._implementation,
        )

    def list(self, *exprs: SparkLikeExpr) -> SparkLikeExpr:
        version = self._version

        def func(df: SparkLikeLazyFrame) -> list[Column]:
            cols = [native_expr for expr in exprs for native_expr in expr(df)]
            return [self._F.array(*cols)]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=version,
            implementation=self._implementation,
        )
