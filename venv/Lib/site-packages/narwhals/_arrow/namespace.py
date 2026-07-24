from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, Literal

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.dataframe import ArrowDataFrame
from narwhals._arrow.expr import ArrowExpr
from narwhals._arrow.selectors import ArrowSelectorNamespace
from narwhals._arrow.series import ArrowSeries
from narwhals._arrow.utils import build_list_array, cast_to_comparable_string_types
from narwhals._compliant import EagerNamespace
from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
)
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from narwhals._arrow.typing import ChunkedArrayAny, Incomplete, ScalarAny
    from narwhals._utils import Version
    from narwhals.typing import CorrelationMethod, IntoDType, PythonLiteral


class ArrowNamespace(
    EagerNamespace[ArrowDataFrame, ArrowSeries, ArrowExpr, pa.Table, "ChunkedArrayAny"]
):
    _implementation = Implementation.PYARROW

    @property
    def _dataframe(self) -> type[ArrowDataFrame]:
        return ArrowDataFrame

    @property
    def _expr(self) -> type[ArrowExpr]:
        return ArrowExpr

    @property
    def _series(self) -> type[ArrowSeries]:
        return ArrowSeries

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def extract_native(
        self, *series: ArrowSeries
    ) -> Iterator[ChunkedArrayAny | ScalarAny]:
        return (s.native[0] if s._broadcast else s.native for s in series)

    def len(self) -> ArrowExpr:
        # coverage bug? this is definitely hit
        return self._expr(  # pragma: no cover
            lambda df: [
                ArrowSeries.from_iterable([len(df.native)], name="len", context=self)
            ],
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            version=self._version,
        )

    def lit(self, value: PythonLiteral, dtype: IntoDType | None) -> ArrowExpr:
        def _lit_arrow_series(_: ArrowDataFrame) -> ArrowSeries:
            arrow_series = ArrowSeries.from_iterable(
                data=[value], name="literal", context=self
            )
            if dtype:
                return arrow_series.cast(dtype)
            return arrow_series

        return self._expr(
            lambda df: [_lit_arrow_series(df)],
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
        )

    def all_horizontal(self, *exprs: ArrowExpr, ignore_nulls: bool) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            series: Iterator[ArrowSeries] = chain.from_iterable(e(df) for e in exprs)
            if ignore_nulls:
                series = (s.fill_null(True, None, None) for s in series)
            return [reduce(operator.and_, series)]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def any_horizontal(self, *exprs: ArrowExpr, ignore_nulls: bool) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            series: Iterator[ArrowSeries] = chain.from_iterable(e(df) for e in exprs)
            if ignore_nulls:
                series = (s.fill_null(False, None, None) for s in series)
            return [reduce(operator.or_, series)]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def sum_horizontal(self, *exprs: ArrowExpr) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            it = chain.from_iterable(expr(df) for expr in exprs)
            series = (s.fill_null(0, strategy=None, limit=None) for s in it)
            return [reduce(operator.add, series)]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def mean_horizontal(self, *exprs: ArrowExpr) -> ArrowExpr:
        int_64 = self._version.dtypes.Int64()

        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            expr_results = tuple(chain.from_iterable(expr(df) for expr in exprs))
            series = [s.fill_null(0, strategy=None, limit=None) for s in expr_results]
            non_na = [1 - s.is_null().cast(int_64) for s in expr_results]
            return [reduce(operator.add, series) / reduce(operator.add, non_na)]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def min_horizontal(self, *exprs: ArrowExpr) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            init_series, *series = tuple(chain.from_iterable(expr(df) for expr in exprs))
            native_series = reduce(
                pc.min_element_wise, [s.native for s in series], init_series.native
            )
            return [
                ArrowSeries(native_series, name=init_series.name, version=self._version)
            ]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def max_horizontal(self, *exprs: ArrowExpr) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            init_series, *series = tuple(chain.from_iterable(expr(df) for expr in exprs))
            native_series = reduce(
                pc.max_element_wise, [s.native for s in series], init_series.native
            )
            return [
                ArrowSeries(native_series, name=init_series.name, version=self._version)
            ]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def _concat_diagonal(self, dfs: Sequence[pa.Table], /) -> pa.Table:
        if self._backend_version >= (14,):
            return pa.concat_tables(dfs, promote_options="default")
        return pa.concat_tables(dfs, promote=True)  # pragma: no cover

    def _concat_horizontal(self, dfs: Sequence[pa.Table], /) -> pa.Table:
        names = list(chain.from_iterable(df.column_names for df in dfs))
        arrays = tuple(chain.from_iterable(df.itercolumns() for df in dfs))
        return pa.Table.from_arrays(arrays, names=names)

    def _concat_vertical(self, dfs: Sequence[pa.Table], /) -> pa.Table:
        cols_0 = dfs[0].column_names
        for i, df in enumerate(dfs[1:], start=1):
            if (cols_current := df.column_names) != cols_0:
                msg = (
                    "unable to vstack, column names don't match:\n"
                    f"   - dataframe 0: {cols_0}\n"
                    f"   - dataframe {i}: {cols_current}\n"
                )
                raise TypeError(msg)
        return pa.concat_tables(dfs)

    @property
    def selectors(self) -> ArrowSelectorNamespace:
        return ArrowSelectorNamespace.from_namespace(self)

    def concat_str(
        self, *exprs: ArrowExpr, separator: str, ignore_nulls: bool
    ) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            series = list(chain.from_iterable(expr(df) for expr in exprs))
            name = series[0].name
            null_handling: Literal["skip", "emit_null"] = (
                "skip" if ignore_nulls else "emit_null"
            )
            it, separator_scalar = cast_to_comparable_string_types(
                *self.extract_native(*series), separator=separator
            )
            # NOTE: stubs indicate `separator` must also be a `ChunkedArray`
            # Reality: `str` is fine
            concat_str: Incomplete = pc.binary_join_element_wise
            compliant = self._series(
                concat_str(*it, separator_scalar, null_handling=null_handling),
                name=name,
                version=self._version,
            )
            return [compliant]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def coalesce(self, *exprs: ArrowExpr) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            align = self._series._align_full_broadcast
            init_series, *series = align(*chain.from_iterable(expr(df) for expr in exprs))
            return [
                ArrowSeries(
                    pc.coalesce(init_series.native, *(s.native for s in series)),
                    name=init_series.name,
                    version=self._version,
                )
            ]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def struct(self, *exprs: ArrowExpr) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            series = tuple(chain.from_iterable(expr(df) for expr in exprs))
            name = series[0].name

            struct_array = pc.make_struct(
                *(s.native if len(s) > 1 else s.native[0] for s in series),
                field_names=tuple(s.name for s in series),
            )
            result = pa.chunked_array([struct_array])
            return [ArrowSeries(result, name=name, version=self._version)]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def list(self, *exprs: ArrowExpr) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            series = tuple(chain.from_iterable(expr(df) for expr in exprs))
            # pa.concat_arrays requires Array, not ChunkedArray
            result = build_list_array([s.native.combine_chunks() for s in series])
            return [ArrowSeries(result, name=series[0].name, version=self._version)]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def _if_then_else(
        self,
        when: ChunkedArrayAny,
        then: ChunkedArrayAny,
        otherwise: ChunkedArrayAny | None = None,
    ) -> ChunkedArrayAny:
        otherwise = pa.nulls(len(when), then.type) if otherwise is None else otherwise
        return pc.if_else(when.fill_null(pa.scalar(False)), then, otherwise)

    def corr(self, a: ArrowExpr, b: ArrowExpr, *, method: CorrelationMethod) -> ArrowExpr:
        if method != "pearson":
            msg = "Only 'pearson' correlation is supported for Pyarrow."
            raise NotImplementedError(msg)

        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            a_series = df._evaluate_single_output_expr(a)
            arr1 = a_series.native
            arr2 = df._evaluate_single_output_expr(b).native
            valid = pc.and_(pc.is_valid(arr1), pc.is_valid(arr2))
            arr1 = pc.filter(arr1, valid)
            arr2 = pc.filter(arr2, valid)

            dev1 = pc.subtract(arr1, pc.mean(arr1))
            dev2 = pc.subtract(arr2, pc.mean(arr2))

            covariance = pc.mean(pc.multiply(dev1, dev2))

            std1 = pc.stddev(arr1)
            std2 = pc.stddev(arr2)

            correlation = pc.divide(covariance, pc.multiply(std1, std2))
            return [
                ArrowSeries.from_iterable(
                    data=[correlation], name=a_series.name, context=self
                )
            ]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(a, b),
            alias_output_names=combine_alias_output_names(a, b),
            context=self,
        )

    def cov(self, a: ArrowExpr, b: ArrowExpr, *, ddof: int) -> ArrowExpr:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            a_series = df._evaluate_single_output_expr(a)
            arr1 = a_series.native
            arr2 = df._evaluate_single_output_expr(b).native
            valid = pc.and_(pc.is_valid(arr1), pc.is_valid(arr2))
            arr1 = pc.filter(arr1, valid)
            arr2 = pc.filter(arr2, valid)
            n = len(arr1)
            if n - ddof <= 0:
                # Matches pandas/Polars. A single valid pair with ddof=0 falls
                # through and yields 0.0, as the deviations are then zero.
                covariance = pa.scalar(None, type=pa.float64())
            else:
                dev1 = pc.subtract(arr1, pc.mean(arr1))
                dev2 = pc.subtract(arr2, pc.mean(arr2))
                numerator = pc.sum(pc.multiply(dev1, dev2))
                covariance = pc.divide(numerator, pa.scalar(n - ddof))
            return [
                ArrowSeries.from_iterable(
                    data=[covariance], name=a_series.name, context=self
                )
            ]

        return self._expr._from_callable(
            func=func,
            evaluate_output_names=combine_evaluate_output_names(a, b),
            alias_output_names=combine_alias_output_names(a, b),
            context=self,
        )
