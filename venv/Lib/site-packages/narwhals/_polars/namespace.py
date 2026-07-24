from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import polars as pl

from narwhals._polars.expr import PolarsExpr
from narwhals._polars.series import PolarsSeries
from narwhals._polars.utils import (
    BACKEND_VERSION,
    extract_args_kwargs,
    narwhals_to_native_dtype,
)
from narwhals._utils import Implementation, Version, isinstance_or_issubclass, requires
from narwhals.dependencies import is_numpy_array_2d
from narwhals.dtypes import DType

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence
    from datetime import timezone

    from typing_extensions import TypeIs

    from narwhals._compliant import CompliantSelectorNamespace
    from narwhals._polars.dataframe import Method, PolarsDataFrame, PolarsLazyFrame
    from narwhals._polars.typing import FrameT
    from narwhals._utils import _LimitedContext
    from narwhals.typing import Into1DArray, IntoDType, TimeUnit, _2DArray


class PolarsNamespace:
    all: Method[PolarsExpr]
    coalesce: Method[PolarsExpr]
    col: Method[PolarsExpr]
    exclude: Method[PolarsExpr]
    sum_horizontal: Method[PolarsExpr]
    min_horizontal: Method[PolarsExpr]
    max_horizontal: Method[PolarsExpr]
    corr: Method[PolarsExpr]

    _implementation: Implementation = Implementation.POLARS
    _version: Version

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def __getattr__(self, attr: str) -> Any:
        def func(*args: Any, **kwargs: Any) -> Any:
            pos, kwds = extract_args_kwargs(args, kwargs)
            return self._expr(getattr(pl, attr)(*pos, **kwds), version=self._version)

        return func

    @property
    def _dataframe(self) -> type[PolarsDataFrame]:
        from narwhals._polars.dataframe import PolarsDataFrame

        return PolarsDataFrame

    @property
    def _lazyframe(self) -> type[PolarsLazyFrame]:
        from narwhals._polars.dataframe import PolarsLazyFrame

        return PolarsLazyFrame

    @property
    def _expr(self) -> type[PolarsExpr]:
        return PolarsExpr

    @property
    def _series(self) -> type[PolarsSeries]:
        return PolarsSeries

    def is_native(self, obj: Any) -> TypeIs[pl.DataFrame | pl.LazyFrame | pl.Series]:
        return isinstance(obj, (pl.DataFrame, pl.LazyFrame, pl.Series))

    @overload
    def from_native(self, data: pl.DataFrame, /) -> PolarsDataFrame: ...
    @overload
    def from_native(self, data: pl.LazyFrame, /) -> PolarsLazyFrame: ...
    @overload
    def from_native(self, data: pl.Series, /) -> PolarsSeries: ...
    def from_native(
        self, data: pl.DataFrame | pl.LazyFrame | pl.Series | Any, /
    ) -> PolarsDataFrame | PolarsLazyFrame | PolarsSeries:
        if self._dataframe._is_native(data):
            return self._dataframe.from_native(data, context=self)
        if self._series._is_native(data):
            return self._series.from_native(data, context=self)
        if self._lazyframe._is_native(data):
            return self._lazyframe.from_native(data, context=self)
        msg = f"Unsupported type: {type(data).__name__!r}"  # pragma: no cover
        raise TypeError(msg)  # pragma: no cover

    @overload
    def from_numpy(self, data: Into1DArray, /, schema: None = ...) -> PolarsSeries: ...

    @overload
    def from_numpy(
        self, data: _2DArray, /, schema: Mapping[str, IntoDType] | Sequence[str] | None
    ) -> PolarsDataFrame: ...

    def from_numpy(
        self,
        data: Into1DArray | _2DArray,
        /,
        schema: Mapping[str, IntoDType] | Sequence[str] | None = None,
    ) -> PolarsDataFrame | PolarsSeries:
        if is_numpy_array_2d(data):
            return self._dataframe.from_numpy(data, schema=schema, context=self)
        return self._series.from_numpy(data, context=self)  # pragma: no cover

    @requires.backend_version(
        (1, 0, 0), "Please use `col` for columns selection instead."
    )
    def nth(self, indices: Sequence[int]) -> PolarsExpr:
        return self._expr(pl.nth(*indices), version=self._version)

    def len(self) -> PolarsExpr:
        if self._backend_version < (0, 20, 5):
            return self._expr(pl.count().alias("len"), self._version)
        return self._expr(pl.len(), self._version)

    def cov(self, a: PolarsExpr, b: PolarsExpr, *, ddof: int) -> PolarsExpr:
        native = pl.cov(a.native, b.native, ddof=ddof)
        n_samples = (a.native.is_not_null() & b.native.is_not_null()).sum()
        result = pl.when(n_samples > pl.lit(ddof)).then(native)
        return self._expr(result, self._version)

    def all_horizontal(self, *exprs: PolarsExpr, ignore_nulls: bool) -> PolarsExpr:
        it = (expr.fill_null(True) for expr in exprs) if ignore_nulls else iter(exprs)
        return self._expr(pl.all_horizontal(*(expr.native for expr in it)), self._version)

    def any_horizontal(self, *exprs: PolarsExpr, ignore_nulls: bool) -> PolarsExpr:
        it = (expr.fill_null(False) for expr in exprs) if ignore_nulls else iter(exprs)
        return self._expr(pl.any_horizontal(*(expr.native for expr in it)), self._version)

    def concat(
        self,
        items: Iterable[FrameT],
        *,
        how: Literal["vertical", "horizontal", "diagonal"],
    ) -> PolarsDataFrame | PolarsLazyFrame:
        _how = (
            "horizontal_extend"
            if how == "horizontal" and self._backend_version >= (1, 42, 1)
            else how
        )
        result = pl.concat((item.native for item in items), how=_how)
        if isinstance(result, pl.DataFrame):
            return self._dataframe(result, version=self._version)
        return self._lazyframe.from_native(result, context=self)

    def lit(self, value: Any, dtype: IntoDType | None) -> PolarsExpr:
        if dtype is not None:
            return self._expr(
                pl.lit(value, dtype=narwhals_to_native_dtype(dtype, self._version)),
                version=self._version,
            )
        return self._expr(pl.lit(value), version=self._version)

    def mean_horizontal(self, *exprs: PolarsExpr) -> PolarsExpr:
        if self._backend_version < (0, 20, 8):
            return self._expr(
                pl.sum_horizontal(e._native_expr for e in exprs)
                / pl.sum_horizontal(1 - e.is_null()._native_expr for e in exprs),
                version=self._version,
            )

        return self._expr(
            pl.mean_horizontal(e._native_expr for e in exprs), version=self._version
        )

    def concat_str(
        self, *exprs: PolarsExpr, separator: str, ignore_nulls: bool
    ) -> PolarsExpr:
        pl_exprs: list[pl.Expr] = [expr._native_expr for expr in exprs]

        if self._backend_version < (0, 20, 6):
            null_mask = [expr.is_null() for expr in pl_exprs]
            sep = pl.lit(separator)

            if not ignore_nulls:
                null_mask_result = pl.any_horizontal(*null_mask)
                output_expr = pl.reduce(
                    lambda x, y: x.cast(pl.String()) + sep + y.cast(pl.String()),  # type: ignore[arg-type,return-value]
                    pl_exprs,
                )
                result = pl.when(~null_mask_result).then(output_expr)
            else:
                init_value, *values = [
                    pl.when(nm).then(pl.lit("")).otherwise(expr.cast(pl.String()))
                    for expr, nm in zip(pl_exprs, null_mask, strict=True)
                ]
                separators = [
                    pl.when(~nm).then(sep).otherwise(pl.lit("")) for nm in null_mask[:-1]
                ]

                result = pl.fold(
                    acc=init_value,
                    function=operator.add,
                    exprs=[s + v for s, v in zip(separators, values, strict=True)],
                )

            return self._expr(result, version=self._version)

        return self._expr(
            pl.concat_str(pl_exprs, separator=separator, ignore_nulls=ignore_nulls),
            version=self._version,
        )

    def struct(self, *exprs: PolarsExpr) -> PolarsExpr:
        pl_exprs: list[pl.Expr] = [expr._native_expr for expr in exprs]
        return self._expr(pl.struct(pl_exprs), version=self._version)

    def list(self, *exprs: PolarsExpr) -> PolarsExpr:
        pl_exprs: list[pl.Expr] = [expr._native_expr for expr in exprs]

        if self._backend_version < (2,):
            length = pl.count() if self._backend_version < (0, 20, 5) else pl.len()
            to_concat = [e.implode().over(pl.int_range(length)) for e in pl_exprs]
            expr = pl.concat_list(to_concat)
        else:
            expr = pl.list(pl_exprs)  # type: ignore[attr-defined]  # pragma: no cover

        return self._expr(expr, version=self._version)

    def when_then(
        self, when: PolarsExpr, then: PolarsExpr, otherwise: PolarsExpr | None = None
    ) -> PolarsExpr:
        if otherwise is None:
            (when_native, then_native), _ = extract_args_kwargs((when, then), {})
            return self._expr(
                pl.when(when_native).then(then_native), version=self._version
            )
        (when_native, then_native, otherwise_native), _ = extract_args_kwargs(
            (when, then, otherwise), {}
        )
        return self._expr(
            pl.when(when_native).then(then_native).otherwise(otherwise_native),
            version=self._version,
        )

    # NOTE: Implementation is too different to annotate correctly (vs other `*SelectorNamespace`)
    # 1. Others have lots of private stuff for code reuse
    #    i. None of that is useful here
    # 2. We don't have a `PolarsSelector` abstraction, and just use `PolarsExpr`
    @property
    def selectors(self) -> CompliantSelectorNamespace[PolarsDataFrame, PolarsSeries]:
        return cast(
            "CompliantSelectorNamespace[PolarsDataFrame, PolarsSeries]",
            PolarsSelectorNamespace(self),
        )


# NOTE: Polars accepts the bare `pl.Enum` class (as opposed to an instantiated `pl.Enum([...])`)
# in `by_dtype` from 0.20.6 onwards, so selecting "any enum" is not expressible on older versions.
_ENUM_SELECTOR_UNSUPPORTED = (
    "Selecting `Enum` columns is only supported for 'polars>=0.20.6', "
    f"found version {BACKEND_VERSION}."
)


class PolarsSelectorNamespace:
    _implementation = Implementation.POLARS

    def __init__(self, context: _LimitedContext, /) -> None:
        self._version = context._version

    def by_dtype(self, dtypes: Iterable[IntoDType]) -> PolarsExpr:
        native_dtypes = [self._to_native_dtype(dtype) for dtype in dtypes]
        return PolarsExpr(pl.selectors.by_dtype(native_dtypes), version=self._version)

    def _to_native_dtype(self, dtype: IntoDType) -> pl.DataType | type[pl.DataType]:
        # `Enum` can only be turned into a concrete `pl.Enum` instance when its
        # categories are known. They aren't when the bare `Enum` class is passed,
        # nor in `stable.v1` (where categories are never tracked). In those cases
        # we match on the `pl.Enum` class, which selects any enum column.
        version = self._version
        dtypes = version.dtypes
        if isinstance_or_issubclass(dtype, dtypes.Enum):
            if version is not Version.V1 and isinstance(dtype, dtypes.Enum):
                return narwhals_to_native_dtype(dtype, version)
            if BACKEND_VERSION < (0, 20, 6):
                raise NotImplementedError(_ENUM_SELECTOR_UNSUPPORTED)
            return pl.Enum
        if isinstance(dtype, type) and issubclass(dtype, DType):
            return narwhals_to_native_dtype(dtype, version).__class__
        return narwhals_to_native_dtype(dtype, version)

    def matches(self, pattern: str) -> PolarsExpr:
        return PolarsExpr(pl.selectors.matches(pattern=pattern), version=self._version)

    def numeric(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.numeric(), version=self._version)

    def boolean(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.boolean(), version=self._version)

    def string(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.string(), version=self._version)

    def categorical(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.categorical(), version=self._version)

    def enum(self) -> PolarsExpr:
        if BACKEND_VERSION < (0, 20, 6):
            raise NotImplementedError(_ENUM_SELECTOR_UNSUPPORTED)
        expr = (
            pl.selectors.by_dtype(pl.Enum)
            if BACKEND_VERSION < (1, 32)
            else pl.selectors.enum()  # Added in v1.32.0 (pola-rs/polars#23351)
        )
        return PolarsExpr(expr, version=self._version)

    def all(self) -> PolarsExpr:
        return PolarsExpr(pl.selectors.all(), version=self._version)

    def datetime(
        self,
        time_unit: TimeUnit | Iterable[TimeUnit] | None,
        time_zone: str | timezone | Iterable[str | timezone | None] | None,
    ) -> PolarsExpr:
        return PolarsExpr(
            pl.selectors.datetime(time_unit=time_unit, time_zone=time_zone),  # type: ignore[arg-type]
            version=self._version,
        )
