from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, cast

from narwhals._compliant import EagerExpr
from narwhals._expression_parsing import evaluate_nodes, evaluate_output_names_and_aliases
from narwhals._pandas_like.group_by import _REMAP_ORDERED_INDEX, PandasLikeGroupBy
from narwhals._pandas_like.series import PandasLikeSeries
from narwhals._pandas_like.utils import make_group_by_kwargs
from narwhals._utils import generate_temporary_column_name

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

    from narwhals._compliant.typing import (
        AliasNames,
        EvalNames,
        EvalSeries,
        NarwhalsAggregation,
    )
    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals._pandas_like.namespace import PandasLikeNamespace
    from narwhals._utils import Implementation, Version, _LimitedContext
    from narwhals.typing import PythonLiteral

WINDOW_FUNCTIONS_TO_PANDAS_EQUIVALENT = {
    "cum_sum": "cumsum",
    "cum_min": "cummin",
    "cum_max": "cummax",
    "cum_prod": "cumprod",
    # Pandas cumcount starts counting from 0 while Polars starts from 1
    # Pandas cumcount counts nulls while Polars does not
    # So, instead of using "cumcount" we use "cumsum" on notna() to get the same result
    "cum_count": "cumsum",
    "rolling_sum": "sum",
    "rolling_mean": "mean",
    "rolling_std": "std",
    "rolling_var": "var",
    "shift": "shift",
    "rank": "rank",
    "diff": "diff",
    "fill_null": "fillna",
    "quantile": "quantile",
    "ewm_mean": "mean",
}


def window_kwargs_to_pandas_equivalent(  # noqa: C901
    function_name: str, kwargs: dict[str, Any]
) -> dict[str, PythonLiteral]:
    if function_name == "shift":
        assert "n" in kwargs  # noqa: S101
        pandas_kwargs: dict[str, PythonLiteral] = {"periods": kwargs["n"]}
    elif function_name == "rank":
        assert "method" in kwargs  # noqa: S101
        assert "descending" in kwargs  # noqa: S101
        _method = kwargs["method"]
        pandas_kwargs = {
            "method": "first" if _method == "ordinal" else _method,
            "ascending": not kwargs["descending"],
            "na_option": "keep",
            "pct": False,
        }
    elif function_name.startswith("cum_"):  # Cumulative operation
        pandas_kwargs = {"skipna": True}
    elif function_name == "n_unique":
        pandas_kwargs = {"dropna": False}
    elif function_name.startswith("rolling_"):  # Rolling operation
        assert "min_samples" in kwargs  # noqa: S101
        assert "window_size" in kwargs  # noqa: S101
        assert "center" in kwargs  # noqa: S101
        pandas_kwargs = {
            "min_periods": kwargs["min_samples"],
            "window": kwargs["window_size"],
            "center": kwargs["center"],
        }
    elif function_name in {"std", "var"}:
        assert "ddof" in kwargs  # noqa: S101
        pandas_kwargs = {"ddof": kwargs["ddof"]}
    elif function_name == "fill_null":
        assert "strategy" in kwargs  # noqa: S101
        assert "limit" in kwargs  # noqa: S101
        pandas_kwargs = {"strategy": kwargs["strategy"], "limit": kwargs["limit"]}
    elif function_name == "quantile":
        assert "quantile" in kwargs  # noqa: S101
        assert "interpolation" in kwargs  # noqa: S101
        pandas_kwargs = {
            "q": kwargs["quantile"],
            "interpolation": kwargs["interpolation"],
        }
    elif function_name.startswith("ewm_"):
        assert "com" in kwargs  # noqa: S101
        assert "span" in kwargs  # noqa: S101
        assert "half_life" in kwargs  # noqa: S101
        assert "alpha" in kwargs  # noqa: S101
        assert "adjust" in kwargs  # noqa: S101
        assert "min_samples" in kwargs  # noqa: S101
        assert "ignore_nulls" in kwargs  # noqa: S101

        pandas_kwargs = {
            "com": kwargs["com"],
            "span": kwargs["span"],
            "halflife": kwargs["half_life"],
            "alpha": kwargs["alpha"],
            "adjust": kwargs["adjust"],
            "min_periods": kwargs["min_samples"],
            "ignore_na": kwargs["ignore_nulls"],
        }
    elif function_name in {"first", "last", "any_value"}:
        if kwargs.get("ignore_nulls"):
            msg = (
                "`Expr.any_value(ignore_nulls=True)` is not supported in a `over` "
                "context for pandas-like backend."
            )
            raise NotImplementedError(msg)
        pandas_kwargs = {
            "n": _REMAP_ORDERED_INDEX[cast("NarwhalsAggregation", function_name)]
        }
    else:  # sum, len, ...
        pandas_kwargs = {}
    return pandas_kwargs


class PandasLikeExpr(EagerExpr["PandasLikeDataFrame", PandasLikeSeries]):
    def __init__(
        self,
        call: EvalSeries[PandasLikeDataFrame, PandasLikeSeries],
        *,
        evaluate_output_names: EvalNames[PandasLikeDataFrame],
        alias_output_names: AliasNames | None,
        implementation: Implementation,
        version: Version,
    ) -> None:
        self._call = call
        self._evaluate_output_names = evaluate_output_names
        self._alias_output_names = alias_output_names
        self._implementation = implementation
        self._version = version

    def __narwhals_namespace__(self) -> PandasLikeNamespace:
        from narwhals._pandas_like.namespace import PandasLikeNamespace

        return PandasLikeNamespace(self._implementation, version=self._version)

    @classmethod
    def from_column_names(
        cls: type[Self],
        evaluate_column_names: EvalNames[PandasLikeDataFrame],
        /,
        *,
        context: _LimitedContext,
    ) -> Self:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            try:
                return [
                    PandasLikeSeries(
                        df._native_frame[column_name],
                        implementation=df._implementation,
                        version=df._version,
                    )
                    for column_name in evaluate_column_names(df)
                ]
            except KeyError as e:
                if error := df._check_columns_exist(evaluate_column_names(df)):
                    raise error from e
                raise

        return cls(
            func,
            evaluate_output_names=evaluate_column_names,
            alias_output_names=None,
            implementation=context._implementation,
            version=context._version,
        )

    @classmethod
    def from_column_indices(cls, *column_indices: int, context: _LimitedContext) -> Self:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            native = df.native
            return [
                PandasLikeSeries.from_native(native.iloc[:, i], context=df)
                for i in column_indices
            ]

        return cls(
            func,
            evaluate_output_names=cls._eval_names_indices(column_indices),
            alias_output_names=None,
            implementation=context._implementation,
            version=context._version,
        )

    def ewm_mean(
        self,
        *,
        com: float | None,
        span: float | None,
        half_life: float | None,
        alpha: float | None,
        adjust: bool,
        min_samples: int,
        ignore_nulls: bool,
    ) -> Self:
        return self._reuse_series(
            "ewm_mean",
            com=com,
            span=span,
            half_life=half_life,
            alpha=alpha,
            adjust=adjust,
            min_samples=min_samples,
            ignore_nulls=ignore_nulls,
        )

    def _over_without_partition_by(self, order_by: Sequence[str]) -> Self:
        # e.g. `nw.col('a').cum_sum().order_by(key)`
        # We can always easily support this as it doesn't require grouping.

        def func(df: PandasLikeDataFrame) -> Sequence[PandasLikeSeries]:
            token = generate_temporary_column_name(8, df.columns)
            df = df.with_row_index(token, order_by=None).sort(
                *order_by, descending=False, nulls_last=False
            )
            results = self(df.drop([token], strict=True))
            meta = self._metadata
            if meta is not None and meta.is_scalar_like:
                # We need to broadcast the result to the original size, since
                # `over` is a length-preserving operation.
                index = df.native.index
                ns = self._implementation.to_native_namespace()
                return [
                    s._with_native(ns.Series(s.item(), index=index, name=s.name))
                    for s in results
                ]

            sorting_indices = df.get_column(token)
            for s in results:
                s._scatter_in_place(sorting_indices, s)
            return results

        return self.__class__(
            func,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            implementation=self._implementation,
            version=self._version,
        )

    def over(  # noqa: C901, PLR0915
        self, partition_by: Sequence[str], order_by: Sequence[str]
    ) -> Self:
        if not partition_by:
            assert order_by  # noqa: S101
            return self._over_without_partition_by(order_by)

        # We have something like prev.leaf().over(...) (e.g. `nw.col('a').sum().over('b')`), where:
        # - `prev` must be elementwise (in the example: `nw.col('a')`)
        # - `leaf` must be a "simple" function, i.e. one that pandas supports in `transform`
        #   (in the example: `sum`)
        #
        # We first evaluate `prev` as-is, and then evaluate `leaf().over(...)`` by using `transform`
        # or other DataFrameGroupBy methods.
        meta = self._metadata
        if partition_by and (meta.prev is not None and not meta.prev.is_elementwise):
            msg = (
                "Only elementary expressions are supported for `.over` in pandas-like backends "
                "when `partition_by` is specified.\n\n"
                "Please see: "
                "https://narwhals-dev.github.io/narwhals/concepts/improve_group_by_operation/"
            )
            raise NotImplementedError(msg)
        nodes = list(reversed(list(self._metadata.iter_nodes_reversed())))

        leaf_node = nodes[-1]
        function_name = leaf_node.name
        pandas_agg = PandasLikeGroupBy._REMAP_AGGS.get(
            cast("NarwhalsAggregation", function_name)
        )
        pandas_function_name = WINDOW_FUNCTIONS_TO_PANDAS_EQUIVALENT.get(
            function_name, pandas_agg
        )
        if pandas_function_name is None:
            msg = (
                f"Unsupported function: {function_name} in `over` context.\n\n"
                f"Supported functions are {', '.join(WINDOW_FUNCTIONS_TO_PANDAS_EQUIVALENT)}\n"
                f"and {', '.join(PandasLikeGroupBy._REMAP_AGGS)}."
            )
            raise NotImplementedError(msg)
        scalar_kwargs = leaf_node.kwargs
        pandas_kwargs = window_kwargs_to_pandas_equivalent(function_name, scalar_kwargs)

        def func(df: PandasLikeDataFrame) -> Sequence[PandasLikeSeries]:  # noqa: C901, PLR0912, PLR0914, PLR0915
            assert pandas_function_name is not None  # help mypy  # noqa: S101
            plx = self.__narwhals_namespace__()
            if meta.prev is not None:
                df = df.with_columns(
                    cast("PandasLikeExpr", evaluate_nodes(nodes[:-1], plx))
                )
            _, aliases = evaluate_output_names_and_aliases(self, df, [])
            if function_name == "cum_count":
                df = df.with_columns(~plx.col(*aliases).is_null())

            if function_name.startswith("cum_"):
                assert "reverse" in scalar_kwargs  # noqa: S101
                reverse = scalar_kwargs["reverse"]
            else:
                assert "reverse" not in scalar_kwargs  # noqa: S101
                reverse = False

            if order_by:
                columns = list(set(partition_by).union(aliases).union(order_by))
                token = generate_temporary_column_name(8, columns)
                df = (
                    df.simple_select(*columns)
                    .with_row_index(token, order_by=None)
                    .sort(*order_by, descending=reverse, nulls_last=reverse)
                )
                sorting_indices = df.get_column(token)
            elif reverse:
                columns = list(set(partition_by).union(aliases))
                df = df.simple_select(*columns)._gather_slice(slice(None, None, -1))
            group_by_kwargs = make_group_by_kwargs(drop_null_keys=False)
            grouped = df._native_frame.groupby(partition_by, **group_by_kwargs)
            if function_name.startswith("rolling"):
                rolling = grouped[list(aliases)].rolling(**pandas_kwargs)
                if pandas_function_name in {"std", "var"}:
                    assert "ddof" in scalar_kwargs  # noqa: S101
                    res_native = getattr(rolling, pandas_function_name)(
                        ddof=scalar_kwargs["ddof"]
                    )
                else:
                    res_native = getattr(rolling, pandas_function_name)()
            elif function_name.startswith("ewm"):
                if self._implementation.is_pandas() and (
                    self._implementation._backend_version()
                ) < (1, 2):  # pragma: no cover
                    msg = (
                        "Exponentially weighted calculation is not available in over "
                        f"context for pandas versions older than 1.2.0, found {self._implementation._backend_version()}."
                    )
                    raise NotImplementedError(msg)
                ewm = grouped[list(aliases)].ewm(**pandas_kwargs)
                assert pandas_function_name is not None  # help mypy  # noqa: S101
                res_native = getattr(ewm, pandas_function_name)()
            elif function_name == "fill_null":
                assert "strategy" in scalar_kwargs  # noqa: S101
                assert "limit" in scalar_kwargs  # noqa: S101
                df_grouped = grouped[list(aliases)]
                if scalar_kwargs["strategy"] == "forward":
                    res_native = df_grouped.ffill(limit=scalar_kwargs["limit"])
                elif scalar_kwargs["strategy"] == "backward":
                    res_native = df_grouped.bfill(limit=scalar_kwargs["limit"])
                else:  # pragma: no cover
                    # This is deprecated in pandas. Indeed, `nw.col('a').fill_null(3).over('b')`
                    # does not seem very useful, and DuckDB doesn't support it either.
                    msg = "`fill_null` with `over` without `strategy` specified is not supported."
                    raise NotImplementedError(msg)
            elif function_name == "len":
                if len(aliases) != 1:  # pragma: no cover
                    msg = "Safety check failed, please report a bug."
                    raise AssertionError(msg)
                res_native = grouped.transform("size").to_frame(aliases[0])
            elif function_name in {"first", "last", "any_value"}:
                with warnings.catch_warnings():
                    # Ignore settingwithcopy warnings/errors, they're false-positives here.
                    warnings.filterwarnings("ignore", message="\n.*copy of a slice")
                    _agg = getattr(
                        grouped[[*partition_by, *aliases]], pandas_function_name
                    )(**pandas_kwargs)
                _agg.reset_index(drop=True, inplace=True)
                keys = list(partition_by)
                res_native = df.native[keys].merge(_agg, on=keys)[list(aliases)]
            else:
                res_native = grouped[list(aliases)].transform(
                    pandas_function_name, **pandas_kwargs
                )
            result_frame = df._with_native(res_native)
            results = [result_frame.get_column(name) for name in aliases]
            if order_by:
                with warnings.catch_warnings():
                    # Ignore settingwithcopy warnings/errors, they're false-positives here.
                    warnings.filterwarnings("ignore", message="\n.*copy of a slice")
                    for s in results:
                        s._scatter_in_place(sorting_indices, s)
                    return results
            if reverse:
                return [s._gather_slice(slice(None, None, -1)) for s in results]
            return results

        return self.__class__(
            func,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            implementation=self._implementation,
            version=self._version,
        )
