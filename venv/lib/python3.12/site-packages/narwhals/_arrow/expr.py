from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.series import ArrowSeries
from narwhals._compliant import EagerExpr
from narwhals._expression_parsing import evaluate_nodes, evaluate_output_names_and_aliases
from narwhals._utils import (
    Implementation,
    generate_temporary_column_name,
    not_implemented,
)
from narwhals.functions import col as nw_col

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

    from narwhals._arrow.dataframe import ArrowDataFrame
    from narwhals._arrow.namespace import ArrowNamespace
    from narwhals._compliant.typing import AliasNames, EvalNames, EvalSeries
    from narwhals._utils import Version, _LimitedContext


class ArrowExpr(EagerExpr["ArrowDataFrame", ArrowSeries]):
    _implementation: Implementation = Implementation.PYARROW

    def __init__(
        self,
        call: EvalSeries[ArrowDataFrame, ArrowSeries],
        *,
        evaluate_output_names: EvalNames[ArrowDataFrame],
        alias_output_names: AliasNames | None,
        version: Version,
        implementation: Implementation = Implementation.PYARROW,
    ) -> None:
        self._call = call
        self._evaluate_output_names = evaluate_output_names
        self._alias_output_names = alias_output_names
        self._version = version

    @classmethod
    def from_column_names(
        cls: type[Self],
        evaluate_column_names: EvalNames[ArrowDataFrame],
        /,
        *,
        context: _LimitedContext,
    ) -> Self:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            try:
                return [
                    ArrowSeries(
                        df.native[column_name], name=column_name, version=df._version
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
            version=context._version,
        )

    @classmethod
    def from_column_indices(cls, *column_indices: int, context: _LimitedContext) -> Self:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            tbl = df.native
            cols = df.columns
            return [
                ArrowSeries.from_native(tbl[i], name=cols[i], context=df)
                for i in column_indices
            ]

        return cls(
            func,
            evaluate_output_names=cls._eval_names_indices(column_indices),
            alias_output_names=None,
            version=context._version,
        )

    def __narwhals_namespace__(self) -> ArrowNamespace:
        from narwhals._arrow.namespace import ArrowNamespace

        return ArrowNamespace(version=self._version)

    def _reuse_series_extra_kwargs(
        self, *, returns_scalar: bool = False
    ) -> dict[str, Any]:
        return {"_return_py_scalar": False} if returns_scalar else {}

    def _over_without_partition_by(self, order_by: Sequence[str]) -> Self:
        # e.g. `nw.col('a').cum_sum().order_by(key)`
        # which we can always easily support, as it doesn't require grouping.
        assert order_by  # noqa: S101
        meta = self._metadata

        def func(df: ArrowDataFrame) -> Sequence[ArrowSeries]:
            token = generate_temporary_column_name(8, df.columns)
            df = df.with_row_index(token, order_by=None).sort(
                *order_by, descending=False, nulls_last=False
            )
            results = self(df.drop([token], strict=True))
            if meta is not None and meta.is_scalar_like:
                # We need to broadcast the results to the original size, since
                # `over` is a length-preserving operation.
                size = len(df)
                return [s._with_native(pa.repeat(s.item(), size)) for s in results]

            # TODO(marco): is there a way to do this efficiently without
            # doing 2 sorts? Here we're sorting the dataframe and then
            # again calling `sort_indices`. `ArrowSeries.scatter` would also sort.
            sorting_indices = pc.sort_indices(df.get_column(token).native)
            return [s._with_native(s.native.take(sorting_indices)) for s in results]

        return self.__class__(
            func,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )

    def over(self, partition_by: Sequence[str], order_by: Sequence[str]) -> Self:
        if not partition_by:
            assert order_by  # noqa: S101
            return self._over_without_partition_by(order_by)

        # We have something like prev.leaf().over(...) (e.g. `nw.col('a').sum().over('b')`), where:
        # - `prev` must be elementwise (in the example: `nw.col('a')`)
        # - `leaf` must be an aggregation (in the example: `sum`)
        #
        # We first evaluate `prev` as-is, and then evaluate `leaf().over(...)`` by doing a `group_by`.
        meta = self._metadata
        if partition_by and (
            not meta.current_node.kind.is_scalar_like
            or (meta.prev is not None and not meta.prev.is_elementwise)
        ):
            msg = (
                "Only elementary aggregations are supported for `.over` in PyArrow backend "
                "when `partition_by` is specified.\n\n"
                "Please see: "
                "https://narwhals-dev.github.io/narwhals/concepts/improve_group_by_operation/"
            )
            raise NotImplementedError(msg)

        nodes = list(reversed(list(self._metadata.iter_nodes_reversed())))

        def func(df: ArrowDataFrame) -> Sequence[ArrowSeries]:  # noqa: PLR0914
            plx = self.__narwhals_namespace__()
            if meta.prev is not None:
                df = df.with_columns(cast("ArrowExpr", evaluate_nodes(nodes[:-1], plx)))
                _, aliases = evaluate_output_names_and_aliases(self, df, [])
                leaf_ce = cast(
                    "ArrowExpr",
                    nw_col(*aliases)._append_node(nodes[-1])._to_compliant_expr(plx),
                )
            else:
                _, aliases = evaluate_output_names_and_aliases(self, df, [])
                leaf_ce = self
            if order_by:
                df = df.sort(*order_by, descending=False, nulls_last=False)

            if overlap := set(aliases).intersection(partition_by):
                # E.g. `df.select(nw.all().sum().over('a'))`. This is well-defined,
                # we just don't support it yet.
                msg = (
                    f"Column names {overlap} appear in both expression output names and in `over` keys.\n"
                    "This is not yet supported."
                )
                raise NotImplementedError(msg)

            partition_tbl = df.simple_select(*partition_by)
            has_nulls = [ca.null_count > 0 for ca in partition_tbl.native.columns]
            if not any(has_nulls):
                tmp = df.group_by(partition_by, drop_null_keys=False).agg(leaf_ce)
                tmp = partition_tbl.join(
                    tmp,
                    how="left",
                    left_on=partition_by,
                    right_on=partition_by,
                    suffix="_right",
                )
                return [tmp.get_column(alias) for alias in aliases]

            tbl_native, current_cols = df.native, df.columns
            plx = self.__narwhals_namespace__()

            group_keys: list[str] = []
            encoded_cols: list[ArrowExpr] = []
            for col_name, has_null in zip(partition_tbl.columns, has_nulls):
                if not has_null:
                    group_keys.append(col_name)
                else:
                    tmp_name = generate_temporary_column_name(8, current_cols)

                    group_keys.append(tmp_name)
                    current_cols.append(tmp_name)

                    indices = plx._series.from_native(
                        (
                            tbl_native.column(col_name)
                            .dictionary_encode("encode")
                            .combine_chunks()
                            .indices  # type: ignore[attr-defined]
                        ),
                        context=plx,
                        name=tmp_name,
                    )

                    encoded_cols.append(plx._expr._from_series(indices))

            tbl_encoded = df.with_columns(*encoded_cols)
            windowed = tbl_encoded.group_by(group_keys, drop_null_keys=False).agg(leaf_ce)
            ret = tbl_encoded.simple_select(*group_keys).join(
                windowed,
                left_on=group_keys,
                right_on=group_keys,
                how="inner",
                suffix="_right",
            )
            return [ret.get_column(alias) for alias in aliases]

        return self.__class__(
            func,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )

    ewm_mean = not_implemented()
