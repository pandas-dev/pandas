from __future__ import annotations

import collections
from typing import TYPE_CHECKING, Any, ClassVar

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.utils import (
    BACKEND_VERSION,
    cast_to_comparable_string_types,
    extract_py_scalar,
)
from narwhals._compliant import EagerGroupBy
from narwhals._expression_parsing import evaluate_output_names_and_aliases
from narwhals._utils import generate_temporary_column_name, requires

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping, Sequence

    from narwhals._arrow.dataframe import ArrowDataFrame
    from narwhals._arrow.expr import ArrowExpr
    from narwhals._arrow.typing import (  # type: ignore[attr-defined]
        AggregateOptions,
        Aggregation,
        Incomplete,
    )
    from narwhals._compliant.typing import NarwhalsAggregation
    from narwhals.typing import UniqueKeepStrategy


class ArrowGroupBy(EagerGroupBy["ArrowDataFrame", "ArrowExpr", "Aggregation"]):
    _REMAP_AGGS: ClassVar[Mapping[NarwhalsAggregation, Aggregation]] = {
        "sum": "sum",
        "mean": "mean",
        "median": "approximate_median",
        "max": "max",
        "min": "min",
        "std": "stddev",
        "var": "variance",
        "len": "count",
        "n_unique": "count_distinct",
        "count": "count",
        "all": "all",
        "any": "any",
        "first": "first",
        "last": "last",
        "any_value": "first",
    }
    _REMAP_UNIQUE: ClassVar[Mapping[UniqueKeepStrategy, Aggregation]] = {
        "any": "min",
        "first": "min",
        "last": "max",
    }
    _OPTION_COUNT_ALL: ClassVar[frozenset[NarwhalsAggregation]] = frozenset(
        ("len", "n_unique")
    )
    _OPTION_COUNT_VALID: ClassVar[frozenset[NarwhalsAggregation]] = frozenset(("count",))
    _OPTION_ORDERED: ClassVar[frozenset[NarwhalsAggregation]] = frozenset(
        ("first", "last", "any_value")
    )
    _OPTION_VARIANCE: ClassVar[frozenset[NarwhalsAggregation]] = frozenset(("std", "var"))
    _OPTION_SCALAR: ClassVar[frozenset[NarwhalsAggregation]] = frozenset(
        ("any", "all", "sum")
    )

    def __init__(
        self,
        df: ArrowDataFrame,
        keys: Sequence[ArrowExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._df = df
        frame, self._keys, self._output_key_names = self._parse_keys(df, keys=keys)
        self._compliant_frame = frame.drop_nulls(self._keys) if drop_null_keys else frame
        self._drop_null_keys = drop_null_keys

    def _configure_agg(
        self, expr: ArrowExpr, /
    ) -> tuple[Aggregation, AggregateOptions | None]:
        option: AggregateOptions | None = None
        function_name = self._leaf_name(expr)
        kwargs = self._kwargs(expr)
        if function_name in self._OPTION_VARIANCE:
            ddof = kwargs["ddof"]
            option = pc.VarianceOptions(ddof=ddof)
        elif function_name in self._OPTION_COUNT_ALL:
            option = pc.CountOptions(mode="all")
        elif function_name in self._OPTION_COUNT_VALID:
            option = pc.CountOptions(mode="only_valid")
        elif function_name in self._OPTION_SCALAR:
            option = pc.ScalarAggregateOptions(min_count=0)
        elif function_name in self._OPTION_ORDERED:
            ignore_nulls = kwargs.get("ignore_nulls", False)
            option = pc.ScalarAggregateOptions(skip_nulls=ignore_nulls)
        return self._remap_expr_name(function_name), option

    def _configure_grouped(self, *exprs: ArrowExpr) -> pa.TableGroupBy:
        order_by = ()
        use_threads = True
        for expr in exprs:
            md = next(expr._metadata.op_nodes_reversed())
            if md.name not in self._OPTION_ORDERED:
                continue
            # [pyarrow-36709]: https://github.com/apache/arrow/issues/36709
            use_threads = False
            if _current_order_by := md.kwargs.get("order_by", ()):
                if order_by and _current_order_by != order_by:
                    msg = f"Only one `order_by` can be specified in `group_by`. Found both {order_by} and {_current_order_by}."
                    raise NotImplementedError(msg)
                order_by = _current_order_by
        if not use_threads and BACKEND_VERSION < (14,):  # pragma: no cover
            msg = (
                f"Using `first/last` in a `group_by().agg(...)` context is only available in 'pyarrow>=14.0.0', "
                f"found version {requires._unparse_version(BACKEND_VERSION)!r}.\n\n"
                f"See https://github.com/apache/arrow/issues/36709"
            )
            raise NotImplementedError(msg)
        if order_by:
            return pa.TableGroupBy(
                self.compliant.sort(*order_by, descending=False, nulls_last=False).native,
                self._keys,
                use_threads=use_threads,
            )
        if not use_threads:
            return pa.TableGroupBy(self.compliant.native, self._keys, use_threads=False)
        # TODO(unassigned): combine with `return` above once PyArrow 15 is the minimum.
        return pa.TableGroupBy(self.compliant.native, self._keys)

    def agg(self, *exprs: ArrowExpr) -> ArrowDataFrame:
        self._ensure_all_simple(exprs)
        grouped = self._configure_grouped(*exprs)

        aggs: list[tuple[str, Aggregation, AggregateOptions | None]] = []
        expected_pyarrow_column_names: list[str] = self._keys.copy()
        new_column_names: list[str] = self._keys.copy()
        exclude = (*self._keys, *self._output_key_names)

        for expr in exprs:
            output_names, aliases = evaluate_output_names_and_aliases(
                expr, self.compliant, exclude
            )
            md = expr._metadata
            op_nodes_reversed = list(md.op_nodes_reversed())
            if len(op_nodes_reversed) == 1:
                # e.g. `agg(nw.len())`
                if op_nodes_reversed[0].name != "len":  # pragma: no cover
                    msg = "Safety assertion failed, please report a bug to https://github.com/narwhals-dev/narwhals/issues"
                    raise AssertionError(msg)

                new_column_names.append(aliases[0])
                expected_pyarrow_column_names.append(f"{self._keys[0]}_count")
                aggs.append((self._keys[0], "count", pc.CountOptions(mode="all")))
                continue

            function_name, option = self._configure_agg(expr)
            new_column_names.extend(aliases)
            expected_pyarrow_column_names.extend(
                [f"{output_name}_{function_name}" for output_name in output_names]
            )
            aggs.extend(
                [(output_name, function_name, option) for output_name in output_names]
            )

        result_simple = grouped.aggregate(aggs)

        # Rename columns, being very careful
        expected_old_names_indices: dict[str, list[int]] = collections.defaultdict(list)
        for idx, item in enumerate(expected_pyarrow_column_names):
            expected_old_names_indices[item].append(idx)
        if not (
            set(result_simple.column_names) == set(expected_pyarrow_column_names)
            and len(result_simple.column_names) == len(expected_pyarrow_column_names)
        ):  # pragma: no cover
            msg = (
                f"Safety assertion failed, expected {expected_pyarrow_column_names} "
                f"got {result_simple.column_names}, "
                "please report a bug at https://github.com/narwhals-dev/narwhals/issues"
            )
            raise AssertionError(msg)
        index_map: list[int] = [
            expected_old_names_indices[item].pop(0) for item in result_simple.column_names
        ]
        new_column_names = [new_column_names[i] for i in index_map]
        result_simple = result_simple.rename_columns(new_column_names)
        return self.compliant._with_native(result_simple).rename(
            dict(zip(self._keys, self._output_key_names))
        )

    def __iter__(self) -> Iterator[tuple[Any, ArrowDataFrame]]:
        col_token = generate_temporary_column_name(
            n_bytes=8, columns=self.compliant.columns
        )
        null_token: str = "__null_token_value__"  # noqa: S105

        table = self.compliant.native
        it, separator_scalar = cast_to_comparable_string_types(
            *(table[key] for key in self._keys), separator=""
        )
        # NOTE: stubs indicate `separator` must also be a `ChunkedArray`
        # Reality: `str` is fine
        concat_str: Incomplete = pc.binary_join_element_wise
        key_values = concat_str(
            *it, separator_scalar, null_handling="replace", null_replacement=null_token
        )
        table = table.add_column(i=0, field_=col_token, column=key_values)

        for v in pc.unique(key_values):
            t = self.compliant._with_native(
                table.filter(pc.equal(table[col_token], v)).drop([col_token])
            )
            row = t.simple_select(*self._keys).row(0)
            yield (
                tuple(extract_py_scalar(el) for el in row),
                t.simple_select(*self._df.columns),
            )
