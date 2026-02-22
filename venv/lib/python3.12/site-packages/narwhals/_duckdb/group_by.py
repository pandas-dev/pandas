from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING

from narwhals._duckdb.utils import join_column_names
from narwhals._sql.group_by import SQLGroupBy

if TYPE_CHECKING:
    from collections.abc import Sequence

    from duckdb import Expression  # noqa: F401

    from narwhals._duckdb.dataframe import DuckDBLazyFrame
    from narwhals._duckdb.expr import DuckDBExpr


class DuckDBGroupBy(SQLGroupBy["DuckDBLazyFrame", "DuckDBExpr", "Expression"]):
    def __init__(
        self,
        df: DuckDBLazyFrame,
        keys: Sequence[DuckDBExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        frame, self._keys, self._output_key_names = self._parse_keys(df, keys=keys)
        self._compliant_frame = frame.drop_nulls(self._keys) if drop_null_keys else frame

    def agg(self, *exprs: DuckDBExpr) -> DuckDBLazyFrame:
        agg_columns = tuple(self._evaluate_exprs(exprs))
        result = self.compliant.native.aggregate(
            tuple(chain(self._keys, agg_columns)),  # type: ignore[arg-type]
            join_column_names(*self._keys),
        )

        return self.compliant._with_native(result).rename(
            dict(zip(self._keys, self._output_key_names))
        )
