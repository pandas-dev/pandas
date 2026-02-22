from __future__ import annotations

from typing import TYPE_CHECKING, Generic, TypeVar

from narwhals._expression_parsing import ExprKind, ExprNode

if TYPE_CHECKING:
    from narwhals.expr import Expr

ExprT = TypeVar("ExprT", bound="Expr")


class ExprCatNamespace(Generic[ExprT]):
    def __init__(self, expr: ExprT) -> None:
        self._expr = expr

    def get_categories(self) -> ExprT:
        """Get unique categories from column.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {"fruits": ["apple", "mango", "mango"]},
            ...     schema={"fruits": pl.Categorical},
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("fruits").cat.get_categories()).to_native()
            shape: (2, 1)
            ┌────────┐
            │ fruits │
            │ ---    │
            │ str    │
            ╞════════╡
            │ apple  │
            │ mango  │
            └────────┘
        """
        return self._expr._append_node(
            ExprNode(ExprKind.FILTRATION, "cat.get_categories")
        )
