from __future__ import annotations

from typing import TYPE_CHECKING, Generic, TypeVar

from narwhals._expression_parsing import ExprKind, ExprNode

if TYPE_CHECKING:
    from narwhals.expr import Expr
    from narwhals.typing import NonNestedLiteral

ExprT = TypeVar("ExprT", bound="Expr")


class ExprListNamespace(Generic[ExprT]):
    def __init__(self, expr: ExprT) -> None:
        self._expr = expr

    def len(self) -> ExprT:
        """Return the number of elements in each list.

        Null values count towards the total.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 2], [3, 4, None], None, []]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_len=nw.col("a").list.len())
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |shape: (4, 2)           |
            |┌──────────────┬───────┐|
            |│ a            ┆ a_len │|
            |│ ---          ┆ ---   │|
            |│ list[i64]    ┆ u32   │|
            |╞══════════════╪═══════╡|
            |│ [1, 2]       ┆ 2     │|
            |│ [3, 4, null] ┆ 3     │|
            |│ null         ┆ null  │|
            |│ []           ┆ 0     │|
            |└──────────────┴───────┘|
            └────────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.len"))

    def unique(self) -> ExprT:
        """Get the unique/distinct values in the list.

        Null values are included in the result. The order of unique values is not guaranteed.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 1, 2], [3, 3, None], None, []]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_unique=nw.col("a").list.unique())
            ┌────────────────────────────┐
            |     Narwhals DataFrame     |
            |----------------------------|
            |shape: (4, 2)               |
            |┌──────────────┬───────────┐|
            |│ a            ┆ a_unique  │|
            |│ ---          ┆ ---       │|
            |│ list[i64]    ┆ list[i64] │|
            |╞══════════════╪═══════════╡|
            |│ [1, 1, 2]    ┆ [1, 2]    │|
            |│ [3, 3, null] ┆ [null, 3] │|
            |│ null         ┆ null      │|
            |│ []           ┆ []        │|
            |└──────────────┴───────────┘|
            └────────────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.unique"))

    def contains(self, item: NonNestedLiteral) -> ExprT:
        """Check if sublists contain the given item.

        Arguments:
            item: Item that will be checked for membership.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 2], None, []]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_contains_1=nw.col("a").list.contains(1))
            ┌────────────────────────────┐
            |     Narwhals DataFrame     |
            |----------------------------|
            |shape: (3, 2)               |
            |┌───────────┬──────────────┐|
            |│ a         ┆ a_contains_1 │|
            |│ ---       ┆ ---          │|
            |│ list[i64] ┆ bool         │|
            |╞═══════════╪══════════════╡|
            |│ [1, 2]    ┆ true         │|
            |│ null      ┆ null         │|
            |│ []        ┆ false        │|
            |└───────────┴──────────────┘|
            └────────────────────────────┘
        """
        return self._expr._append_node(
            ExprNode(ExprKind.ELEMENTWISE, "list.contains", item=item)
        )

    def get(self, index: int) -> ExprT:
        """Return the value by index in each list.

        Negative indices are not accepted.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1, 2], [3, 4, None], [None, 5]]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_first=nw.col("a").list.get(0))
            ┌──────────────────────────┐
            |    Narwhals DataFrame    |
            |--------------------------|
            |shape: (3, 2)             |
            |┌──────────────┬─────────┐|
            |│ a            ┆ a_first │|
            |│ ---          ┆ ---     │|
            |│ list[i64]    ┆ i64     │|
            |╞══════════════╪═════════╡|
            |│ [1, 2]       ┆ 1       │|
            |│ [3, 4, null] ┆ 3       │|
            |│ [null, 5]    ┆ null    │|
            |└──────────────┴─────────┘|
            └──────────────────────────┘
        """
        if not isinstance(index, int):
            msg = (
                f"Index must be of type 'int'. Got type '{type(index).__name__}' instead."
            )
            raise TypeError(msg)

        if index < 0:
            msg = f"Index {index} is out of bounds: should be greater than or equal to 0."
            raise ValueError(msg)

        return self._expr._append_node(
            ExprNode(ExprKind.ELEMENTWISE, "list.get", index=index)
        )

    def min(self) -> ExprT:
        """Compute the min value of the lists in the array.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql("SELECT * FROM VALUES ([1]), ([3, 4, NULL]) df(a)")
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_min=nw.col("a").list.min())
            ┌────────────────────────┐
            |   Narwhals LazyFrame   |
            |------------------------|
            |┌──────────────┬───────┐|
            |│      a       │ a_min │|
            |│   int32[]    │ int32 │|
            |├──────────────┼───────┤|
            |│ [1]          │     1 │|
            |│ [3, 4, NULL] │     3 │|
            |└──────────────┴───────┘|
            └────────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.min"))

    def max(self) -> ExprT:
        """Compute the max value of the lists in the array.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1], [3, 4, None]]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_max=nw.col("a").list.max())
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |shape: (2, 2)           |
            |┌──────────────┬───────┐|
            |│ a            ┆ a_max │|
            |│ ---          ┆ ---   │|
            |│ list[i64]    ┆ i64   │|
            |╞══════════════╪═══════╡|
            |│ [1]          ┆ 1     │|
            |│ [3, 4, null] ┆ 4     │|
            |└──────────────┴───────┘|
            └────────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.max"))

    def mean(self) -> ExprT:
        """Compute the mean value of the lists in the array.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"a": [[1], [3, 4, None]]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_mean=nw.col("a").list.mean())
            ┌──────────────────────┐
            |  Narwhals DataFrame  |
            |----------------------|
            |pyarrow.Table         |
            |a: list<item: int64>  |
            |  child 0, item: int64|
            |a_mean: double        |
            |----                  |
            |a: [[[1],[3,4,null]]] |
            |a_mean: [[1,3.5]]     |
            └──────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.mean"))

    def median(self) -> ExprT:
        """Compute the median value of the lists in the array.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql("SELECT * FROM VALUES ([1]), ([3, 4, NULL]) df(a)")
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_median=nw.col("a").list.median())
            ┌───────────────────────────┐
            |    Narwhals LazyFrame     |
            |---------------------------|
            |┌──────────────┬──────────┐|
            |│      a       │ a_median │|
            |│   int32[]    │  double  │|
            |├──────────────┼──────────┤|
            |│ [1]          │      1.0 │|
            |│ [3, 4, NULL] │      3.5 │|
            |└──────────────┴──────────┘|
            └───────────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.median"))

    def sum(self) -> ExprT:
        """Compute the sum value of the lists in the array.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [[1], [3, 4, None]]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_sum=nw.col("a").list.sum())
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |shape: (2, 2)           |
            |┌──────────────┬───────┐|
            |│ a            ┆ a_sum │|
            |│ ---          ┆ ---   │|
            |│ list[i64]    ┆ i64   │|
            |╞══════════════╪═══════╡|
            |│ [1]          ┆ 1     │|
            |│ [3, 4, null] ┆ 7     │|
            |└──────────────┴───────┘|
            └────────────────────────┘
        """
        return self._expr._append_node(ExprNode(ExprKind.ELEMENTWISE, "list.sum"))

    def sort(self, *, descending: bool = False, nulls_last: bool = False) -> ExprT:
        """Sort the lists of the expression.

        Arguments:
            descending: Sort in descending order.
            nulls_last: Place null values last.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES ([2, -1, 1]), ([3, -4, NULL]) df(a)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_sorted=nw.col("a").list.sort())
            ┌─────────────────────────────────┐
            |       Narwhals LazyFrame        |
            |---------------------------------|
            |┌───────────────┬───────────────┐|
            |│       a       │   a_sorted    │|
            |│    int32[]    │    int32[]    │|
            |├───────────────┼───────────────┤|
            |│ [2, -1, 1]    │ [-1, 1, 2]    │|
            |│ [3, -4, NULL] │ [NULL, -4, 3] │|
            |└───────────────┴───────────────┘|
            └─────────────────────────────────┘
        """
        return self._expr._append_node(
            ExprNode(
                ExprKind.ELEMENTWISE,
                "list.sort",
                descending=descending,
                nulls_last=nulls_last,
            )
        )
