from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from narwhals._duckdb.utils import DeferredTimeZone, narwhals_to_native_dtype
from narwhals.dataframe import LazyFrame
from narwhals.translate import from_native
from narwhals.utils import Version

if TYPE_CHECKING:
    from narwhals._compliant.typing import CompliantLazyFrameAny
    from narwhals.typing import IntoSchema

try:
    import duckdb  # ignore-banned-import
except ImportError as _exc:  # pragma: no cover
    msg = (
        "`narwhals.sql` requires DuckDB to be installed.\n\n"
        "Hint: run `pip install -U narwhals[sql]`"
    )
    raise ModuleNotFoundError(msg) from _exc

CONN = duckdb.connect()
TZ = DeferredTimeZone(
    CONN.sql("select value from duckdb_settings() where name = 'TimeZone'")
)


class SQLTable(LazyFrame[duckdb.DuckDBPyRelation]):
    """A LazyFrame with an additional `to_sql` method."""

    def __init__(
        self, df: CompliantLazyFrameAny, level: Literal["full", "interchange", "lazy"]
    ) -> None:
        super().__init__(df, level=level)

    def to_sql(self, *, pretty: bool = False) -> str:
        """Convert to SQL query.

        Arguments:
            pretty: Whether to pretty-print SQL query. If `True`, requires `sqlparse`
                to be installed.

        Examples:
            >>> import narwhals as nw
            >>> from narwhals.sql import table
            >>> schema = {"date": nw.Date, "price": nw.Int64, "symbol": nw.String}
            >>> assets = table("assets", schema)
            >>> result = assets.filter(nw.col("price") > 100)
            >>> print(result.to_sql())
            SELECT * FROM main.assets WHERE (price > 100)
        """
        sql_query = self.to_native().sql_query()
        if not pretty:
            return sql_query
        try:
            import sqlparse
        except ImportError as _exc:  # pragma: no cover
            msg = (
                "`SQLTable.to_sql` with `pretty=True`"
                "requires `sqlparse` to be installed.\n\n"
                "Hint: run `pip install -U narwhals[sql]`"
            )
            raise ModuleNotFoundError(msg) from _exc
        return sqlparse.format(sql_query, reindent=True, keyword_case="upper")  # type: ignore[no-any-return]


def table(name: str, schema: IntoSchema) -> SQLTable:
    """Generate standalone LazyFrame which you can use to generate SQL.

    Note that this requires DuckDB to be installed.

    Parameters:
        name: Table name.
        schema: Table schema.

    Examples:
        >>> import narwhals as nw
        >>> from narwhals.sql import table
        >>> schema = {"date": nw.Date, "price": nw.List(nw.Int64), "symbol": nw.String}
        >>> table("t", schema)
        ┌────────────────────────────┐
        |     Narwhals LazyFrame     |
        |----------------------------|
        |┌──────┬─────────┬─────────┐|
        |│ date │  price  │ symbol  │|
        |│ date │ int64[] │ varchar │|
        |├──────┴─────────┴─────────┤|
        |│          0 rows          │|
        |└──────────────────────────┘|
        └────────────────────────────┘
    """
    column_mapping = {
        col: narwhals_to_native_dtype(dtype, Version.MAIN, TZ)
        for col, dtype in schema.items()
    }
    dtypes = ", ".join(f'"{col}" {dtype}' for col, dtype in column_mapping.items())
    CONN.sql(f"""
        CREATE TABLE "{name}"
        ({dtypes});
        """)
    lf = from_native(CONN.table(name))
    return SQLTable(lf._compliant_frame, level=lf._level)


__all__ = ["table"]
