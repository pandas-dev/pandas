"""
ADBC (Arrow Database Connectivity) backend for ``read_sql`` / ``to_sql``.
"""

from __future__ import annotations

from contextlib import (
    contextmanager,
)
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
)

from pandas._libs import lib
from pandas.compat._optional import (
    import_optional_dependency,
)
from pandas.errors import (
    DatabaseError,
)

from pandas.core.common import maybe_make_list

from pandas.io._util import (
    arrow_table_to_pandas,
    suppress_pyarrow_values_warning,
)
from pandas.io.sql.base import (
    PandasSQL,
    get_unicode_name,
    parse_date_columns,
)

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Iterator,
    )

    from sqlalchemy.sql.expression import (
        Select,
        TextClause,
    )

    from pandas._typing import (
        DtypeArg,
        DtypeBackend,
    )

    from pandas.core.api import (
        DataFrame,
    )


def _wrap_result_adbc(
    df: DataFrame,
    *,
    index_col=None,
    parse_dates=None,
    dtype: DtypeArg | None = None,
    dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
) -> DataFrame:
    """Wrap result set of a SQLAlchemy query in a DataFrame."""
    if dtype:
        df = df.astype(dtype)

    df = parse_date_columns(df, parse_dates)

    if index_col is not None:
        df = df.set_index(index_col)

    return df


def _quote_identifier(name: object) -> str:
    """
    Escape a SQL identifier (table name or schema name) for safe use in
    ADBC-generated SQL statements.

    Uses ANSI SQL double-quote escaping, which is supported by PostgreSQL,
    DuckDB, SQLite, and other databases targeted by ADBC drivers.

    Parameters
    ----------
    name : object
        Identifier to escape.

    Returns
    -------
    str
        The identifier wrapped in double quotes with any internal double
        quotes doubled, e.g. ``my"table`` becomes ``"my""table"``.

    Raises
    ------
    ValueError
        If *name* is empty or contains a NUL character.
    """
    uname = get_unicode_name(name)
    if not len(uname):
        raise ValueError("Empty table or column name specified")

    nul_index = uname.find("\x00")
    if nul_index >= 0:
        raise ValueError("SQL identifier cannot contain NUL characters")

    return '"' + uname.replace('"', '""') + '"'


class ADBCDatabase(PandasSQL):
    """
    This class enables conversion between DataFrame and SQL databases
    using ADBC to handle DataBase abstraction.

    Parameters
    ----------
    con : adbc_driver_manager.dbapi.Connection
    """

    def __init__(self, con) -> None:
        self.con = con

    @contextmanager
    def run_transaction(self):
        with self.con.cursor() as cur:
            try:
                yield cur
            except Exception:
                self.con.rollback()
                raise
            self.con.commit()

    def execute(self, sql: str | Select[Any] | TextClause, params=None):
        from adbc_driver_manager import Error

        if not isinstance(sql, str):
            raise TypeError("Query must be a string unless using sqlalchemy.")
        args = [] if params is None else [params]
        cur = self.con.cursor()
        try:
            cur.execute(sql, *args)
            return cur
        except Error as exc:
            try:
                self.con.rollback()
            except Error as inner_exc:  # pragma: no cover
                ex = DatabaseError(
                    f"Execution failed on sql: {sql}\n{exc}\nunable to rollback"
                )
                raise ex from inner_exc

            ex = DatabaseError(f"Execution failed on sql '{sql}': {exc}")
            raise ex from exc

    def read_table(
        self,
        table_name: str,
        index_col: str | list[str] | None = None,
        coerce_float: bool = True,
        parse_dates=None,
        columns=None,
        schema: str | None = None,
        chunksize: int | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> DataFrame | Iterator[DataFrame]:
        """
        Read SQL database table into a DataFrame.

        Parameters
        ----------
        table_name : str
            Name of SQL table in database.
        coerce_float : bool, default True
            Raises NotImplementedError
        parse_dates : list or dict, default: None
            - List of column names to parse as dates.
            - Dict of ``{column_name: format string}`` where format string is
              strftime compatible in case of parsing string times, or is one of
              (D, s, ns, ms, us) in case of parsing integer timestamps.
            - Dict of ``{column_name: arg}``, where the arg corresponds
              to the keyword arguments of :func:`pandas.to_datetime`.
              Especially useful with databases without native Datetime support,
              such as SQLite.
        columns : list, default: None
            List of column names to select from SQL table.
        schema : string, default None
            Name of SQL schema in database to query (if database flavor
            supports this).  If specified, this overwrites the default
            schema of the SQL database object.
        chunksize : int, default None
            Raises NotImplementedError
        dtype_backend : {'numpy_nullable', 'pyarrow'}
            Back-end data type applied to the resultant :class:`DataFrame`
            (still experimental). If not specified, the default behavior
            is to not use nullable data types. If specified, the behavior
            is as follows:

            * ``"numpy_nullable"``: returns nullable-dtype-backed :class:`DataFrame`
            * ``"pyarrow"``: returns pyarrow-backed nullable
              :class:`ArrowDtype` :class:`DataFrame`

            .. versionadded:: 2.0

        Returns
        -------
        DataFrame

        See Also
        --------
        pandas.read_sql_table
        SQLDatabase.read_query

        """
        if coerce_float is not True:
            raise NotImplementedError(
                "'coerce_float' is not implemented for ADBC drivers"
            )
        if chunksize:
            raise NotImplementedError("'chunksize' is not implemented for ADBC drivers")

        if columns:
            if index_col:
                index_select = maybe_make_list(index_col)
            else:
                index_select = []
            to_select = index_select + columns
            select_list = ", ".join(_quote_identifier(x) for x in to_select)
        else:
            select_list = "*"
        quoted_table_name = _quote_identifier(table_name)
        if schema:
            quoted_schema = _quote_identifier(schema)
            stmt = f"SELECT {select_list} FROM {quoted_schema}.{quoted_table_name}"
        else:
            stmt = f"SELECT {select_list} FROM {quoted_table_name}"

        with self.execute(stmt) as cur:
            pa_table = cur.fetch_arrow_table()
            df = arrow_table_to_pandas(pa_table, dtype_backend=dtype_backend)

        return _wrap_result_adbc(
            df,
            index_col=index_col,
            parse_dates=parse_dates,
        )

    def read_query(
        self,
        sql: str,
        index_col: str | list[str] | None = None,
        coerce_float: bool = True,
        parse_dates=None,
        params=None,
        chunksize: int | None = None,
        dtype: DtypeArg | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> DataFrame | Iterator[DataFrame]:
        """
        Read SQL query into a DataFrame.

        Parameters
        ----------
        sql : str
            SQL query to be executed.
        index_col : string, optional, default: None
            Column name to use as index for the returned DataFrame object.
        coerce_float : bool, default True
            Raises NotImplementedError
        params : list, tuple or dict, optional, default: None
            Raises NotImplementedError
        parse_dates : list or dict, default: None
            - List of column names to parse as dates.
            - Dict of ``{column_name: format string}`` where format string is
              strftime compatible in case of parsing string times, or is one of
              (D, s, ns, ms, us) in case of parsing integer timestamps.
            - Dict of ``{column_name: arg dict}``, where the arg dict
              corresponds to the keyword arguments of
              :func:`pandas.to_datetime` Especially useful with databases
              without native Datetime support, such as SQLite.
        chunksize : int, default None
            Raises NotImplementedError
        dtype : Type name or dict of columns
            Data type for data or columns. E.g. np.float64 or
            {'a': np.float64, 'b': np.int32, 'c': 'Int64'}

        Returns
        -------
        DataFrame

        See Also
        --------
        read_sql_table : Read SQL database table into a DataFrame.
        read_sql

        """
        if coerce_float is not True:
            raise NotImplementedError(
                "'coerce_float' is not implemented for ADBC drivers"
            )
        if params:
            raise NotImplementedError("'params' is not implemented for ADBC drivers")
        if chunksize:
            raise NotImplementedError("'chunksize' is not implemented for ADBC drivers")

        with self.execute(sql) as cur:
            pa_table = cur.fetch_arrow_table()
            df = arrow_table_to_pandas(pa_table, dtype_backend=dtype_backend)

        return _wrap_result_adbc(
            df,
            index_col=index_col,
            parse_dates=parse_dates,
            dtype=dtype,
        )

    read_sql = read_query

    def to_sql(
        self,
        frame,
        name: str,
        if_exists: Literal["fail", "replace", "append", "delete_rows"] = "fail",
        index: bool = True,
        index_label=None,
        schema: str | None = None,
        chunksize: int | None = None,
        dtype: DtypeArg | None = None,
        method: Literal["multi"] | Callable[..., Any] | None = None,
        engine: str = "auto",
        **engine_kwargs,
    ) -> int | None:
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        frame : DataFrame
        name : string
            Name of SQL table.
        if_exists : {'fail', 'replace', 'append'}, default 'fail'
            - fail: If table exists, do nothing.
            - replace: If table exists, drop it, recreate it, and insert data.
            - append: If table exists, insert data. Create if does not exist.
            - delete_rows: If a table exists, delete all records and insert data.
        index : boolean, default True
            Write DataFrame index as a column.
        index_label : string or sequence, default None
            Raises NotImplementedError
        schema : string, default None
            Name of SQL schema in database to write to (if database flavor
            supports this). If specified, this overwrites the default
            schema of the SQLDatabase object.
        chunksize : int, default None
            Raises NotImplementedError
        dtype : single type or dict of column name to SQL type, default None
            Raises NotImplementedError
        method : {None', 'multi', callable}, default None
            Raises NotImplementedError
        engine : {'auto', 'sqlalchemy'}, default 'auto'
            Raises NotImplementedError if not set to 'auto'
        """
        pa = import_optional_dependency("pyarrow")
        from adbc_driver_manager import Error

        if index_label:
            raise NotImplementedError(
                "'index_label' is not implemented for ADBC drivers"
            )
        if chunksize:
            raise NotImplementedError("'chunksize' is not implemented for ADBC drivers")
        if dtype:
            raise NotImplementedError("'dtype' is not implemented for ADBC drivers")
        if method:
            raise NotImplementedError("'method' is not implemented for ADBC drivers")
        if engine != "auto":
            raise NotImplementedError(
                "engine != 'auto' not implemented for ADBC drivers"
            )

        quoted_name = _quote_identifier(name)
        if schema:
            quoted_schema = _quote_identifier(schema)
            quoted_table_name = f"{quoted_schema}.{quoted_name}"
        else:
            quoted_table_name = quoted_name

        # pandas if_exists="append" will still create the
        # table if it does not exist; ADBC is more explicit with append/create
        # as applicable modes, so the semantics get blurred across
        # the libraries
        mode = "create"
        if self.has_table(name, schema):
            if if_exists == "fail":
                raise ValueError(f"Table '{name}' already exists.")
            elif if_exists == "replace":
                sql_statement = f"DROP TABLE {quoted_table_name}"
                self.execute(sql_statement).close()
            elif if_exists == "append":
                mode = "append"
            elif if_exists == "delete_rows":
                mode = "append"
                self.delete_rows(name, schema)

        try:
            with suppress_pyarrow_values_warning():
                tbl = pa.Table.from_pandas(frame, preserve_index=index)
        except pa.ArrowNotImplementedError as exc:
            raise ValueError("datatypes not supported") from exc

        with self.con.cursor() as cur:
            try:
                total_inserted = cur.adbc_ingest(
                    table_name=name, data=tbl, mode=mode, db_schema_name=schema
                )
            except Error as exc:
                raise DatabaseError(
                    f"Failed to insert records on table={name} with {mode=}"
                ) from exc

        self.con.commit()
        # adbc_ingest returns -1 when the driver can't count the rows
        if lib.is_integer(total_inserted) and total_inserted < 0:
            return None
        return total_inserted

    def has_table(self, name: str, schema: str | None = None) -> bool:
        meta = self.con.adbc_get_objects(
            depth="tables", db_schema_filter=schema, table_name_filter=name
        ).read_all()

        for catalog_schema in meta["catalog_db_schemas"].to_pylist():
            if not catalog_schema:
                continue
            for schema_record in catalog_schema:
                if not schema_record:
                    continue

                for table_record in schema_record["db_schema_tables"]:
                    if table_record["table_name"] == name:
                        return True

        return False

    def delete_rows(self, name: str, schema: str | None = None) -> None:
        quoted_table_name = _quote_identifier(name)
        if schema:
            quoted_schema = _quote_identifier(schema)
            quoted_table_name = f"{quoted_schema}.{quoted_table_name}"

        if self.has_table(name, schema):
            self.execute(f"DELETE FROM {quoted_table_name}").close()

    def _create_sql_schema(
        self,
        frame: DataFrame,
        table_name: str,
        keys: list[str] | None = None,
        dtype: DtypeArg | None = None,
        schema: str | None = None,
    ) -> str:
        raise NotImplementedError("not implemented for adbc")
