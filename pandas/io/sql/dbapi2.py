"""
Generic DB-API 2.0 (PEP 249) backend for ``read_sql`` / ``to_sql``.

Builds SQL statements as strings and runs them through the connection's
cursor, without SQLAlchemy. Driver- and dialect-specific details are left to
subclasses (see ``pandas.io.sql.sqlite``).
"""

from __future__ import annotations

from abc import abstractmethod
from contextlib import (
    contextmanager,
)
import re
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    cast,
)
import warnings

from pandas._libs import lib
from pandas.errors import (
    DatabaseError,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_dict_like,
    is_list_like,
)

from pandas.core.api import (
    DataFrame,
)

from pandas.io.sql.base import (
    PandasSQL,
    SQLTable,
    wrap_result,
)

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Generator,
        Hashable,
        Iterator,
        Mapping,
    )

    from sqlalchemy.sql.expression import (
        Select,
        TextClause,
    )

    from pandas._typing import (
        DtypeArg,
        DtypeBackend,
    )


class DBAPI2Table(SQLTable):
    """
    SQLTable for DB-API 2.0 connections, without SQLAlchemy.

    Instead of a SQLAlchemy ``Table``, ``self.table`` holds the
    ``CREATE TABLE`` / ``CREATE INDEX`` statements as strings. Identifier
    quoting, placeholders and type names come from the owning
    ``DBAPI2Database``.
    """

    def sql_schema(self) -> str:
        return str(";\n".join(self.table))

    def _execute_create(self) -> None:
        with self.pd_sql.run_transaction() as cur:
            for stmt in self.table:
                cur.execute(stmt)

    def insert_statement(self, *, num_rows: int) -> str:
        names = list(map(str, self.frame.columns))
        escape = self.pd_sql._escape_identifier

        if self.index is not None:
            for idx in self.index[::-1]:
                names.insert(0, idx)

        bracketed_names = [escape(column) for column in names]
        col_names = ",".join(bracketed_names)

        num_cols = len(names)
        placeholders = self.pd_sql._placeholders(num_rows * num_cols)
        wildcards = ",".join(
            "(" + ",".join(placeholders[row * num_cols : (row + 1) * num_cols]) + ")"
            for row in range(num_rows)
        )
        insert_statement = (
            f"INSERT INTO {escape(self.name)} ({col_names}) VALUES {wildcards}"
        )
        return insert_statement

    def _execute_insert(self, conn, keys, data_iter) -> int:
        driver_error = self.pd_sql._driver_error()

        data_list = list(data_iter)
        try:
            conn.executemany(self.insert_statement(num_rows=1), data_list)
        except driver_error as exc:
            raise DatabaseError("Execution failed") from exc
        return conn.rowcount

    def _execute_insert_multi(self, conn, keys, data_iter) -> int:
        data_list = list(data_iter)
        flattened_data = [x for row in data_list for x in row]
        conn.execute(self.insert_statement(num_rows=len(data_list)), flattened_data)
        return conn.rowcount

    def _create_table_setup(self):
        """
        Return a list of SQL statements that creates a table reflecting the
        structure of a DataFrame.  The first entry will be a CREATE TABLE
        statement while the rest will be CREATE INDEX statements.
        """
        column_names_and_types = self._get_column_names_and_types(self._sql_type_name)
        escape = self.pd_sql._escape_identifier

        create_tbl_stmts = [
            escape(cname) + " " + ctype for cname, ctype, _ in column_names_and_types
        ]

        if self.keys is not None and len(self.keys):
            if not is_list_like(self.keys):
                keys = [self.keys]
            else:
                keys = self.keys
            cnames_br = ", ".join([escape(c) for c in keys])
            create_tbl_stmts.append(
                f"CONSTRAINT {self.name}_pk PRIMARY KEY ({cnames_br})"
            )
        if self.schema:
            schema_name = self.schema + "."
        else:
            schema_name = ""
        create_stmts = [
            "CREATE TABLE "
            + schema_name
            + escape(self.name)
            + " (\n"
            + ",\n  ".join(create_tbl_stmts)
            + "\n)"
        ]

        ix_cols = [cname for cname, _, is_index in column_names_and_types if is_index]
        if ix_cols:
            cnames = "_".join(ix_cols)
            cnames_br = ",".join([escape(c) for c in ix_cols])
            create_stmts.append(
                "CREATE INDEX "
                + escape("ix_" + self.name + "_" + cnames)
                + "ON "
                + escape(self.name)
                + " ("
                + cnames_br
                + ")"
            )

        return create_stmts

    def _sql_type_name(self, col):
        dtype: DtypeArg = self.dtype or {}
        if is_dict_like(dtype):
            dtype = cast("Mapping[Hashable, Any]", dtype)
            if col.name in dtype:
                return dtype[col.name]

        # Infer type of column, while ignoring missing values.
        # Needed for inserting typed data containing NULLs, GH 8778.
        col_type = lib.infer_dtype(col, skipna=True)

        if col_type == "timedelta64":
            warnings.warn(
                "the 'timedelta' type is not supported, and will be "
                "written as integer values (ns frequency) to the database.",
                UserWarning,
                stacklevel=find_stack_level(),
            )
            col_type = "integer"

        elif col_type == "datetime64":
            col_type = "datetime"

        elif col_type == "empty":
            col_type = "string"

        elif col_type == "complex":
            raise ValueError("Complex datatypes not supported")

        sql_types = self.pd_sql._sql_type_map()
        if col_type not in sql_types:
            col_type = "string"

        return sql_types[col_type]


class DBAPI2Database(PandasSQL):
    """
    PandasSQL for DB-API 2.0 connections, without SQLAlchemy.

    Subclasses provide the driver- and dialect-specific pieces: the driver's
    exception class, identifier quoting, the parameter placeholder, SQL type
    names and ``has_table``.

    Parameters
    ----------
    con : DB-API 2.0 connection object
    """

    _table_class: type[DBAPI2Table] = DBAPI2Table

    def __init__(self, con) -> None:
        self.con = con

    @abstractmethod
    def _driver_error(self) -> type[Exception]:
        """
        Return the driver's base exception class (PEP 249 ``Error``).
        """

    @abstractmethod
    def _escape_identifier(self, name: object) -> str:
        """
        Quote a table or column name for use in generated SQL.
        """

    @abstractmethod
    def _placeholders(self, count: int) -> list[str]:
        """
        Return the parameter placeholders for ``count`` positional values.
        """

    @abstractmethod
    def _sql_type_map(self) -> Mapping[str, str]:
        """
        Map ``lib.infer_dtype`` results to SQL column type names.
        """

    @contextmanager
    def run_transaction(self):
        cur = self.con.cursor()
        try:
            yield cur
            self.con.commit()
        except Exception:
            self.con.rollback()
            raise
        finally:
            cur.close()

    def execute(self, sql: str | Select[Any] | TextClause, params=None):
        driver_error = self._driver_error()

        if not isinstance(sql, str):
            raise TypeError("Query must be a string unless using sqlalchemy.")
        args = [] if params is None else [params]
        cur = self.con.cursor()
        try:
            cur.execute(sql, *args)
            return cur
        except driver_error as exc:
            try:
                self.con.rollback()
            except driver_error as inner_exc:  # pragma: no cover
                ex = DatabaseError(
                    f"Execution failed on sql: {sql}\n{exc}\nunable to rollback"
                )
                raise ex from inner_exc

            if re.fullmatch(r"[A-Za-z_]\w*", sql.strip()):
                # GH#54233 a bare identifier is almost certainly an attempt to
                #  read a table by name.
                ex = DatabaseError(
                    f"Execution failed on sql '{sql}': {exc}\n\n"
                    "Reading a table by name is only supported when using a "
                    "SQLAlchemy connection. With a DBAPI connection, pass a SQL "
                    f"query instead, e.g. 'SELECT * FROM {sql.strip()}'."
                )
                raise ex from exc

            ex = DatabaseError(f"Execution failed on sql '{sql}': {exc}")
            raise ex from exc

    @staticmethod
    def _query_iterator(
        cursor,
        chunksize: int,
        columns,
        index_col=None,
        coerce_float: bool = True,
        parse_dates=None,
        dtype: DtypeArg | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> Generator[DataFrame]:
        """Return generator through chunked result set"""
        has_read_data = False
        while True:
            data = cursor.fetchmany(chunksize)
            if type(data) == tuple:
                data = list(data)
            if not data:
                cursor.close()
                if not has_read_data:
                    result = DataFrame.from_records(
                        [], columns=columns, coerce_float=coerce_float
                    )
                    if dtype:
                        result = result.astype(dtype)
                    yield result
                break

            has_read_data = True
            yield wrap_result(
                data,
                columns,
                index_col=index_col,
                coerce_float=coerce_float,
                parse_dates=parse_dates,
                dtype=dtype,
                dtype_backend=dtype_backend,
            )

    def read_query(
        self,
        sql,
        index_col=None,
        coerce_float: bool = True,
        parse_dates=None,
        params=None,
        chunksize: int | None = None,
        dtype: DtypeArg | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> DataFrame | Iterator[DataFrame]:
        cursor = self.execute(sql, params)
        columns = [col_desc[0] for col_desc in cursor.description]

        if chunksize is not None:
            return self._query_iterator(
                cursor,
                chunksize,
                columns,
                index_col=index_col,
                coerce_float=coerce_float,
                parse_dates=parse_dates,
                dtype=dtype,
                dtype_backend=dtype_backend,
            )
        else:
            data = self._fetchall_as_list(cursor)
            cursor.close()

            frame = wrap_result(
                data,
                columns,
                index_col=index_col,
                coerce_float=coerce_float,
                parse_dates=parse_dates,
                dtype=dtype,
                dtype_backend=dtype_backend,
            )
            return frame

    def _fetchall_as_list(self, cur):
        result = cur.fetchall()
        if not isinstance(result, list):
            result = list(result)
        return result

    def to_sql(
        self,
        frame,
        name: str,
        if_exists: Literal["fail", "replace", "append", "delete_rows"] = "fail",
        index: bool = True,
        index_label=None,
        schema=None,
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
        frame: DataFrame
        name: string
            Name of SQL table.
        if_exists: {'fail', 'replace', 'append', 'delete_rows'}, default 'fail'
            fail: If table exists, do nothing.
            replace: If table exists, drop it, recreate it, and insert data.
            append: If table exists, insert data. Create if it does not exist.
            delete_rows: If a table exists, delete all records and insert data.
        index : bool, default True
            Write DataFrame index as a column
        index_label : string or sequence, default None
            Column label for index column(s). If None is given (default) and
            `index` is True, then the index names are used.
            A sequence should be given if the DataFrame uses MultiIndex.
        schema : string, default None
            Ignored parameter included for compatibility with SQLAlchemy
            version of ``to_sql``.
        chunksize : int, default None
            If not None, then rows will be written in batches of this
            size at a time. If None, all rows will be written at once.
        dtype : single type or dict of column name to SQL type, default None
            Optional specifying the datatype for columns. The SQL type should
            be a string. If all columns are of the same type, one single value
            can be used.
        method : {None, 'multi', callable}, default None
            Controls the SQL insertion clause used:

            * None : Uses standard SQL ``INSERT`` clause (one per row).
            * 'multi': Pass multiple values in a single ``INSERT`` clause.
            * callable with signature ``(pd_table, conn, keys, data_iter)``.

            Details and a sample callable implementation can be found in the
            section :ref:`insert method <io.sql.method>`.
        """
        if dtype:
            if not is_dict_like(dtype):
                # error: Value expression in dictionary comprehension has incompatible
                # type "Union[ExtensionDtype, str, dtype[Any], Type[object],
                # Dict[Hashable, Union[ExtensionDtype, Union[str, dtype[Any]],
                # Type[str], Type[float], Type[int], Type[complex], Type[bool],
                # Type[object]]]]"; expected type "Union[ExtensionDtype, str,
                # dtype[Any], Type[object]]"
                dtype = dict.fromkeys(frame, dtype)  # type: ignore[arg-type]
            else:
                dtype = cast("Mapping[Hashable, Any]", dtype)

            for col, my_type in dtype.items():
                if not isinstance(my_type, str):
                    raise ValueError(
                        f"Invalid type '{my_type}' for dtype of column '{col}': "
                        "expected a string. When using a DB-API connection "
                        "(e.g. sqlite3), dtype values must be SQL type "
                        "strings (e.g. 'TEXT', 'FLOAT')."
                    )

        table = self._table_class(
            name,
            self,
            frame=frame,
            index=index,
            if_exists=if_exists,
            index_label=index_label,
            dtype=dtype,
        )
        table.create()
        return table.insert(chunksize, method)

    def get_table(self, table_name: str, schema: str | None = None) -> None:
        return None  # not supported in fallback mode

    def drop_table(self, name: str, schema: str | None = None) -> None:
        drop_sql = f"DROP TABLE {self._escape_identifier(name)}"
        self.execute(drop_sql).close()

    def delete_rows(self, name: str, schema: str | None = None) -> None:
        delete_sql = f"DELETE FROM {self._escape_identifier(name)}"
        if self.has_table(name, schema):
            self.execute(delete_sql).close()

    def _create_sql_schema(
        self,
        frame,
        table_name: str,
        keys=None,
        dtype: DtypeArg | None = None,
        schema: str | None = None,
    ) -> str:
        table = self._table_class(
            table_name,
            self,
            frame=frame,
            index=False,
            keys=keys,
            dtype=dtype,
            schema=schema,
        )
        return str(table.sql_schema())
