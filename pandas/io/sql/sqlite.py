"""
sqlite3 backend for ``read_sql`` / ``to_sql`` (used without SQLAlchemy).
"""

from __future__ import annotations

from datetime import (
    date,
    datetime,
    time,
)
from typing import (
    TYPE_CHECKING,
)

from pandas.io.sql.base import get_unicode_name
from pandas.io.sql.dbapi2 import (
    DBAPI2Database,
    DBAPI2Table,
)

if TYPE_CHECKING:
    from collections.abc import (
        Mapping,
    )


# sqlite-specific sql strings and handler class
# dictionary used for readability purposes
_SQL_TYPES = {
    "string": "TEXT",
    "floating": "REAL",
    "integer": "INTEGER",
    "datetime": "TIMESTAMP",
    "date": "DATE",
    "time": "TIME",
    "boolean": "INTEGER",
}


def _get_valid_sqlite_name(name: object) -> str:
    # See https://stackoverflow.com/questions/6514274/how-do-you-escape-strings\
    # -for-sqlite-table-column-names-in-python
    # Ensure the string can be encoded as UTF-8.
    # Ensure the string does not include any NUL characters.
    # Replace all " with "".
    # Wrap the entire thing in double quotes.

    uname = get_unicode_name(name)
    if not len(uname):
        raise ValueError("Empty table or column name specified")

    nul_index = uname.find("\x00")
    if nul_index >= 0:
        raise ValueError("SQLite identifier cannot contain NULs")
    return '"' + uname.replace('"', '""') + '"'


class SQLiteTable(DBAPI2Table):
    """
    DBAPI2Table for sqlite3: registers the sqlite3 date/time adapters.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        self._register_date_adapters()

    def _register_date_adapters(self) -> None:
        # GH 8341
        # register an adapter callable for datetime.time object
        import sqlite3

        # this will transform time(12,34,56,789) into '12:34:56.000789'
        # (this is what sqlalchemy does)
        def _adapt_time(t) -> str:
            # This is faster than strftime
            return f"{t.hour:02d}:{t.minute:02d}:{t.second:02d}.{t.microsecond:06d}"

        # Also register adapters for date/datetime and co
        # xref https://docs.python.org/3.12/library/sqlite3.html#adapter-and-converter-recipes
        # Python 3.12+ doesn't auto-register adapters for us anymore

        adapt_date_iso = lambda val: val.isoformat()
        adapt_datetime_iso = lambda val: val.isoformat(" ")

        sqlite3.register_adapter(time, _adapt_time)

        sqlite3.register_adapter(date, adapt_date_iso)
        sqlite3.register_adapter(datetime, adapt_datetime_iso)

        convert_date = lambda val: date.fromisoformat(val.decode())
        convert_timestamp = lambda val: datetime.fromisoformat(val.decode())

        sqlite3.register_converter("date", convert_date)
        sqlite3.register_converter("timestamp", convert_timestamp)


class SQLiteDatabase(DBAPI2Database):
    """
    DBAPI2Database for sqlite3 connections (fallback without SQLAlchemy).
    This should only be used internally.

    Parameters
    ----------
    con : sqlite connection object

    """

    _table_class = SQLiteTable

    def _driver_error(self) -> type[Exception]:
        from sqlite3 import Error

        return Error

    def _escape_identifier(self, name: object) -> str:
        return _get_valid_sqlite_name(name)

    def _placeholders(self, count: int) -> list[str]:
        return ["?"] * count

    def _sql_type_map(self) -> Mapping[str, str]:
        return _SQL_TYPES

    def has_table(self, name: str, schema: str | None = None) -> bool:
        wld = "?"
        query = f"""
        SELECT
            name
        FROM
            sqlite_master
        WHERE
            type IN ('table', 'view')
            AND name={wld};
        """

        return len(self.execute(query, [name]).fetchall()) > 0
