"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""

from pandas.errors import DatabaseError

from pandas.io.sql.adbc import ADBCDatabase
from pandas.io.sql.api import (
    get_schema,
    has_table,
    pandasSQL_builder,
    read_sql,
    read_sql_query,
    read_sql_table,
    table_exists,
    to_sql,
)
from pandas.io.sql.base import (
    PandasSQL,
    SQLTable,
)
from pandas.io.sql.dbapi2 import (
    DBAPI2Database,
    DBAPI2Table,
)
from pandas.io.sql.sqla import (
    BaseEngine,
    SQLAlchemyEngine,
    SQLDatabase,
    get_engine,
)
from pandas.io.sql.sqlite import (
    SQLiteDatabase,
    SQLiteTable,
)

__all__ = [
    "ADBCDatabase",
    "BaseEngine",
    "DBAPI2Database",
    "DBAPI2Table",
    "DatabaseError",
    "PandasSQL",
    "SQLAlchemyEngine",
    "SQLDatabase",
    "SQLTable",
    "SQLiteDatabase",
    "SQLiteTable",
    "get_engine",
    "get_schema",
    "has_table",
    "pandasSQL_builder",
    "read_sql",
    "read_sql_query",
    "read_sql_table",
    "table_exists",
    "to_sql",
]
