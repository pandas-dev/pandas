# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

"""PEP 249 (DB-API 2.0) API wrapper for the ADBC Driver Manager.

Resource Management
===================

You must ``close()`` Connection and Cursor objects, or else driver
resources may be leaked.  A ``__del__`` is implemented as a fallback,
but Python does not guarantee the timing of when this is called.  For
development, ``__del__`` will raise a ResourceWarning when running
under pytest, or when the environment variable
``_ADBC_DRIVER_MANAGER_WARN_UNCLOSED_RESOURCE`` is set to ``1``.

"""

import abc
import datetime
import os
import threading
import time
import typing
import warnings
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

try:
    import pyarrow
except ImportError as e:
    raise ImportError("PyArrow is required for the DBAPI-compatible interface") from e

try:
    import pyarrow.dataset
except ImportError:
    _pya_dataset = ()
    _pya_scanner = ()
else:
    _pya_dataset = (pyarrow.dataset.Dataset,)
    _pya_scanner = (pyarrow.dataset.Scanner,)

import adbc_driver_manager

from . import _lib, _reader
from ._lib import _blocking_call

if typing.TYPE_CHECKING:
    import pandas
    from typing_extensions import Self

# ----------------------------------------------------------
# Globals

#: The DB-API API level (2.0).
apilevel = "2.0"
#: The thread safety level (connections may not be shared).
threadsafety = 1
#: The parameter style (qmark). This is hardcoded, but actually
#: depends on the driver.
paramstyle = "qmark"

Warning = _lib.Warning
Error = _lib.Error
InterfaceError = _lib.InterfaceError
DatabaseError = _lib.DatabaseError
DataError = _lib.DataError
OperationalError = _lib.OperationalError
IntegrityError = _lib.IntegrityError
InternalError = _lib.InternalError
ProgrammingError = _lib.ProgrammingError
NotSupportedError = _lib.NotSupportedError

_KNOWN_INFO_VALUES = {
    0: "vendor_name",
    1: "vendor_version",
    2: "vendor_arrow_version",
    100: "driver_name",
    101: "driver_version",
    102: "driver_arrow_version",
    103: "driver_adbc_version",
}

# ----------------------------------------------------------
# Types

#: The type for date values.
Date = datetime.date
#: The type for time values.
Time = datetime.time
#: The type for timestamp values.
Timestamp = datetime.datetime


def DateFromTicks(ticks: int) -> Date:
    """Construct a date value from a count of seconds."""
    # Standard implementations from PEP 249 itself
    return Date(*time.localtime(ticks)[:3])


def TimeFromTicks(ticks: int) -> Time:
    """Construct a time value from a count of seconds."""
    return Time(*time.localtime(ticks)[3:6])


def TimestampFromTicks(ticks: int) -> Timestamp:
    """Construct a timestamp value from a count of seconds."""
    return Timestamp(*time.localtime(ticks)[:6])


class _TypeSet(frozenset):
    """A set of PyArrow type IDs that compares equal to subsets of self."""

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, _TypeSet):
            return not (other - self)
        elif isinstance(other, pyarrow.DataType):
            return other.id in self
        return False


#: The type of binary columns.
BINARY = _TypeSet({pyarrow.binary().id, pyarrow.large_binary().id})
#: The type of datetime columns.
DATETIME = _TypeSet(
    [
        pyarrow.date32().id,
        pyarrow.date64().id,
        pyarrow.time32("s").id,
        pyarrow.time64("ns").id,
        pyarrow.timestamp("s").id,
    ]
)
#: The type of numeric columns.
NUMBER = _TypeSet(
    [
        pyarrow.int8().id,
        pyarrow.int16().id,
        pyarrow.int32().id,
        pyarrow.int64().id,
        pyarrow.uint8().id,
        pyarrow.uint16().id,
        pyarrow.uint32().id,
        pyarrow.uint64().id,
        pyarrow.float32().id,
        pyarrow.float64().id,
    ]
)
#: The type of "row ID" columns.
ROWID = _TypeSet([pyarrow.int64().id])
#: The type of string columns.
STRING = _TypeSet([pyarrow.string().id, pyarrow.large_string().id])

# ----------------------------------------------------------
# Functions


def connect(
    *,
    driver: str,
    entrypoint: Optional[str] = None,
    db_kwargs: Optional[Dict[str, str]] = None,
    conn_kwargs: Optional[Dict[str, str]] = None,
    autocommit=False,
) -> "Connection":
    """
    Connect to a database via ADBC.

    Parameters
    ----------
    driver
        The driver name. For example, "adbc_driver_sqlite" will
        attempt to load libadbc_driver_sqlite.so on Linux systems,
        libadbc_driver_sqlite.dylib on MacOS, and
        adbc_driver_sqlite.dll on Windows. This may also be a path to
        the library to load.
    entrypoint
        The driver-specific entrypoint, if different than the default.
    db_kwargs
        Key-value parameters to pass to the driver to initialize the
        database.
    conn_kwargs
        Key-value parameters to pass to the driver to initialize the
        connection.
    autocommit
        Whether to enable autocommit.  For compliance with DB-API,
        this is disabled by default.  A warning will be emitted if it
        cannot be disabled.
    """
    db = None
    conn = None

    db_kwargs = dict(db_kwargs or {})
    db_kwargs["driver"] = driver
    if entrypoint:
        db_kwargs["entrypoint"] = entrypoint
    if conn_kwargs is None:
        conn_kwargs = {}

    try:
        db = _lib.AdbcDatabase(**db_kwargs)
        conn = _lib.AdbcConnection(db, **conn_kwargs)
        return Connection(db, conn, conn_kwargs, autocommit=autocommit)
    except Exception:
        if conn:
            conn.close()
        if db:
            db.close()
        raise


# ----------------------------------------------------------
# Classes


class _Closeable(abc.ABC):
    """Base class providing context manager interface."""

    def __enter__(self) -> "Self":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    @abc.abstractmethod
    def close(self) -> None: ...


class _SharedDatabase(_Closeable):
    """A holder for a shared AdbcDatabase."""

    def __init__(self, db: _lib.AdbcDatabase) -> None:
        self._db = db
        self._lock = threading.Lock()
        self._refcount = 1

    def _inc(self) -> None:
        with self._lock:
            self._refcount += 1

    def _dec(self) -> int:
        with self._lock:
            self._refcount -= 1
            return self._refcount

    def clone(self) -> "Self":
        self._inc()
        return self

    def close(self) -> None:
        if self._dec() == 0:
            self._db.close()


class Connection(_Closeable):
    """
    A DB-API 2.0 (PEP 249) connection.

    Do not create this object directly; use connect().
    """

    # Optional extension: expose exception classes on Connection
    Warning = _lib.Warning
    Error = _lib.Error
    InterfaceError = _lib.InterfaceError
    DatabaseError = _lib.DatabaseError
    DataError = _lib.DataError
    OperationalError = _lib.OperationalError
    IntegrityError = _lib.IntegrityError
    InternalError = _lib.InternalError
    ProgrammingError = _lib.ProgrammingError
    NotSupportedError = _lib.NotSupportedError

    def __init__(
        self,
        db: Union[_lib.AdbcDatabase, _SharedDatabase],
        conn: _lib.AdbcConnection,
        conn_kwargs: Optional[Dict[str, str]] = None,
        *,
        autocommit=False,
    ) -> None:
        self._closed = False
        if isinstance(db, _SharedDatabase):
            self._db = db.clone()
        else:
            self._db = _SharedDatabase(db)
        self._conn = conn
        self._conn_kwargs = conn_kwargs

        try:
            self._conn.set_autocommit(False)
        except _lib.NotSupportedError:
            self._commit_supported = False
            if not autocommit:
                warnings.warn(
                    "Cannot disable autocommit; conn will not be DB-API 2.0 compliant",
                    category=Warning,
                )
            self._autocommit = True
        else:
            self._autocommit = False
            self._commit_supported = True

        if autocommit and self._commit_supported:
            self._conn.set_autocommit(True)
            self._autocommit = True

    def close(self) -> None:
        """
        Close the connection.

        Warnings
        --------
        Failure to close a connection may leak memory or database
        connections.
        """
        if self._closed:
            return

        self._conn.close()
        self._db.close()
        self._closed = True

    def commit(self) -> None:
        """Explicitly commit."""
        if self._commit_supported:
            self._conn.commit()

    def cursor(self) -> "Cursor":
        """Create a new cursor for querying the database."""
        return Cursor(self)

    def rollback(self) -> None:
        """Explicitly rollback."""
        if self._commit_supported:
            self._conn.rollback()

    def __del__(self) -> None:
        if not self._closed:
            self.close()
            _warn_unclosed("adbc_driver_manager.dbapi.Connection")

    # ------------------------------------------------------------
    # API Extensions
    # ------------------------------------------------------------

    def adbc_cancel(self) -> None:
        """
        Cancel any ongoing operations on this connection.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        self._conn.cancel()

    def adbc_clone(self) -> "Connection":
        """
        Create a new Connection sharing the same underlying database.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        conn = _lib.AdbcConnection(self._db._db, **(self._conn_kwargs or {}))
        return Connection(self._db, conn)

    def adbc_get_info(self) -> Dict[Union[str, int], Any]:
        """
        Get metadata about the database and driver.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        handle = _blocking_call(self._conn.get_info, (), {}, self._conn.cancel)
        reader = pyarrow.RecordBatchReader._import_from_c(handle.address)
        table = _blocking_call(reader.read_all, (), {}, self._conn.cancel)
        info = table.to_pylist()
        return dict(
            {
                _KNOWN_INFO_VALUES.get(row["info_name"], row["info_name"]): row[
                    "info_value"
                ]
                for row in info
            }
        )

    def adbc_get_objects(
        self,
        *,
        depth: Literal["all", "catalogs", "db_schemas", "tables", "columns"] = "all",
        catalog_filter: Optional[str] = None,
        db_schema_filter: Optional[str] = None,
        table_name_filter: Optional[str] = None,
        table_types_filter: Optional[List[str]] = None,
        column_name_filter: Optional[str] = None,
    ) -> pyarrow.RecordBatchReader:
        """
        List catalogs, schemas, tables, etc. in the database.

        Parameters
        ----------
        depth
            What objects to return info on.
        catalog_filter
            An optional filter on the catalog names returned.
        db_schema_filter
            An optional filter on the database schema names returned.
        table_name_filter
            An optional filter on the table names returned.
        table_types_filter
            An optional list of types of tables returned.
        column_name_filter
            An optional filter on the column names returned.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        if depth in ("all", "columns"):
            c_depth = _lib.GetObjectsDepth.ALL
        elif depth == "catalogs":
            c_depth = _lib.GetObjectsDepth.CATALOGS
        elif depth == "db_schemas":
            c_depth = _lib.GetObjectsDepth.DB_SCHEMAS
        elif depth == "tables":
            c_depth = _lib.GetObjectsDepth.TABLES
        else:
            raise ValueError(f"Invalid value for 'depth': {depth}")
        handle = _blocking_call(
            self._conn.get_objects,
            (c_depth,),
            dict(
                catalog=catalog_filter,
                db_schema=db_schema_filter,
                table_name=table_name_filter,
                table_types=table_types_filter,
                column_name=column_name_filter,
            ),
            self._conn.cancel,
        )
        return pyarrow.RecordBatchReader._import_from_c(handle.address)

    def adbc_get_table_schema(
        self,
        table_name: str,
        *,
        catalog_filter: Optional[str] = None,
        db_schema_filter: Optional[str] = None,
    ) -> pyarrow.Schema:
        """
        Get the Arrow schema of a table by name.

        Parameters
        ----------
        table_name
            The table to get the schema of.
        catalog_filter
            An optional filter on the catalog name of the table.
        db_schema_filter
            An optional filter on the database schema name of the table.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        handle = _blocking_call(
            self._conn.get_table_schema,
            (
                catalog_filter,
                db_schema_filter,
                table_name,
            ),
            {},
            self._conn.cancel,
        )
        return pyarrow.Schema._import_from_c(handle.address)

    def adbc_get_table_types(self) -> List[str]:
        """
        List the types of tables that the server knows about.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        handle = _blocking_call(
            self._conn.get_table_types,
            (),
            {},
            self._conn.cancel,
        )
        reader = pyarrow.RecordBatchReader._import_from_c(handle.address)
        table = reader.read_all()
        return table[0].to_pylist()

    @property
    def adbc_connection(self) -> _lib.AdbcConnection:
        """
        Get the underlying ADBC connection.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        return self._conn

    @property
    def adbc_current_catalog(self) -> str:
        """
        The name of the current catalog.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        key = adbc_driver_manager.ConnectionOptions.CURRENT_CATALOG.value
        return self._conn.get_option(key)

    @adbc_current_catalog.setter
    def adbc_current_catalog(self, catalog: str) -> None:
        key = adbc_driver_manager.ConnectionOptions.CURRENT_CATALOG.value
        self._conn.set_options(**{key: catalog})

    @property
    def adbc_current_db_schema(self) -> str:
        """
        The name of the current schema.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        key = adbc_driver_manager.ConnectionOptions.CURRENT_DB_SCHEMA.value
        return self._conn.get_option(key)

    @adbc_current_db_schema.setter
    def adbc_current_db_schema(self, db_schema: str) -> None:
        key = adbc_driver_manager.ConnectionOptions.CURRENT_DB_SCHEMA.value
        self._conn.set_options(**{key: db_schema})

    @property
    def adbc_database(self) -> _lib.AdbcDatabase:
        """
        Get the underlying ADBC database.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        return self._db._db


class Cursor(_Closeable):
    """
    A DB-API 2.0 (PEP 249) cursor.

    Do not create this object directly; use Connection.cursor().
    """

    def __init__(self, conn: Connection) -> None:
        # Must be at top in case __init__ is interrupted and then __del__ is called
        self._closed = True
        self._conn = conn
        self._stmt = _lib.AdbcStatement(conn._conn)
        self._closed = False

        self._last_query: Optional[Union[str, bytes]] = None
        self._results: Optional["_RowIterator"] = None
        self._arraysize = 1
        self._rowcount = -1

    @property
    def arraysize(self) -> int:
        """The number of rows to fetch at a time with fetchmany()."""
        return self._arraysize

    @arraysize.setter
    def arraysize(self, size: int) -> None:
        self._arraysize = size

    @property
    def connection(self) -> Connection:
        """
        Get the connection associated with this cursor.

        This is an optional DB-API extension.
        """
        return self._conn

    @property
    def description(self) -> Optional[List[tuple]]:
        """The schema of the result set."""
        if self._results is None:
            return None
        return self._results.description

    @property
    def rowcount(self) -> int:
        """
        Get the row count of the result set, or -1 if not known.
        """
        return self._rowcount

    @property
    def rownumber(self) -> Optional[int]:
        """Get the current row number, or None if not applicable."""
        if self._results is not None:
            return self._results.rownumber
        return None

    def callproc(self, procname, parameters):
        """Call a stored procedure (not supported)."""
        raise NotSupportedError("Cursor.callproc")

    def close(self):
        """Close the cursor and free resources."""
        if self._closed:
            return

        if self._results is not None:
            self._results.close()
            self._results = None
        self._stmt.close()
        self._closed = True

    def _bind(self, parameters) -> None:
        if hasattr(parameters, "__arrow_c_array__"):
            self._stmt.bind(parameters)
        elif hasattr(parameters, "__arrow_c_stream__"):
            self._stmt.bind_stream(parameters)
        elif isinstance(parameters, pyarrow.RecordBatch):
            arr_handle = _lib.ArrowArrayHandle()
            sch_handle = _lib.ArrowSchemaHandle()
            parameters._export_to_c(arr_handle.address, sch_handle.address)
            self._stmt.bind(arr_handle, sch_handle)
        else:
            if isinstance(parameters, pyarrow.Table):
                parameters = parameters.to_reader()
            stream_handle = _lib.ArrowArrayStreamHandle()
            parameters._export_to_c(stream_handle.address)
            self._stmt.bind_stream(stream_handle)

    def _prepare_execute(self, operation, parameters=None) -> None:
        self._results = None
        if operation != self._last_query:
            self._last_query = operation
            if isinstance(operation, bytes):
                # Serialized Substrait plan
                self._stmt.set_substrait_plan(operation)
            else:
                self._stmt.set_sql_query(operation)
            try:
                _blocking_call(self._stmt.prepare, (), {}, self._stmt.cancel)
            except NotSupportedError:
                # Not all drivers support it
                pass

        if _is_arrow_data(parameters):
            self._bind(parameters)
        elif parameters:
            rb = pyarrow.record_batch(
                [[param_value] for param_value in parameters],
                names=[str(i) for i in range(len(parameters))],
            )
            self._bind(rb)

    def execute(self, operation: Union[bytes, str], parameters=None) -> None:
        """
        Execute a query.

        Parameters
        ----------
        operation : bytes or str
            The query to execute.  Pass SQL queries as strings,
            (serialized) Substrait plans as bytes.
        parameters
            Parameters to bind.  Can be a Python sequence (to provide
            a single set of parameters), or an Arrow record batch,
            table, or record batch reader (to provide multiple
            parameters, which will each be bound in turn).
        """
        self._prepare_execute(operation, parameters)

        handle, self._rowcount = _blocking_call(
            self._stmt.execute_query, (), {}, self._stmt.cancel
        )
        self._results = _RowIterator(
            self._stmt, _reader.AdbcRecordBatchReader._import_from_c(handle.address)
        )

    def executemany(self, operation: Union[bytes, str], seq_of_parameters) -> None:
        """
        Execute a query with multiple parameter sets.

        This method does not generate a result set.

        Parameters
        ----------
        operation : bytes or str
            The query to execute.  Pass SQL queries as strings,
            (serialized) Substrait plans as bytes.
        seq_of_parameters
            Parameters to bind.  Can be a list of Python sequences, or
            an Arrow record batch, table, or record batch reader.  If
            None, then the query will be executed once, else it will
            be executed once per row.
        """
        self._results = None
        if operation != self._last_query:
            self._last_query = operation
            self._stmt.set_sql_query(operation)
            self._stmt.prepare()

        if _is_arrow_data(seq_of_parameters):
            arrow_parameters = seq_of_parameters
        elif seq_of_parameters:
            arrow_parameters = pyarrow.RecordBatch.from_pydict(
                {
                    str(col_idx): pyarrow.array(x)
                    for col_idx, x in enumerate(map(list, zip(*seq_of_parameters)))
                },
            )
        else:
            arrow_parameters = pyarrow.record_batch([])

        self._bind(arrow_parameters)
        self._rowcount = _blocking_call(
            self._stmt.execute_update, (), {}, self._stmt.cancel
        )

    def fetchone(self) -> Optional[tuple]:
        """Fetch one row of the result."""
        if self._results is None:
            raise ProgrammingError(
                "Cannot fetchone() before execute()",
                status_code=_lib.AdbcStatusCode.INVALID_STATE,
            )
        return self._results.fetchone()

    def fetchmany(self, size: Optional[int] = None) -> List[tuple]:
        """Fetch some rows of the result."""
        if self._results is None:
            raise ProgrammingError(
                "Cannot fetchmany() before execute()",
                status_code=_lib.AdbcStatusCode.INVALID_STATE,
            )
        if size is None:
            size = self.arraysize
        return self._results.fetchmany(size)

    def fetchall(self) -> List[tuple]:
        """Fetch all rows of the result."""
        if self._results is None:
            raise ProgrammingError(
                "Cannot fetchall() before execute()",
                status_code=_lib.AdbcStatusCode.INVALID_STATE,
            )
        return self._results.fetchall()

    def next(self):
        """Fetch the next row, or raise StopIteration."""
        row = self.fetchone()
        if row is None:
            raise StopIteration
        return row

    def nextset(self):
        """Move to the next available result set (not supported)."""
        raise NotSupportedError("Cursor.nextset")

    def setinputsizes(self, sizes):
        """Preallocate memory for the parameters (no-op)."""
        pass

    def setoutputsize(self, size, column=None):
        """Preallocate memory for the result set (no-op)."""
        pass

    def __del__(self) -> None:
        if not self._closed:
            self.close()
            _warn_unclosed("adbc_driver_manager.dbapi.Cursor")

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    # ------------------------------------------------------------
    # API Extensions
    # ------------------------------------------------------------

    def adbc_cancel(self) -> None:
        """
        Cancel any ongoing operations on this statement.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        self._stmt.cancel()

    def adbc_ingest(
        self,
        table_name: str,
        data: Union[pyarrow.RecordBatch, pyarrow.Table, pyarrow.RecordBatchReader],
        mode: Literal["append", "create", "replace", "create_append"] = "create",
        *,
        catalog_name: Optional[str] = None,
        db_schema_name: Optional[str] = None,
        temporary: bool = False,
    ) -> int:
        """Ingest Arrow data into a database table.

        Depending on the driver, this can avoid per-row overhead that
        would result from a typical prepare-bind-insert loop.

        Parameters
        ----------
        table_name
            The table to insert into.
        data
            The Arrow data to insert. This can be a pyarrow RecordBatch, Table
            or RecordBatchReader, or any Arrow-compatible data that implements
            the Arrow PyCapsule Protocol (i.e. has an ``__arrow_c_array__``
            or ``__arrow_c_stream__`` method).
        mode
            How to deal with existing data:

            - 'append': append to a table (error if table does not exist)
            - 'create': create a table and insert (error if table exists)
            - 'create_append': create a table (if not exists) and insert
            - 'replace': drop existing table (if any), then same as 'create'
        catalog_name
            If given, the catalog to create/locate the table in.
            **This API is EXPERIMENTAL.**
        db_schema_name
            If given, the schema to create/locate the table in.
            **This API is EXPERIMENTAL.**
        temporary
            Whether to ingest to a temporary table or not.  Most drivers will
            not support setting this along with catalog_name and/or
            db_schema_name.
            **This API is EXPERIMENTAL.**

        Returns
        -------
        int
            The number of rows inserted, or -1 if the driver cannot
            provide this information.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.

        """
        if mode == "append":
            c_mode = _lib.INGEST_OPTION_MODE_APPEND
        elif mode == "create":
            c_mode = _lib.INGEST_OPTION_MODE_CREATE
        elif mode == "create_append":
            c_mode = _lib.INGEST_OPTION_MODE_CREATE_APPEND
        elif mode == "replace":
            c_mode = _lib.INGEST_OPTION_MODE_REPLACE
        else:
            raise ValueError(f"Invalid value for 'mode': {mode}")

        options = {
            _lib.INGEST_OPTION_TARGET_TABLE: table_name,
            _lib.INGEST_OPTION_MODE: c_mode,
        }
        if catalog_name is not None:
            options[
                adbc_driver_manager.StatementOptions.INGEST_TARGET_CATALOG.value
            ] = catalog_name
        if db_schema_name is not None:
            options[
                adbc_driver_manager.StatementOptions.INGEST_TARGET_DB_SCHEMA.value
            ] = db_schema_name
        self._stmt.set_options(**options)

        if temporary:
            self._stmt.set_options(
                **{
                    adbc_driver_manager.StatementOptions.INGEST_TEMPORARY.value: "true",
                }
            )
        else:
            # Need to explicitly clear it, but not all drivers support this
            options = {
                adbc_driver_manager.StatementOptions.INGEST_TEMPORARY.value: "false",
            }
            try:
                self._stmt.set_options(**options)
            except NotSupportedError:
                pass

        if hasattr(data, "__arrow_c_array__"):
            self._stmt.bind(data)
        elif hasattr(data, "__arrow_c_stream__"):
            self._stmt.bind_stream(data)
        elif isinstance(data, pyarrow.RecordBatch):
            array = _lib.ArrowArrayHandle()
            schema = _lib.ArrowSchemaHandle()
            data._export_to_c(array.address, schema.address)
            self._stmt.bind(array, schema)
        else:
            if isinstance(data, pyarrow.Table):
                data = data.to_reader()
            elif isinstance(data, pyarrow.dataset.Dataset):
                data = data.scanner().to_reader()
            elif isinstance(data, pyarrow.dataset.Scanner):
                data = data.to_reader()
            elif not hasattr(data, "_export_to_c"):
                data = pyarrow.Table.from_batches(data)
                data = data.to_reader()
            handle = _lib.ArrowArrayStreamHandle()
            data._export_to_c(handle.address)
            self._stmt.bind_stream(handle)

        self._last_query = None
        return _blocking_call(self._stmt.execute_update, (), {}, self._stmt.cancel)

    def adbc_execute_partitions(
        self,
        operation,
        parameters=None,
    ) -> Tuple[List[bytes], pyarrow.Schema]:
        """
        Execute a query and get the partitions of a distributed result set.

        Returns
        -------
        partitions : list of byte
            A list of partition descriptors, which can be read with
            read_partition.
        schema : pyarrow.Schema or None
            The schema of the result set.  May be None if incremental query
            execution is enabled and the server has not returned a schema.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        self._prepare_execute(operation, parameters)
        partitions, schema_handle, self._rowcount = _blocking_call(
            self._stmt.execute_partitions, (), {}, self._stmt.cancel
        )
        if schema_handle and schema_handle.address:
            schema = pyarrow.Schema._import_from_c(schema_handle.address)
        else:
            schema = None
        return partitions, schema

    def adbc_execute_schema(self, operation, parameters=None) -> pyarrow.Schema:
        """
        Get the schema of the result set of a query without executing it.

        Returns
        -------
        pyarrow.Schema
            The schema of the result set.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        self._prepare_execute(operation, parameters)
        schema = _blocking_call(self._stmt.execute_schema, (), {}, self._stmt.cancel)
        return pyarrow.Schema._import_from_c(schema.address)

    def adbc_prepare(self, operation: Union[bytes, str]) -> Optional[pyarrow.Schema]:
        """
        Prepare a query without executing it.

        To execute the query afterwards, call :meth:`execute` or
        :meth:`executemany` with the same query.  This will not
        prepare the query a second time.

        Returns
        -------
        pyarrow.Schema or None
            The schema of the bind parameters, or None if the schema
            could not be determined.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        self._prepare_execute(operation)

        try:
            handle = _blocking_call(
                self._stmt.get_parameter_schema, (), {}, self._stmt.cancel
            )
        except NotSupportedError:
            return None
        return pyarrow.Schema._import_from_c(handle.address)

    def adbc_read_partition(self, partition: bytes) -> None:
        """
        Read a partition of a distributed result set.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        self._results = None
        handle = _blocking_call(
            self._conn._conn.read_partition, (partition,), {}, self._stmt.cancel
        )
        self._rowcount = -1
        self._results = _RowIterator(
            self._stmt, pyarrow.RecordBatchReader._import_from_c(handle.address)
        )

    @property
    def adbc_statement(self) -> _lib.AdbcStatement:
        """
        Get the underlying ADBC statement.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        return self._stmt

    def executescript(self, operation: str) -> None:
        """
        Execute multiple statements.

        If there is a pending transaction, commits first.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        if not self._conn._autocommit:
            self._conn.commit()

        self._last_query = None
        self._results = None
        self._stmt.set_sql_query(operation)
        _blocking_call(self._stmt.execute_update, (), {}, self._stmt.cancel)

    def fetchallarrow(self) -> pyarrow.Table:
        """
        Fetch all rows of the result as a PyArrow Table.

        This implements a similar API as turbodbc.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        return self.fetch_arrow_table()

    def fetch_arrow_table(self) -> pyarrow.Table:
        """
        Fetch all rows of the result as a PyArrow Table.

        This implements a similar API as DuckDB.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        if self._results is None:
            raise ProgrammingError(
                "Cannot fetch_arrow_table() before execute()",
                status_code=_lib.AdbcStatusCode.INVALID_STATE,
            )
        return self._results.fetch_arrow_table()

    def fetch_df(self) -> "pandas.DataFrame":
        """
        Fetch all rows of the result as a Pandas DataFrame.

        This implements a similar API as DuckDB.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        if self._results is None:
            raise ProgrammingError(
                "Cannot fetch_df() before execute()",
                status_code=_lib.AdbcStatusCode.INVALID_STATE,
            )
        return self._results.fetch_df()

    def fetch_record_batch(self) -> pyarrow.RecordBatchReader:
        """
        Fetch the result as a PyArrow RecordBatchReader.

        This implements a similar API as DuckDB:
        https://duckdb.org/docs/guides/python/export_arrow.html#export-as-a-recordbatchreader

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        if self._results is None:
            raise ProgrammingError(
                "Cannot fetch_record_batch() before execute()",
                status_code=_lib.AdbcStatusCode.INVALID_STATE,
            )
        # XXX(https://github.com/apache/arrow-adbc/issues/1523): return the
        # "real" PyArrow reader since PyArrow may try to poke the internal C++
        # reader pointer
        return self._results._reader._reader


# ----------------------------------------------------------
# Utilities


class _RowIterator(_Closeable):
    """Track state needed to iterate over the result set."""

    def __init__(self, stmt, reader: pyarrow.RecordBatchReader) -> None:
        self._stmt = stmt
        self._reader = reader
        self._current_batch = None
        self._next_row = 0
        self._finished = False
        self.rownumber = 0

    def close(self) -> None:
        if hasattr(self._reader, "close"):
            # Only in recent PyArrow
            self._reader.close()

    @property
    def description(self) -> List[tuple]:
        return [
            (field.name, field.type, None, None, None, None, None)
            for field in self._reader.schema
        ]

    def fetchone(self) -> Optional[tuple]:
        if self._current_batch is None or self._next_row >= len(self._current_batch):
            try:
                while True:
                    self._current_batch = _blocking_call(
                        self._reader.read_next_batch, (), {}, self._stmt.cancel
                    )
                    if self._current_batch.num_rows > 0:
                        break
                self._next_row = 0
            except StopIteration:
                self._current_batch = None
                self._finished = True

        if self._finished or self._current_batch is None:
            return None

        row = tuple(arr[self._next_row].as_py() for arr in self._current_batch.columns)
        self._next_row += 1
        self.rownumber += 1
        return row

    def fetchmany(self, size: int) -> List[tuple]:
        rows = []
        for _ in range(size):
            row = self.fetchone()
            if row is None:
                break
            rows.append(row)
        return rows

    def fetchall(self) -> List[tuple]:
        rows = []
        while True:
            row = self.fetchone()
            if row is None:
                break
            rows.append(row)
        return rows

    def fetch_arrow_table(self) -> pyarrow.Table:
        return _blocking_call(self._reader.read_all, (), {}, self._stmt.cancel)

    def fetch_df(self) -> "pandas.DataFrame":
        return _blocking_call(self._reader.read_pandas, (), {}, self._stmt.cancel)


_PYTEST_ENV_VAR = "PYTEST_CURRENT_TEST"
_ADBC_ENV_VAR = "_ADBC_DRIVER_MANAGER_WARN_UNCLOSED_RESOURCE"


def _warn_unclosed(name):
    if _PYTEST_ENV_VAR in os.environ or os.environ.get(_ADBC_ENV_VAR) == "1":
        warnings.warn(
            f"A {name} was not explicitly close()d, which may leak "
            f"driver resources. (This warning is only emitted if "
            f"{_PYTEST_ENV_VAR} or {_ADBC_ENV_VAR} are set.)",
            category=ResourceWarning,
            stacklevel=2,
        )


def _is_arrow_data(data):
    return (
        hasattr(data, "__arrow_c_array__")
        or hasattr(data, "__arrow_c_stream__")
        or isinstance(
            data, (pyarrow.RecordBatch, pyarrow.Table, pyarrow.RecordBatchReader)
        )
    )
