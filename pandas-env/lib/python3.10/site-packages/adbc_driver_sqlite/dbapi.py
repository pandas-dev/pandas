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

"""
DBAPI 2.0-compatible facade for the ADBC SQLite driver.
"""

import typing

import adbc_driver_manager
import adbc_driver_manager.dbapi
import adbc_driver_sqlite

__all__ = [
    "BINARY",
    "DATETIME",
    "NUMBER",
    "ROWID",
    "STRING",
    "Connection",
    "Cursor",
    "DataError",
    "DatabaseError",
    "Date",
    "DateFromTicks",
    "Error",
    "IntegrityError",
    "InterfaceError",
    "InternalError",
    "NotSupportedError",
    "OperationalError",
    "ProgrammingError",
    "Time",
    "TimeFromTicks",
    "Timestamp",
    "TimestampFromTicks",
    "Warning",
    "apilevel",
    "connect",
    "paramstyle",
    "threadsafety",
]

# ----------------------------------------------------------
# Globals

apilevel = adbc_driver_manager.dbapi.apilevel
threadsafety = adbc_driver_manager.dbapi.threadsafety
paramstyle = "qmark"

Warning = adbc_driver_manager.dbapi.Warning
Error = adbc_driver_manager.dbapi.Error
InterfaceError = adbc_driver_manager.dbapi.InterfaceError
DatabaseError = adbc_driver_manager.dbapi.DatabaseError
DataError = adbc_driver_manager.dbapi.DataError
OperationalError = adbc_driver_manager.dbapi.OperationalError
IntegrityError = adbc_driver_manager.dbapi.IntegrityError
InternalError = adbc_driver_manager.dbapi.InternalError
ProgrammingError = adbc_driver_manager.dbapi.ProgrammingError
NotSupportedError = adbc_driver_manager.dbapi.NotSupportedError

# ----------------------------------------------------------
# Types

Date = adbc_driver_manager.dbapi.Date
Time = adbc_driver_manager.dbapi.Time
Timestamp = adbc_driver_manager.dbapi.Timestamp
DateFromTicks = adbc_driver_manager.dbapi.DateFromTicks
TimeFromTicks = adbc_driver_manager.dbapi.TimeFromTicks
TimestampFromTicks = adbc_driver_manager.dbapi.TimestampFromTicks
STRING = adbc_driver_manager.dbapi.STRING
BINARY = adbc_driver_manager.dbapi.BINARY
NUMBER = adbc_driver_manager.dbapi.NUMBER
DATETIME = adbc_driver_manager.dbapi.DATETIME
ROWID = adbc_driver_manager.dbapi.ROWID

# ----------------------------------------------------------
# Functions


def connect(uri: typing.Optional[str] = None, **kwargs) -> "Connection":
    """Connect to SQLite via ADBC."""
    db = None
    conn = None

    try:
        db = adbc_driver_sqlite.connect(uri)
        conn = adbc_driver_manager.AdbcConnection(db)
        return Connection(db, conn, **kwargs)
    except Exception:
        if conn:
            conn.close()
        if db:
            db.close()
        raise


# ----------------------------------------------------------
# Classes


class AdbcSqliteConnection(adbc_driver_manager.dbapi.Connection):
    """
    A connection to an SQLite 3 database.

    This adds SQLite-specific functionality to the base ADBC-DBAPI bindings in
    the adbc_driver_manager.dbapi module.
    """

    def enable_load_extension(self, enabled: bool) -> None:
        """
        Toggle whether extension loading is allowed.

        Parameters
        ----------
        enabled
            Whether extension loading is allowed or not.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.
        """
        flag = adbc_driver_sqlite.ConnectionOptions.LOAD_EXTENSION_ENABLED.value
        self.adbc_connection.set_options(**{flag: "true" if enabled else "false"})

    def load_extension(
        self, path: str, *, entrypoint: typing.Optional[str] = None
    ) -> None:
        """
        Load an extension into the current connection.

        Parameters
        ----------
        path
            The path to the extension to load.
        entrypoint
            The entrypoint to the extension.  If not provided or None, then
            SQLite will derive its own entrypoint name.

        Notes
        -----
        This is an extension and not part of the DBAPI standard.

        See the SQLite documentation for general information on extensions:
        https://www.sqlite.org/loadext.html
        """
        flag = adbc_driver_sqlite.ConnectionOptions.LOAD_EXTENSION_PATH.value
        self.adbc_connection.set_options(**{flag: path})
        flag = adbc_driver_sqlite.ConnectionOptions.LOAD_EXTENSION_ENTRYPOINT.value
        self.adbc_connection.set_options(**{flag: entrypoint})


Connection = AdbcSqliteConnection
Cursor = adbc_driver_manager.dbapi.Cursor
