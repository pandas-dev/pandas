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
DBAPI 2.0-compatible facade for the ADBC libpq driver.
"""

import typing

import adbc_driver_manager
import adbc_driver_manager.dbapi
import adbc_driver_postgresql

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
# XXX: PostgreSQL doesn't fit any of the param styles
# We'll need some Python-side wrangling specific to this driver
paramstyle = "pyformat"

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


def connect(
    uri: str,
    db_kwargs: typing.Optional[typing.Dict[str, str]] = None,
    conn_kwargs: typing.Optional[typing.Dict[str, str]] = None,
    **kwargs,
) -> "Connection":
    """
    Connect to PostgreSQL via ADBC.

    Parameters
    ----------
    uri : str
        The URI to connect to.
    db_kwargs : dict, optional
        Initial database connection parameters.
    conn_kwargs : dict, optional
        Connection-specific parameters.  (ADBC differentiates between
        a 'database' object shared between multiple 'connection'
        objects.)
    """
    db = None
    conn = None

    try:
        db = adbc_driver_postgresql.connect(uri, db_kwargs=db_kwargs)
        conn = adbc_driver_manager.AdbcConnection(db, **(conn_kwargs or {}))
        return adbc_driver_manager.dbapi.Connection(
            db, conn, conn_kwargs=conn_kwargs, **kwargs
        )
    except Exception:
        if conn:
            conn.close()
        if db:
            db.close()
        raise


# ----------------------------------------------------------
# Classes

Connection = adbc_driver_manager.dbapi.Connection
Cursor = adbc_driver_manager.dbapi.Cursor
