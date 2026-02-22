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

"""Low-level ADBC bindings for Python.

The root module provides a fairly direct, 1:1 mapping to the C API
definitions in Python.  For a higher-level interface, use
:mod:`adbc_driver_manager.dbapi`.  (This requires PyArrow.)
"""

import enum

from ._lib import (
    INGEST_OPTION_MODE,
    INGEST_OPTION_MODE_APPEND,
    INGEST_OPTION_MODE_CREATE,
    INGEST_OPTION_TARGET_TABLE,
    AdbcConnection,
    AdbcDatabase,
    AdbcInfoCode,
    AdbcStatement,
    AdbcStatusCode,
    ArrowArrayHandle,
    ArrowArrayStreamHandle,
    ArrowSchemaHandle,
    DatabaseError,
    DataError,
    Error,
    GetObjectsDepth,
    IntegrityError,
    InterfaceError,
    InternalError,
    NotSupportedError,
    OperationalError,
    ProgrammingError,
    Warning,
)
from ._version import __version__  # noqa:F401

__all__ = [
    "INGEST_OPTION_MODE",
    "INGEST_OPTION_MODE_APPEND",
    "INGEST_OPTION_MODE_CREATE",
    "INGEST_OPTION_TARGET_TABLE",
    "AdbcConnection",
    "AdbcDatabase",
    "AdbcInfoCode",
    "AdbcStatement",
    "AdbcStatusCode",
    "ArrowArrayHandle",
    "ArrowArrayStreamHandle",
    "ArrowSchemaHandle",
    "ConnectionOptions",
    "DatabaseError",
    "DatabaseOptions",
    "DataError",
    "Error",
    "GetObjectsDepth",
    "IntegrityError",
    "InterfaceError",
    "InternalError",
    "NotSupportedError",
    "OperationalError",
    "ProgrammingError",
    "StatementOptions",
    "Warning",
]


class DatabaseOptions(enum.Enum):
    """
    Database options that are standardized between drivers.

    Not all drivers support all options.
    """

    #: Set the password to use for username-password authentication.
    PASSWORD = "password"
    #: The URI to connect to.
    URI = "uri"
    #: Set the username to use for username-password authentication.
    USERNAME = "username"


class ConnectionOptions(enum.Enum):
    """Connection options that are standardized between drivers.

    Not all drivers support all options.
    """

    #: Get/set the current catalog.
    CURRENT_CATALOG = "adbc.connection.catalog"
    #: Get/set the current schema.
    CURRENT_DB_SCHEMA = "adbc.connection.db_schema"
    #: Set the transaction isolation level.
    ISOLATION_LEVEL = "adbc.connection.transaction.isolation_level"


class StatementOptions(enum.Enum):
    """Statement options that are standardized between drivers.

    Not all drivers support all options.
    """

    #: Bind parameters by name instead of by position.
    BIND_BY_NAME = "adbc.statement.bind_by_name"
    #: Enable incremental execution on ExecutePartitions.
    INCREMENTAL = "adbc.statement.exec.incremental"
    #: For bulk ingestion, whether to create or append to the table.
    INGEST_MODE = INGEST_OPTION_MODE
    #: For bulk ingestion, the table to ingest into.
    INGEST_TARGET_TABLE = INGEST_OPTION_TARGET_TABLE
    #: For bulk ingestion, the catalog to create/locate the table in.
    #: **This API is EXPERIMENTAL.**
    INGEST_TARGET_CATALOG = "adbc.ingest.target_catalog"
    #: For bulk ingestion, the schema to create/locate the table in.
    #: **This API is EXPERIMENTAL.**
    INGEST_TARGET_DB_SCHEMA = "adbc.ingest.target_db_schema"
    #: For bulk ingestion, use a temporary table.
    #: **This API is EXPERIMENTAL.**
    INGEST_TEMPORARY = "adbc.ingest.temporary"
    #: Get progress of a query.
    PROGRESS = "adbc.statement.exec.progress"
