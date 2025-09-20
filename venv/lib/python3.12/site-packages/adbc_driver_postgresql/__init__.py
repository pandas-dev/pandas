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

import enum
import functools
import typing

import adbc_driver_manager

from ._version import __version__

__all__ = [
    "ConnectionOptions",
    "StatementOptions",
    "connect",
    "__version__",
]


class ConnectionOptions(enum.Enum):
    """Connection options specific to the PostgreSQL driver."""

    #: Get the transaction status.
    TRANSACTION_STATUS = "adbc.postgresql.transaction_status"


class StatementOptions(enum.Enum):
    """Statement options specific to the PostgreSQL driver."""

    #: Try to limit returned batches to this size (in bytes).
    #:
    #: This is merely a hint, and because the size is estimated, the
    #: actual size may differ.
    BATCH_SIZE_HINT_BYTES = "adbc.postgresql.batch_size_hint_bytes"

    #: Enable or disable the ``COPY`` optimization (default: enabled).
    #:
    #: This is necessary for some queries since PostgreSQL does not support
    #: ``COPY`` with those queries, e.g. queries using ``SHOW``.
    USE_COPY = "adbc.postgresql.use_copy"


def connect(
    uri: str,
    db_kwargs: typing.Optional[typing.Dict[str, str]] = None,
) -> adbc_driver_manager.AdbcDatabase:
    """Create a low level ADBC connection to PostgreSQL."""
    db_options = dict(db_kwargs or {})
    db_options["driver"] = _driver_path()
    db_options["uri"] = uri
    return adbc_driver_manager.AdbcDatabase(**db_options)


@functools.lru_cache
def _driver_path() -> str:
    import pathlib
    import sys

    import importlib_resources

    driver = "adbc_driver_postgresql"

    # Wheels bundle the shared library
    root = importlib_resources.files(driver)
    # The filename is always the same regardless of platform
    entrypoint = root.joinpath(f"lib{driver}.so")
    if entrypoint.is_file():
        return str(entrypoint)

    # Search sys.prefix + '/lib' (Unix, Conda on Unix)
    root = pathlib.Path(sys.prefix)
    for filename in (f"lib{driver}.so", f"lib{driver}.dylib"):
        entrypoint = root.joinpath("lib", filename)
        if entrypoint.is_file():
            return str(entrypoint)

    # Conda on Windows
    entrypoint = root.joinpath("bin", f"{driver}.dll")
    if entrypoint.is_file():
        return str(entrypoint)

    # Let the driver manager fall back to (DY)LD_LIBRARY_PATH/PATH
    # (It will insert 'lib', 'so', etc. as needed)
    return driver
