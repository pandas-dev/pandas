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

"""Low-level ADBC bindings for the SQLite driver."""

import enum
import functools
import typing

import adbc_driver_manager

from ._version import __version__  # noqa:F401

__all__ = ["ConnectionOptions", "StatementOptions", "connect"]


class ConnectionOptions(enum.Enum):
    """Connection options specific to the SQLite driver."""

    #: Whether to enable ("true") or disable ("false") extension loading.
    #: Default is disabled.
    LOAD_EXTENSION_ENABLED = "adbc.sqlite.load_extension.enabled"

    #: The path to an extension to load.
    #: Set this option after LOAD_EXTENSION_PATH.  This will actually
    #: load the extension.
    LOAD_EXTENSION_ENTRYPOINT = "adbc.sqlite.load_extension.entrypoint"

    #: The path to an extension to load.
    #: First set this option, then LOAD_EXTENSION_ENTRYPOINT.  The second
    #: call will actually load the extension.
    LOAD_EXTENSION_PATH = "adbc.sqlite.load_extension.path"


class StatementOptions(enum.Enum):
    """Statement options specific to the SQLite driver."""

    #: The number of rows per batch. Defaults to 1024.
    BATCH_ROWS = "adbc.sqlite.query.batch_rows"


def connect(uri: typing.Optional[str] = None) -> adbc_driver_manager.AdbcDatabase:
    """Create a low level ADBC connection to SQLite."""
    if uri is None:
        return adbc_driver_manager.AdbcDatabase(driver=_driver_path())
    return adbc_driver_manager.AdbcDatabase(driver=_driver_path(), uri=uri)


@functools.lru_cache
def _driver_path() -> str:
    import pathlib
    import sys

    import importlib_resources

    driver = "adbc_driver_sqlite"

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
