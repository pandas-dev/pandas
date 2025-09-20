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
"""Base FileIO classes for implementing reading and writing table files.

The FileIO abstraction includes a subset of full filesystem implementations. Specifically,
Iceberg needs to read or write a file at a given location (as a seekable stream), as well
as check if a file exists. An implementation of the FileIO abstract base class is responsible
for returning an InputFile instance, an OutputFile instance, and deleting a file given
its location.
"""

from __future__ import annotations

import importlib
import logging
import warnings
from abc import ABC, abstractmethod
from io import SEEK_SET
from types import TracebackType
from typing import (
    Dict,
    List,
    Optional,
    Protocol,
    Type,
    Union,
    runtime_checkable,
)
from urllib.parse import urlparse

from pyiceberg.typedef import EMPTY_DICT, Properties

logger = logging.getLogger(__name__)

AWS_REGION = "client.region"
AWS_ACCESS_KEY_ID = "client.access-key-id"
AWS_SECRET_ACCESS_KEY = "client.secret-access-key"
AWS_SESSION_TOKEN = "client.session-token"
AWS_ROLE_ARN = "client.role-arn"
AWS_ROLE_SESSION_NAME = "client.role-session-name"
S3_ANONYMOUS = "s3.anonymous"
S3_ENDPOINT = "s3.endpoint"
S3_ACCESS_KEY_ID = "s3.access-key-id"
S3_SECRET_ACCESS_KEY = "s3.secret-access-key"
S3_SESSION_TOKEN = "s3.session-token"
S3_REGION = "s3.region"
S3_RESOLVE_REGION = "s3.resolve-region"
S3_PROXY_URI = "s3.proxy-uri"
S3_CONNECT_TIMEOUT = "s3.connect-timeout"
S3_REQUEST_TIMEOUT = "s3.request-timeout"
S3_SIGNER = "s3.signer"
S3_SIGNER_URI = "s3.signer.uri"
S3_SIGNER_ENDPOINT = "s3.signer.endpoint"
S3_SIGNER_ENDPOINT_DEFAULT = "v1/aws/s3/sign"
S3_ROLE_ARN = "s3.role-arn"
S3_ROLE_SESSION_NAME = "s3.role-session-name"
S3_FORCE_VIRTUAL_ADDRESSING = "s3.force-virtual-addressing"
S3_RETRY_STRATEGY_IMPL = "s3.retry-strategy-impl"
HDFS_HOST = "hdfs.host"
HDFS_PORT = "hdfs.port"
HDFS_USER = "hdfs.user"
HDFS_KERB_TICKET = "hdfs.kerberos_ticket"
ADLS_CONNECTION_STRING = "adls.connection-string"
ADLS_CREDENTIAL = "adls.credential"
ADLS_ACCOUNT_NAME = "adls.account-name"
ADLS_ACCOUNT_KEY = "adls.account-key"
ADLS_SAS_TOKEN = "adls.sas-token"
ADLS_TENANT_ID = "adls.tenant-id"
ADLS_CLIENT_ID = "adls.client-id"
ADLS_CLIENT_SECRET = "adls.client-secret"
ADLS_ACCOUNT_HOST = "adls.account-host"
ADLS_BLOB_STORAGE_AUTHORITY = "adls.blob-storage-authority"
ADLS_DFS_STORAGE_AUTHORITY = "adls.dfs-storage-authority"
ADLS_BLOB_STORAGE_SCHEME = "adls.blob-storage-scheme"
ADLS_DFS_STORAGE_SCHEME = "adls.dfs-storage-scheme"
ADLS_TOKEN = "adls.token"
GCS_TOKEN = "gcs.oauth2.token"
GCS_TOKEN_EXPIRES_AT_MS = "gcs.oauth2.token-expires-at"
GCS_PROJECT_ID = "gcs.project-id"
GCS_ACCESS = "gcs.access"
GCS_CONSISTENCY = "gcs.consistency"
GCS_CACHE_TIMEOUT = "gcs.cache-timeout"
GCS_REQUESTER_PAYS = "gcs.requester-pays"
GCS_SESSION_KWARGS = "gcs.session-kwargs"
GCS_SERVICE_HOST = "gcs.service.host"
GCS_DEFAULT_LOCATION = "gcs.default-bucket-location"
GCS_VERSION_AWARE = "gcs.version-aware"
HF_ENDPOINT = "hf.endpoint"
HF_TOKEN = "hf.token"
PYARROW_USE_LARGE_TYPES_ON_READ = "pyarrow.use-large-types-on-read"


@runtime_checkable
class InputStream(Protocol):
    """A protocol for the file-like object returned by InputFile.open(...).

    This outlines the minimally required methods for a seekable input stream returned from an InputFile
    implementation's `open(...)` method. These methods are a subset of IOBase/RawIOBase.
    """

    @abstractmethod
    def read(self, size: int = 0) -> bytes: ...

    @abstractmethod
    def seek(self, offset: int, whence: int = SEEK_SET) -> int: ...

    @abstractmethod
    def tell(self) -> int: ...

    @abstractmethod
    def close(self) -> None: ...

    def __enter__(self) -> InputStream:
        """Provide setup when opening an InputStream using a 'with' statement."""

    @abstractmethod
    def __exit__(
        self, exctype: Optional[Type[BaseException]], excinst: Optional[BaseException], exctb: Optional[TracebackType]
    ) -> None:
        """Perform cleanup when exiting the scope of a 'with' statement."""


@runtime_checkable
class OutputStream(Protocol):  # pragma: no cover
    """A protocol for the file-like object returned by OutputFile.create(...).

    This outlines the minimally required methods for a writable output stream returned from an OutputFile
    implementation's `create(...)` method. These methods are a subset of IOBase/RawIOBase.
    """

    @abstractmethod
    def write(self, b: bytes) -> int: ...

    @abstractmethod
    def close(self) -> None: ...

    @abstractmethod
    def __enter__(self) -> OutputStream:
        """Provide setup when opening an OutputStream using a 'with' statement."""

    @abstractmethod
    def __exit__(
        self, exctype: Optional[Type[BaseException]], excinst: Optional[BaseException], exctb: Optional[TracebackType]
    ) -> None:
        """Perform cleanup when exiting the scope of a 'with' statement."""


class InputFile(ABC):
    """A base class for InputFile implementations.

    Args:
        location (str): A URI or a path to a local file.

    Attributes:
        location (str): The URI or path to a local file for an InputFile instance.
        exists (bool): Whether the file exists or not.
    """

    def __init__(self, location: str):
        self._location = location

    @abstractmethod
    def __len__(self) -> int:
        """Return the total length of the file, in bytes."""

    @property
    def location(self) -> str:
        """The fully-qualified location of the input file."""
        return self._location

    @abstractmethod
    def exists(self) -> bool:
        """Check whether the location exists.

        Raises:
            PermissionError: If the file at self.location cannot be accessed due to a permission error.
        """

    @abstractmethod
    def open(self, seekable: bool = True) -> InputStream:
        """Return an object that matches the InputStream protocol.

        Args:
            seekable: If the stream should support seek, or if it is consumed sequential.

        Returns:
            InputStream: An object that matches the InputStream protocol.

        Raises:
            PermissionError: If the file at self.location cannot be accessed due to a permission error.
            FileNotFoundError: If the file at self.location does not exist.
        """


class OutputFile(ABC):
    """A base class for OutputFile implementations.

    Args:
        location (str): A URI or a path to a local file.

    Attributes:
        location (str): The URI or path to a local file for an OutputFile instance.
        exists (bool): Whether the file exists or not.
    """

    def __init__(self, location: str):
        self._location = location

    @abstractmethod
    def __len__(self) -> int:
        """Return the total length of the file, in bytes."""

    @property
    def location(self) -> str:
        """The fully-qualified location of the output file."""
        return self._location

    @abstractmethod
    def exists(self) -> bool:
        """Check whether the location exists.

        Raises:
            PermissionError: If the file at self.location cannot be accessed due to a permission error.
        """

    @abstractmethod
    def to_input_file(self) -> InputFile:
        """Return an InputFile for the location of this output file."""

    @abstractmethod
    def create(self, overwrite: bool = False) -> OutputStream:
        """Return an object that matches the OutputStream protocol.

        Args:
            overwrite (bool): If the file already exists at `self.location`
                and `overwrite` is False a FileExistsError should be raised.

        Returns:
            OutputStream: An object that matches the OutputStream protocol.

        Raises:
            PermissionError: If the file at self.location cannot be accessed due to a permission error.
            FileExistsError: If the file at self.location already exists and `overwrite=False`.
        """


class FileIO(ABC):
    """A base class for FileIO implementations."""

    properties: Properties

    def __init__(self, properties: Properties = EMPTY_DICT):
        self.properties = properties

    @abstractmethod
    def new_input(self, location: str) -> InputFile:
        """Get an InputFile instance to read bytes from the file at the given location.

        Args:
            location (str): A URI or a path to a local file.
        """

    @abstractmethod
    def new_output(self, location: str) -> OutputFile:
        """Get an OutputFile instance to write bytes to the file at the given location.

        Args:
            location (str): A URI or a path to a local file.
        """

    @abstractmethod
    def delete(self, location: Union[str, InputFile, OutputFile]) -> None:
        """Delete the file at the given path.

        Args:
            location (Union[str, InputFile, OutputFile]): A URI or a path to a local file--if an InputFile instance or
                an OutputFile instance is provided, the location attribute for that instance is used as the URI to delete.

        Raises:
            PermissionError: If the file at location cannot be accessed due to a permission error.
            FileNotFoundError: When the file at the provided location does not exist.
        """


LOCATION = "location"
WAREHOUSE = "warehouse"

ARROW_FILE_IO = "pyiceberg.io.pyarrow.PyArrowFileIO"
FSSPEC_FILE_IO = "pyiceberg.io.fsspec.FsspecFileIO"

# Mappings from the Java FileIO impl to a Python one. The list is ordered by preference.
# If an implementation isn't installed, it will fall back to the next one.
SCHEMA_TO_FILE_IO: Dict[str, List[str]] = {
    "s3": [ARROW_FILE_IO, FSSPEC_FILE_IO],
    "s3a": [ARROW_FILE_IO, FSSPEC_FILE_IO],
    "s3n": [ARROW_FILE_IO, FSSPEC_FILE_IO],
    "oss": [ARROW_FILE_IO],
    "gs": [ARROW_FILE_IO],
    "file": [ARROW_FILE_IO, FSSPEC_FILE_IO],
    "hdfs": [ARROW_FILE_IO],
    "viewfs": [ARROW_FILE_IO],
    "abfs": [FSSPEC_FILE_IO, ARROW_FILE_IO],
    "abfss": [FSSPEC_FILE_IO, ARROW_FILE_IO],
    "wasb": [FSSPEC_FILE_IO, ARROW_FILE_IO],
    "wasbs": [FSSPEC_FILE_IO, ARROW_FILE_IO],
    "hf": [FSSPEC_FILE_IO],
}


def _import_file_io(io_impl: str, properties: Properties) -> Optional[FileIO]:
    try:
        path_parts = io_impl.split(".")
        if len(path_parts) < 2:
            raise ValueError(f"py-io-impl should be full path (module.CustomFileIO), got: {io_impl}")
        module_name, class_name = ".".join(path_parts[:-1]), path_parts[-1]
        module = importlib.import_module(module_name)
        class_ = getattr(module, class_name)
        return class_(properties)
    except ModuleNotFoundError as exc:
        logger.warning(f"Could not initialize FileIO: {io_impl}", exc_info=exc)
        return None


PY_IO_IMPL = "py-io-impl"


def _infer_file_io_from_scheme(path: str, properties: Properties) -> Optional[FileIO]:
    parsed_url = urlparse(path)
    if parsed_url.scheme:
        if file_ios := SCHEMA_TO_FILE_IO.get(parsed_url.scheme):
            for file_io_path in file_ios:
                if file_io := _import_file_io(file_io_path, properties):
                    return file_io
        else:
            warnings.warn(f"No preferred file implementation for scheme: {parsed_url.scheme}")
    return None


def load_file_io(properties: Properties = EMPTY_DICT, location: Optional[str] = None) -> FileIO:
    # First look for the py-io-impl property to directly load the class
    if io_impl := properties.get(PY_IO_IMPL):
        if file_io := _import_file_io(io_impl, properties):
            logger.info("Loaded FileIO: %s", io_impl)
            return file_io
        else:
            raise ValueError(f"Could not initialize FileIO: {io_impl}")

    # Check the table location
    if location:
        if file_io := _infer_file_io_from_scheme(location, properties):
            return file_io

    # Look at the schema of the warehouse
    if warehouse_location := properties.get(WAREHOUSE):
        if file_io := _infer_file_io_from_scheme(warehouse_location, properties):
            return file_io

    try:
        # Default to PyArrow
        logger.info("Defaulting to PyArrow FileIO")
        from pyiceberg.io.pyarrow import PyArrowFileIO

        return PyArrowFileIO(properties)
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            'Could not load a FileIO, please consider installing one: pip3 install "pyiceberg[pyarrow]", for more options refer to the docs.'
        ) from e
