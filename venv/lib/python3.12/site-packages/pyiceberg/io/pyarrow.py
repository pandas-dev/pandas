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
# pylint: disable=redefined-outer-name,arguments-renamed,fixme
"""FileIO implementation for reading and writing table files that uses pyarrow.fs.

This file contains a FileIO implementation that relies on the filesystem interface provided
by PyArrow. It relies on PyArrow's `from_uri` method that infers the correct filesystem
type to use. Theoretically, this allows the supported storage types to grow naturally
with the pyarrow library.
"""

from __future__ import annotations

import fnmatch
import functools
import importlib
import itertools
import logging
import operator
import os
import re
import uuid
import warnings
from abc import ABC, abstractmethod
from copy import copy
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache, singledispatch
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Generic,
    Iterable,
    Iterator,
    List,
    Optional,
    Set,
    Tuple,
    TypeVar,
    Union,
    cast,
)
from urllib.parse import urlparse

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.lib
import pyarrow.parquet as pq
from pyarrow import ChunkedArray
from pyarrow._s3fs import S3RetryStrategy
from pyarrow.fs import (
    FileInfo,
    FileSystem,
    FileType,
)

from pyiceberg.conversions import to_bytes
from pyiceberg.exceptions import ResolveError
from pyiceberg.expressions import AlwaysTrue, BooleanExpression, BoundIsNaN, BoundIsNull, BoundTerm, Not, Or
from pyiceberg.expressions.literals import Literal
from pyiceberg.expressions.visitors import (
    BoundBooleanExpressionVisitor,
    bind,
    extract_field_ids,
    translate_column_names,
)
from pyiceberg.expressions.visitors import visit as boolean_expression_visit
from pyiceberg.io import (
    ADLS_ACCOUNT_KEY,
    ADLS_ACCOUNT_NAME,
    ADLS_BLOB_STORAGE_AUTHORITY,
    ADLS_BLOB_STORAGE_SCHEME,
    ADLS_CLIENT_ID,
    ADLS_CLIENT_SECRET,
    ADLS_DFS_STORAGE_AUTHORITY,
    ADLS_DFS_STORAGE_SCHEME,
    ADLS_SAS_TOKEN,
    ADLS_TENANT_ID,
    AWS_ACCESS_KEY_ID,
    AWS_REGION,
    AWS_ROLE_ARN,
    AWS_ROLE_SESSION_NAME,
    AWS_SECRET_ACCESS_KEY,
    AWS_SESSION_TOKEN,
    GCS_DEFAULT_LOCATION,
    GCS_SERVICE_HOST,
    GCS_TOKEN,
    GCS_TOKEN_EXPIRES_AT_MS,
    HDFS_HOST,
    HDFS_KERB_TICKET,
    HDFS_PORT,
    HDFS_USER,
    PYARROW_USE_LARGE_TYPES_ON_READ,
    S3_ACCESS_KEY_ID,
    S3_ANONYMOUS,
    S3_CONNECT_TIMEOUT,
    S3_ENDPOINT,
    S3_FORCE_VIRTUAL_ADDRESSING,
    S3_PROXY_URI,
    S3_REGION,
    S3_REQUEST_TIMEOUT,
    S3_RESOLVE_REGION,
    S3_RETRY_STRATEGY_IMPL,
    S3_ROLE_ARN,
    S3_ROLE_SESSION_NAME,
    S3_SECRET_ACCESS_KEY,
    S3_SESSION_TOKEN,
    FileIO,
    InputFile,
    InputStream,
    OutputFile,
    OutputStream,
)
from pyiceberg.manifest import (
    DataFile,
    DataFileContent,
    FileFormat,
)
from pyiceberg.partitioning import PartitionField, PartitionFieldValue, PartitionKey, PartitionSpec, partition_record_value
from pyiceberg.schema import (
    PartnerAccessor,
    PreOrderSchemaVisitor,
    Schema,
    SchemaVisitorPerPrimitiveType,
    SchemaWithPartnerVisitor,
    _check_schema_compatible,
    build_position_accessors,
    pre_order_visit,
    promote,
    prune_columns,
    sanitize_column_names,
    visit,
    visit_with_partner,
)
from pyiceberg.table import DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE, TableProperties
from pyiceberg.table.locations import load_location_provider
from pyiceberg.table.metadata import TableMetadata
from pyiceberg.table.name_mapping import NameMapping, apply_name_mapping
from pyiceberg.table.puffin import PuffinFile
from pyiceberg.transforms import IdentityTransform, TruncateTransform
from pyiceberg.typedef import EMPTY_DICT, Properties, Record, TableVersion
from pyiceberg.types import (
    BinaryType,
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    FixedType,
    FloatType,
    IcebergType,
    IntegerType,
    ListType,
    LongType,
    MapType,
    NestedField,
    PrimitiveType,
    StringType,
    StructType,
    TimestampNanoType,
    TimestampType,
    TimestamptzNanoType,
    TimestamptzType,
    TimeType,
    UnknownType,
    UUIDType,
    strtobool,
)
from pyiceberg.utils.concurrent import ExecutorFactory
from pyiceberg.utils.config import Config
from pyiceberg.utils.datetime import millis_to_datetime
from pyiceberg.utils.decimal import unscaled_to_decimal
from pyiceberg.utils.deprecated import deprecation_message
from pyiceberg.utils.properties import get_first_property_value, property_as_bool, property_as_int
from pyiceberg.utils.singleton import Singleton
from pyiceberg.utils.truncate import truncate_upper_bound_binary_string, truncate_upper_bound_text_string

if TYPE_CHECKING:
    from pyiceberg.table import FileScanTask, WriteTask

logger = logging.getLogger(__name__)

ONE_MEGABYTE = 1024 * 1024
BUFFER_SIZE = "buffer-size"
ICEBERG_SCHEMA = b"iceberg.schema"
# The PARQUET: in front means that it is Parquet specific, in this case the field_id
PYARROW_PARQUET_FIELD_ID_KEY = b"PARQUET:field_id"
PYARROW_FIELD_DOC_KEY = b"doc"
LIST_ELEMENT_NAME = "element"
MAP_KEY_NAME = "key"
MAP_VALUE_NAME = "value"
DOC = "doc"
UTC_ALIASES = {"UTC", "+00:00", "Etc/UTC", "Z"}

T = TypeVar("T")


@lru_cache
def _cached_resolve_s3_region(bucket: str) -> Optional[str]:
    from pyarrow.fs import resolve_s3_region

    try:
        return resolve_s3_region(bucket=bucket)
    except (OSError, TypeError):
        logger.warning(f"Unable to resolve region for bucket {bucket}")
        return None


def _import_retry_strategy(impl: str) -> Optional[S3RetryStrategy]:
    try:
        path_parts = impl.split(".")
        if len(path_parts) < 2:
            raise ValueError(f"retry-strategy-impl should be full path (module.CustomS3RetryStrategy), got: {impl}")
        module_name, class_name = ".".join(path_parts[:-1]), path_parts[-1]
        module = importlib.import_module(module_name)
        class_ = getattr(module, class_name)
        return class_()
    except (ModuleNotFoundError, AttributeError):
        warnings.warn(f"Could not initialize S3 retry strategy: {impl}")
        return None


class UnsupportedPyArrowTypeException(Exception):
    """Cannot convert PyArrow type to corresponding Iceberg type."""

    def __init__(self, field: pa.Field, *args: Any):
        self.field = field
        super().__init__(*args)


class PyArrowLocalFileSystem(pyarrow.fs.LocalFileSystem):
    def open_output_stream(self, path: str, *args: Any, **kwargs: Any) -> pyarrow.NativeFile:
        # In LocalFileSystem, parent directories must be first created before opening an output stream
        self.create_dir(os.path.dirname(path), recursive=True)
        return super().open_output_stream(path, *args, **kwargs)


class PyArrowFile(InputFile, OutputFile):
    """A combined InputFile and OutputFile implementation that uses a pyarrow filesystem to generate pyarrow.lib.NativeFile instances.

    Args:
        location (str): A URI or a path to a local file.

    Attributes:
        location(str): The URI or path to a local file for a PyArrowFile instance.

    Examples:
        >>> from pyiceberg.io.pyarrow import PyArrowFile
        >>> # input_file = PyArrowFile("s3://foo/bar.txt")
        >>> # Read the contents of the PyArrowFile instance
        >>> # Make sure that you have permissions to read/write
        >>> # file_content = input_file.open().read()

        >>> # output_file = PyArrowFile("s3://baz/qux.txt")
        >>> # Write bytes to a file
        >>> # Make sure that you have permissions to read/write
        >>> # output_file.create().write(b'foobytes')
    """

    _filesystem: FileSystem
    _path: str
    _buffer_size: int

    def __init__(self, location: str, path: str, fs: FileSystem, buffer_size: int = ONE_MEGABYTE):
        self._filesystem = fs
        self._path = path
        self._buffer_size = buffer_size
        super().__init__(location=location)

    def _file_info(self) -> FileInfo:
        """Retrieve a pyarrow.fs.FileInfo object for the location.

        Raises:
            PermissionError: If the file at self.location cannot be accessed due to a permission error such as
                an AWS error code 15.
        """
        try:
            file_info = self._filesystem.get_file_info(self._path)
        except OSError as e:
            if e.errno == 13 or "AWS Error [code 15]" in str(e):
                raise PermissionError(f"Cannot get file info, access denied: {self.location}") from e
            raise  # pragma: no cover - If some other kind of OSError, raise the raw error

        if file_info.type == FileType.NotFound:
            raise FileNotFoundError(f"Cannot get file info, file not found: {self.location}")
        return file_info

    def __len__(self) -> int:
        """Return the total length of the file, in bytes."""
        file_info = self._file_info()
        return file_info.size

    def exists(self) -> bool:
        """Check whether the location exists."""
        try:
            self._file_info()  # raises FileNotFoundError if it does not exist
            return True
        except FileNotFoundError:
            return False

    def open(self, seekable: bool = True) -> InputStream:
        """Open the location using a PyArrow FileSystem inferred from the location.

        Args:
            seekable: If the stream should support seek, or if it is consumed sequential.

        Returns:
            pyarrow.lib.NativeFile: A NativeFile instance for the file located at `self.location`.

        Raises:
            FileNotFoundError: If the file at self.location does not exist.
            PermissionError: If the file at self.location cannot be accessed due to a permission error such as
                an AWS error code 15.
        """
        try:
            if seekable:
                input_file = self._filesystem.open_input_file(self._path)
            else:
                input_file = self._filesystem.open_input_stream(self._path, buffer_size=self._buffer_size)
        except (FileNotFoundError, PermissionError):
            raise
        except OSError as e:
            if e.errno == 2 or "Path does not exist" in str(e):
                raise FileNotFoundError(f"Cannot open file, does not exist: {self.location}") from e
            elif e.errno == 13 or "AWS Error [code 15]" in str(e):
                raise PermissionError(f"Cannot open file, access denied: {self.location}") from e
            raise  # pragma: no cover - If some other kind of OSError, raise the raw error
        return input_file

    def create(self, overwrite: bool = False) -> OutputStream:
        """Create a writable pyarrow.lib.NativeFile for this PyArrowFile's location.

        Args:
            overwrite (bool): Whether to overwrite the file if it already exists.

        Returns:
            pyarrow.lib.NativeFile: A NativeFile instance for the file located at self.location.

        Raises:
            FileExistsError: If the file already exists at `self.location` and `overwrite` is False.

        Note:
            This retrieves a pyarrow NativeFile by opening an output stream. If overwrite is set to False,
            a check is first performed to verify that the file does not exist. This is not thread-safe and
            a possibility does exist that the file can be created by a concurrent process after the existence
            check yet before the output stream is created. In such a case, the default pyarrow behavior will
            truncate the contents of the existing file when opening the output stream.
        """
        try:
            if not overwrite and self.exists() is True:
                raise FileExistsError(f"Cannot create file, already exists: {self.location}")
            output_file = self._filesystem.open_output_stream(self._path, buffer_size=self._buffer_size)
        except PermissionError:
            raise
        except OSError as e:
            if e.errno == 13 or "AWS Error [code 15]" in str(e):
                raise PermissionError(f"Cannot create file, access denied: {self.location}") from e
            raise  # pragma: no cover - If some other kind of OSError, raise the raw error
        return output_file

    def to_input_file(self) -> PyArrowFile:
        """Return a new PyArrowFile for the location of an existing PyArrowFile instance.

        This method is included to abide by the OutputFile abstract base class. Since this implementation uses a single
        PyArrowFile class (as opposed to separate InputFile and OutputFile implementations), this method effectively returns
        a copy of the same instance.
        """
        return self


class PyArrowFileIO(FileIO):
    fs_by_scheme: Callable[[str, Optional[str]], FileSystem]

    def __init__(self, properties: Properties = EMPTY_DICT):
        self.fs_by_scheme: Callable[[str, Optional[str]], FileSystem] = lru_cache(self._initialize_fs)
        super().__init__(properties=properties)

    @staticmethod
    def parse_location(location: str, properties: Properties = EMPTY_DICT) -> Tuple[str, str, str]:
        """Return (scheme, netloc, path) for the given location.

        Uses DEFAULT_SCHEME and DEFAULT_NETLOC if scheme/netloc are missing.
        """
        uri = urlparse(location)

        if not uri.scheme:
            default_scheme = properties.get("DEFAULT_SCHEME", "file")
            default_netloc = properties.get("DEFAULT_NETLOC", "")
            return default_scheme, default_netloc, os.path.abspath(location)
        elif uri.scheme in ("hdfs", "viewfs"):
            return uri.scheme, uri.netloc, uri.path
        else:
            return uri.scheme, uri.netloc, f"{uri.netloc}{uri.path}"

    def _initialize_fs(self, scheme: str, netloc: Optional[str] = None) -> FileSystem:
        """Initialize FileSystem for different scheme."""
        if scheme in {"oss"}:
            return self._initialize_oss_fs()

        elif scheme in {"s3", "s3a", "s3n"}:
            return self._initialize_s3_fs(netloc)

        elif scheme in {"hdfs", "viewfs"}:
            return self._initialize_hdfs_fs(scheme, netloc)

        elif scheme in {"gs", "gcs"}:
            return self._initialize_gcs_fs()

        elif scheme in {"abfs", "abfss", "wasb", "wasbs"}:
            return self._initialize_azure_fs()

        elif scheme in {"file"}:
            return self._initialize_local_fs()

        else:
            raise ValueError(f"Unrecognized filesystem type in URI: {scheme}")

    def _initialize_oss_fs(self) -> FileSystem:
        from pyarrow.fs import S3FileSystem

        client_kwargs: Dict[str, Any] = {
            "endpoint_override": self.properties.get(S3_ENDPOINT),
            "access_key": get_first_property_value(self.properties, S3_ACCESS_KEY_ID, AWS_ACCESS_KEY_ID),
            "secret_key": get_first_property_value(self.properties, S3_SECRET_ACCESS_KEY, AWS_SECRET_ACCESS_KEY),
            "session_token": get_first_property_value(self.properties, S3_SESSION_TOKEN, AWS_SESSION_TOKEN),
            "region": get_first_property_value(self.properties, S3_REGION, AWS_REGION),
            "force_virtual_addressing": property_as_bool(self.properties, S3_FORCE_VIRTUAL_ADDRESSING, True),
        }

        if proxy_uri := self.properties.get(S3_PROXY_URI):
            client_kwargs["proxy_options"] = proxy_uri

        if connect_timeout := self.properties.get(S3_CONNECT_TIMEOUT):
            client_kwargs["connect_timeout"] = float(connect_timeout)

        if request_timeout := self.properties.get(S3_REQUEST_TIMEOUT):
            client_kwargs["request_timeout"] = float(request_timeout)

        if role_arn := get_first_property_value(self.properties, S3_ROLE_ARN, AWS_ROLE_ARN):
            client_kwargs["role_arn"] = role_arn

        if session_name := get_first_property_value(self.properties, S3_ROLE_SESSION_NAME, AWS_ROLE_SESSION_NAME):
            client_kwargs["session_name"] = session_name

        if s3_anonymous := self.properties.get(S3_ANONYMOUS):
            client_kwargs["anonymous"] = strtobool(s3_anonymous)

        return S3FileSystem(**client_kwargs)

    def _initialize_s3_fs(self, netloc: Optional[str]) -> FileSystem:
        from pyarrow.fs import S3FileSystem

        provided_region = get_first_property_value(self.properties, S3_REGION, AWS_REGION)

        # Do this when we don't provide the region at all, or when we explicitly enable it
        if provided_region is None or property_as_bool(self.properties, S3_RESOLVE_REGION, False) is True:
            # Resolve region from netloc(bucket), fallback to user-provided region
            # Only supported by buckets hosted by S3
            bucket_region = _cached_resolve_s3_region(bucket=netloc) or provided_region
            if provided_region is not None and bucket_region != provided_region:
                logger.warning(
                    f"PyArrow FileIO overriding S3 bucket region for bucket {netloc}: "
                    f"provided region {provided_region}, actual region {bucket_region}"
                )
        else:
            bucket_region = provided_region

        client_kwargs: Dict[str, Any] = {
            "endpoint_override": self.properties.get(S3_ENDPOINT),
            "access_key": get_first_property_value(self.properties, S3_ACCESS_KEY_ID, AWS_ACCESS_KEY_ID),
            "secret_key": get_first_property_value(self.properties, S3_SECRET_ACCESS_KEY, AWS_SECRET_ACCESS_KEY),
            "session_token": get_first_property_value(self.properties, S3_SESSION_TOKEN, AWS_SESSION_TOKEN),
            "region": bucket_region,
        }

        if proxy_uri := self.properties.get(S3_PROXY_URI):
            client_kwargs["proxy_options"] = proxy_uri

        if connect_timeout := self.properties.get(S3_CONNECT_TIMEOUT):
            client_kwargs["connect_timeout"] = float(connect_timeout)

        if request_timeout := self.properties.get(S3_REQUEST_TIMEOUT):
            client_kwargs["request_timeout"] = float(request_timeout)

        if role_arn := get_first_property_value(self.properties, S3_ROLE_ARN, AWS_ROLE_ARN):
            client_kwargs["role_arn"] = role_arn

        if session_name := get_first_property_value(self.properties, S3_ROLE_SESSION_NAME, AWS_ROLE_SESSION_NAME):
            client_kwargs["session_name"] = session_name

        if self.properties.get(S3_FORCE_VIRTUAL_ADDRESSING) is not None:
            client_kwargs["force_virtual_addressing"] = property_as_bool(self.properties, S3_FORCE_VIRTUAL_ADDRESSING, False)

        if (retry_strategy_impl := self.properties.get(S3_RETRY_STRATEGY_IMPL)) and (
            retry_instance := _import_retry_strategy(retry_strategy_impl)
        ):
            client_kwargs["retry_strategy"] = retry_instance

        if s3_anonymous := self.properties.get(S3_ANONYMOUS):
            client_kwargs["anonymous"] = strtobool(s3_anonymous)

        return S3FileSystem(**client_kwargs)

    def _initialize_azure_fs(self) -> FileSystem:
        # https://arrow.apache.org/docs/python/generated/pyarrow.fs.AzureFileSystem.html
        from packaging import version

        MIN_PYARROW_VERSION_SUPPORTING_AZURE_FS = "20.0.0"
        if version.parse(pyarrow.__version__) < version.parse(MIN_PYARROW_VERSION_SUPPORTING_AZURE_FS):
            raise ImportError(
                f"pyarrow version >= {MIN_PYARROW_VERSION_SUPPORTING_AZURE_FS} required for AzureFileSystem support, "
                f"but found version {pyarrow.__version__}."
            )

        from pyarrow.fs import AzureFileSystem

        client_kwargs: Dict[str, str] = {}

        if account_name := self.properties.get(ADLS_ACCOUNT_NAME):
            client_kwargs["account_name"] = account_name

        if account_key := self.properties.get(ADLS_ACCOUNT_KEY):
            client_kwargs["account_key"] = account_key

        if blob_storage_authority := self.properties.get(ADLS_BLOB_STORAGE_AUTHORITY):
            client_kwargs["blob_storage_authority"] = blob_storage_authority

        if dfs_storage_authority := self.properties.get(ADLS_DFS_STORAGE_AUTHORITY):
            client_kwargs["dfs_storage_authority"] = dfs_storage_authority

        if blob_storage_scheme := self.properties.get(ADLS_BLOB_STORAGE_SCHEME):
            client_kwargs["blob_storage_scheme"] = blob_storage_scheme

        if dfs_storage_scheme := self.properties.get(ADLS_DFS_STORAGE_SCHEME):
            client_kwargs["dfs_storage_scheme"] = dfs_storage_scheme

        if sas_token := self.properties.get(ADLS_SAS_TOKEN):
            client_kwargs["sas_token"] = sas_token

        if client_id := self.properties.get(ADLS_CLIENT_ID):
            client_kwargs["client_id"] = client_id
        if client_secret := self.properties.get(ADLS_CLIENT_SECRET):
            client_kwargs["client_secret"] = client_secret
        if tenant_id := self.properties.get(ADLS_TENANT_ID):
            client_kwargs["tenant_id"] = tenant_id

        # Validate that all three are provided together for ClientSecretCredential
        credential_keys = ["client_id", "client_secret", "tenant_id"]
        provided_keys = [key for key in credential_keys if key in client_kwargs]
        if provided_keys and len(provided_keys) != len(credential_keys):
            missing_keys = [key for key in credential_keys if key not in client_kwargs]
            raise ValueError(
                f"client_id, client_secret, and tenant_id must all be provided together "
                f"to use ClientSecretCredential for Azure authentication. "
                f"Provided: {provided_keys}, Missing: {missing_keys}"
            )

        return AzureFileSystem(**client_kwargs)

    def _initialize_hdfs_fs(self, scheme: str, netloc: Optional[str]) -> FileSystem:
        from pyarrow.fs import HadoopFileSystem

        hdfs_kwargs: Dict[str, Any] = {}
        if netloc:
            return HadoopFileSystem.from_uri(f"{scheme}://{netloc}")
        if host := self.properties.get(HDFS_HOST):
            hdfs_kwargs["host"] = host
        if port := self.properties.get(HDFS_PORT):
            # port should be an integer type
            hdfs_kwargs["port"] = int(port)
        if user := self.properties.get(HDFS_USER):
            hdfs_kwargs["user"] = user
        if kerb_ticket := self.properties.get(HDFS_KERB_TICKET):
            hdfs_kwargs["kerb_ticket"] = kerb_ticket

        return HadoopFileSystem(**hdfs_kwargs)

    def _initialize_gcs_fs(self) -> FileSystem:
        from pyarrow.fs import GcsFileSystem

        gcs_kwargs: Dict[str, Any] = {}
        if access_token := self.properties.get(GCS_TOKEN):
            gcs_kwargs["access_token"] = access_token
        if expiration := self.properties.get(GCS_TOKEN_EXPIRES_AT_MS):
            gcs_kwargs["credential_token_expiration"] = millis_to_datetime(int(expiration))
        if bucket_location := self.properties.get(GCS_DEFAULT_LOCATION):
            gcs_kwargs["default_bucket_location"] = bucket_location
        if endpoint := self.properties.get(GCS_SERVICE_HOST):
            url_parts = urlparse(endpoint)
            gcs_kwargs["scheme"] = url_parts.scheme
            gcs_kwargs["endpoint_override"] = url_parts.netloc

        return GcsFileSystem(**gcs_kwargs)

    def _initialize_local_fs(self) -> FileSystem:
        return PyArrowLocalFileSystem()

    def new_input(self, location: str) -> PyArrowFile:
        """Get a PyArrowFile instance to read bytes from the file at the given location.

        Args:
            location (str): A URI or a path to a local file.

        Returns:
            PyArrowFile: A PyArrowFile instance for the given location.
        """
        scheme, netloc, path = self.parse_location(location, self.properties)
        return PyArrowFile(
            fs=self.fs_by_scheme(scheme, netloc),
            location=location,
            path=path,
            buffer_size=int(self.properties.get(BUFFER_SIZE, ONE_MEGABYTE)),
        )

    def new_output(self, location: str) -> PyArrowFile:
        """Get a PyArrowFile instance to write bytes to the file at the given location.

        Args:
            location (str): A URI or a path to a local file.

        Returns:
            PyArrowFile: A PyArrowFile instance for the given location.
        """
        scheme, netloc, path = self.parse_location(location, self.properties)
        return PyArrowFile(
            fs=self.fs_by_scheme(scheme, netloc),
            location=location,
            path=path,
            buffer_size=int(self.properties.get(BUFFER_SIZE, ONE_MEGABYTE)),
        )

    def delete(self, location: Union[str, InputFile, OutputFile]) -> None:
        """Delete the file at the given location.

        Args:
            location (Union[str, InputFile, OutputFile]): The URI to the file--if an InputFile instance or an OutputFile instance is provided,
                the location attribute for that instance is used as the location to delete.

        Raises:
            FileNotFoundError: When the file at the provided location does not exist.
            PermissionError: If the file at the provided location cannot be accessed due to a permission error such as
                an AWS error code 15.
        """
        str_location = location.location if isinstance(location, (InputFile, OutputFile)) else location
        scheme, netloc, path = self.parse_location(str_location, self.properties)
        fs = self.fs_by_scheme(scheme, netloc)

        try:
            fs.delete_file(path)
        except FileNotFoundError:
            raise
        except PermissionError:
            raise
        except OSError as e:
            if e.errno == 2 or "Path does not exist" in str(e):
                raise FileNotFoundError(f"Cannot delete file, does not exist: {location}") from e
            elif e.errno == 13 or "AWS Error [code 15]" in str(e):
                raise PermissionError(f"Cannot delete file, access denied: {location}") from e
            raise  # pragma: no cover - If some other kind of OSError, raise the raw error

    def __getstate__(self) -> Dict[str, Any]:
        """Create a dictionary of the PyArrowFileIO fields used when pickling."""
        fileio_copy = copy(self.__dict__)
        fileio_copy["fs_by_scheme"] = None
        return fileio_copy

    def __setstate__(self, state: Dict[str, Any]) -> None:
        """Deserialize the state into a PyArrowFileIO instance."""
        self.__dict__ = state
        self.fs_by_scheme = lru_cache(self._initialize_fs)


def schema_to_pyarrow(
    schema: Union[Schema, IcebergType],
    metadata: Dict[bytes, bytes] = EMPTY_DICT,
    include_field_ids: bool = True,
) -> pa.schema:
    return visit(schema, _ConvertToArrowSchema(metadata, include_field_ids))


class _ConvertToArrowSchema(SchemaVisitorPerPrimitiveType[pa.DataType]):
    _metadata: Dict[bytes, bytes]

    def __init__(self, metadata: Dict[bytes, bytes] = EMPTY_DICT, include_field_ids: bool = True) -> None:
        self._metadata = metadata
        self._include_field_ids = include_field_ids

    def schema(self, _: Schema, struct_result: pa.StructType) -> pa.schema:
        return pa.schema(list(struct_result), metadata=self._metadata)

    def struct(self, _: StructType, field_results: List[pa.DataType]) -> pa.DataType:
        return pa.struct(field_results)

    def field(self, field: NestedField, field_result: pa.DataType) -> pa.Field:
        metadata = {}
        if field.doc:
            metadata[PYARROW_FIELD_DOC_KEY] = field.doc
        if self._include_field_ids:
            metadata[PYARROW_PARQUET_FIELD_ID_KEY] = str(field.field_id)

        return pa.field(
            name=field.name,
            type=field_result,
            nullable=field.optional,
            metadata=metadata,
        )

    def list(self, list_type: ListType, element_result: pa.DataType) -> pa.DataType:
        element_field = self.field(list_type.element_field, element_result)
        return pa.large_list(value_type=element_field)

    def map(self, map_type: MapType, key_result: pa.DataType, value_result: pa.DataType) -> pa.DataType:
        key_field = self.field(map_type.key_field, key_result)
        value_field = self.field(map_type.value_field, value_result)
        return pa.map_(key_type=key_field, item_type=value_field)

    def visit_fixed(self, fixed_type: FixedType) -> pa.DataType:
        return pa.binary(len(fixed_type))

    def visit_decimal(self, decimal_type: DecimalType) -> pa.DataType:
        # It looks like decimal{32,64} is not fully implemented:
        # https://github.com/apache/arrow/issues/25483
        # https://github.com/apache/arrow/issues/43956
        # However, if we keep it as 128 in memory, and based on the
        # precision/scale Arrow will map it to INT{32,64}
        # https://github.com/apache/arrow/blob/598938711a8376cbfdceaf5c77ab0fd5057e6c02/cpp/src/parquet/arrow/schema.cc#L380-L392
        return pa.decimal128(decimal_type.precision, decimal_type.scale)

    def visit_boolean(self, _: BooleanType) -> pa.DataType:
        return pa.bool_()

    def visit_integer(self, _: IntegerType) -> pa.DataType:
        return pa.int32()

    def visit_long(self, _: LongType) -> pa.DataType:
        return pa.int64()

    def visit_float(self, _: FloatType) -> pa.DataType:
        # 32-bit IEEE 754 floating point
        return pa.float32()

    def visit_double(self, _: DoubleType) -> pa.DataType:
        # 64-bit IEEE 754 floating point
        return pa.float64()

    def visit_date(self, _: DateType) -> pa.DataType:
        # Date encoded as an int
        return pa.date32()

    def visit_time(self, _: TimeType) -> pa.DataType:
        return pa.time64("us")

    def visit_timestamp(self, _: TimestampType) -> pa.DataType:
        return pa.timestamp(unit="us")

    def visit_timestamp_ns(self, _: TimestampNanoType) -> pa.DataType:
        return pa.timestamp(unit="ns")

    def visit_timestamptz(self, _: TimestamptzType) -> pa.DataType:
        return pa.timestamp(unit="us", tz="UTC")

    def visit_timestamptz_ns(self, _: TimestamptzNanoType) -> pa.DataType:
        return pa.timestamp(unit="ns", tz="UTC")

    def visit_string(self, _: StringType) -> pa.DataType:
        return pa.large_string()

    def visit_uuid(self, _: UUIDType) -> pa.DataType:
        return pa.uuid()

    def visit_unknown(self, _: UnknownType) -> pa.DataType:
        """Type `UnknownType` can be promoted to any primitive type in V3+ tables per the Iceberg spec."""
        return pa.null()

    def visit_binary(self, _: BinaryType) -> pa.DataType:
        return pa.large_binary()


def _convert_scalar(value: Any, iceberg_type: IcebergType) -> pa.scalar:
    if not isinstance(iceberg_type, PrimitiveType):
        raise ValueError(f"Expected primitive type, got: {iceberg_type}")
    return pa.scalar(value=value, type=schema_to_pyarrow(iceberg_type))


class _ConvertToArrowExpression(BoundBooleanExpressionVisitor[pc.Expression]):
    def visit_in(self, term: BoundTerm[Any], literals: Set[Any]) -> pc.Expression:
        pyarrow_literals = pa.array(literals, type=schema_to_pyarrow(term.ref().field.field_type))
        return pc.field(term.ref().field.name).isin(pyarrow_literals)

    def visit_not_in(self, term: BoundTerm[Any], literals: Set[Any]) -> pc.Expression:
        pyarrow_literals = pa.array(literals, type=schema_to_pyarrow(term.ref().field.field_type))
        return ~pc.field(term.ref().field.name).isin(pyarrow_literals)

    def visit_is_nan(self, term: BoundTerm[Any]) -> pc.Expression:
        ref = pc.field(term.ref().field.name)
        return pc.is_nan(ref)

    def visit_not_nan(self, term: BoundTerm[Any]) -> pc.Expression:
        ref = pc.field(term.ref().field.name)
        return ~pc.is_nan(ref)

    def visit_is_null(self, term: BoundTerm[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name).is_null(nan_is_null=False)

    def visit_not_null(self, term: BoundTerm[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name).is_valid()

    def visit_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name) == _convert_scalar(literal.value, term.ref().field.field_type)

    def visit_not_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name) != _convert_scalar(literal.value, term.ref().field.field_type)

    def visit_greater_than_or_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name) >= _convert_scalar(literal.value, term.ref().field.field_type)

    def visit_greater_than(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name) > _convert_scalar(literal.value, term.ref().field.field_type)

    def visit_less_than(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name) < _convert_scalar(literal.value, term.ref().field.field_type)

    def visit_less_than_or_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.field(term.ref().field.name) <= _convert_scalar(literal.value, term.ref().field.field_type)

    def visit_starts_with(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return pc.starts_with(pc.field(term.ref().field.name), literal.value)

    def visit_not_starts_with(self, term: BoundTerm[Any], literal: Literal[Any]) -> pc.Expression:
        return ~pc.starts_with(pc.field(term.ref().field.name), literal.value)

    def visit_true(self) -> pc.Expression:
        return pc.scalar(True)

    def visit_false(self) -> pc.Expression:
        return pc.scalar(False)

    def visit_not(self, child_result: pc.Expression) -> pc.Expression:
        return ~child_result

    def visit_and(self, left_result: pc.Expression, right_result: pc.Expression) -> pc.Expression:
        return left_result & right_result

    def visit_or(self, left_result: pc.Expression, right_result: pc.Expression) -> pc.Expression:
        return left_result | right_result


class _NullNaNUnmentionedTermsCollector(BoundBooleanExpressionVisitor[None]):
    # BoundTerms which have either is_null or is_not_null appearing at least once in the boolean expr.
    is_null_or_not_bound_terms: set[BoundTerm[Any]]
    # The remaining BoundTerms appearing in the boolean expr.
    null_unmentioned_bound_terms: set[BoundTerm[Any]]
    # BoundTerms which have either is_nan or is_not_nan appearing at least once in the boolean expr.
    is_nan_or_not_bound_terms: set[BoundTerm[Any]]
    # The remaining BoundTerms appearing in the boolean expr.
    nan_unmentioned_bound_terms: set[BoundTerm[Any]]

    def __init__(self) -> None:
        super().__init__()
        self.is_null_or_not_bound_terms = set()
        self.null_unmentioned_bound_terms = set()
        self.is_nan_or_not_bound_terms = set()
        self.nan_unmentioned_bound_terms = set()

    def _handle_explicit_is_null_or_not(self, term: BoundTerm[Any]) -> None:
        """Handle the predicate case where either is_null or is_not_null is included."""
        if term in self.null_unmentioned_bound_terms:
            self.null_unmentioned_bound_terms.remove(term)
        self.is_null_or_not_bound_terms.add(term)

    def _handle_null_unmentioned(self, term: BoundTerm[Any]) -> None:
        """Handle the predicate case where neither is_null or is_not_null is included."""
        if term not in self.is_null_or_not_bound_terms:
            self.null_unmentioned_bound_terms.add(term)

    def _handle_explicit_is_nan_or_not(self, term: BoundTerm[Any]) -> None:
        """Handle the predicate case where either is_nan or is_not_nan is included."""
        if term in self.nan_unmentioned_bound_terms:
            self.nan_unmentioned_bound_terms.remove(term)
        self.is_nan_or_not_bound_terms.add(term)

    def _handle_nan_unmentioned(self, term: BoundTerm[Any]) -> None:
        """Handle the predicate case where neither is_nan or is_not_nan is included."""
        if term not in self.is_nan_or_not_bound_terms:
            self.nan_unmentioned_bound_terms.add(term)

    def visit_in(self, term: BoundTerm[Any], literals: Set[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_not_in(self, term: BoundTerm[Any], literals: Set[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_is_nan(self, term: BoundTerm[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_explicit_is_nan_or_not(term)

    def visit_not_nan(self, term: BoundTerm[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_explicit_is_nan_or_not(term)

    def visit_is_null(self, term: BoundTerm[Any]) -> None:
        self._handle_explicit_is_null_or_not(term)
        self._handle_nan_unmentioned(term)

    def visit_not_null(self, term: BoundTerm[Any]) -> None:
        self._handle_explicit_is_null_or_not(term)
        self._handle_nan_unmentioned(term)

    def visit_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_not_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_greater_than_or_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_greater_than(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_less_than(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_less_than_or_equal(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_starts_with(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_not_starts_with(self, term: BoundTerm[Any], literal: Literal[Any]) -> None:
        self._handle_null_unmentioned(term)
        self._handle_nan_unmentioned(term)

    def visit_true(self) -> None:
        return

    def visit_false(self) -> None:
        return

    def visit_not(self, child_result: None) -> None:
        return

    def visit_and(self, left_result: None, right_result: None) -> None:
        return

    def visit_or(self, left_result: None, right_result: None) -> None:
        return

    def collect(
        self,
        expr: BooleanExpression,
    ) -> None:
        """Collect the bound references categorized by having at least one is_null or is_not_null in the expr and the remaining."""
        boolean_expression_visit(expr, self)


def expression_to_pyarrow(expr: BooleanExpression) -> pc.Expression:
    return boolean_expression_visit(expr, _ConvertToArrowExpression())


def _expression_to_complementary_pyarrow(expr: BooleanExpression) -> pc.Expression:
    """Complementary filter conversion function of expression_to_pyarrow.

    Could not use expression_to_pyarrow(Not(expr)) to achieve this complementary effect because ~ in pyarrow.compute.Expression does not handle null.
    """
    collector = _NullNaNUnmentionedTermsCollector()
    collector.collect(expr)

    # Convert the set of terms to a sorted list so that layout of the expression to build is deterministic.
    null_unmentioned_bound_terms: List[BoundTerm[Any]] = sorted(
        collector.null_unmentioned_bound_terms, key=lambda term: term.ref().field.name
    )
    nan_unmentioned_bound_terms: List[BoundTerm[Any]] = sorted(
        collector.nan_unmentioned_bound_terms, key=lambda term: term.ref().field.name
    )

    preserve_expr: BooleanExpression = Not(expr)
    for term in null_unmentioned_bound_terms:
        preserve_expr = Or(preserve_expr, BoundIsNull(term=term))
    for term in nan_unmentioned_bound_terms:
        preserve_expr = Or(preserve_expr, BoundIsNaN(term=term))
    return expression_to_pyarrow(preserve_expr)


@lru_cache
def _get_file_format(file_format: FileFormat, **kwargs: Dict[str, Any]) -> ds.FileFormat:
    if file_format == FileFormat.PARQUET:
        return ds.ParquetFileFormat(**kwargs)
    else:
        raise ValueError(f"Unsupported file format: {file_format}")


def _read_deletes(io: FileIO, data_file: DataFile) -> Dict[str, pa.ChunkedArray]:
    if data_file.file_format == FileFormat.PARQUET:
        with io.new_input(data_file.file_path).open() as fi:
            delete_fragment = _get_file_format(
                data_file.file_format, dictionary_columns=("file_path",), pre_buffer=True, buffer_size=ONE_MEGABYTE
            ).make_fragment(fi)
            table = ds.Scanner.from_fragment(fragment=delete_fragment).to_table()
        table = table.unify_dictionaries()
        return {
            file.as_py(): table.filter(pc.field("file_path") == file).column("pos")
            for file in table.column("file_path").chunks[0].dictionary
        }
    elif data_file.file_format == FileFormat.PUFFIN:
        with io.new_input(data_file.file_path).open() as fi:
            payload = fi.read()

        return PuffinFile(payload).to_vector()
    else:
        raise ValueError(f"Delete file format not supported: {data_file.file_format}")


def _combine_positional_deletes(positional_deletes: List[pa.ChunkedArray], start_index: int, end_index: int) -> pa.Array:
    if len(positional_deletes) == 1:
        all_chunks = positional_deletes[0]
    else:
        all_chunks = pa.chunked_array(itertools.chain(*[arr.chunks for arr in positional_deletes]))

    # Create the full range array with pyarrow
    full_range = pa.array(range(start_index, end_index))
    # When available, replace with Arrow generator to improve performance
    # See https://github.com/apache/iceberg-python/issues/1271 for details

    # Filter out values in all_chunks from full_range
    result = pc.filter(full_range, pc.invert(pc.is_in(full_range, value_set=all_chunks)))

    # Subtract the start_index from each element in the result
    return pc.subtract(result, pa.scalar(start_index))


def pyarrow_to_schema(
    schema: pa.Schema,
    name_mapping: Optional[NameMapping] = None,
    downcast_ns_timestamp_to_us: bool = False,
    format_version: TableVersion = TableProperties.DEFAULT_FORMAT_VERSION,
) -> Schema:
    has_ids = visit_pyarrow(schema, _HasIds())
    if has_ids:
        return visit_pyarrow(
            schema, _ConvertToIceberg(downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us, format_version=format_version)
        )
    elif name_mapping is not None:
        schema_without_ids = _pyarrow_to_schema_without_ids(
            schema, downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us, format_version=format_version
        )
        return apply_name_mapping(schema_without_ids, name_mapping)
    else:
        raise ValueError(
            "Parquet file does not have field-ids and the Iceberg table does not have 'schema.name-mapping.default' defined"
        )


def _pyarrow_to_schema_without_ids(
    schema: pa.Schema,
    downcast_ns_timestamp_to_us: bool = False,
    format_version: TableVersion = TableProperties.DEFAULT_FORMAT_VERSION,
) -> Schema:
    return visit_pyarrow(
        schema,
        _ConvertToIcebergWithoutIDs(downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us, format_version=format_version),
    )


def _pyarrow_schema_ensure_large_types(schema: pa.Schema) -> pa.Schema:
    return visit_pyarrow(schema, _ConvertToLargeTypes())


def _pyarrow_schema_ensure_small_types(schema: pa.Schema) -> pa.Schema:
    return visit_pyarrow(schema, _ConvertToSmallTypes())


@singledispatch
def visit_pyarrow(obj: Union[pa.DataType, pa.Schema], visitor: PyArrowSchemaVisitor[T]) -> T:
    """Apply a pyarrow schema visitor to any point within a schema.

    The function traverses the schema in post-order fashion.

    Args:
        obj (Union[pa.DataType, pa.Schema]): An instance of a Schema or an IcebergType.
        visitor (PyArrowSchemaVisitor[T]): An instance of an implementation of the generic PyarrowSchemaVisitor base class.

    Raises:
        NotImplementedError: If attempting to visit an unrecognized object type.
    """
    raise NotImplementedError(f"Cannot visit non-type: {obj}")


@visit_pyarrow.register(pa.Schema)
def _(obj: pa.Schema, visitor: PyArrowSchemaVisitor[T]) -> T:
    return visitor.schema(obj, visit_pyarrow(pa.struct(obj), visitor))


@visit_pyarrow.register(pa.StructType)
def _(obj: pa.StructType, visitor: PyArrowSchemaVisitor[T]) -> T:
    results = [visit_pyarrow(field, visitor) for field in obj]

    return visitor.struct(obj, results)


@visit_pyarrow.register(pa.ListType)
@visit_pyarrow.register(pa.FixedSizeListType)
@visit_pyarrow.register(pa.LargeListType)
def _(obj: Union[pa.ListType, pa.LargeListType, pa.FixedSizeListType], visitor: PyArrowSchemaVisitor[T]) -> T:
    visitor.before_list_element(obj.value_field)
    result = visit_pyarrow(obj.value_type, visitor)
    visitor.after_list_element(obj.value_field)
    return visitor.list(obj, result)


@visit_pyarrow.register(pa.MapType)
def _(obj: pa.MapType, visitor: PyArrowSchemaVisitor[T]) -> T:
    visitor.before_map_key(obj.key_field)
    key_result = visit_pyarrow(obj.key_type, visitor)
    visitor.after_map_key(obj.key_field)

    visitor.before_map_value(obj.item_field)
    value_result = visit_pyarrow(obj.item_type, visitor)
    visitor.after_map_value(obj.item_field)

    return visitor.map(obj, key_result, value_result)


@visit_pyarrow.register(pa.DictionaryType)
def _(obj: pa.DictionaryType, visitor: PyArrowSchemaVisitor[T]) -> T:
    # Parquet has no dictionary type. dictionary-encoding is handled
    # as an encoding detail, not as a separate type.
    # We will follow this approach in determining the Iceberg Type,
    # as we only support parquet in PyIceberg for now.
    logger.warning(f"Iceberg does not have a dictionary type. {type(obj)} will be inferred as {obj.value_type} on read.")
    return visit_pyarrow(obj.value_type, visitor)


@visit_pyarrow.register(pa.Field)
def _(obj: pa.Field, visitor: PyArrowSchemaVisitor[T]) -> T:
    field_type = obj.type

    visitor.before_field(obj)
    try:
        result = visit_pyarrow(field_type, visitor)
    except TypeError as e:
        raise UnsupportedPyArrowTypeException(obj, f"Column '{obj.name}' has an unsupported type: {field_type}") from e
    visitor.after_field(obj)

    return visitor.field(obj, result)


@visit_pyarrow.register(pa.DataType)
def _(obj: pa.DataType, visitor: PyArrowSchemaVisitor[T]) -> T:
    if pa.types.is_nested(obj):
        raise TypeError(f"Expected primitive type, got: {type(obj)}")
    return visitor.primitive(obj)


class PyArrowSchemaVisitor(Generic[T], ABC):
    def before_field(self, field: pa.Field) -> None:
        """Override this method to perform an action immediately before visiting a field."""

    def after_field(self, field: pa.Field) -> None:
        """Override this method to perform an action immediately after visiting a field."""

    def before_list_element(self, element: pa.Field) -> None:
        """Override this method to perform an action immediately before visiting an element within a ListType."""

    def after_list_element(self, element: pa.Field) -> None:
        """Override this method to perform an action immediately after visiting an element within a ListType."""

    def before_map_key(self, key: pa.Field) -> None:
        """Override this method to perform an action immediately before visiting a key within a MapType."""

    def after_map_key(self, key: pa.Field) -> None:
        """Override this method to perform an action immediately after visiting a key within a MapType."""

    def before_map_value(self, value: pa.Field) -> None:
        """Override this method to perform an action immediately before visiting a value within a MapType."""

    def after_map_value(self, value: pa.Field) -> None:
        """Override this method to perform an action immediately after visiting a value within a MapType."""

    @abstractmethod
    def schema(self, schema: pa.Schema, struct_result: T) -> T:
        """Visit a schema."""

    @abstractmethod
    def struct(self, struct: pa.StructType, field_results: List[T]) -> T:
        """Visit a struct."""

    @abstractmethod
    def field(self, field: pa.Field, field_result: T) -> T:
        """Visit a field."""

    @abstractmethod
    def list(self, list_type: pa.ListType, element_result: T) -> T:
        """Visit a list."""

    @abstractmethod
    def map(self, map_type: pa.MapType, key_result: T, value_result: T) -> T:
        """Visit a map."""

    @abstractmethod
    def primitive(self, primitive: pa.DataType) -> T:
        """Visit a primitive type."""


def _get_field_id(field: pa.Field) -> Optional[int]:
    return (
        int(field_id_str.decode())
        if (field.metadata and (field_id_str := field.metadata.get(PYARROW_PARQUET_FIELD_ID_KEY)))
        else None
    )


class _HasIds(PyArrowSchemaVisitor[bool]):
    def schema(self, schema: pa.Schema, struct_result: bool) -> bool:
        return struct_result

    def struct(self, struct: pa.StructType, field_results: List[bool]) -> bool:
        return all(field_results)

    def field(self, field: pa.Field, field_result: bool) -> bool:
        return all([_get_field_id(field) is not None, field_result])

    def list(self, list_type: pa.ListType, element_result: bool) -> bool:
        element_field = list_type.value_field
        element_id = _get_field_id(element_field)
        return element_result and element_id is not None

    def map(self, map_type: pa.MapType, key_result: bool, value_result: bool) -> bool:
        key_field = map_type.key_field
        key_id = _get_field_id(key_field)
        value_field = map_type.item_field
        value_id = _get_field_id(value_field)
        return all([key_id is not None, value_id is not None, key_result, value_result])

    def primitive(self, primitive: pa.DataType) -> bool:
        return True


class _ConvertToIceberg(PyArrowSchemaVisitor[Union[IcebergType, Schema]]):
    """Converts PyArrowSchema to Iceberg Schema. Applies the IDs from name_mapping if provided."""

    _field_names: List[str]

    def __init__(
        self, downcast_ns_timestamp_to_us: bool = False, format_version: TableVersion = TableProperties.DEFAULT_FORMAT_VERSION
    ) -> None:  # noqa: F821
        self._field_names = []
        self._downcast_ns_timestamp_to_us = downcast_ns_timestamp_to_us
        self._format_version = format_version

    def _field_id(self, field: pa.Field) -> int:
        if (field_id := _get_field_id(field)) is not None:
            return field_id
        else:
            raise ValueError(f"Cannot convert {field} to Iceberg Field as field_id is empty.")

    def schema(self, schema: pa.Schema, struct_result: StructType) -> Schema:
        return Schema(*struct_result.fields)

    def struct(self, struct: pa.StructType, field_results: List[NestedField]) -> StructType:
        return StructType(*field_results)

    def field(self, field: pa.Field, field_result: IcebergType) -> NestedField:
        field_id = self._field_id(field)
        field_doc = doc_str.decode() if (field.metadata and (doc_str := field.metadata.get(PYARROW_FIELD_DOC_KEY))) else None
        field_type = field_result
        return NestedField(field_id, field.name, field_type, required=not field.nullable, doc=field_doc)

    def list(self, list_type: pa.ListType, element_result: IcebergType) -> ListType:
        element_field = list_type.value_field
        self._field_names.append(LIST_ELEMENT_NAME)
        element_id = self._field_id(element_field)
        self._field_names.pop()
        return ListType(element_id, element_result, element_required=not element_field.nullable)

    def map(self, map_type: pa.MapType, key_result: IcebergType, value_result: IcebergType) -> MapType:
        key_field = map_type.key_field
        self._field_names.append(MAP_KEY_NAME)
        key_id = self._field_id(key_field)
        self._field_names.pop()
        value_field = map_type.item_field
        self._field_names.append(MAP_VALUE_NAME)
        value_id = self._field_id(value_field)
        self._field_names.pop()
        return MapType(key_id, key_result, value_id, value_result, value_required=not value_field.nullable)

    def primitive(self, primitive: pa.DataType) -> PrimitiveType:
        if pa.types.is_boolean(primitive):
            return BooleanType()
        elif pa.types.is_integer(primitive):
            width = primitive.bit_width
            if width <= 32:
                return IntegerType()
            elif width <= 64:
                return LongType()
            else:
                # Does not exist (yet)
                raise TypeError(f"Unsupported integer type: {primitive}")
        elif pa.types.is_float32(primitive):
            return FloatType()
        elif pa.types.is_float64(primitive):
            return DoubleType()
        elif isinstance(primitive, pa.Decimal128Type):
            primitive = cast(pa.Decimal128Type, primitive)
            return DecimalType(primitive.precision, primitive.scale)
        elif pa.types.is_string(primitive) or pa.types.is_large_string(primitive) or pa.types.is_string_view(primitive):
            return StringType()
        elif pa.types.is_date32(primitive):
            return DateType()
        elif isinstance(primitive, pa.Time64Type) and primitive.unit == "us":
            return TimeType()
        elif pa.types.is_timestamp(primitive):
            primitive = cast(pa.TimestampType, primitive)
            if primitive.unit in ("s", "ms", "us"):
                # Supported types, will be upcast automatically to 'us'
                pass
            elif primitive.unit == "ns":
                if self._downcast_ns_timestamp_to_us:
                    logger.warning("Iceberg does not yet support 'ns' timestamp precision. Downcasting to 'us'.")
                elif self._format_version >= 3:
                    if primitive.tz in UTC_ALIASES:
                        return TimestamptzNanoType()
                    else:
                        return TimestampNanoType()
                else:
                    raise TypeError(
                        "Iceberg does not yet support 'ns' timestamp precision. Use 'downcast-ns-timestamp-to-us-on-write' configuration property to automatically downcast 'ns' to 'us' on write.",
                    )
            else:
                raise TypeError(f"Unsupported precision for timestamp type: {primitive.unit}")

            if primitive.tz in UTC_ALIASES:
                return TimestamptzType()
            elif primitive.tz is None:
                return TimestampType()

        elif pa.types.is_binary(primitive) or pa.types.is_large_binary(primitive) or pa.types.is_binary_view(primitive):
            return BinaryType()
        elif pa.types.is_fixed_size_binary(primitive):
            primitive = cast(pa.FixedSizeBinaryType, primitive)
            return FixedType(primitive.byte_width)
        elif pa.types.is_null(primitive):
            # PyArrow null type (pa.null()) is converted to Iceberg UnknownType
            # UnknownType can be promoted to any primitive type in V3+ tables per the Iceberg spec
            return UnknownType()
        elif isinstance(primitive, pa.UuidType):
            return UUIDType()

        raise TypeError(f"Unsupported type: {primitive}")

    def before_field(self, field: pa.Field) -> None:
        self._field_names.append(field.name)

    def after_field(self, field: pa.Field) -> None:
        self._field_names.pop()

    def before_list_element(self, element: pa.Field) -> None:
        self._field_names.append(LIST_ELEMENT_NAME)

    def after_list_element(self, element: pa.Field) -> None:
        self._field_names.pop()

    def before_map_key(self, key: pa.Field) -> None:
        self._field_names.append(MAP_KEY_NAME)

    def after_map_key(self, element: pa.Field) -> None:
        self._field_names.pop()

    def before_map_value(self, value: pa.Field) -> None:
        self._field_names.append(MAP_VALUE_NAME)

    def after_map_value(self, element: pa.Field) -> None:
        self._field_names.pop()


class _ConvertToLargeTypes(PyArrowSchemaVisitor[Union[pa.DataType, pa.Schema]]):
    def schema(self, schema: pa.Schema, struct_result: pa.StructType) -> pa.Schema:
        return pa.schema(struct_result)

    def struct(self, struct: pa.StructType, field_results: List[pa.Field]) -> pa.StructType:
        return pa.struct(field_results)

    def field(self, field: pa.Field, field_result: pa.DataType) -> pa.Field:
        return field.with_type(field_result)

    def list(self, list_type: pa.ListType, element_result: pa.DataType) -> pa.DataType:
        return pa.large_list(element_result)

    def map(self, map_type: pa.MapType, key_result: pa.DataType, value_result: pa.DataType) -> pa.DataType:
        return pa.map_(key_result, value_result)

    def primitive(self, primitive: pa.DataType) -> pa.DataType:
        if primitive == pa.string():
            return pa.large_string()
        elif primitive == pa.binary():
            return pa.large_binary()
        return primitive


class _ConvertToSmallTypes(PyArrowSchemaVisitor[Union[pa.DataType, pa.Schema]]):
    def schema(self, schema: pa.Schema, struct_result: pa.StructType) -> pa.Schema:
        return pa.schema(struct_result)

    def struct(self, struct: pa.StructType, field_results: List[pa.Field]) -> pa.StructType:
        return pa.struct(field_results)

    def field(self, field: pa.Field, field_result: pa.DataType) -> pa.Field:
        return field.with_type(field_result)

    def list(self, list_type: pa.ListType, element_result: pa.DataType) -> pa.DataType:
        return pa.list_(element_result)

    def map(self, map_type: pa.MapType, key_result: pa.DataType, value_result: pa.DataType) -> pa.DataType:
        return pa.map_(key_result, value_result)

    def primitive(self, primitive: pa.DataType) -> pa.DataType:
        if primitive == pa.large_string():
            return pa.string()
        elif primitive == pa.large_binary():
            return pa.binary()
        return primitive


class _ConvertToIcebergWithoutIDs(_ConvertToIceberg):
    """
    Converts PyArrowSchema to Iceberg Schema with all -1 ids.

    The schema generated through this visitor should always be
    used in conjunction with `new_table_metadata` function to
    assign new field ids in order. This is currently used only
    when creating an Iceberg Schema from a PyArrow schema when
    creating a new Iceberg table.
    """

    def _field_id(self, field: pa.Field) -> int:
        return -1


def _get_column_projection_values(
    file: DataFile, projected_schema: Schema, partition_spec: Optional[PartitionSpec], file_project_field_ids: Set[int]
) -> Dict[int, Any]:
    """Apply Column Projection rules to File Schema."""
    project_schema_diff = projected_schema.field_ids.difference(file_project_field_ids)
    if len(project_schema_diff) == 0 or partition_spec is None:
        return EMPTY_DICT

    partition_schema = partition_spec.partition_type(projected_schema)
    accessors = build_position_accessors(partition_schema)

    projected_missing_fields = {}
    for field_id in project_schema_diff:
        for partition_field in partition_spec.fields_by_source_id(field_id):
            if isinstance(partition_field.transform, IdentityTransform):
                if partition_value := accessors[partition_field.field_id].get(file.partition):
                    projected_missing_fields[field_id] = partition_value

    return projected_missing_fields


def _task_to_record_batches(
    io: FileIO,
    task: FileScanTask,
    bound_row_filter: BooleanExpression,
    projected_schema: Schema,
    projected_field_ids: Set[int],
    positional_deletes: Optional[List[ChunkedArray]],
    case_sensitive: bool,
    name_mapping: Optional[NameMapping] = None,
    partition_spec: Optional[PartitionSpec] = None,
    format_version: TableVersion = TableProperties.DEFAULT_FORMAT_VERSION,
    downcast_ns_timestamp_to_us: Optional[bool] = None,
) -> Iterator[pa.RecordBatch]:
    arrow_format = ds.ParquetFileFormat(pre_buffer=True, buffer_size=(ONE_MEGABYTE * 8))
    with io.new_input(task.file.file_path).open() as fin:
        fragment = arrow_format.make_fragment(fin)
        physical_schema = fragment.physical_schema

        # For V1 and V2, we only support Timestamp 'us' in Iceberg Schema, therefore it is reasonable to always cast 'ns' timestamp to 'us' on read.
        # For V3 this has to set explicitly to avoid nanosecond timestamp to be down-casted by default
        downcast_ns_timestamp_to_us = (
            downcast_ns_timestamp_to_us if downcast_ns_timestamp_to_us is not None else format_version <= 2
        )
        file_schema = pyarrow_to_schema(
            physical_schema, name_mapping, downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us, format_version=format_version
        )

        # Apply column projection rules: https://iceberg.apache.org/spec/#column-projection
        projected_missing_fields = _get_column_projection_values(
            task.file, projected_schema, partition_spec, file_schema.field_ids
        )

        pyarrow_filter = None
        if bound_row_filter is not AlwaysTrue():
            translated_row_filter = translate_column_names(
                bound_row_filter, file_schema, case_sensitive=case_sensitive, projected_field_values=projected_missing_fields
            )
            bound_file_filter = bind(file_schema, translated_row_filter, case_sensitive=case_sensitive)
            pyarrow_filter = expression_to_pyarrow(bound_file_filter)

        file_project_schema = prune_columns(file_schema, projected_field_ids, select_full_types=False)

        fragment_scanner = ds.Scanner.from_fragment(
            fragment=fragment,
            schema=physical_schema,
            # This will push down the query to Arrow.
            # But in case there are positional deletes, we have to apply them first
            filter=pyarrow_filter if not positional_deletes else None,
            columns=[col.name for col in file_project_schema.columns],
        )

        next_index = 0
        batches = fragment_scanner.to_batches()
        for batch in batches:
            next_index = next_index + len(batch)
            current_index = next_index - len(batch)
            current_batch = batch

            if positional_deletes:
                # Create the mask of indices that we're interested in
                indices = _combine_positional_deletes(positional_deletes, current_index, current_index + len(batch))
                current_batch = current_batch.take(indices)

            # skip empty batches
            if current_batch.num_rows == 0:
                continue

            # Apply the user filter
            if pyarrow_filter is not None:
                # Temporary fix until PyArrow 21 is released ( https://github.com/apache/arrow/pull/46057 )
                table = pa.Table.from_batches([current_batch])
                table = table.filter(pyarrow_filter)
                # skip empty batches
                if table.num_rows == 0:
                    continue

                current_batch = table.combine_chunks().to_batches()[0]

            yield _to_requested_schema(
                projected_schema,
                file_project_schema,
                current_batch,
                downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
                projected_missing_fields=projected_missing_fields,
            )


def _read_all_delete_files(io: FileIO, tasks: Iterable[FileScanTask]) -> Dict[str, List[ChunkedArray]]:
    deletes_per_file: Dict[str, List[ChunkedArray]] = {}
    unique_deletes = set(itertools.chain.from_iterable([task.delete_files for task in tasks]))
    if len(unique_deletes) > 0:
        executor = ExecutorFactory.get_or_create()
        deletes_per_files: Iterator[Dict[str, ChunkedArray]] = executor.map(
            lambda args: _read_deletes(*args),
            [(io, delete_file) for delete_file in unique_deletes],
        )
        for delete in deletes_per_files:
            for file, arr in delete.items():
                if file in deletes_per_file:
                    deletes_per_file[file].append(arr)
                else:
                    deletes_per_file[file] = [arr]

    return deletes_per_file


class ArrowScan:
    _table_metadata: TableMetadata
    _io: FileIO
    _projected_schema: Schema
    _bound_row_filter: BooleanExpression
    _case_sensitive: bool
    _limit: Optional[int]
    _downcast_ns_timestamp_to_us: Optional[bool]
    """Scan the Iceberg Table and create an Arrow construct.

    Attributes:
        _table_metadata: Current table metadata of the Iceberg table
        _io: PyIceberg FileIO implementation from which to fetch the io properties
        _projected_schema: Iceberg Schema to project onto the data files
        _bound_row_filter: Schema bound row expression to filter the data with
        _case_sensitive: Case sensitivity when looking up column names
        _limit: Limit the number of records.
    """

    def __init__(
        self,
        table_metadata: TableMetadata,
        io: FileIO,
        projected_schema: Schema,
        row_filter: BooleanExpression,
        case_sensitive: bool = True,
        limit: Optional[int] = None,
    ) -> None:
        self._table_metadata = table_metadata
        self._io = io
        self._projected_schema = projected_schema
        self._bound_row_filter = bind(table_metadata.schema(), row_filter, case_sensitive=case_sensitive)
        self._case_sensitive = case_sensitive
        self._limit = limit
        self._downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE)

    @property
    def _projected_field_ids(self) -> Set[int]:
        """Set of field IDs that should be projected from the data files."""
        return {
            id
            for id in self._projected_schema.field_ids
            if not isinstance(self._projected_schema.find_type(id), (MapType, ListType))
        }.union(extract_field_ids(self._bound_row_filter))

    def to_table(self, tasks: Iterable[FileScanTask]) -> pa.Table:
        """Scan the Iceberg table and return a pa.Table.

        Returns a pa.Table with data from the Iceberg table by resolving the
        right columns that match the current table schema. Only data that
        matches the provided row_filter expression is returned.

        Args:
            tasks: FileScanTasks representing the data files and delete files to read from.

        Returns:
            A PyArrow table. Total number of rows will be capped if specified.

        Raises:
            ResolveError: When a required field cannot be found in the file
            ValueError: When a field type in the file cannot be projected to the schema type
        """
        arrow_schema = schema_to_pyarrow(self._projected_schema, include_field_ids=False)

        batches = self.to_record_batches(tasks)
        try:
            first_batch = next(batches)
        except StopIteration:
            # Empty
            return arrow_schema.empty_table()

        # Note: cannot use pa.Table.from_batches(itertools.chain([first_batch], batches)))
        #       as different batches can use different schema's (due to large_ types)
        result = pa.concat_tables(
            (pa.Table.from_batches([batch]) for batch in itertools.chain([first_batch], batches)), promote_options="permissive"
        )

        if property_as_bool(self._io.properties, PYARROW_USE_LARGE_TYPES_ON_READ, False):
            deprecation_message(
                deprecated_in="0.10.0",
                removed_in="0.11.0",
                help_message=f"Property `{PYARROW_USE_LARGE_TYPES_ON_READ}` will be removed.",
            )
            result = result.cast(arrow_schema)

        return result

    def to_record_batches(self, tasks: Iterable[FileScanTask]) -> Iterator[pa.RecordBatch]:
        """Scan the Iceberg table and return an Iterator[pa.RecordBatch].

        Returns an Iterator of pa.RecordBatch with data from the Iceberg table
        by resolving the right columns that match the current table schema.
        Only data that matches the provided row_filter expression is returned.

        Args:
            tasks: FileScanTasks representing the data files and delete files to read from.

        Returns:
            An Iterator of PyArrow RecordBatches.
            Total number of rows will be capped if specified.

        Raises:
            ResolveError: When a required field cannot be found in the file
            ValueError: When a field type in the file cannot be projected to the schema type
        """
        deletes_per_file = _read_all_delete_files(self._io, tasks)

        total_row_count = 0
        executor = ExecutorFactory.get_or_create()

        def batches_for_task(task: FileScanTask) -> List[pa.RecordBatch]:
            # Materialize the iterator here to ensure execution happens within the executor.
            # Otherwise, the iterator would be lazily consumed later (in the main thread),
            # defeating the purpose of using executor.map.
            return list(self._record_batches_from_scan_tasks_and_deletes([task], deletes_per_file))

        limit_reached = False
        for batches in executor.map(batches_for_task, tasks):
            for batch in batches:
                current_batch_size = len(batch)
                if self._limit is not None and total_row_count + current_batch_size >= self._limit:
                    yield batch.slice(0, self._limit - total_row_count)

                    limit_reached = True
                    break
                else:
                    yield batch
                    total_row_count += current_batch_size

            if limit_reached:
                # This break will also cancel all running tasks in the executor
                break

    def _record_batches_from_scan_tasks_and_deletes(
        self, tasks: Iterable[FileScanTask], deletes_per_file: Dict[str, List[ChunkedArray]]
    ) -> Iterator[pa.RecordBatch]:
        total_row_count = 0
        for task in tasks:
            if self._limit is not None and total_row_count >= self._limit:
                break
            batches = _task_to_record_batches(
                self._io,
                task,
                self._bound_row_filter,
                self._projected_schema,
                self._projected_field_ids,
                deletes_per_file.get(task.file.file_path),
                self._case_sensitive,
                self._table_metadata.name_mapping(),
                self._table_metadata.specs().get(task.file.spec_id),
                self._table_metadata.format_version,
                self._downcast_ns_timestamp_to_us,
            )
            for batch in batches:
                if self._limit is not None:
                    if total_row_count >= self._limit:
                        break
                    elif total_row_count + len(batch) >= self._limit:
                        batch = batch.slice(0, self._limit - total_row_count)
                yield batch
                total_row_count += len(batch)


def _to_requested_schema(
    requested_schema: Schema,
    file_schema: Schema,
    batch: pa.RecordBatch,
    downcast_ns_timestamp_to_us: bool = False,
    include_field_ids: bool = False,
    projected_missing_fields: Dict[int, Any] = EMPTY_DICT,
) -> pa.RecordBatch:
    # We could reuse some of these visitors
    struct_array = visit_with_partner(
        requested_schema,
        batch,
        ArrowProjectionVisitor(
            file_schema, downcast_ns_timestamp_to_us, include_field_ids, projected_missing_fields=projected_missing_fields
        ),
        ArrowAccessor(file_schema),
    )
    return pa.RecordBatch.from_struct_array(struct_array)


class ArrowProjectionVisitor(SchemaWithPartnerVisitor[pa.Array, Optional[pa.Array]]):
    _file_schema: Schema
    _include_field_ids: bool
    _downcast_ns_timestamp_to_us: bool
    _use_large_types: Optional[bool]
    _projected_missing_fields: Dict[int, Any]

    def __init__(
        self,
        file_schema: Schema,
        downcast_ns_timestamp_to_us: bool = False,
        include_field_ids: bool = False,
        use_large_types: Optional[bool] = None,
        projected_missing_fields: Dict[int, Any] = EMPTY_DICT,
    ) -> None:
        self._file_schema = file_schema
        self._include_field_ids = include_field_ids
        self._downcast_ns_timestamp_to_us = downcast_ns_timestamp_to_us
        self._use_large_types = use_large_types
        self._projected_missing_fields = projected_missing_fields

        if use_large_types is not None:
            deprecation_message(
                deprecated_in="0.10.0",
                removed_in="0.11.0",
                help_message="Argument `use_large_types` will be removed from ArrowProjectionVisitor",
            )

    def _cast_if_needed(self, field: NestedField, values: pa.Array) -> pa.Array:
        file_field = self._file_schema.find_field(field.field_id)

        if field.field_type.is_primitive:
            if (target_type := schema_to_pyarrow(field.field_type, include_field_ids=self._include_field_ids)) != values.type:
                if field.field_type == TimestampType():
                    # Downcasting of nanoseconds to microseconds
                    if (
                        pa.types.is_timestamp(target_type)
                        and not target_type.tz
                        and pa.types.is_timestamp(values.type)
                        and not values.type.tz
                    ):
                        if target_type.unit == "us" and values.type.unit == "ns" and self._downcast_ns_timestamp_to_us:
                            return values.cast(target_type, safe=False)
                        elif target_type.unit == "us" and values.type.unit in {"s", "ms"}:
                            return values.cast(target_type)
                    raise ValueError(f"Unsupported schema projection from {values.type} to {target_type}")
                elif field.field_type == TimestamptzType():
                    if (
                        pa.types.is_timestamp(target_type)
                        and target_type.tz == "UTC"
                        and pa.types.is_timestamp(values.type)
                        and (values.type.tz in UTC_ALIASES or values.type.tz is None)
                    ):
                        if target_type.unit == "us" and values.type.unit == "ns" and self._downcast_ns_timestamp_to_us:
                            return values.cast(target_type, safe=False)
                        elif target_type.unit == "us" and values.type.unit in {"s", "ms", "us"}:
                            return values.cast(target_type)
                    raise ValueError(f"Unsupported schema projection from {values.type} to {target_type}")

            if field.field_type != file_field.field_type:
                target_schema = schema_to_pyarrow(
                    promote(file_field.field_type, field.field_type), include_field_ids=self._include_field_ids
                )
                if self._use_large_types is False:
                    target_schema = _pyarrow_schema_ensure_small_types(target_schema)
                return values.cast(target_schema)

        return values

    def _construct_field(self, field: NestedField, arrow_type: pa.DataType) -> pa.Field:
        metadata = {}
        if field.doc:
            metadata[PYARROW_FIELD_DOC_KEY] = field.doc
        if self._include_field_ids:
            metadata[PYARROW_PARQUET_FIELD_ID_KEY] = str(field.field_id)

        return pa.field(
            name=field.name,
            type=arrow_type,
            nullable=field.optional,
            metadata=metadata,
        )

    def schema(self, schema: Schema, schema_partner: Optional[pa.Array], struct_result: Optional[pa.Array]) -> Optional[pa.Array]:
        return struct_result

    def struct(
        self, struct: StructType, struct_array: Optional[pa.Array], field_results: List[Optional[pa.Array]]
    ) -> Optional[pa.Array]:
        if struct_array is None:
            return None
        field_arrays: List[pa.Array] = []
        fields: List[pa.Field] = []
        for field, field_array in zip(struct.fields, field_results):
            if field_array is not None:
                array = self._cast_if_needed(field, field_array)
                field_arrays.append(array)
                fields.append(self._construct_field(field, array.type))
            elif field.optional or field.initial_default is not None:
                # When an optional field is added, or when a required field with a non-null initial default is added
                arrow_type = schema_to_pyarrow(field.field_type, include_field_ids=self._include_field_ids)
                if projected_value := self._projected_missing_fields.get(field.field_id):
                    field_arrays.append(pa.repeat(pa.scalar(projected_value, type=arrow_type), len(struct_array)))
                elif field.initial_default is None:
                    field_arrays.append(pa.nulls(len(struct_array), type=arrow_type))
                else:
                    field_arrays.append(pa.repeat(pa.scalar(field.initial_default, type=arrow_type), len(struct_array)))
                fields.append(self._construct_field(field, arrow_type))
            else:
                raise ResolveError(f"Field is required, and could not be found in the file: {field}")

        return pa.StructArray.from_arrays(
            arrays=field_arrays,
            fields=pa.struct(fields),
            mask=struct_array.is_null() if isinstance(struct_array, pa.StructArray) else None,
        )

    def field(self, field: NestedField, _: Optional[pa.Array], field_array: Optional[pa.Array]) -> Optional[pa.Array]:
        return field_array

    def list(self, list_type: ListType, list_array: Optional[pa.Array], value_array: Optional[pa.Array]) -> Optional[pa.Array]:
        if isinstance(list_array, (pa.ListArray, pa.LargeListArray, pa.FixedSizeListArray)) and value_array is not None:
            list_initializer = pa.large_list if isinstance(list_array, pa.LargeListArray) else pa.list_
            if isinstance(value_array, pa.StructArray):
                # This can be removed once this has been fixed:
                # https://github.com/apache/arrow/issues/38809
                list_array = pa.LargeListArray.from_arrays(list_array.offsets, value_array)
            value_array = self._cast_if_needed(list_type.element_field, value_array)
            arrow_field = list_initializer(self._construct_field(list_type.element_field, value_array.type))
            return list_array.cast(arrow_field)
        else:
            return None

    def map(
        self, map_type: MapType, map_array: Optional[pa.Array], key_result: Optional[pa.Array], value_result: Optional[pa.Array]
    ) -> Optional[pa.Array]:
        if isinstance(map_array, pa.MapArray) and key_result is not None and value_result is not None:
            key_result = self._cast_if_needed(map_type.key_field, key_result)
            value_result = self._cast_if_needed(map_type.value_field, value_result)
            arrow_field = pa.map_(
                self._construct_field(map_type.key_field, key_result.type),
                self._construct_field(map_type.value_field, value_result.type),
            )
            if isinstance(value_result, pa.StructArray):
                # Arrow does not allow reordering of fields, therefore we have to copy the array :(
                return pa.MapArray.from_arrays(map_array.offsets, key_result, value_result, arrow_field)
            else:
                return map_array.cast(arrow_field)
        else:
            return None

    def primitive(self, _: PrimitiveType, array: Optional[pa.Array]) -> Optional[pa.Array]:
        return array


class ArrowAccessor(PartnerAccessor[pa.Array]):
    file_schema: Schema

    def __init__(self, file_schema: Schema):
        self.file_schema = file_schema

    def schema_partner(self, partner: Optional[pa.Array]) -> Optional[pa.Array]:
        return partner

    def field_partner(self, partner_struct: Optional[pa.Array], field_id: int, _: str) -> Optional[pa.Array]:
        if partner_struct is not None:
            # use the field name from the file schema
            try:
                name = self.file_schema.find_field(field_id).name
            except ValueError:
                return None

            if isinstance(partner_struct, pa.StructArray):
                return partner_struct.field(name)
            elif isinstance(partner_struct, pa.Table):
                return partner_struct.column(name).combine_chunks()
            elif isinstance(partner_struct, pa.RecordBatch):
                return partner_struct.column(name)
            else:
                raise ValueError(f"Cannot find {name} in expected partner_struct type {type(partner_struct)}")

        return None

    def list_element_partner(self, partner_list: Optional[pa.Array]) -> Optional[pa.Array]:
        return partner_list.values if isinstance(partner_list, (pa.ListArray, pa.LargeListArray, pa.FixedSizeListArray)) else None

    def map_key_partner(self, partner_map: Optional[pa.Array]) -> Optional[pa.Array]:
        return partner_map.keys if isinstance(partner_map, pa.MapArray) else None

    def map_value_partner(self, partner_map: Optional[pa.Array]) -> Optional[pa.Array]:
        return partner_map.items if isinstance(partner_map, pa.MapArray) else None


def _primitive_to_physical(iceberg_type: PrimitiveType) -> str:
    return visit(iceberg_type, _PRIMITIVE_TO_PHYSICAL_TYPE_VISITOR)


class PrimitiveToPhysicalType(SchemaVisitorPerPrimitiveType[str]):
    def schema(self, schema: Schema, struct_result: str) -> str:
        raise ValueError(f"Expected primitive-type, got: {schema}")

    def struct(self, struct: StructType, field_results: List[str]) -> str:
        raise ValueError(f"Expected primitive-type, got: {struct}")

    def field(self, field: NestedField, field_result: str) -> str:
        raise ValueError(f"Expected primitive-type, got: {field}")

    def list(self, list_type: ListType, element_result: str) -> str:
        raise ValueError(f"Expected primitive-type, got: {list_type}")

    def map(self, map_type: MapType, key_result: str, value_result: str) -> str:
        raise ValueError(f"Expected primitive-type, got: {map_type}")

    def visit_fixed(self, fixed_type: FixedType) -> str:
        return "FIXED_LEN_BYTE_ARRAY"

    def visit_decimal(self, decimal_type: DecimalType) -> str:
        return "INT32" if decimal_type.precision <= 9 else "INT64" if decimal_type.precision <= 18 else "FIXED_LEN_BYTE_ARRAY"

    def visit_boolean(self, boolean_type: BooleanType) -> str:
        return "BOOLEAN"

    def visit_integer(self, integer_type: IntegerType) -> str:
        return "INT32"

    def visit_long(self, long_type: LongType) -> str:
        return "INT64"

    def visit_float(self, float_type: FloatType) -> str:
        return "FLOAT"

    def visit_double(self, double_type: DoubleType) -> str:
        return "DOUBLE"

    def visit_date(self, date_type: DateType) -> str:
        return "INT32"

    def visit_time(self, time_type: TimeType) -> str:
        return "INT64"

    def visit_timestamp(self, timestamp_type: TimestampType) -> str:
        return "INT64"

    def visit_timestamp_ns(self, timestamp_type: TimestampNanoType) -> str:
        return "INT64"

    def visit_timestamptz(self, timestamptz_type: TimestamptzType) -> str:
        return "INT64"

    def visit_timestamptz_ns(self, timestamptz_ns_type: TimestamptzNanoType) -> str:
        return "INT64"

    def visit_string(self, string_type: StringType) -> str:
        return "BYTE_ARRAY"

    def visit_uuid(self, uuid_type: UUIDType) -> str:
        return "FIXED_LEN_BYTE_ARRAY"

    def visit_binary(self, binary_type: BinaryType) -> str:
        return "BYTE_ARRAY"

    def visit_unknown(self, unknown_type: UnknownType) -> str:
        return "UNKNOWN"


_PRIMITIVE_TO_PHYSICAL_TYPE_VISITOR = PrimitiveToPhysicalType()


class StatsAggregator:
    current_min: Any
    current_max: Any
    trunc_length: Optional[int]

    def __init__(self, iceberg_type: PrimitiveType, physical_type_string: str, trunc_length: Optional[int] = None) -> None:
        self.current_min = None
        self.current_max = None
        self.trunc_length = trunc_length

        expected_physical_type = _primitive_to_physical(iceberg_type)
        if expected_physical_type != physical_type_string:
            # Allow promotable physical types
            # INT32 -> INT64 and FLOAT -> DOUBLE are safe type casts
            if (physical_type_string == "INT32" and expected_physical_type == "INT64") or (
                physical_type_string == "FLOAT" and expected_physical_type == "DOUBLE"
            ):
                pass
            else:
                raise ValueError(
                    f"Unexpected physical type {physical_type_string} for {iceberg_type}, expected {expected_physical_type}"
                )

        self.primitive_type = iceberg_type

    def serialize(self, value: Any) -> bytes:
        return to_bytes(self.primitive_type, value)

    def update_min(self, val: Optional[Any]) -> None:
        if self.current_min is None:
            self.current_min = val
        elif val is not None:
            self.current_min = min(val, self.current_min)

    def update_max(self, val: Optional[Any]) -> None:
        if self.current_max is None:
            self.current_max = val
        elif val is not None:
            self.current_max = max(val, self.current_max)

    def min_as_bytes(self) -> Optional[bytes]:
        if self.current_min is None:
            return None

        return self.serialize(
            self.current_min
            if self.trunc_length is None
            else TruncateTransform(width=self.trunc_length).transform(self.primitive_type)(self.current_min)
        )

    def max_as_bytes(self) -> Optional[bytes]:
        if self.current_max is None:
            return None

        if self.primitive_type == StringType():
            if not isinstance(self.current_max, str):
                raise ValueError("Expected the current_max to be a string")
            s_result = truncate_upper_bound_text_string(self.current_max, self.trunc_length)
            return self.serialize(s_result) if s_result is not None else None
        elif self.primitive_type == BinaryType():
            if not isinstance(self.current_max, bytes):
                raise ValueError("Expected the current_max to be bytes")
            b_result = truncate_upper_bound_binary_string(self.current_max, self.trunc_length)
            return self.serialize(b_result) if b_result is not None else None
        else:
            if self.trunc_length is not None:
                raise ValueError(f"{self.primitive_type} cannot be truncated")
            return self.serialize(self.current_max)


DEFAULT_TRUNCATION_LENGTH = 16
TRUNCATION_EXPR = r"^truncate\((\d+)\)$"


class MetricModeTypes(Enum):
    TRUNCATE = "truncate"
    NONE = "none"
    COUNTS = "counts"
    FULL = "full"


@dataclass(frozen=True)
class MetricsMode(Singleton):
    type: MetricModeTypes
    length: Optional[int] = None


def match_metrics_mode(mode: str) -> MetricsMode:
    sanitized_mode = mode.strip().lower()
    if sanitized_mode.startswith("truncate"):
        m = re.match(TRUNCATION_EXPR, sanitized_mode)
        if m:
            length = int(m[1])
            if length < 1:
                raise ValueError("Truncation length must be larger than 0")
            return MetricsMode(MetricModeTypes.TRUNCATE, int(m[1]))
        else:
            raise ValueError(f"Malformed truncate: {mode}")
    elif sanitized_mode == "none":
        return MetricsMode(MetricModeTypes.NONE)
    elif sanitized_mode == "counts":
        return MetricsMode(MetricModeTypes.COUNTS)
    elif sanitized_mode == "full":
        return MetricsMode(MetricModeTypes.FULL)
    else:
        raise ValueError(f"Unsupported metrics mode: {mode}")


@dataclass(frozen=True)
class StatisticsCollector:
    field_id: int
    iceberg_type: PrimitiveType
    mode: MetricsMode
    column_name: str


class PyArrowStatisticsCollector(PreOrderSchemaVisitor[List[StatisticsCollector]]):
    _field_id: int = 0
    _schema: Schema
    _properties: Dict[str, str]
    _default_mode: str

    def __init__(self, schema: Schema, properties: Dict[str, str]):
        from pyiceberg.table import TableProperties

        self._schema = schema
        self._properties = properties
        self._default_mode = self._properties.get(
            TableProperties.DEFAULT_WRITE_METRICS_MODE, TableProperties.DEFAULT_WRITE_METRICS_MODE_DEFAULT
        )

    def schema(self, schema: Schema, struct_result: Callable[[], List[StatisticsCollector]]) -> List[StatisticsCollector]:
        return struct_result()

    def struct(
        self, struct: StructType, field_results: List[Callable[[], List[StatisticsCollector]]]
    ) -> List[StatisticsCollector]:
        return list(itertools.chain(*[result() for result in field_results]))

    def field(self, field: NestedField, field_result: Callable[[], List[StatisticsCollector]]) -> List[StatisticsCollector]:
        self._field_id = field.field_id
        return field_result()

    def list(self, list_type: ListType, element_result: Callable[[], List[StatisticsCollector]]) -> List[StatisticsCollector]:
        self._field_id = list_type.element_id
        return element_result()

    def map(
        self,
        map_type: MapType,
        key_result: Callable[[], List[StatisticsCollector]],
        value_result: Callable[[], List[StatisticsCollector]],
    ) -> List[StatisticsCollector]:
        self._field_id = map_type.key_id
        k = key_result()
        self._field_id = map_type.value_id
        v = value_result()
        return k + v

    def primitive(self, primitive: PrimitiveType) -> List[StatisticsCollector]:
        from pyiceberg.table import TableProperties

        column_name = self._schema.find_column_name(self._field_id)
        if column_name is None:
            return []

        metrics_mode = match_metrics_mode(self._default_mode)

        col_mode = self._properties.get(f"{TableProperties.METRICS_MODE_COLUMN_CONF_PREFIX}.{column_name}")
        if col_mode:
            metrics_mode = match_metrics_mode(col_mode)

        if (
            not (isinstance(primitive, StringType) or isinstance(primitive, BinaryType))
            and metrics_mode.type == MetricModeTypes.TRUNCATE
        ):
            metrics_mode = MetricsMode(MetricModeTypes.FULL)

        is_nested = column_name.find(".") >= 0

        if is_nested and metrics_mode.type in [MetricModeTypes.TRUNCATE, MetricModeTypes.FULL]:
            metrics_mode = MetricsMode(MetricModeTypes.COUNTS)

        return [StatisticsCollector(field_id=self._field_id, iceberg_type=primitive, mode=metrics_mode, column_name=column_name)]


def compute_statistics_plan(
    schema: Schema,
    table_properties: Dict[str, str],
) -> Dict[int, StatisticsCollector]:
    """
    Compute the statistics plan for all columns.

    The resulting list is assumed to have the same length and same order as the columns in the pyarrow table.
    This allows the list to map from the column index to the Iceberg column ID.
    For each element, the desired metrics collection that was provided by the user in the configuration
    is computed and then adjusted according to the data type of the column. For nested columns the minimum
    and maximum values are not computed. And truncation is only applied to text of binary strings.

    Args:
        table_properties (from pyiceberg.table.metadata.TableMetadata): The Iceberg table metadata properties.
            They are required to compute the mapping of column position to iceberg schema type id. It's also
            used to set the mode for column metrics collection
    """
    stats_cols = pre_order_visit(schema, PyArrowStatisticsCollector(schema, table_properties))
    result: Dict[int, StatisticsCollector] = {}
    for stats_col in stats_cols:
        result[stats_col.field_id] = stats_col
    return result


@dataclass(frozen=True)
class ID2ParquetPath:
    field_id: int
    parquet_path: str


class ID2ParquetPathVisitor(PreOrderSchemaVisitor[List[ID2ParquetPath]]):
    _field_id: int = 0
    _path: List[str]

    def __init__(self) -> None:
        self._path = []

    def schema(self, schema: Schema, struct_result: Callable[[], List[ID2ParquetPath]]) -> List[ID2ParquetPath]:
        return struct_result()

    def struct(self, struct: StructType, field_results: List[Callable[[], List[ID2ParquetPath]]]) -> List[ID2ParquetPath]:
        return list(itertools.chain(*[result() for result in field_results]))

    def field(self, field: NestedField, field_result: Callable[[], List[ID2ParquetPath]]) -> List[ID2ParquetPath]:
        self._field_id = field.field_id
        self._path.append(field.name)
        result = field_result()
        self._path.pop()
        return result

    def list(self, list_type: ListType, element_result: Callable[[], List[ID2ParquetPath]]) -> List[ID2ParquetPath]:
        self._field_id = list_type.element_id
        self._path.append("list.element")
        result = element_result()
        self._path.pop()
        return result

    def map(
        self,
        map_type: MapType,
        key_result: Callable[[], List[ID2ParquetPath]],
        value_result: Callable[[], List[ID2ParquetPath]],
    ) -> List[ID2ParquetPath]:
        self._field_id = map_type.key_id
        self._path.append("key_value.key")
        k = key_result()
        self._path.pop()
        self._field_id = map_type.value_id
        self._path.append("key_value.value")
        v = value_result()
        self._path.pop()
        return k + v

    def primitive(self, primitive: PrimitiveType) -> List[ID2ParquetPath]:
        return [ID2ParquetPath(field_id=self._field_id, parquet_path=".".join(self._path))]


def parquet_path_to_id_mapping(
    schema: Schema,
) -> Dict[str, int]:
    """
    Compute the mapping of parquet column path to Iceberg ID.

    For each column, the parquet file metadata has a path_in_schema attribute that follows
    a specific naming scheme for nested columns. This function computes a mapping of
    the full paths to the corresponding Iceberg IDs.

    Args:
        schema (pyiceberg.schema.Schema): The current table schema.
    """
    result: Dict[str, int] = {}
    for pair in pre_order_visit(schema, ID2ParquetPathVisitor()):
        result[pair.parquet_path] = pair.field_id
    return result


@dataclass(frozen=True)
class DataFileStatistics:
    record_count: int
    column_sizes: Dict[int, int]
    value_counts: Dict[int, int]
    null_value_counts: Dict[int, int]
    nan_value_counts: Dict[int, int]
    column_aggregates: Dict[int, StatsAggregator]
    split_offsets: List[int]

    def _partition_value(self, partition_field: PartitionField, schema: Schema) -> Any:
        if partition_field.source_id not in self.column_aggregates:
            return None

        source_field = schema.find_field(partition_field.source_id)
        iceberg_transform = partition_field.transform

        if not iceberg_transform.preserves_order:
            raise ValueError(
                f"Cannot infer partition value from parquet metadata for a non-linear Partition Field: {partition_field.name} with transform {partition_field.transform}"
            )

        transform_func = iceberg_transform.transform(source_field.field_type)

        lower_value = transform_func(
            partition_record_value(
                partition_field=partition_field,
                value=self.column_aggregates[partition_field.source_id].current_min,
                schema=schema,
            )
        )
        upper_value = transform_func(
            partition_record_value(
                partition_field=partition_field,
                value=self.column_aggregates[partition_field.source_id].current_max,
                schema=schema,
            )
        )
        if lower_value != upper_value:
            raise ValueError(
                f"Cannot infer partition value from parquet metadata as there are more than one partition values for Partition Field: {partition_field.name}. {lower_value=}, {upper_value=}"
            )

        return lower_value

    def partition(self, partition_spec: PartitionSpec, schema: Schema) -> Record:
        return Record(*[self._partition_value(field, schema) for field in partition_spec.fields])

    def to_serialized_dict(self) -> Dict[str, Any]:
        lower_bounds = {}
        upper_bounds = {}

        for k, agg in self.column_aggregates.items():
            _min = agg.min_as_bytes()
            if _min is not None:
                lower_bounds[k] = _min
            _max = agg.max_as_bytes()
            if _max is not None:
                upper_bounds[k] = _max
        return {
            "record_count": self.record_count,
            "column_sizes": self.column_sizes,
            "value_counts": self.value_counts,
            "null_value_counts": self.null_value_counts,
            "nan_value_counts": self.nan_value_counts,
            "lower_bounds": lower_bounds,
            "upper_bounds": upper_bounds,
            "split_offsets": self.split_offsets,
        }


def data_file_statistics_from_parquet_metadata(
    parquet_metadata: pq.FileMetaData,
    stats_columns: Dict[int, StatisticsCollector],
    parquet_column_mapping: Dict[str, int],
) -> DataFileStatistics:
    """
    Compute and return DataFileStatistics that includes the following.

    - record_count
    - column_sizes
    - value_counts
    - null_value_counts
    - nan_value_counts
    - column_aggregates
    - split_offsets

    Args:
        parquet_metadata (pyarrow.parquet.FileMetaData): A pyarrow metadata object.
        stats_columns (Dict[int, StatisticsCollector]): The statistics gathering plan. It is required to
            set the mode for column metrics collection
        parquet_column_mapping (Dict[str, int]): The mapping of the parquet file name to the field ID
    """
    column_sizes: Dict[int, int] = {}
    value_counts: Dict[int, int] = {}
    split_offsets: List[int] = []

    null_value_counts: Dict[int, int] = {}
    nan_value_counts: Dict[int, int] = {}

    col_aggs = {}

    invalidate_col: Set[int] = set()
    for r in range(parquet_metadata.num_row_groups):
        # References:
        # https://github.com/apache/iceberg/blob/fc381a81a1fdb8f51a0637ca27cd30673bd7aad3/parquet/src/main/java/org/apache/iceberg/parquet/ParquetUtil.java#L232
        # https://github.com/apache/parquet-mr/blob/ac29db4611f86a07cc6877b416aa4b183e09b353/parquet-hadoop/src/main/java/org/apache/parquet/hadoop/metadata/ColumnChunkMetaData.java#L184

        row_group = parquet_metadata.row_group(r)

        data_offset = row_group.column(0).data_page_offset
        dictionary_offset = row_group.column(0).dictionary_page_offset

        if row_group.column(0).has_dictionary_page and dictionary_offset < data_offset:
            split_offsets.append(dictionary_offset)
        else:
            split_offsets.append(data_offset)

        for pos in range(parquet_metadata.num_columns):
            column = row_group.column(pos)
            field_id = parquet_column_mapping[column.path_in_schema]

            stats_col = stats_columns[field_id]

            column_sizes.setdefault(field_id, 0)
            column_sizes[field_id] += column.total_compressed_size

            if stats_col.mode == MetricsMode(MetricModeTypes.NONE):
                continue

            value_counts[field_id] = value_counts.get(field_id, 0) + column.num_values

            if column.is_stats_set:
                try:
                    statistics = column.statistics

                    if statistics.has_null_count:
                        null_value_counts[field_id] = null_value_counts.get(field_id, 0) + statistics.null_count

                    if stats_col.mode == MetricsMode(MetricModeTypes.COUNTS):
                        continue

                    if field_id not in col_aggs:
                        try:
                            col_aggs[field_id] = StatsAggregator(
                                stats_col.iceberg_type, statistics.physical_type, stats_col.mode.length
                            )
                        except ValueError as e:
                            raise ValueError(f"{e} for column '{stats_col.column_name}'") from e

                    if isinstance(stats_col.iceberg_type, DecimalType) and statistics.physical_type != "FIXED_LEN_BYTE_ARRAY":
                        scale = stats_col.iceberg_type.scale
                        col_aggs[field_id].update_min(
                            unscaled_to_decimal(statistics.min_raw, scale)
                        ) if statistics.min_raw is not None else None
                        col_aggs[field_id].update_max(
                            unscaled_to_decimal(statistics.max_raw, scale)
                        ) if statistics.max_raw is not None else None
                    else:
                        col_aggs[field_id].update_min(statistics.min)
                        col_aggs[field_id].update_max(statistics.max)

                except pyarrow.lib.ArrowNotImplementedError as e:
                    invalidate_col.add(field_id)
                    logger.warning(e)
            else:
                invalidate_col.add(field_id)
                logger.warning("PyArrow statistics missing for column %d when writing file", pos)

    split_offsets.sort()

    for field_id in invalidate_col:
        col_aggs.pop(field_id, None)
        null_value_counts.pop(field_id, None)

    return DataFileStatistics(
        record_count=parquet_metadata.num_rows,
        column_sizes=column_sizes,
        value_counts=value_counts,
        null_value_counts=null_value_counts,
        nan_value_counts=nan_value_counts,
        column_aggregates=col_aggs,
        split_offsets=split_offsets,
    )


def write_file(io: FileIO, table_metadata: TableMetadata, tasks: Iterator[WriteTask]) -> Iterator[DataFile]:
    from pyiceberg.table import DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE, TableProperties

    parquet_writer_kwargs = _get_parquet_writer_kwargs(table_metadata.properties)
    row_group_size = property_as_int(
        properties=table_metadata.properties,
        property_name=TableProperties.PARQUET_ROW_GROUP_LIMIT,
        default=TableProperties.PARQUET_ROW_GROUP_LIMIT_DEFAULT,
    )
    location_provider = load_location_provider(table_location=table_metadata.location, table_properties=table_metadata.properties)

    def write_parquet(task: WriteTask) -> DataFile:
        table_schema = table_metadata.schema()
        # if schema needs to be transformed, use the transformed schema and adjust the arrow table accordingly
        # otherwise use the original schema
        if (sanitized_schema := sanitize_column_names(table_schema)) != table_schema:
            file_schema = sanitized_schema
        else:
            file_schema = table_schema

        downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
        batches = [
            _to_requested_schema(
                requested_schema=file_schema,
                file_schema=task.schema,
                batch=batch,
                downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
                include_field_ids=True,
            )
            for batch in task.record_batches
        ]
        arrow_table = pa.Table.from_batches(batches)
        file_path = location_provider.new_data_location(
            data_file_name=task.generate_data_file_filename("parquet"),
            partition_key=task.partition_key,
        )
        fo = io.new_output(file_path)
        with fo.create(overwrite=True) as fos:
            with pq.ParquetWriter(
                fos, schema=arrow_table.schema, store_decimal_as_integer=True, **parquet_writer_kwargs
            ) as writer:
                writer.write(arrow_table, row_group_size=row_group_size)
        statistics = data_file_statistics_from_parquet_metadata(
            parquet_metadata=writer.writer.metadata,
            stats_columns=compute_statistics_plan(file_schema, table_metadata.properties),
            parquet_column_mapping=parquet_path_to_id_mapping(file_schema),
        )
        data_file = DataFile.from_args(
            content=DataFileContent.DATA,
            file_path=file_path,
            file_format=FileFormat.PARQUET,
            partition=task.partition_key.partition if task.partition_key else Record(),
            file_size_in_bytes=len(fo),
            # After this has been fixed:
            # https://github.com/apache/iceberg-python/issues/271
            # sort_order_id=task.sort_order_id,
            sort_order_id=None,
            # Just copy these from the table for now
            spec_id=table_metadata.default_spec_id,
            equality_ids=None,
            key_metadata=None,
            **statistics.to_serialized_dict(),
        )

        return data_file

    executor = ExecutorFactory.get_or_create()
    data_files = executor.map(write_parquet, tasks)

    return iter(data_files)


def bin_pack_arrow_table(tbl: pa.Table, target_file_size: int) -> Iterator[List[pa.RecordBatch]]:
    from pyiceberg.utils.bin_packing import PackingIterator

    avg_row_size_bytes = tbl.nbytes / tbl.num_rows
    target_rows_per_file = target_file_size // avg_row_size_bytes
    batches = tbl.to_batches(max_chunksize=target_rows_per_file)
    bin_packed_record_batches = PackingIterator(
        items=batches,
        target_weight=target_file_size,
        lookback=len(batches),  # ignore lookback
        weight_func=lambda x: x.nbytes,
        largest_bin_first=False,
    )
    return bin_packed_record_batches


def _check_pyarrow_schema_compatible(
    requested_schema: Schema,
    provided_schema: pa.Schema,
    downcast_ns_timestamp_to_us: bool = False,
    format_version: TableVersion = TableProperties.DEFAULT_FORMAT_VERSION,
) -> None:
    """
    Check if the `requested_schema` is compatible with `provided_schema`.

    Two schemas are considered compatible when they are equal in terms of the Iceberg Schema type.

    Raises:
        ValueError: If the schemas are not compatible.
    """
    name_mapping = requested_schema.name_mapping
    try:
        provided_schema = pyarrow_to_schema(
            provided_schema,
            name_mapping=name_mapping,
            downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
            format_version=format_version,
        )
    except ValueError as e:
        provided_schema = _pyarrow_to_schema_without_ids(
            provided_schema, downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us, format_version=format_version
        )
        additional_names = set(provided_schema._name_to_id.keys()) - set(requested_schema._name_to_id.keys())
        raise ValueError(
            f"PyArrow table contains more columns: {', '.join(sorted(additional_names))}. Update the schema first (hint, use union_by_name)."
        ) from e
    _check_schema_compatible(requested_schema, provided_schema)


def parquet_files_to_data_files(io: FileIO, table_metadata: TableMetadata, file_paths: Iterator[str]) -> Iterator[DataFile]:
    for file_path in file_paths:
        data_file = parquet_file_to_data_file(io=io, table_metadata=table_metadata, file_path=file_path)
        yield data_file


def parquet_file_to_data_file(io: FileIO, table_metadata: TableMetadata, file_path: str) -> DataFile:
    input_file = io.new_input(file_path)
    with input_file.open() as input_stream:
        parquet_metadata = pq.read_metadata(input_stream)

    arrow_schema = parquet_metadata.schema.to_arrow_schema()
    if visit_pyarrow(arrow_schema, _HasIds()):
        raise NotImplementedError(
            f"Cannot add file {file_path} because it has field IDs. `add_files` only supports addition of files without field_ids"
        )

    schema = table_metadata.schema()
    _check_pyarrow_schema_compatible(schema, arrow_schema, format_version=table_metadata.format_version)

    statistics = data_file_statistics_from_parquet_metadata(
        parquet_metadata=parquet_metadata,
        stats_columns=compute_statistics_plan(schema, table_metadata.properties),
        parquet_column_mapping=parquet_path_to_id_mapping(schema),
    )
    data_file = DataFile.from_args(
        content=DataFileContent.DATA,
        file_path=file_path,
        file_format=FileFormat.PARQUET,
        partition=statistics.partition(table_metadata.spec(), table_metadata.schema()),
        file_size_in_bytes=len(input_file),
        sort_order_id=None,
        spec_id=table_metadata.default_spec_id,
        equality_ids=None,
        key_metadata=None,
        **statistics.to_serialized_dict(),
    )

    return data_file


ICEBERG_UNCOMPRESSED_CODEC = "uncompressed"
PYARROW_UNCOMPRESSED_CODEC = "none"


def _get_parquet_writer_kwargs(table_properties: Properties) -> Dict[str, Any]:
    from pyiceberg.table import TableProperties

    for key_pattern in [
        TableProperties.PARQUET_ROW_GROUP_SIZE_BYTES,
        TableProperties.PARQUET_BLOOM_FILTER_MAX_BYTES,
        f"{TableProperties.PARQUET_BLOOM_FILTER_COLUMN_ENABLED_PREFIX}.*",
    ]:
        if unsupported_keys := fnmatch.filter(table_properties, key_pattern):
            warnings.warn(f"Parquet writer option(s) {unsupported_keys} not implemented")

    compression_codec = table_properties.get(TableProperties.PARQUET_COMPRESSION, TableProperties.PARQUET_COMPRESSION_DEFAULT)
    compression_level = property_as_int(
        properties=table_properties,
        property_name=TableProperties.PARQUET_COMPRESSION_LEVEL,
        default=TableProperties.PARQUET_COMPRESSION_LEVEL_DEFAULT,
    )
    if compression_codec == ICEBERG_UNCOMPRESSED_CODEC:
        compression_codec = PYARROW_UNCOMPRESSED_CODEC

    return {
        "compression": compression_codec,
        "compression_level": compression_level,
        "data_page_size": property_as_int(
            properties=table_properties,
            property_name=TableProperties.PARQUET_PAGE_SIZE_BYTES,
            default=TableProperties.PARQUET_PAGE_SIZE_BYTES_DEFAULT,
        ),
        "dictionary_pagesize_limit": property_as_int(
            properties=table_properties,
            property_name=TableProperties.PARQUET_DICT_SIZE_BYTES,
            default=TableProperties.PARQUET_DICT_SIZE_BYTES_DEFAULT,
        ),
        "write_batch_size": property_as_int(
            properties=table_properties,
            property_name=TableProperties.PARQUET_PAGE_ROW_LIMIT,
            default=TableProperties.PARQUET_PAGE_ROW_LIMIT_DEFAULT,
        ),
    }


def _dataframe_to_data_files(
    table_metadata: TableMetadata,
    df: pa.Table,
    io: FileIO,
    write_uuid: Optional[uuid.UUID] = None,
    counter: Optional[itertools.count[int]] = None,
) -> Iterable[DataFile]:
    """Convert a PyArrow table into a DataFile.

    Returns:
        An iterable that supplies datafiles that represent the table.
    """
    from pyiceberg.table import DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE, TableProperties, WriteTask

    counter = counter or itertools.count(0)
    write_uuid = write_uuid or uuid.uuid4()
    target_file_size: int = property_as_int(  # type: ignore  # The property is set with non-None value.
        properties=table_metadata.properties,
        property_name=TableProperties.WRITE_TARGET_FILE_SIZE_BYTES,
        default=TableProperties.WRITE_TARGET_FILE_SIZE_BYTES_DEFAULT,
    )
    name_mapping = table_metadata.schema().name_mapping
    downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
    task_schema = pyarrow_to_schema(
        df.schema,
        name_mapping=name_mapping,
        downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
        format_version=table_metadata.format_version,
    )

    if table_metadata.spec().is_unpartitioned():
        yield from write_file(
            io=io,
            table_metadata=table_metadata,
            tasks=iter(
                [
                    WriteTask(write_uuid=write_uuid, task_id=next(counter), record_batches=batches, schema=task_schema)
                    for batches in bin_pack_arrow_table(df, target_file_size)
                ]
            ),
        )
    else:
        partitions = _determine_partitions(spec=table_metadata.spec(), schema=table_metadata.schema(), arrow_table=df)
        yield from write_file(
            io=io,
            table_metadata=table_metadata,
            tasks=iter(
                [
                    WriteTask(
                        write_uuid=write_uuid,
                        task_id=next(counter),
                        record_batches=batches,
                        partition_key=partition.partition_key,
                        schema=task_schema,
                    )
                    for partition in partitions
                    for batches in bin_pack_arrow_table(partition.arrow_table_partition, target_file_size)
                ]
            ),
        )


@dataclass(frozen=True)
class _TablePartition:
    partition_key: PartitionKey
    arrow_table_partition: pa.Table


def _determine_partitions(spec: PartitionSpec, schema: Schema, arrow_table: pa.Table) -> List[_TablePartition]:
    """Based on the iceberg table partition spec, filter the arrow table into partitions with their keys.

    Example:
    Input:
    An arrow table with partition key of ['n_legs', 'year'] and with data of
    {'year': [2020, 2022, 2022, 2021, 2022, 2022, 2022, 2019, 2021],
     'n_legs': [2, 2, 2, 4, 4, 4, 4, 5, 100],
     'animal': ["Flamingo", "Parrot", "Parrot", "Dog", "Horse", "Horse", "Horse","Brittle stars", "Centipede"]}.
    The algorithm:
    - We determine the set of unique partition keys
    - Then we produce a set of partitions by filtering on each of the combinations
    - We combine the chunks to create a copy to avoid GIL congestion on the original table
    """
    # Assign unique names to columns where the partition transform has been applied
    # to avoid conflicts
    partition_fields = [f"_partition_{field.name}" for field in spec.fields]

    for partition, name in zip(spec.fields, partition_fields):
        source_field = schema.find_field(partition.source_id)
        full_field_name = schema.find_column_name(partition.source_id)
        if full_field_name is None:
            raise ValueError(f"Could not find column name for field ID: {partition.source_id}")
        field_array = _get_field_from_arrow_table(arrow_table, full_field_name)
        arrow_table = arrow_table.append_column(name, partition.transform.pyarrow_transform(source_field.field_type)(field_array))

    unique_partition_fields = arrow_table.select(partition_fields).group_by(partition_fields).aggregate([])

    table_partitions = []
    # TODO: As a next step, we could also play around with yielding instead of materializing the full list
    for unique_partition in unique_partition_fields.to_pylist():
        partition_key = PartitionKey(
            field_values=[
                PartitionFieldValue(field=field, value=unique_partition[name])
                for field, name in zip(spec.fields, partition_fields)
            ],
            partition_spec=spec,
            schema=schema,
        )
        filtered_table = arrow_table.filter(
            functools.reduce(
                operator.and_,
                [
                    (
                        pc.field(partition_field_name) == unique_partition[partition_field_name]
                        if unique_partition[partition_field_name] is not None
                        else pc.field(partition_field_name).is_null()
                    )
                    for field, partition_field_name in zip(spec.fields, partition_fields)
                ],
            )
        )
        filtered_table = filtered_table.drop_columns(partition_fields)

        # The combine_chunks seems to be counter-intuitive to do, but it actually returns
        # fresh buffers that don't interfere with each other when it is written out to file
        table_partitions.append(
            _TablePartition(partition_key=partition_key, arrow_table_partition=filtered_table.combine_chunks())
        )

    return table_partitions


def _get_field_from_arrow_table(arrow_table: pa.Table, field_path: str) -> pa.Array:
    """Get a field from an Arrow table, supporting both literal field names and nested field paths.

    This function handles two cases:
    1. Literal field names that may contain dots (e.g., "some.id")
    2. Nested field paths using dot notation (e.g., "bar.baz" for nested access)

    Args:
        arrow_table: The Arrow table containing the field
        field_path: Field name or dot-separated path

    Returns:
        The field as a PyArrow Array

    Raises:
        KeyError: If the field path cannot be resolved
    """
    # Try exact column name match (handles field names containing literal dots)
    if field_path in arrow_table.column_names:
        return arrow_table[field_path]

    # If not found as exact name, treat as nested field path
    path_parts = field_path.split(".")
    # Get the struct column from the table (e.g., "bar" from "bar.baz")
    field_array = arrow_table[path_parts[0]]
    # Navigate into the struct using the remaining path parts
    return pc.struct_field(field_array, path_parts[1:])
