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
"""FileIO implementation for reading and writing table files that uses fsspec compatible filesystems."""

import abc
import errno
import json
import logging
import os
import threading
from collections.abc import Callable
from copy import copy
from functools import lru_cache
from typing import (
    TYPE_CHECKING,
    Any,
)
from urllib.parse import urlparse

import requests
from fsspec import AbstractFileSystem
from fsspec.implementations.local import LocalFileSystem
from requests import HTTPError

from pyiceberg.catalog import TOKEN, URI
from pyiceberg.catalog.rest.auth import AUTH_MANAGER
from pyiceberg.exceptions import SignError
from pyiceberg.io import (
    ADLS_ACCOUNT_HOST,
    ADLS_ACCOUNT_KEY,
    ADLS_ACCOUNT_NAME,
    ADLS_ANON,
    ADLS_CLIENT_ID,
    ADLS_CLIENT_SECRET,
    ADLS_CONNECTION_STRING,
    ADLS_CREDENTIAL,
    ADLS_SAS_TOKEN,
    ADLS_TENANT_ID,
    ADLS_TOKEN,
    AWS_ACCESS_KEY_ID,
    AWS_PROFILE_NAME,
    AWS_REGION,
    AWS_SECRET_ACCESS_KEY,
    AWS_SESSION_TOKEN,
    GCS_ACCESS,
    GCS_CACHE_TIMEOUT,
    GCS_CONSISTENCY,
    GCS_DEFAULT_LOCATION,
    GCS_PROJECT_ID,
    GCS_REQUESTER_PAYS,
    GCS_SERVICE_HOST,
    GCS_SESSION_KWARGS,
    GCS_TOKEN,
    GCS_VERSION_AWARE,
    HF_ENDPOINT,
    HF_TOKEN,
    S3_ACCESS_KEY_ID,
    S3_ANONYMOUS,
    S3_CONNECT_TIMEOUT,
    S3_ENDPOINT,
    S3_FORCE_VIRTUAL_ADDRESSING,
    S3_PROFILE_NAME,
    S3_PROXY_URI,
    S3_REGION,
    S3_REQUEST_TIMEOUT,
    S3_SECRET_ACCESS_KEY,
    S3_SESSION_TOKEN,
    S3_SIGNER,
    S3_SIGNER_ENDPOINT,
    S3_SIGNER_ENDPOINT_DEFAULT,
    S3_SIGNER_URI,
    FileIO,
    InputFile,
    InputStream,
    OutputFile,
    OutputStream,
)
from pyiceberg.typedef import Properties
from pyiceberg.types import strtobool
from pyiceberg.utils.properties import get_first_property_value, get_header_properties, property_as_bool

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from botocore.awsrequest import AWSRequest


class S3RequestSigner(abc.ABC):
    """Abstract base class for S3 request signers."""

    properties: Properties

    def __init__(self, properties: Properties) -> None:
        self.properties = properties

    @abc.abstractmethod
    def __call__(self, request: "AWSRequest", **_: Any) -> None:
        pass


class S3V4RestSigner(S3RequestSigner):
    """An S3 request signer that uses an external REST signing service to sign requests."""

    _session: requests.Session

    def __init__(self, properties: Properties) -> None:
        super().__init__(properties)
        self._session = requests.Session()

    def __call__(self, request: "AWSRequest", **_: Any) -> None:
        signer_url = self.properties.get(S3_SIGNER_URI, self.properties[URI]).rstrip("/")  # type: ignore
        signer_endpoint = self.properties.get(S3_SIGNER_ENDPOINT, S3_SIGNER_ENDPOINT_DEFAULT)

        signer_headers: dict[str, str] = {}

        auth_header: str | None = None
        if auth_manager := self.properties.get(AUTH_MANAGER):
            auth_header = auth_manager.auth_header()
        elif token := self.properties.get(TOKEN):
            auth_header = f"Bearer {token}"

        if auth_header:
            signer_headers["Authorization"] = auth_header

        signer_headers.update(get_header_properties(self.properties))

        signer_body = {
            "method": request.method,
            "region": request.context["client_region"],
            "uri": request.url,
            "headers": {key: [val] for key, val in request.headers.items()},
        }

        response = self._session.post(f"{signer_url}/{signer_endpoint.strip()}", headers=signer_headers, json=signer_body)
        try:
            response.raise_for_status()
            response_json = response.json()
        except HTTPError as e:
            raise SignError(f"Failed to sign request {response.status_code}: {signer_body}") from e

        for key, value in response_json["headers"].items():
            request.headers.add_header(key, ", ".join(value))

        request.url = response_json["uri"]


SIGNERS: dict[str, type[S3RequestSigner]] = {"S3V4RestSigner": S3V4RestSigner}


def _file(_: Properties) -> LocalFileSystem:
    return LocalFileSystem(auto_mkdir=True)


def _s3(properties: Properties) -> AbstractFileSystem:
    from s3fs import S3FileSystem

    client_kwargs = {
        "endpoint_url": properties.get(S3_ENDPOINT),
        "aws_access_key_id": get_first_property_value(properties, S3_ACCESS_KEY_ID, AWS_ACCESS_KEY_ID),
        "aws_secret_access_key": get_first_property_value(properties, S3_SECRET_ACCESS_KEY, AWS_SECRET_ACCESS_KEY),
        "aws_session_token": get_first_property_value(properties, S3_SESSION_TOKEN, AWS_SESSION_TOKEN),
        "region_name": get_first_property_value(properties, S3_REGION, AWS_REGION),
    }
    config_kwargs = {}
    register_events: dict[str, Callable[[AWSRequest], None]] = {}

    if signer := properties.get(S3_SIGNER):
        logger.info("Loading signer %s", signer)
        if signer_cls := SIGNERS.get(signer):
            signer = signer_cls(properties)
            register_events["before-sign.s3"] = signer

            # Disable the AWS Signer
            from botocore import UNSIGNED

            config_kwargs["signature_version"] = UNSIGNED
        else:
            raise ValueError(f"Signer not available: {signer}")

    if proxy_uri := properties.get(S3_PROXY_URI):
        config_kwargs["proxies"] = {"http": proxy_uri, "https": proxy_uri}

    if connect_timeout := properties.get(S3_CONNECT_TIMEOUT):
        config_kwargs["connect_timeout"] = float(connect_timeout)

    if request_timeout := properties.get(S3_REQUEST_TIMEOUT):
        config_kwargs["read_timeout"] = float(request_timeout)

    if _force_virtual_addressing := properties.get(S3_FORCE_VIRTUAL_ADDRESSING):
        config_kwargs["s3"] = {"addressing_style": "virtual"}

    if s3_anonymous := properties.get(S3_ANONYMOUS):
        anon = strtobool(s3_anonymous)
    else:
        anon = False

    s3_fs_kwargs = {
        "anon": anon,
        "client_kwargs": client_kwargs,
        "config_kwargs": config_kwargs,
    }

    if profile_name := get_first_property_value(properties, S3_PROFILE_NAME, AWS_PROFILE_NAME):
        s3_fs_kwargs["profile"] = profile_name

    fs = S3FileSystem(**s3_fs_kwargs)

    for event_name, event_function in register_events.items():
        fs.s3.meta.events.unregister(event_name, unique_id=1925)
        fs.s3.meta.events.register_last(event_name, event_function, unique_id=1925)

    return fs


def _gs(properties: Properties) -> AbstractFileSystem:
    # https://gcsfs.readthedocs.io/en/latest/api.html#gcsfs.core.GCSFileSystem
    from gcsfs import GCSFileSystem

    return GCSFileSystem(
        project=properties.get(GCS_PROJECT_ID),
        access=properties.get(GCS_ACCESS, "full_control"),
        token=properties.get(GCS_TOKEN),
        consistency=properties.get(GCS_CONSISTENCY, "none"),
        cache_timeout=properties.get(GCS_CACHE_TIMEOUT),
        requester_pays=property_as_bool(properties, GCS_REQUESTER_PAYS, False),
        session_kwargs=json.loads(properties.get(GCS_SESSION_KWARGS, "{}")),
        endpoint_url=properties.get(GCS_SERVICE_HOST),
        default_location=properties.get(GCS_DEFAULT_LOCATION),
        version_aware=property_as_bool(properties, GCS_VERSION_AWARE, False),
    )


def _adls(properties: Properties) -> AbstractFileSystem:
    # https://fsspec.github.io/adlfs/api/

    from adlfs import AzureBlobFileSystem
    from azure.core.credentials import AccessToken
    from azure.core.credentials_async import AsyncTokenCredential

    for key, sas_token in {
        key.replace(f"{ADLS_SAS_TOKEN}.", ""): value for key, value in properties.items() if key.startswith(f"{ADLS_SAS_TOKEN}.")
    }.items():
        if ADLS_ACCOUNT_NAME not in properties:
            properties[ADLS_ACCOUNT_NAME] = key.split(".")[0]
        if ADLS_SAS_TOKEN not in properties:
            properties[ADLS_SAS_TOKEN] = sas_token

    class StaticTokenCredential(AsyncTokenCredential):
        _DEFAULT_EXPIRY_SECONDS = 3600

        def __init__(self, token_string: str) -> None:
            self._token = token_string

        async def get_token(self, *scopes: str, **kwargs: Any) -> AccessToken:
            import time

            # Set expiration 1 hour from now
            expires_on = int(time.time()) + self._DEFAULT_EXPIRY_SECONDS
            return AccessToken(self._token, expires_on)

    if token := properties.get(ADLS_TOKEN):
        credential = StaticTokenCredential(token)
    else:
        credential = properties.get(ADLS_CREDENTIAL)  # type: ignore

    return AzureBlobFileSystem(
        connection_string=properties.get(ADLS_CONNECTION_STRING),
        credential=credential,
        account_name=properties.get(ADLS_ACCOUNT_NAME),
        account_key=properties.get(ADLS_ACCOUNT_KEY),
        sas_token=properties.get(ADLS_SAS_TOKEN),
        tenant_id=properties.get(ADLS_TENANT_ID),
        client_id=properties.get(ADLS_CLIENT_ID),
        client_secret=properties.get(ADLS_CLIENT_SECRET),
        account_host=properties.get(ADLS_ACCOUNT_HOST),
        anon=properties.get(ADLS_ANON),
    )


def _hf(properties: Properties) -> AbstractFileSystem:
    from huggingface_hub import HfFileSystem

    return HfFileSystem(
        endpoint=properties.get(HF_ENDPOINT),
        token=properties.get(HF_TOKEN),
    )


SCHEME_TO_FS = {
    "": _file,
    "file": _file,
    "s3": _s3,
    "s3a": _s3,
    "s3n": _s3,
    "abfs": _adls,
    "abfss": _adls,
    "gs": _gs,
    "gcs": _gs,
    "hf": _hf,
}


class FsspecInputFile(InputFile):
    """An input file implementation for the FsspecFileIO.

    Args:
        location (str): A URI to a file location.
        fs (AbstractFileSystem): An fsspec filesystem instance.
    """

    def __init__(self, location: str, fs: AbstractFileSystem):
        self._fs = fs
        super().__init__(location=location)

    def __len__(self) -> int:
        """Return the total length of the file, in bytes."""
        object_info = self._fs.info(self.location)
        if size := object_info.get("Size"):
            return size
        elif size := object_info.get("size"):
            return size
        raise RuntimeError(f"Cannot retrieve object info: {self.location}")

    def exists(self) -> bool:
        """Check whether the location exists."""
        return self._fs.lexists(self.location)

    def open(self, seekable: bool = True) -> InputStream:
        """Create an input stream for reading the contents of the file.

        Args:
            seekable: If the stream should support seek, or if it is consumed sequential.

        Returns:
            OpenFile: An fsspec compliant file-like object.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        try:
            return self._fs.open(self.location, "rb")
        except FileNotFoundError as e:
            # To have a consistent error handling experience, make sure exception contains missing file location.
            raise e if e.filename else FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.location) from e


class FsspecOutputFile(OutputFile):
    """An output file implementation for the FsspecFileIO.

    Args:
        location (str): A URI to a file location.
        fs (AbstractFileSystem): An fsspec filesystem instance.
    """

    def __init__(self, location: str, fs: AbstractFileSystem):
        self._fs = fs
        super().__init__(location=location)

    def __len__(self) -> int:
        """Return the total length of the file, in bytes."""
        object_info = self._fs.info(self.location)
        if size := object_info.get("Size"):
            return size
        elif size := object_info.get("size"):
            return size
        raise RuntimeError(f"Cannot retrieve object info: {self.location}")

    def exists(self) -> bool:
        """Check whether the location exists."""
        return self._fs.lexists(self.location)

    def create(self, overwrite: bool = False) -> OutputStream:
        """Create an output stream for reading the contents of the file.

        Args:
            overwrite (bool): Whether to overwrite the file if it already exists.

        Returns:
            OpenFile: An fsspec compliant file-like object.

        Raises:
            FileExistsError: If the file already exists at the location and overwrite is set to False.

        Note:
            If overwrite is set to False, a check is first performed to verify that the file does not exist.
            This is not thread-safe and a possibility does exist that the file can be created by a concurrent
            process after the existence check yet before the output stream is created. In such a case, the default
            behavior will truncate the contents of the existing file when opening the output stream.
        """
        if not overwrite and self.exists():
            raise FileExistsError(f"Cannot create file, file already exists: {self.location}")
        return self._fs.open(self.location, "wb")

    def to_input_file(self) -> FsspecInputFile:
        """Return a new FsspecInputFile for the location at `self.location`."""
        return FsspecInputFile(location=self.location, fs=self._fs)


class FsspecFileIO(FileIO):
    """A FileIO implementation that uses fsspec."""

    def __init__(self, properties: Properties):
        self._scheme_to_fs = {}
        self._scheme_to_fs.update(SCHEME_TO_FS)
        self._thread_locals = threading.local()
        super().__init__(properties=properties)

    def new_input(self, location: str) -> FsspecInputFile:
        """Get an FsspecInputFile instance to read bytes from the file at the given location.

        Args:
            location (str): A URI or a path to a local file.

        Returns:
            FsspecInputFile: An FsspecInputFile instance for the given location.
        """
        uri = urlparse(location)
        fs = self.get_fs(uri.scheme)
        return FsspecInputFile(location=location, fs=fs)

    def new_output(self, location: str) -> FsspecOutputFile:
        """Get an FsspecOutputFile instance to write bytes to the file at the given location.

        Args:
            location (str): A URI or a path to a local file.

        Returns:
            FsspecOutputFile: An FsspecOutputFile instance for the given location.
        """
        uri = urlparse(location)
        fs = self.get_fs(uri.scheme)
        return FsspecOutputFile(location=location, fs=fs)

    def delete(self, location: str | InputFile | OutputFile) -> None:
        """Delete the file at the given location.

        Args:
            location (Union[str, InputFile, OutputFile]): The URI to the file--if an InputFile instance or an
                OutputFile instance is provided, the location attribute for that instance is used as the location
                to delete.
        """
        if isinstance(location, (InputFile, OutputFile)):
            str_location = location.location  # Use InputFile or OutputFile location
        else:
            str_location = location

        uri = urlparse(str_location)
        fs = self.get_fs(uri.scheme)
        fs.rm(str_location)

    def get_fs(self, scheme: str) -> AbstractFileSystem:
        """Get a filesystem for a specific scheme, cached per thread."""
        if not hasattr(self._thread_locals, "get_fs_cached"):
            self._thread_locals.get_fs_cached = lru_cache(self._get_fs)

        return self._thread_locals.get_fs_cached(scheme)

    def _get_fs(self, scheme: str) -> AbstractFileSystem:
        """Get a filesystem for a specific scheme."""
        if scheme not in self._scheme_to_fs:
            raise ValueError(f"No registered filesystem for scheme: {scheme}")
        return self._scheme_to_fs[scheme](self.properties)

    def __getstate__(self) -> dict[str, Any]:
        """Create a dictionary of the FsSpecFileIO fields used when pickling."""
        fileio_copy = copy(self.__dict__)
        del fileio_copy["_thread_locals"]
        return fileio_copy

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Deserialize the state into a FsSpecFileIO instance."""
        self.__dict__ = state
        self._thread_locals = threading.local()
