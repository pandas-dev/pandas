# -*- coding: utf-8 -*-
# Copyright 2025 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import logging as std_logging
from collections import OrderedDict
import re
from typing import (
    Dict,
    Callable,
    Mapping,
    MutableMapping,
    MutableSequence,
    Optional,
    AsyncIterable,
    Awaitable,
    AsyncIterator,
    Sequence,
    Tuple,
    Type,
    Union,
)

from google.cloud._storage_v2 import gapic_version as package_version

from google.api_core.client_options import ClientOptions
from google.api_core import exceptions as core_exceptions
from google.api_core import gapic_v1
from google.api_core import retry_async as retries
from google.auth import credentials as ga_credentials  # type: ignore
from google.oauth2 import service_account  # type: ignore
import google.protobuf


try:
    OptionalRetry = Union[retries.AsyncRetry, gapic_v1.method._MethodDefault, None]
except AttributeError:  # pragma: NO COVER
    OptionalRetry = Union[retries.AsyncRetry, object, None]  # type: ignore

from google.cloud._storage_v2.services.storage import pagers
from google.cloud._storage_v2.types import storage
from google.iam.v1 import iam_policy_pb2  # type: ignore
from google.iam.v1 import policy_pb2  # type: ignore
from google.longrunning import operations_pb2  # type: ignore
from google.protobuf import field_mask_pb2  # type: ignore
from google.protobuf import timestamp_pb2  # type: ignore
from .transports.base import StorageTransport, DEFAULT_CLIENT_INFO
from .transports.grpc_asyncio import StorageGrpcAsyncIOTransport
from .client import StorageClient

try:
    from google.api_core import client_logging  # type: ignore

    CLIENT_LOGGING_SUPPORTED = True  # pragma: NO COVER
except ImportError:  # pragma: NO COVER
    CLIENT_LOGGING_SUPPORTED = False

_LOGGER = std_logging.getLogger(__name__)


class StorageAsyncClient:
    """API Overview and Naming Syntax
    ------------------------------

    The Cloud Storage gRPC API allows applications to read and write
    data through the abstractions of buckets and objects. For a
    description of these abstractions please see
    https://cloud.google.com/storage/docs.

    Resources are named as follows:

    - Projects are referred to as they are defined by the Resource
      Manager API, using strings like ``projects/123456`` or
      ``projects/my-string-id``.

    - Buckets are named using string names of the form:
      ``projects/{project}/buckets/{bucket}`` For globally unique
      buckets, ``_`` may be substituted for the project.

    - Objects are uniquely identified by their name along with the name
      of the bucket they belong to, as separate strings in this API. For
      example:

      ReadObjectRequest { bucket: 'projects/\_/buckets/my-bucket'
      object: 'my-object' } Note that object names can contain ``/``
      characters, which are treated as any other character (no special
      directory semantics).
    """

    _client: StorageClient

    # Copy defaults from the synchronous client for use here.
    # Note: DEFAULT_ENDPOINT is deprecated. Use _DEFAULT_ENDPOINT_TEMPLATE instead.
    DEFAULT_ENDPOINT = StorageClient.DEFAULT_ENDPOINT
    DEFAULT_MTLS_ENDPOINT = StorageClient.DEFAULT_MTLS_ENDPOINT
    _DEFAULT_ENDPOINT_TEMPLATE = StorageClient._DEFAULT_ENDPOINT_TEMPLATE
    _DEFAULT_UNIVERSE = StorageClient._DEFAULT_UNIVERSE

    bucket_path = staticmethod(StorageClient.bucket_path)
    parse_bucket_path = staticmethod(StorageClient.parse_bucket_path)
    crypto_key_path = staticmethod(StorageClient.crypto_key_path)
    parse_crypto_key_path = staticmethod(StorageClient.parse_crypto_key_path)
    common_billing_account_path = staticmethod(
        StorageClient.common_billing_account_path
    )
    parse_common_billing_account_path = staticmethod(
        StorageClient.parse_common_billing_account_path
    )
    common_folder_path = staticmethod(StorageClient.common_folder_path)
    parse_common_folder_path = staticmethod(StorageClient.parse_common_folder_path)
    common_organization_path = staticmethod(StorageClient.common_organization_path)
    parse_common_organization_path = staticmethod(
        StorageClient.parse_common_organization_path
    )
    common_project_path = staticmethod(StorageClient.common_project_path)
    parse_common_project_path = staticmethod(StorageClient.parse_common_project_path)
    common_location_path = staticmethod(StorageClient.common_location_path)
    parse_common_location_path = staticmethod(StorageClient.parse_common_location_path)

    @classmethod
    def from_service_account_info(cls, info: dict, *args, **kwargs):
        """Creates an instance of this client using the provided credentials
            info.

        Args:
            info (dict): The service account private key info.
            args: Additional arguments to pass to the constructor.
            kwargs: Additional arguments to pass to the constructor.

        Returns:
            StorageAsyncClient: The constructed client.
        """
        return StorageClient.from_service_account_info.__func__(StorageAsyncClient, info, *args, **kwargs)  # type: ignore

    @classmethod
    def from_service_account_file(cls, filename: str, *args, **kwargs):
        """Creates an instance of this client using the provided credentials
            file.

        Args:
            filename (str): The path to the service account private key json
                file.
            args: Additional arguments to pass to the constructor.
            kwargs: Additional arguments to pass to the constructor.

        Returns:
            StorageAsyncClient: The constructed client.
        """
        return StorageClient.from_service_account_file.__func__(StorageAsyncClient, filename, *args, **kwargs)  # type: ignore

    from_service_account_json = from_service_account_file

    @classmethod
    def get_mtls_endpoint_and_cert_source(
        cls, client_options: Optional[ClientOptions] = None
    ):
        """Return the API endpoint and client cert source for mutual TLS.

        The client cert source is determined in the following order:
        (1) if `GOOGLE_API_USE_CLIENT_CERTIFICATE` environment variable is not "true", the
        client cert source is None.
        (2) if `client_options.client_cert_source` is provided, use the provided one; if the
        default client cert source exists, use the default one; otherwise the client cert
        source is None.

        The API endpoint is determined in the following order:
        (1) if `client_options.api_endpoint` if provided, use the provided one.
        (2) if `GOOGLE_API_USE_CLIENT_CERTIFICATE` environment variable is "always", use the
        default mTLS endpoint; if the environment variable is "never", use the default API
        endpoint; otherwise if client cert source exists, use the default mTLS endpoint, otherwise
        use the default API endpoint.

        More details can be found at https://google.aip.dev/auth/4114.

        Args:
            client_options (google.api_core.client_options.ClientOptions): Custom options for the
                client. Only the `api_endpoint` and `client_cert_source` properties may be used
                in this method.

        Returns:
            Tuple[str, Callable[[], Tuple[bytes, bytes]]]: returns the API endpoint and the
                client cert source to use.

        Raises:
            google.auth.exceptions.MutualTLSChannelError: If any errors happen.
        """
        return StorageClient.get_mtls_endpoint_and_cert_source(client_options)  # type: ignore

    @property
    def transport(self) -> StorageTransport:
        """Returns the transport used by the client instance.

        Returns:
            StorageTransport: The transport used by the client instance.
        """
        return self._client.transport

    @property
    def api_endpoint(self):
        """Return the API endpoint used by the client instance.

        Returns:
            str: The API endpoint used by the client instance.
        """
        return self._client._api_endpoint

    @property
    def universe_domain(self) -> str:
        """Return the universe domain used by the client instance.

        Returns:
            str: The universe domain used
                by the client instance.
        """
        return self._client._universe_domain

    get_transport_class = StorageClient.get_transport_class

    def __init__(
        self,
        *,
        credentials: Optional[ga_credentials.Credentials] = None,
        transport: Optional[
            Union[str, StorageTransport, Callable[..., StorageTransport]]
        ] = "grpc_asyncio",
        client_options: Optional[ClientOptions] = None,
        client_info: gapic_v1.client_info.ClientInfo = DEFAULT_CLIENT_INFO,
    ) -> None:
        """Instantiates the storage async client.

        Args:
            credentials (Optional[google.auth.credentials.Credentials]): The
                authorization credentials to attach to requests. These
                credentials identify the application to the service; if none
                are specified, the client will attempt to ascertain the
                credentials from the environment.
            transport (Optional[Union[str,StorageTransport,Callable[..., StorageTransport]]]):
                The transport to use, or a Callable that constructs and returns a new transport to use.
                If a Callable is given, it will be called with the same set of initialization
                arguments as used in the StorageTransport constructor.
                If set to None, a transport is chosen automatically.
            client_options (Optional[Union[google.api_core.client_options.ClientOptions, dict]]):
                Custom options for the client.

                1. The ``api_endpoint`` property can be used to override the
                default endpoint provided by the client when ``transport`` is
                not explicitly provided. Only if this property is not set and
                ``transport`` was not explicitly provided, the endpoint is
                determined by the GOOGLE_API_USE_MTLS_ENDPOINT environment
                variable, which have one of the following values:
                "always" (always use the default mTLS endpoint), "never" (always
                use the default regular endpoint) and "auto" (auto-switch to the
                default mTLS endpoint if client certificate is present; this is
                the default value).

                2. If the GOOGLE_API_USE_CLIENT_CERTIFICATE environment variable
                is "true", then the ``client_cert_source`` property can be used
                to provide a client certificate for mTLS transport. If
                not provided, the default SSL client certificate will be used if
                present. If GOOGLE_API_USE_CLIENT_CERTIFICATE is "false" or not
                set, no client certificate will be used.

                3. The ``universe_domain`` property can be used to override the
                default "googleapis.com" universe. Note that ``api_endpoint``
                property still takes precedence; and ``universe_domain`` is
                currently not supported for mTLS.

            client_info (google.api_core.gapic_v1.client_info.ClientInfo):
                The client info used to send a user-agent string along with
                API requests. If ``None``, then default info will be used.
                Generally, you only need to set this if you're developing
                your own client library.

        Raises:
            google.auth.exceptions.MutualTlsChannelError: If mutual TLS transport
                creation failed for any reason.
        """
        self._client = StorageClient(
            credentials=credentials,
            transport=transport,
            client_options=client_options,
            client_info=client_info,
        )

        if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
            std_logging.DEBUG
        ):  # pragma: NO COVER
            _LOGGER.debug(
                "Created client `google.storage_v2.StorageAsyncClient`.",
                extra={
                    "serviceName": "google.storage.v2.Storage",
                    "universeDomain": getattr(
                        self._client._transport._credentials, "universe_domain", ""
                    ),
                    "credentialsType": f"{type(self._client._transport._credentials).__module__}.{type(self._client._transport._credentials).__qualname__}",
                    "credentialsInfo": getattr(
                        self.transport._credentials, "get_cred_info", lambda: None
                    )(),
                }
                if hasattr(self._client._transport, "_credentials")
                else {
                    "serviceName": "google.storage.v2.Storage",
                    "credentialsType": None,
                },
            )

    async def delete_bucket(
        self,
        request: Optional[Union[storage.DeleteBucketRequest, dict]] = None,
        *,
        name: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> None:
        r"""Permanently deletes an empty bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_delete_bucket():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.DeleteBucketRequest(
                    name="name_value",
                )

                # Make the request
                await client.delete_bucket(request=request)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.DeleteBucketRequest, dict]]):
                The request object. Request message for DeleteBucket.
            name (:class:`str`):
                Required. Name of a bucket to delete.
                This corresponds to the ``name`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [name]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.DeleteBucketRequest):
            request = storage.DeleteBucketRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if name is not None:
            request.name = name

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.delete_bucket
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.name)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

    async def get_bucket(
        self,
        request: Optional[Union[storage.GetBucketRequest, dict]] = None,
        *,
        name: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Bucket:
        r"""Returns metadata for the specified bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_get_bucket():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.GetBucketRequest(
                    name="name_value",
                )

                # Make the request
                response = await client.get_bucket(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.GetBucketRequest, dict]]):
                The request object. Request message for GetBucket.
            name (:class:`str`):
                Required. Name of a bucket.
                This corresponds to the ``name`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Bucket:
                A bucket.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [name]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.GetBucketRequest):
            request = storage.GetBucketRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if name is not None:
            request.name = name

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.get_bucket
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.name)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def create_bucket(
        self,
        request: Optional[Union[storage.CreateBucketRequest, dict]] = None,
        *,
        parent: Optional[str] = None,
        bucket: Optional[storage.Bucket] = None,
        bucket_id: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Bucket:
        r"""Creates a new bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_create_bucket():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.CreateBucketRequest(
                    parent="parent_value",
                    bucket_id="bucket_id_value",
                )

                # Make the request
                response = await client.create_bucket(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.CreateBucketRequest, dict]]):
                The request object. Request message for CreateBucket.
            parent (:class:`str`):
                Required. The project to which this bucket will belong.
                This field must either be empty or ``projects/_``. The
                project ID that owns this bucket should be specified in
                the ``bucket.project`` field.

                This corresponds to the ``parent`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            bucket (:class:`google.cloud._storage_v2.types.Bucket`):
                Optional. Properties of the new bucket being inserted.
                The name of the bucket is specified in the ``bucket_id``
                field. Populating ``bucket.name`` field will result in
                an error. The project of the bucket must be specified in
                the ``bucket.project`` field. This field must be in
                ``projects/{projectIdentifier}`` format,
                {projectIdentifier} can be the project ID or project
                number. The ``parent`` field must be either empty or
                ``projects/_``.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            bucket_id (:class:`str`):
                Required. The ID to use for this bucket, which will
                become the final component of the bucket's resource
                name. For example, the value ``foo`` might result in a
                bucket with the name ``projects/123456/buckets/foo``.

                This corresponds to the ``bucket_id`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Bucket:
                A bucket.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [parent, bucket, bucket_id]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.CreateBucketRequest):
            request = storage.CreateBucketRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if parent is not None:
            request.parent = parent
        if bucket is not None:
            request.bucket = bucket
        if bucket_id is not None:
            request.bucket_id = bucket_id

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.create_bucket
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<project>.*)$")
        regex_match = routing_param_regex.match(request.parent)
        if regex_match and regex_match.group("project"):
            header_params["project"] = regex_match.group("project")

        routing_param_regex = re.compile("^(?P<project>.*)$")
        regex_match = routing_param_regex.match(request.bucket.project)
        if regex_match and regex_match.group("project"):
            header_params["project"] = regex_match.group("project")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def list_buckets(
        self,
        request: Optional[Union[storage.ListBucketsRequest, dict]] = None,
        *,
        parent: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> pagers.ListBucketsAsyncPager:
        r"""Retrieves a list of buckets for a given project.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_list_buckets():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.ListBucketsRequest(
                    parent="parent_value",
                )

                # Make the request
                page_result = client.list_buckets(request=request)

                # Handle the response
                async for response in page_result:
                    print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.ListBucketsRequest, dict]]):
                The request object. Request message for ListBuckets.
            parent (:class:`str`):
                Required. The project whose buckets
                we are listing.

                This corresponds to the ``parent`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.services.storage.pagers.ListBucketsAsyncPager:
                The result of a call to
                Buckets.ListBuckets
                Iterating over this object will yield
                results and resolve additional pages
                automatically.

        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [parent]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.ListBucketsRequest):
            request = storage.ListBucketsRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if parent is not None:
            request.parent = parent

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.list_buckets
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<project>.*)$")
        regex_match = routing_param_regex.match(request.parent)
        if regex_match and regex_match.group("project"):
            header_params["project"] = regex_match.group("project")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # This method is paged; wrap the response in a pager, which provides
        # an `__aiter__` convenience method.
        response = pagers.ListBucketsAsyncPager(
            method=rpc,
            request=request,
            response=response,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def lock_bucket_retention_policy(
        self,
        request: Optional[Union[storage.LockBucketRetentionPolicyRequest, dict]] = None,
        *,
        bucket: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Bucket:
        r"""Locks retention policy on a bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_lock_bucket_retention_policy():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.LockBucketRetentionPolicyRequest(
                    bucket="bucket_value",
                    if_metageneration_match=2413,
                )

                # Make the request
                response = await client.lock_bucket_retention_policy(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.LockBucketRetentionPolicyRequest, dict]]):
                The request object. Request message for
                LockBucketRetentionPolicyRequest.
            bucket (:class:`str`):
                Required. Name of a bucket.
                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Bucket:
                A bucket.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.LockBucketRetentionPolicyRequest):
            request = storage.LockBucketRetentionPolicyRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.lock_bucket_retention_policy
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def get_iam_policy(
        self,
        request: Optional[Union[iam_policy_pb2.GetIamPolicyRequest, dict]] = None,
        *,
        resource: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> policy_pb2.Policy:
        r"""Gets the IAM policy for a specified bucket. The ``resource``
        field in the request should be ``projects/_/buckets/{bucket}``
        for a bucket, or
        ``projects/_/buckets/{bucket}/managedFolders/{managedFolder}``
        for a managed folder.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2
            from google.iam.v1 import iam_policy_pb2  # type: ignore

            async def sample_get_iam_policy():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = iam_policy_pb2.GetIamPolicyRequest(
                    resource="resource_value",
                )

                # Make the request
                response = await client.get_iam_policy(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.iam.v1.iam_policy_pb2.GetIamPolicyRequest, dict]]):
                The request object. Request message for ``GetIamPolicy`` method.
            resource (:class:`str`):
                REQUIRED: The resource for which the
                policy is being requested. See the
                operation documentation for the
                appropriate value for this field.

                This corresponds to the ``resource`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.iam.v1.policy_pb2.Policy:
                An Identity and Access Management (IAM) policy, which specifies access
                   controls for Google Cloud resources.

                   A Policy is a collection of bindings. A binding binds
                   one or more members, or principals, to a single role.
                   Principals can be user accounts, service accounts,
                   Google groups, and domains (such as G Suite). A role
                   is a named list of permissions; each role can be an
                   IAM predefined role or a user-created custom role.

                   For some types of Google Cloud resources, a binding
                   can also specify a condition, which is a logical
                   expression that allows access to a resource only if
                   the expression evaluates to true. A condition can add
                   constraints based on attributes of the request, the
                   resource, or both. To learn which resources support
                   conditions in their IAM policies, see the [IAM
                   documentation](https://cloud.google.com/iam/help/conditions/resource-policies).

                   **JSON example:**

                   :literal:``     {       "bindings": [         {           "role": "roles/resourcemanager.organizationAdmin",           "members": [             "user:mike@example.com",             "group:admins@example.com",             "domain:google.com",             "serviceAccount:my-project-id@appspot.gserviceaccount.com"           ]         },         {           "role": "roles/resourcemanager.organizationViewer",           "members": [             "user:eve@example.com"           ],           "condition": {             "title": "expirable access",             "description": "Does not grant access after Sep 2020",             "expression": "request.time <             timestamp('2020-10-01T00:00:00.000Z')",           }         }       ],       "etag": "BwWWja0YfJA=",       "version": 3     }`\ \`

                   **YAML example:**

                   :literal:``     bindings:     - members:       - user:mike@example.com       - group:admins@example.com       - domain:google.com       - serviceAccount:my-project-id@appspot.gserviceaccount.com       role: roles/resourcemanager.organizationAdmin     - members:       - user:eve@example.com       role: roles/resourcemanager.organizationViewer       condition:         title: expirable access         description: Does not grant access after Sep 2020         expression: request.time < timestamp('2020-10-01T00:00:00.000Z')     etag: BwWWja0YfJA=     version: 3`\ \`

                   For a description of IAM and its features, see the
                   [IAM
                   documentation](https://cloud.google.com/iam/docs/).

        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [resource]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - The request isn't a proto-plus wrapped type,
        #   so it must be constructed via keyword expansion.
        if isinstance(request, dict):
            request = iam_policy_pb2.GetIamPolicyRequest(**request)
        elif not request:
            request = iam_policy_pb2.GetIamPolicyRequest(resource=resource)

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.get_iam_policy
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        routing_param_regex = re.compile(
            "^(?P<bucket>projects/[^/]+/buckets/[^/]+)(?:/.*)?$"
        )
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def set_iam_policy(
        self,
        request: Optional[Union[iam_policy_pb2.SetIamPolicyRequest, dict]] = None,
        *,
        resource: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> policy_pb2.Policy:
        r"""Updates an IAM policy for the specified bucket. The ``resource``
        field in the request should be ``projects/_/buckets/{bucket}``
        for a bucket, or
        ``projects/_/buckets/{bucket}/managedFolders/{managedFolder}``
        for a managed folder.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2
            from google.iam.v1 import iam_policy_pb2  # type: ignore

            async def sample_set_iam_policy():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = iam_policy_pb2.SetIamPolicyRequest(
                    resource="resource_value",
                )

                # Make the request
                response = await client.set_iam_policy(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.iam.v1.iam_policy_pb2.SetIamPolicyRequest, dict]]):
                The request object. Request message for ``SetIamPolicy`` method.
            resource (:class:`str`):
                REQUIRED: The resource for which the
                policy is being specified. See the
                operation documentation for the
                appropriate value for this field.

                This corresponds to the ``resource`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.iam.v1.policy_pb2.Policy:
                An Identity and Access Management (IAM) policy, which specifies access
                   controls for Google Cloud resources.

                   A Policy is a collection of bindings. A binding binds
                   one or more members, or principals, to a single role.
                   Principals can be user accounts, service accounts,
                   Google groups, and domains (such as G Suite). A role
                   is a named list of permissions; each role can be an
                   IAM predefined role or a user-created custom role.

                   For some types of Google Cloud resources, a binding
                   can also specify a condition, which is a logical
                   expression that allows access to a resource only if
                   the expression evaluates to true. A condition can add
                   constraints based on attributes of the request, the
                   resource, or both. To learn which resources support
                   conditions in their IAM policies, see the [IAM
                   documentation](https://cloud.google.com/iam/help/conditions/resource-policies).

                   **JSON example:**

                   :literal:``     {       "bindings": [         {           "role": "roles/resourcemanager.organizationAdmin",           "members": [             "user:mike@example.com",             "group:admins@example.com",             "domain:google.com",             "serviceAccount:my-project-id@appspot.gserviceaccount.com"           ]         },         {           "role": "roles/resourcemanager.organizationViewer",           "members": [             "user:eve@example.com"           ],           "condition": {             "title": "expirable access",             "description": "Does not grant access after Sep 2020",             "expression": "request.time <             timestamp('2020-10-01T00:00:00.000Z')",           }         }       ],       "etag": "BwWWja0YfJA=",       "version": 3     }`\ \`

                   **YAML example:**

                   :literal:``     bindings:     - members:       - user:mike@example.com       - group:admins@example.com       - domain:google.com       - serviceAccount:my-project-id@appspot.gserviceaccount.com       role: roles/resourcemanager.organizationAdmin     - members:       - user:eve@example.com       role: roles/resourcemanager.organizationViewer       condition:         title: expirable access         description: Does not grant access after Sep 2020         expression: request.time < timestamp('2020-10-01T00:00:00.000Z')     etag: BwWWja0YfJA=     version: 3`\ \`

                   For a description of IAM and its features, see the
                   [IAM
                   documentation](https://cloud.google.com/iam/docs/).

        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [resource]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - The request isn't a proto-plus wrapped type,
        #   so it must be constructed via keyword expansion.
        if isinstance(request, dict):
            request = iam_policy_pb2.SetIamPolicyRequest(**request)
        elif not request:
            request = iam_policy_pb2.SetIamPolicyRequest(resource=resource)

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.set_iam_policy
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        routing_param_regex = re.compile(
            "^(?P<bucket>projects/[^/]+/buckets/[^/]+)(?:/.*)?$"
        )
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def test_iam_permissions(
        self,
        request: Optional[Union[iam_policy_pb2.TestIamPermissionsRequest, dict]] = None,
        *,
        resource: Optional[str] = None,
        permissions: Optional[MutableSequence[str]] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> iam_policy_pb2.TestIamPermissionsResponse:
        r"""Tests a set of permissions on the given bucket, object, or
        managed folder to see which, if any, are held by the caller. The
        ``resource`` field in the request should be
        ``projects/_/buckets/{bucket}`` for a bucket,
        ``projects/_/buckets/{bucket}/objects/{object}`` for an object,
        or
        ``projects/_/buckets/{bucket}/managedFolders/{managedFolder}``
        for a managed folder.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2
            from google.iam.v1 import iam_policy_pb2  # type: ignore

            async def sample_test_iam_permissions():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = iam_policy_pb2.TestIamPermissionsRequest(
                    resource="resource_value",
                    permissions=['permissions_value1', 'permissions_value2'],
                )

                # Make the request
                response = await client.test_iam_permissions(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.iam.v1.iam_policy_pb2.TestIamPermissionsRequest, dict]]):
                The request object. Request message for ``TestIamPermissions`` method.
            resource (:class:`str`):
                REQUIRED: The resource for which the
                policy detail is being requested. See
                the operation documentation for the
                appropriate value for this field.

                This corresponds to the ``resource`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            permissions (:class:`MutableSequence[str]`):
                The set of permissions to check for the ``resource``.
                Permissions with wildcards (such as '*' or 'storage.*')
                are not allowed. For more information see `IAM
                Overview <https://cloud.google.com/iam/docs/overview#permissions>`__.

                This corresponds to the ``permissions`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.iam.v1.iam_policy_pb2.TestIamPermissionsResponse:
                Response message for TestIamPermissions method.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [resource, permissions]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - The request isn't a proto-plus wrapped type,
        #   so it must be constructed via keyword expansion.
        if isinstance(request, dict):
            request = iam_policy_pb2.TestIamPermissionsRequest(**request)
        elif not request:
            request = iam_policy_pb2.TestIamPermissionsRequest(
                resource=resource, permissions=permissions
            )

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.test_iam_permissions
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        routing_param_regex = re.compile(
            "^(?P<bucket>projects/[^/]+/buckets/[^/]+)/objects(?:/.*)?$"
        )
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        routing_param_regex = re.compile(
            "^(?P<bucket>projects/[^/]+/buckets/[^/]+)/managedFolders(?:/.*)?$"
        )
        regex_match = routing_param_regex.match(request.resource)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def update_bucket(
        self,
        request: Optional[Union[storage.UpdateBucketRequest, dict]] = None,
        *,
        bucket: Optional[storage.Bucket] = None,
        update_mask: Optional[field_mask_pb2.FieldMask] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Bucket:
        r"""Updates a bucket. Equivalent to JSON API's
        storage.buckets.patch method.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_update_bucket():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.UpdateBucketRequest(
                )

                # Make the request
                response = await client.update_bucket(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.UpdateBucketRequest, dict]]):
                The request object. Request for UpdateBucket method.
            bucket (:class:`google.cloud._storage_v2.types.Bucket`):
                Required. The bucket to update. The bucket's ``name``
                field will be used to identify the bucket.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            update_mask (:class:`google.protobuf.field_mask_pb2.FieldMask`):
                Required. List of fields to be updated.

                To specify ALL fields, equivalent to the JSON API's
                "update" function, specify a single field with the value
                ``*``. Note: not recommended. If a new field is
                introduced at a later time, an older client updating
                with the ``*`` may accidentally reset the new field's
                value.

                Not specifying any fields is an error.

                This corresponds to the ``update_mask`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Bucket:
                A bucket.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket, update_mask]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.UpdateBucketRequest):
            request = storage.UpdateBucketRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket
        if update_mask is not None:
            request.update_mask = update_mask

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.update_bucket
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket.name)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def compose_object(
        self,
        request: Optional[Union[storage.ComposeObjectRequest, dict]] = None,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Object:
        r"""Concatenates a list of existing objects into a new
        object in the same bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_compose_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.ComposeObjectRequest(
                )

                # Make the request
                response = await client.compose_object(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.ComposeObjectRequest, dict]]):
                The request object. Request message for ComposeObject.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Object:
                An object.
        """
        # Create or coerce a protobuf request object.
        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.ComposeObjectRequest):
            request = storage.ComposeObjectRequest(request)

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.compose_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.destination.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def delete_object(
        self,
        request: Optional[Union[storage.DeleteObjectRequest, dict]] = None,
        *,
        bucket: Optional[str] = None,
        object_: Optional[str] = None,
        generation: Optional[int] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> None:
        r"""Deletes an object and its metadata. Deletions are permanent if
        versioning is not enabled for the bucket, or if the generation
        parameter is used, or if `soft
        delete <https://cloud.google.com/storage/docs/soft-delete>`__ is
        not enabled for the bucket. When this API is used to delete an
        object from a bucket that has soft delete policy enabled, the
        object becomes soft deleted, and the ``softDeleteTime`` and
        ``hardDeleteTime`` properties are set on the object. This API
        cannot be used to permanently delete soft-deleted objects.
        Soft-deleted objects are permanently deleted according to their
        ``hardDeleteTime``.

        You can use the
        [``RestoreObject``][google.storage.v2.Storage.RestoreObject] API
        to restore soft-deleted objects until the soft delete retention
        period has passed.

        **IAM Permissions**:

        Requires ``storage.objects.delete`` `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        on the bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_delete_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.DeleteObjectRequest(
                    bucket="bucket_value",
                    object_="object__value",
                )

                # Make the request
                await client.delete_object(request=request)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.DeleteObjectRequest, dict]]):
                The request object. Message for deleting an object. ``bucket`` and
                ``object`` **must** be set.
            bucket (:class:`str`):
                Required. Name of the bucket in which
                the object resides.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            object_ (:class:`str`):
                Required. The name of the finalized object to delete.
                Note: If you want to delete an unfinalized resumable
                upload please use ``CancelResumableWrite``.

                This corresponds to the ``object_`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            generation (:class:`int`):
                Optional. If present, permanently
                deletes a specific revision of this
                object (as opposed to the latest
                version, the default).

                This corresponds to the ``generation`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket, object_, generation]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.DeleteObjectRequest):
            request = storage.DeleteObjectRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket
        if object_ is not None:
            request.object_ = object_
        if generation is not None:
            request.generation = generation

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.delete_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

    async def restore_object(
        self,
        request: Optional[Union[storage.RestoreObjectRequest, dict]] = None,
        *,
        bucket: Optional[str] = None,
        object_: Optional[str] = None,
        generation: Optional[int] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Object:
        r"""Restores a soft-deleted object.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_restore_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.RestoreObjectRequest(
                    bucket="bucket_value",
                    object_="object__value",
                    generation=1068,
                )

                # Make the request
                response = await client.restore_object(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.RestoreObjectRequest, dict]]):
                The request object. Message for restoring an object. ``bucket``, ``object``,
                and ``generation`` **must** be set.
            bucket (:class:`str`):
                Required. Name of the bucket in which
                the object resides.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            object_ (:class:`str`):
                Required. The name of the object to
                restore.

                This corresponds to the ``object_`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            generation (:class:`int`):
                Required. The specific revision of
                the object to restore.

                This corresponds to the ``generation`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Object:
                An object.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket, object_, generation]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.RestoreObjectRequest):
            request = storage.RestoreObjectRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket
        if object_ is not None:
            request.object_ = object_
        if generation is not None:
            request.generation = generation

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.restore_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def cancel_resumable_write(
        self,
        request: Optional[Union[storage.CancelResumableWriteRequest, dict]] = None,
        *,
        upload_id: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.CancelResumableWriteResponse:
        r"""Cancels an in-progress resumable upload.

        Any attempts to write to the resumable upload after
        cancelling the upload will fail.

        The behavior for currently in progress write operations
        is not guaranteed - they could either complete before
        the cancellation or fail if the cancellation completes
        first.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_cancel_resumable_write():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.CancelResumableWriteRequest(
                    upload_id="upload_id_value",
                )

                # Make the request
                response = await client.cancel_resumable_write(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.CancelResumableWriteRequest, dict]]):
                The request object. Message for canceling an in-progress resumable upload.
                ``upload_id`` **must** be set.
            upload_id (:class:`str`):
                Required. The upload_id of the resumable upload to
                cancel. This should be copied from the ``upload_id``
                field of ``StartResumableWriteResponse``.

                This corresponds to the ``upload_id`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.CancelResumableWriteResponse:
                Empty response message for canceling
                an in-progress resumable upload, will be
                extended as needed.

        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [upload_id]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.CancelResumableWriteRequest):
            request = storage.CancelResumableWriteRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if upload_id is not None:
            request.upload_id = upload_id

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.cancel_resumable_write
        ]

        header_params = {}

        routing_param_regex = re.compile(
            "^(?P<bucket>projects/[^/]+/buckets/[^/]+)(?:/.*)?$"
        )
        regex_match = routing_param_regex.match(request.upload_id)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def get_object(
        self,
        request: Optional[Union[storage.GetObjectRequest, dict]] = None,
        *,
        bucket: Optional[str] = None,
        object_: Optional[str] = None,
        generation: Optional[int] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Object:
        r"""Retrieves object metadata.

        **IAM Permissions**:

        Requires ``storage.objects.get`` `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        on the bucket. To return object ACLs, the authenticated user
        must also have the ``storage.objects.getIamPolicy`` permission.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_get_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.GetObjectRequest(
                    bucket="bucket_value",
                    object_="object__value",
                )

                # Make the request
                response = await client.get_object(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.GetObjectRequest, dict]]):
                The request object. Request message for GetObject.
            bucket (:class:`str`):
                Required. Name of the bucket in which
                the object resides.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            object_ (:class:`str`):
                Required. Name of the object.
                This corresponds to the ``object_`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            generation (:class:`int`):
                Optional. If present, selects a
                specific revision of this object (as
                opposed to the latest version, the
                default).

                This corresponds to the ``generation`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Object:
                An object.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket, object_, generation]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.GetObjectRequest):
            request = storage.GetObjectRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket
        if object_ is not None:
            request.object_ = object_
        if generation is not None:
            request.generation = generation

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.get_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    def read_object(
        self,
        request: Optional[Union[storage.ReadObjectRequest, dict]] = None,
        *,
        bucket: Optional[str] = None,
        object_: Optional[str] = None,
        generation: Optional[int] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> Awaitable[AsyncIterable[storage.ReadObjectResponse]]:
        r"""Retrieves object data.

        **IAM Permissions**:

        Requires ``storage.objects.get`` `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        on the bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_read_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.ReadObjectRequest(
                    bucket="bucket_value",
                    object_="object__value",
                )

                # Make the request
                stream = await client.read_object(request=request)

                # Handle the response
                async for response in stream:
                    print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.ReadObjectRequest, dict]]):
                The request object. Request message for ReadObject.
            bucket (:class:`str`):
                Required. The name of the bucket
                containing the object to read.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            object_ (:class:`str`):
                Required. The name of the object to
                read.

                This corresponds to the ``object_`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            generation (:class:`int`):
                Optional. If present, selects a
                specific revision of this object (as
                opposed to the latest version, the
                default).

                This corresponds to the ``generation`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            AsyncIterable[google.cloud._storage_v2.types.ReadObjectResponse]:
                Response message for ReadObject.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket, object_, generation]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.ReadObjectRequest):
            request = storage.ReadObjectRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket
        if object_ is not None:
            request.object_ = object_
        if generation is not None:
            request.generation = generation

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.read_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    def bidi_read_object(
        self,
        requests: Optional[AsyncIterator[storage.BidiReadObjectRequest]] = None,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> Awaitable[AsyncIterable[storage.BidiReadObjectResponse]]:
        r"""Reads an object's data.

        This is a bi-directional API with the added support for reading
        multiple ranges within one stream both within and across
        multiple messages. If the server encountered an error for any of
        the inputs, the stream will be closed with the relevant error
        code. Because the API allows for multiple outstanding requests,
        when the stream is closed the error response will contain a
        BidiReadObjectRangesError proto in the error extension
        describing the error for each outstanding read_id.

        **IAM Permissions**:

        Requires ``storage.objects.get``

        `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        on the bucket.

        This API is currently in preview and is not yet available for
        general use.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_bidi_read_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.BidiReadObjectRequest(
                )

                # This method expects an iterator which contains
                # 'storage_v2.BidiReadObjectRequest' objects
                # Here we create a generator that yields a single `request` for
                # demonstrative purposes.
                requests = [request]

                def request_generator():
                    for request in requests:
                        yield request

                # Make the request
                stream = await client.bidi_read_object(requests=request_generator())

                # Handle the response
                async for response in stream:
                    print(response)

        Args:
            requests (AsyncIterator[`google.cloud._storage_v2.types.BidiReadObjectRequest`]):
                The request object AsyncIterator. Request message for BidiReadObject.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            AsyncIterable[google.cloud._storage_v2.types.BidiReadObjectResponse]:
                Response message for BidiReadObject.
        """

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.bidi_read_object
        ]

        header_params = {}

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = rpc(
            requests,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def update_object(
        self,
        request: Optional[Union[storage.UpdateObjectRequest, dict]] = None,
        *,
        object_: Optional[storage.Object] = None,
        update_mask: Optional[field_mask_pb2.FieldMask] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Object:
        r"""Updates an object's metadata.
        Equivalent to JSON API's storage.objects.patch.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_update_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.UpdateObjectRequest(
                )

                # Make the request
                response = await client.update_object(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.UpdateObjectRequest, dict]]):
                The request object. Request message for UpdateObject.
            object_ (:class:`google.cloud._storage_v2.types.Object`):
                Required. The object to update.
                The object's bucket and name fields are
                used to identify the object to update.
                If present, the object's generation
                field selects a specific revision of
                this object whose metadata should be
                updated. Otherwise, assumes the live
                version of the object.

                This corresponds to the ``object_`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            update_mask (:class:`google.protobuf.field_mask_pb2.FieldMask`):
                Required. List of fields to be updated.

                To specify ALL fields, equivalent to the JSON API's
                "update" function, specify a single field with the value
                ``*``. Note: not recommended. If a new field is
                introduced at a later time, an older client updating
                with the ``*`` may accidentally reset the new field's
                value.

                Not specifying any fields is an error.

                This corresponds to the ``update_mask`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Object:
                An object.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [object_, update_mask]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.UpdateObjectRequest):
            request = storage.UpdateObjectRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if object_ is not None:
            request.object_ = object_
        if update_mask is not None:
            request.update_mask = update_mask

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.update_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.object.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def write_object(
        self,
        requests: Optional[AsyncIterator[storage.WriteObjectRequest]] = None,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.WriteObjectResponse:
        r"""Stores a new object and metadata.

        An object can be written either in a single message stream or in
        a resumable sequence of message streams. To write using a single
        stream, the client should include in the first message of the
        stream an ``WriteObjectSpec`` describing the destination bucket,
        object, and any preconditions. Additionally, the final message
        must set 'finish_write' to true, or else it is an error.

        For a resumable write, the client should instead call
        ``StartResumableWrite()``, populating a ``WriteObjectSpec`` into
        that request. They should then attach the returned ``upload_id``
        to the first message of each following call to ``WriteObject``.
        If the stream is closed before finishing the upload (either
        explicitly by the client or due to a network error or an error
        response from the server), the client should do as follows:

        - Check the result Status of the stream, to determine if writing
          can be resumed on this stream or must be restarted from
          scratch (by calling ``StartResumableWrite()``). The resumable
          errors are DEADLINE_EXCEEDED, INTERNAL, and UNAVAILABLE. For
          each case, the client should use binary exponential backoff
          before retrying. Additionally, writes can be resumed after
          RESOURCE_EXHAUSTED errors, but only after taking appropriate
          measures, which may include reducing aggregate send rate
          across clients and/or requesting a quota increase for your
          project.
        - If the call to ``WriteObject`` returns ``ABORTED``, that
          indicates concurrent attempts to update the resumable write,
          caused either by multiple racing clients or by a single client
          where the previous request was timed out on the client side
          but nonetheless reached the server. In this case the client
          should take steps to prevent further concurrent writes (e.g.,
          increase the timeouts, stop using more than one process to
          perform the upload, etc.), and then should follow the steps
          below for resuming the upload.
        - For resumable errors, the client should call
          ``QueryWriteStatus()`` and then continue writing from the
          returned ``persisted_size``. This may be less than the amount
          of data the client previously sent. Note also that it is
          acceptable to send data starting at an offset earlier than the
          returned ``persisted_size``; in this case, the service will
          skip data at offsets that were already persisted (without
          checking that it matches the previously written data), and
          write only the data starting from the persisted offset. Even
          though the data isn't written, it may still incur a
          performance cost over resuming at the correct write offset.
          This behavior can make client-side handling simpler in some
          cases.
        - Clients must only send data that is a multiple of 256 KiB per
          message, unless the object is being finished with
          ``finish_write`` set to ``true``.

        The service will not view the object as complete until the
        client has sent a ``WriteObjectRequest`` with ``finish_write``
        set to ``true``. Sending any requests on a stream after sending
        a request with ``finish_write`` set to ``true`` will cause an
        error. The client **should** check the response it receives to
        determine how much data the service was able to commit and
        whether the service views the object as complete.

        Attempting to resume an already finalized object will result in
        an OK status, with a ``WriteObjectResponse`` containing the
        finalized object's metadata.

        Alternatively, the BidiWriteObject operation may be used to
        write an object with controls over flushing and the ability to
        fetch the ability to determine the current persisted size.

        **IAM Permissions**:

        Requires ``storage.objects.create`` `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        on the bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_write_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.WriteObjectRequest(
                    upload_id="upload_id_value",
                    write_offset=1297,
                )

                # This method expects an iterator which contains
                # 'storage_v2.WriteObjectRequest' objects
                # Here we create a generator that yields a single `request` for
                # demonstrative purposes.
                requests = [request]

                def request_generator():
                    for request in requests:
                        yield request

                # Make the request
                response = await client.write_object(requests=request_generator())

                # Handle the response
                print(response)

        Args:
            requests (AsyncIterator[`google.cloud._storage_v2.types.WriteObjectRequest`]):
                The request object AsyncIterator. Request message for WriteObject.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.WriteObjectResponse:
                Response message for WriteObject.
        """

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.write_object
        ]

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            requests,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    def bidi_write_object(
        self,
        requests: Optional[AsyncIterator[storage.BidiWriteObjectRequest]] = None,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> Awaitable[AsyncIterable[storage.BidiWriteObjectResponse]]:
        r"""Stores a new object and metadata.

        This is similar to the WriteObject call with the added support
        for manual flushing of persisted state, and the ability to
        determine current persisted size without closing the stream.

        The client may specify one or both of the ``state_lookup`` and
        ``flush`` fields in each BidiWriteObjectRequest. If ``flush`` is
        specified, the data written so far will be persisted to storage.
        If ``state_lookup`` is specified, the service will respond with
        a BidiWriteObjectResponse that contains the persisted size. If
        both ``flush`` and ``state_lookup`` are specified, the flush
        will always occur before a ``state_lookup``, so that both may be
        set in the same request and the returned state will be the state
        of the object post-flush. When the stream is closed, a
        BidiWriteObjectResponse will always be sent to the client,
        regardless of the value of ``state_lookup``.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_bidi_write_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.BidiWriteObjectRequest(
                    upload_id="upload_id_value",
                    write_offset=1297,
                )

                # This method expects an iterator which contains
                # 'storage_v2.BidiWriteObjectRequest' objects
                # Here we create a generator that yields a single `request` for
                # demonstrative purposes.
                requests = [request]

                def request_generator():
                    for request in requests:
                        yield request

                # Make the request
                stream = await client.bidi_write_object(requests=request_generator())

                # Handle the response
                async for response in stream:
                    print(response)

        Args:
            requests (AsyncIterator[`google.cloud._storage_v2.types.BidiWriteObjectRequest`]):
                The request object AsyncIterator. Request message for BidiWriteObject.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            AsyncIterable[google.cloud._storage_v2.types.BidiWriteObjectResponse]:
                Response message for BidiWriteObject.
        """

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.bidi_write_object
        ]

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = rpc(
            requests,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def list_objects(
        self,
        request: Optional[Union[storage.ListObjectsRequest, dict]] = None,
        *,
        parent: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> pagers.ListObjectsAsyncPager:
        r"""Retrieves a list of objects matching the criteria.

        **IAM Permissions**:

        The authenticated user requires ``storage.objects.list`` `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        to use this method. To return object ACLs, the authenticated
        user must also have the ``storage.objects.getIamPolicy``
        permission.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_list_objects():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.ListObjectsRequest(
                    parent="parent_value",
                )

                # Make the request
                page_result = client.list_objects(request=request)

                # Handle the response
                async for response in page_result:
                    print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.ListObjectsRequest, dict]]):
                The request object. Request message for ListObjects.
            parent (:class:`str`):
                Required. Name of the bucket in which
                to look for objects.

                This corresponds to the ``parent`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.services.storage.pagers.ListObjectsAsyncPager:
                The result of a call to
                Objects.ListObjects
                Iterating over this object will yield
                results and resolve additional pages
                automatically.

        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [parent]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.ListObjectsRequest):
            request = storage.ListObjectsRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if parent is not None:
            request.parent = parent

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.list_objects
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.parent)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # This method is paged; wrap the response in a pager, which provides
        # an `__aiter__` convenience method.
        response = pagers.ListObjectsAsyncPager(
            method=rpc,
            request=request,
            response=response,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def rewrite_object(
        self,
        request: Optional[Union[storage.RewriteObjectRequest, dict]] = None,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.RewriteResponse:
        r"""Rewrites a source object to a destination object.
        Optionally overrides metadata.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_rewrite_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.RewriteObjectRequest(
                    destination_name="destination_name_value",
                    destination_bucket="destination_bucket_value",
                    source_bucket="source_bucket_value",
                    source_object="source_object_value",
                )

                # Make the request
                response = await client.rewrite_object(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.RewriteObjectRequest, dict]]):
                The request object. Request message for RewriteObject. If the source object
                is encrypted using a Customer-Supplied Encryption Key
                the key information must be provided in the
                copy_source_encryption_algorithm,
                copy_source_encryption_key_bytes, and
                copy_source_encryption_key_sha256_bytes fields. If the
                destination object should be encrypted the keying
                information should be provided in the
                encryption_algorithm, encryption_key_bytes, and
                encryption_key_sha256_bytes fields of the
                common_object_request_params.customer_encryption field.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.RewriteResponse:
                A rewrite response.
        """
        # Create or coerce a protobuf request object.
        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.RewriteObjectRequest):
            request = storage.RewriteObjectRequest(request)

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.rewrite_object
        ]

        header_params = {}

        if request.source_bucket:
            header_params["source_bucket"] = request.source_bucket

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.destination_bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def start_resumable_write(
        self,
        request: Optional[Union[storage.StartResumableWriteRequest, dict]] = None,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.StartResumableWriteResponse:
        r"""Starts a resumable write operation. This method is part of the
        `Resumable
        upload <https://cloud.google.com/storage/docs/resumable-uploads>`__
        feature. This allows you to upload large objects in multiple
        chunks, which is more resilient to network interruptions than a
        single upload. The validity duration of the write operation, and
        the consequences of it becoming invalid, are service-dependent.

        **IAM Permissions**:

        Requires ``storage.objects.create`` `IAM
        permission <https://cloud.google.com/iam/docs/overview#permissions>`__
        on the bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_start_resumable_write():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.StartResumableWriteRequest(
                )

                # Make the request
                response = await client.start_resumable_write(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.StartResumableWriteRequest, dict]]):
                The request object. Request message StartResumableWrite.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.StartResumableWriteResponse:
                Response object for StartResumableWrite.
        """
        # Create or coerce a protobuf request object.
        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.StartResumableWriteRequest):
            request = storage.StartResumableWriteRequest(request)

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.start_resumable_write
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(
            request.write_object_spec.resource.bucket
        )
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def query_write_status(
        self,
        request: Optional[Union[storage.QueryWriteStatusRequest, dict]] = None,
        *,
        upload_id: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.QueryWriteStatusResponse:
        r"""Determines the ``persisted_size`` of an object that is being
        written. This method is part of the `resumable
        upload <https://cloud.google.com/storage/docs/resumable-uploads>`__
        feature. The returned value is the size of the object that has
        been persisted so far. The value can be used as the
        ``write_offset`` for the next ``Write()`` call.

        If the object does not exist, meaning if it was deleted, or the
        first ``Write()`` has not yet reached the service, this method
        returns the error ``NOT_FOUND``.

        This method is useful for clients that buffer data and need to
        know which data can be safely evicted. The client can call
        ``QueryWriteStatus()`` at any time to determine how much data
        has been logged for this object. For any sequence of
        ``QueryWriteStatus()`` calls for a given object name, the
        sequence of returned ``persisted_size`` values are
        non-decreasing.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_query_write_status():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.QueryWriteStatusRequest(
                    upload_id="upload_id_value",
                )

                # Make the request
                response = await client.query_write_status(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.QueryWriteStatusRequest, dict]]):
                The request object. Request object for ``QueryWriteStatus``.
            upload_id (:class:`str`):
                Required. The name of the resume
                token for the object whose write status
                is being requested.

                This corresponds to the ``upload_id`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.QueryWriteStatusResponse:
                Response object for QueryWriteStatus.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [upload_id]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.QueryWriteStatusRequest):
            request = storage.QueryWriteStatusRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if upload_id is not None:
            request.upload_id = upload_id

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.query_write_status
        ]

        header_params = {}

        routing_param_regex = re.compile(
            "^(?P<bucket>projects/[^/]+/buckets/[^/]+)(?:/.*)?$"
        )
        regex_match = routing_param_regex.match(request.upload_id)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def move_object(
        self,
        request: Optional[Union[storage.MoveObjectRequest, dict]] = None,
        *,
        bucket: Optional[str] = None,
        source_object: Optional[str] = None,
        destination_object: Optional[str] = None,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
    ) -> storage.Object:
        r"""Moves the source object to the destination object in
        the same bucket.

        .. code-block:: python

            # This snippet has been automatically generated and should be regarded as a
            # code template only.
            # It will require modifications to work:
            # - It may require correct/in-range values for request initialization.
            # - It may require specifying regional endpoints when creating the service
            #   client as shown in:
            #   https://googleapis.dev/python/google-api-core/latest/client_options.html
            from google.cloud import storage_v2

            async def sample_move_object():
                # Create a client
                client = storage_v2.StorageAsyncClient()

                # Initialize request argument(s)
                request = storage_v2.MoveObjectRequest(
                    bucket="bucket_value",
                    source_object="source_object_value",
                    destination_object="destination_object_value",
                )

                # Make the request
                response = await client.move_object(request=request)

                # Handle the response
                print(response)

        Args:
            request (Optional[Union[google.cloud._storage_v2.types.MoveObjectRequest, dict]]):
                The request object. Request message for MoveObject.
            bucket (:class:`str`):
                Required. Name of the bucket in which
                the object resides.

                This corresponds to the ``bucket`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            source_object (:class:`str`):
                Required. Name of the source object.
                This corresponds to the ``source_object`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            destination_object (:class:`str`):
                Required. Name of the destination
                object.

                This corresponds to the ``destination_object`` field
                on the ``request`` instance; if ``request`` is provided, this
                should not be set.
            retry (google.api_core.retry_async.AsyncRetry): Designation of what errors, if any,
                should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.

        Returns:
            google.cloud._storage_v2.types.Object:
                An object.
        """
        # Create or coerce a protobuf request object.
        # - Quick check: If we got a request object, we should *not* have
        #   gotten any keyword arguments that map to the request.
        flattened_params = [bucket, source_object, destination_object]
        has_flattened_params = (
            len([param for param in flattened_params if param is not None]) > 0
        )
        if request is not None and has_flattened_params:
            raise ValueError(
                "If the `request` argument is set, then none of "
                "the individual field arguments should be set."
            )

        # - Use the request object if provided (there's no risk of modifying the input as
        #   there are no flattened fields), or create one.
        if not isinstance(request, storage.MoveObjectRequest):
            request = storage.MoveObjectRequest(request)

        # If we have keyword arguments corresponding to fields on the
        # request, apply these.
        if bucket is not None:
            request.bucket = bucket
        if source_object is not None:
            request.source_object = source_object
        if destination_object is not None:
            request.destination_object = destination_object

        # Wrap the RPC method; this adds retry and timeout information,
        # and friendly error handling.
        rpc = self._client._transport._wrapped_methods[
            self._client._transport.move_object
        ]

        header_params = {}

        routing_param_regex = re.compile("^(?P<bucket>.*)$")
        regex_match = routing_param_regex.match(request.bucket)
        if regex_match and regex_match.group("bucket"):
            header_params["bucket"] = regex_match.group("bucket")

        if header_params:
            metadata = tuple(metadata) + (
                gapic_v1.routing_header.to_grpc_metadata(header_params),
            )

        # Validate the universe domain.
        self._client._validate_universe_domain()

        # Send the request.
        response = await rpc(
            request,
            retry=retry,
            timeout=timeout,
            metadata=metadata,
        )

        # Done; return the response.
        return response

    async def __aenter__(self) -> "StorageAsyncClient":
        return self

    async def __aexit__(self, exc_type, exc, tb):
        await self.transport.close()


DEFAULT_CLIENT_INFO = gapic_v1.client_info.ClientInfo(
    gapic_version=package_version.__version__
)

if hasattr(DEFAULT_CLIENT_INFO, "protobuf_runtime_version"):  # pragma: NO COVER
    DEFAULT_CLIENT_INFO.protobuf_runtime_version = google.protobuf.__version__


__all__ = ("StorageAsyncClient",)
