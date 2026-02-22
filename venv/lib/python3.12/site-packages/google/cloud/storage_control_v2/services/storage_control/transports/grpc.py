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
import json
import logging as std_logging
import pickle
from typing import Callable, Dict, Optional, Sequence, Tuple, Union
import warnings

from google.api_core import gapic_v1, grpc_helpers, operations_v1
import google.auth  # type: ignore
from google.auth import credentials as ga_credentials  # type: ignore
from google.auth.transport.grpc import SslCredentials  # type: ignore
import google.iam.v1.iam_policy_pb2 as iam_policy_pb2  # type: ignore
import google.iam.v1.policy_pb2 as policy_pb2  # type: ignore
from google.longrunning import operations_pb2  # type: ignore
import google.protobuf.empty_pb2 as empty_pb2  # type: ignore
from google.protobuf.json_format import MessageToJson
import google.protobuf.message
import grpc  # type: ignore
import proto  # type: ignore

from google.cloud.storage_control_v2.types import storage_control

from .base import DEFAULT_CLIENT_INFO, StorageControlTransport

try:
    from google.api_core import client_logging  # type: ignore

    CLIENT_LOGGING_SUPPORTED = True  # pragma: NO COVER
except ImportError:  # pragma: NO COVER
    CLIENT_LOGGING_SUPPORTED = False

_LOGGER = std_logging.getLogger(__name__)


class _LoggingClientInterceptor(grpc.UnaryUnaryClientInterceptor):  # pragma: NO COVER
    def intercept_unary_unary(self, continuation, client_call_details, request):
        logging_enabled = CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
            std_logging.DEBUG
        )
        if logging_enabled:  # pragma: NO COVER
            request_metadata = client_call_details.metadata
            if isinstance(request, proto.Message):
                request_payload = type(request).to_json(request)
            elif isinstance(request, google.protobuf.message.Message):
                request_payload = MessageToJson(request)
            else:
                request_payload = f"{type(request).__name__}: {pickle.dumps(request)}"

            request_metadata = {
                key: value.decode("utf-8") if isinstance(value, bytes) else value
                for key, value in request_metadata
            }
            grpc_request = {
                "payload": request_payload,
                "requestMethod": "grpc",
                "metadata": dict(request_metadata),
            }
            _LOGGER.debug(
                f"Sending request for {client_call_details.method}",
                extra={
                    "serviceName": "google.storage.control.v2.StorageControl",
                    "rpcName": str(client_call_details.method),
                    "request": grpc_request,
                    "metadata": grpc_request["metadata"],
                },
            )
        response = continuation(client_call_details, request)
        if logging_enabled:  # pragma: NO COVER
            response_metadata = response.trailing_metadata()
            # Convert gRPC metadata `<class 'grpc.aio._metadata.Metadata'>` to list of tuples
            metadata = (
                dict([(k, str(v)) for k, v in response_metadata])
                if response_metadata
                else None
            )
            result = response.result()
            if isinstance(result, proto.Message):
                response_payload = type(result).to_json(result)
            elif isinstance(result, google.protobuf.message.Message):
                response_payload = MessageToJson(result)
            else:
                response_payload = f"{type(result).__name__}: {pickle.dumps(result)}"
            grpc_response = {
                "payload": response_payload,
                "metadata": metadata,
                "status": "OK",
            }
            _LOGGER.debug(
                f"Received response for {client_call_details.method}.",
                extra={
                    "serviceName": "google.storage.control.v2.StorageControl",
                    "rpcName": client_call_details.method,
                    "response": grpc_response,
                    "metadata": grpc_response["metadata"],
                },
            )
        return response


class StorageControlGrpcTransport(StorageControlTransport):
    """gRPC backend transport for StorageControl.

    StorageControl service includes selected control plane
    operations.

    This class defines the same methods as the primary client, so the
    primary client can load the underlying transport implementation
    and call it.

    It sends protocol buffers over the wire using gRPC (which is built on
    top of HTTP/2); the ``grpcio`` package must be installed.
    """

    _stubs: Dict[str, Callable]

    def __init__(
        self,
        *,
        host: str = "storage.googleapis.com",
        credentials: Optional[ga_credentials.Credentials] = None,
        credentials_file: Optional[str] = None,
        scopes: Optional[Sequence[str]] = None,
        channel: Optional[Union[grpc.Channel, Callable[..., grpc.Channel]]] = None,
        api_mtls_endpoint: Optional[str] = None,
        client_cert_source: Optional[Callable[[], Tuple[bytes, bytes]]] = None,
        ssl_channel_credentials: Optional[grpc.ChannelCredentials] = None,
        client_cert_source_for_mtls: Optional[Callable[[], Tuple[bytes, bytes]]] = None,
        quota_project_id: Optional[str] = None,
        client_info: gapic_v1.client_info.ClientInfo = DEFAULT_CLIENT_INFO,
        always_use_jwt_access: Optional[bool] = False,
        api_audience: Optional[str] = None,
    ) -> None:
        """Instantiate the transport.

        Args:
            host (Optional[str]):
                 The hostname to connect to (default: 'storage.googleapis.com').
            credentials (Optional[google.auth.credentials.Credentials]): The
                authorization credentials to attach to requests. These
                credentials identify the application to the service; if none
                are specified, the client will attempt to ascertain the
                credentials from the environment.
                This argument is ignored if a ``channel`` instance is provided.
            credentials_file (Optional[str]): Deprecated. A file with credentials that can
                be loaded with :func:`google.auth.load_credentials_from_file`.
                This argument is ignored if a ``channel`` instance is provided.
                This argument will be removed in the next major version of this library.
            scopes (Optional(Sequence[str])): A list of scopes. This argument is
                ignored if a ``channel`` instance is provided.
            channel (Optional[Union[grpc.Channel, Callable[..., grpc.Channel]]]):
                A ``Channel`` instance through which to make calls, or a Callable
                that constructs and returns one. If set to None, ``self.create_channel``
                is used to create the channel. If a Callable is given, it will be called
                with the same arguments as used in ``self.create_channel``.
            api_mtls_endpoint (Optional[str]): Deprecated. The mutual TLS endpoint.
                If provided, it overrides the ``host`` argument and tries to create
                a mutual TLS channel with client SSL credentials from
                ``client_cert_source`` or application default SSL credentials.
            client_cert_source (Optional[Callable[[], Tuple[bytes, bytes]]]):
                Deprecated. A callback to provide client SSL certificate bytes and
                private key bytes, both in PEM format. It is ignored if
                ``api_mtls_endpoint`` is None.
            ssl_channel_credentials (grpc.ChannelCredentials): SSL credentials
                for the grpc channel. It is ignored if a ``channel`` instance is provided.
            client_cert_source_for_mtls (Optional[Callable[[], Tuple[bytes, bytes]]]):
                A callback to provide client certificate bytes and private key bytes,
                both in PEM format. It is used to configure a mutual TLS channel. It is
                ignored if a ``channel`` instance or ``ssl_channel_credentials`` is provided.
            quota_project_id (Optional[str]): An optional project to use for billing
                and quota.
            client_info (google.api_core.gapic_v1.client_info.ClientInfo):
                The client info used to send a user-agent string along with
                API requests. If ``None``, then default info will be used.
                Generally, you only need to set this if you're developing
                your own client library.
            always_use_jwt_access (Optional[bool]): Whether self signed JWT should
                be used for service account credentials.

        Raises:
          google.auth.exceptions.MutualTLSChannelError: If mutual TLS transport
              creation failed for any reason.
          google.api_core.exceptions.DuplicateCredentialArgs: If both ``credentials``
              and ``credentials_file`` are passed.
        """
        self._grpc_channel = None
        self._ssl_channel_credentials = ssl_channel_credentials
        self._stubs: Dict[str, Callable] = {}
        self._operations_client: Optional[operations_v1.OperationsClient] = None

        if api_mtls_endpoint:
            warnings.warn("api_mtls_endpoint is deprecated", DeprecationWarning)
        if client_cert_source:
            warnings.warn("client_cert_source is deprecated", DeprecationWarning)

        if isinstance(channel, grpc.Channel):
            # Ignore credentials if a channel was passed.
            credentials = None
            self._ignore_credentials = True
            # If a channel was explicitly provided, set it.
            self._grpc_channel = channel
            self._ssl_channel_credentials = None

        else:
            if api_mtls_endpoint:
                host = api_mtls_endpoint

                # Create SSL credentials with client_cert_source or application
                # default SSL credentials.
                if client_cert_source:
                    cert, key = client_cert_source()
                    self._ssl_channel_credentials = grpc.ssl_channel_credentials(
                        certificate_chain=cert, private_key=key
                    )
                else:
                    self._ssl_channel_credentials = SslCredentials().ssl_credentials

            else:
                if client_cert_source_for_mtls and not ssl_channel_credentials:
                    cert, key = client_cert_source_for_mtls()
                    self._ssl_channel_credentials = grpc.ssl_channel_credentials(
                        certificate_chain=cert, private_key=key
                    )

        # The base transport sets the host, credentials and scopes
        super().__init__(
            host=host,
            credentials=credentials,
            credentials_file=credentials_file,
            scopes=scopes,
            quota_project_id=quota_project_id,
            client_info=client_info,
            always_use_jwt_access=always_use_jwt_access,
            api_audience=api_audience,
        )

        if not self._grpc_channel:
            # initialize with the provided callable or the default channel
            channel_init = channel or type(self).create_channel
            self._grpc_channel = channel_init(
                self._host,
                # use the credentials which are saved
                credentials=self._credentials,
                # Set ``credentials_file`` to ``None`` here as
                # the credentials that we saved earlier should be used.
                credentials_file=None,
                scopes=self._scopes,
                ssl_credentials=self._ssl_channel_credentials,
                quota_project_id=quota_project_id,
                options=[
                    ("grpc.max_send_message_length", -1),
                    ("grpc.max_receive_message_length", -1),
                ],
            )

        self._interceptor = _LoggingClientInterceptor()
        self._logged_channel = grpc.intercept_channel(
            self._grpc_channel, self._interceptor
        )

        # Wrap messages. This must be done after self._logged_channel exists
        self._prep_wrapped_messages(client_info)

    @classmethod
    def create_channel(
        cls,
        host: str = "storage.googleapis.com",
        credentials: Optional[ga_credentials.Credentials] = None,
        credentials_file: Optional[str] = None,
        scopes: Optional[Sequence[str]] = None,
        quota_project_id: Optional[str] = None,
        **kwargs,
    ) -> grpc.Channel:
        """Create and return a gRPC channel object.
        Args:
            host (Optional[str]): The host for the channel to use.
            credentials (Optional[~.Credentials]): The
                authorization credentials to attach to requests. These
                credentials identify this application to the service. If
                none are specified, the client will attempt to ascertain
                the credentials from the environment.
            credentials_file (Optional[str]): Deprecated. A file with credentials that can
                be loaded with :func:`google.auth.load_credentials_from_file`.
                This argument is mutually exclusive with credentials.  This argument will be
                removed in the next major version of this library.
            scopes (Optional[Sequence[str]]): A optional list of scopes needed for this
                service. These are only used when credentials are not specified and
                are passed to :func:`google.auth.default`.
            quota_project_id (Optional[str]): An optional project to use for billing
                and quota.
            kwargs (Optional[dict]): Keyword arguments, which are passed to the
                channel creation.
        Returns:
            grpc.Channel: A gRPC channel object.

        Raises:
            google.api_core.exceptions.DuplicateCredentialArgs: If both ``credentials``
              and ``credentials_file`` are passed.
        """

        return grpc_helpers.create_channel(
            host,
            credentials=credentials,
            credentials_file=credentials_file,
            quota_project_id=quota_project_id,
            default_scopes=cls.AUTH_SCOPES,
            scopes=scopes,
            default_host=cls.DEFAULT_HOST,
            **kwargs,
        )

    @property
    def grpc_channel(self) -> grpc.Channel:
        """Return the channel designed to connect to this service."""
        return self._grpc_channel

    @property
    def operations_client(self) -> operations_v1.OperationsClient:
        """Create the client designed to process long-running operations.

        This property caches on the instance; repeated calls return the same
        client.
        """
        # Quick check: Only create a new client if we do not already have one.
        if self._operations_client is None:
            self._operations_client = operations_v1.OperationsClient(
                self._logged_channel
            )

        # Return the client from cache.
        return self._operations_client

    @property
    def create_folder(
        self,
    ) -> Callable[[storage_control.CreateFolderRequest], storage_control.Folder]:
        r"""Return a callable for the create folder method over gRPC.

        Creates a new folder. This operation is only
        applicable to a hierarchical namespace enabled bucket.

        Returns:
            Callable[[~.CreateFolderRequest],
                    ~.Folder]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "create_folder" not in self._stubs:
            self._stubs["create_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/CreateFolder",
                request_serializer=storage_control.CreateFolderRequest.serialize,
                response_deserializer=storage_control.Folder.deserialize,
            )
        return self._stubs["create_folder"]

    @property
    def delete_folder(
        self,
    ) -> Callable[[storage_control.DeleteFolderRequest], empty_pb2.Empty]:
        r"""Return a callable for the delete folder method over gRPC.

        Permanently deletes an empty folder. This operation
        is only applicable to a hierarchical namespace enabled
        bucket.

        Returns:
            Callable[[~.DeleteFolderRequest],
                    ~.Empty]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "delete_folder" not in self._stubs:
            self._stubs["delete_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/DeleteFolder",
                request_serializer=storage_control.DeleteFolderRequest.serialize,
                response_deserializer=empty_pb2.Empty.FromString,
            )
        return self._stubs["delete_folder"]

    @property
    def get_folder(
        self,
    ) -> Callable[[storage_control.GetFolderRequest], storage_control.Folder]:
        r"""Return a callable for the get folder method over gRPC.

        Returns metadata for the specified folder. This
        operation is only applicable to a hierarchical namespace
        enabled bucket.

        Returns:
            Callable[[~.GetFolderRequest],
                    ~.Folder]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_folder" not in self._stubs:
            self._stubs["get_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetFolder",
                request_serializer=storage_control.GetFolderRequest.serialize,
                response_deserializer=storage_control.Folder.deserialize,
            )
        return self._stubs["get_folder"]

    @property
    def list_folders(
        self,
    ) -> Callable[
        [storage_control.ListFoldersRequest], storage_control.ListFoldersResponse
    ]:
        r"""Return a callable for the list folders method over gRPC.

        Retrieves a list of folders. This operation is only
        applicable to a hierarchical namespace enabled bucket.

        Returns:
            Callable[[~.ListFoldersRequest],
                    ~.ListFoldersResponse]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "list_folders" not in self._stubs:
            self._stubs["list_folders"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/ListFolders",
                request_serializer=storage_control.ListFoldersRequest.serialize,
                response_deserializer=storage_control.ListFoldersResponse.deserialize,
            )
        return self._stubs["list_folders"]

    @property
    def rename_folder(
        self,
    ) -> Callable[[storage_control.RenameFolderRequest], operations_pb2.Operation]:
        r"""Return a callable for the rename folder method over gRPC.

        Renames a source folder to a destination folder. This
        operation is only applicable to a hierarchical namespace
        enabled bucket. During a rename, the source and
        destination folders are locked until the long running
        operation completes.

        Returns:
            Callable[[~.RenameFolderRequest],
                    ~.Operation]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "rename_folder" not in self._stubs:
            self._stubs["rename_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/RenameFolder",
                request_serializer=storage_control.RenameFolderRequest.serialize,
                response_deserializer=operations_pb2.Operation.FromString,
            )
        return self._stubs["rename_folder"]

    @property
    def delete_folder_recursive(
        self,
    ) -> Callable[
        [storage_control.DeleteFolderRecursiveRequest], operations_pb2.Operation
    ]:
        r"""Return a callable for the delete folder recursive method over gRPC.

        Deletes a folder recursively. This operation is only
        applicable to a hierarchical namespace enabled bucket.

        Returns:
            Callable[[~.DeleteFolderRecursiveRequest],
                    ~.Operation]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "delete_folder_recursive" not in self._stubs:
            self._stubs["delete_folder_recursive"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/DeleteFolderRecursive",
                request_serializer=storage_control.DeleteFolderRecursiveRequest.serialize,
                response_deserializer=operations_pb2.Operation.FromString,
            )
        return self._stubs["delete_folder_recursive"]

    @property
    def get_storage_layout(
        self,
    ) -> Callable[
        [storage_control.GetStorageLayoutRequest], storage_control.StorageLayout
    ]:
        r"""Return a callable for the get storage layout method over gRPC.

        Returns the storage layout configuration for a given
        bucket.

        Returns:
            Callable[[~.GetStorageLayoutRequest],
                    ~.StorageLayout]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_storage_layout" not in self._stubs:
            self._stubs["get_storage_layout"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetStorageLayout",
                request_serializer=storage_control.GetStorageLayoutRequest.serialize,
                response_deserializer=storage_control.StorageLayout.deserialize,
            )
        return self._stubs["get_storage_layout"]

    @property
    def create_managed_folder(
        self,
    ) -> Callable[
        [storage_control.CreateManagedFolderRequest], storage_control.ManagedFolder
    ]:
        r"""Return a callable for the create managed folder method over gRPC.

        Creates a new managed folder.

        Returns:
            Callable[[~.CreateManagedFolderRequest],
                    ~.ManagedFolder]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "create_managed_folder" not in self._stubs:
            self._stubs["create_managed_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/CreateManagedFolder",
                request_serializer=storage_control.CreateManagedFolderRequest.serialize,
                response_deserializer=storage_control.ManagedFolder.deserialize,
            )
        return self._stubs["create_managed_folder"]

    @property
    def delete_managed_folder(
        self,
    ) -> Callable[[storage_control.DeleteManagedFolderRequest], empty_pb2.Empty]:
        r"""Return a callable for the delete managed folder method over gRPC.

        Permanently deletes an empty managed folder.

        Returns:
            Callable[[~.DeleteManagedFolderRequest],
                    ~.Empty]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "delete_managed_folder" not in self._stubs:
            self._stubs["delete_managed_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/DeleteManagedFolder",
                request_serializer=storage_control.DeleteManagedFolderRequest.serialize,
                response_deserializer=empty_pb2.Empty.FromString,
            )
        return self._stubs["delete_managed_folder"]

    @property
    def get_managed_folder(
        self,
    ) -> Callable[
        [storage_control.GetManagedFolderRequest], storage_control.ManagedFolder
    ]:
        r"""Return a callable for the get managed folder method over gRPC.

        Returns metadata for the specified managed folder.

        Returns:
            Callable[[~.GetManagedFolderRequest],
                    ~.ManagedFolder]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_managed_folder" not in self._stubs:
            self._stubs["get_managed_folder"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetManagedFolder",
                request_serializer=storage_control.GetManagedFolderRequest.serialize,
                response_deserializer=storage_control.ManagedFolder.deserialize,
            )
        return self._stubs["get_managed_folder"]

    @property
    def list_managed_folders(
        self,
    ) -> Callable[
        [storage_control.ListManagedFoldersRequest],
        storage_control.ListManagedFoldersResponse,
    ]:
        r"""Return a callable for the list managed folders method over gRPC.

        Retrieves a list of managed folders for a given
        bucket.

        Returns:
            Callable[[~.ListManagedFoldersRequest],
                    ~.ListManagedFoldersResponse]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "list_managed_folders" not in self._stubs:
            self._stubs["list_managed_folders"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/ListManagedFolders",
                request_serializer=storage_control.ListManagedFoldersRequest.serialize,
                response_deserializer=storage_control.ListManagedFoldersResponse.deserialize,
            )
        return self._stubs["list_managed_folders"]

    @property
    def create_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.CreateAnywhereCacheRequest], operations_pb2.Operation
    ]:
        r"""Return a callable for the create anywhere cache method over gRPC.

        Creates an Anywhere Cache instance.

        Returns:
            Callable[[~.CreateAnywhereCacheRequest],
                    ~.Operation]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "create_anywhere_cache" not in self._stubs:
            self._stubs["create_anywhere_cache"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/CreateAnywhereCache",
                request_serializer=storage_control.CreateAnywhereCacheRequest.serialize,
                response_deserializer=operations_pb2.Operation.FromString,
            )
        return self._stubs["create_anywhere_cache"]

    @property
    def update_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.UpdateAnywhereCacheRequest], operations_pb2.Operation
    ]:
        r"""Return a callable for the update anywhere cache method over gRPC.

        Updates an Anywhere Cache instance. Mutable fields include
        ``ttl`` and ``admission_policy``.

        Returns:
            Callable[[~.UpdateAnywhereCacheRequest],
                    ~.Operation]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "update_anywhere_cache" not in self._stubs:
            self._stubs["update_anywhere_cache"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/UpdateAnywhereCache",
                request_serializer=storage_control.UpdateAnywhereCacheRequest.serialize,
                response_deserializer=operations_pb2.Operation.FromString,
            )
        return self._stubs["update_anywhere_cache"]

    @property
    def disable_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.DisableAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        r"""Return a callable for the disable anywhere cache method over gRPC.

        Disables an Anywhere Cache instance. A disabled
        instance is read-only. The disablement could be revoked
        by calling ResumeAnywhereCache. The cache instance will
        be deleted automatically if it remains in the disabled
        state for at least one hour.

        Returns:
            Callable[[~.DisableAnywhereCacheRequest],
                    ~.AnywhereCache]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "disable_anywhere_cache" not in self._stubs:
            self._stubs["disable_anywhere_cache"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/DisableAnywhereCache",
                request_serializer=storage_control.DisableAnywhereCacheRequest.serialize,
                response_deserializer=storage_control.AnywhereCache.deserialize,
            )
        return self._stubs["disable_anywhere_cache"]

    @property
    def pause_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.PauseAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        r"""Return a callable for the pause anywhere cache method over gRPC.

        Pauses an Anywhere Cache instance.

        Returns:
            Callable[[~.PauseAnywhereCacheRequest],
                    ~.AnywhereCache]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "pause_anywhere_cache" not in self._stubs:
            self._stubs["pause_anywhere_cache"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/PauseAnywhereCache",
                request_serializer=storage_control.PauseAnywhereCacheRequest.serialize,
                response_deserializer=storage_control.AnywhereCache.deserialize,
            )
        return self._stubs["pause_anywhere_cache"]

    @property
    def resume_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.ResumeAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        r"""Return a callable for the resume anywhere cache method over gRPC.

        Resumes a disabled or paused Anywhere Cache instance.

        Returns:
            Callable[[~.ResumeAnywhereCacheRequest],
                    ~.AnywhereCache]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "resume_anywhere_cache" not in self._stubs:
            self._stubs["resume_anywhere_cache"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/ResumeAnywhereCache",
                request_serializer=storage_control.ResumeAnywhereCacheRequest.serialize,
                response_deserializer=storage_control.AnywhereCache.deserialize,
            )
        return self._stubs["resume_anywhere_cache"]

    @property
    def get_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.GetAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        r"""Return a callable for the get anywhere cache method over gRPC.

        Gets an Anywhere Cache instance.

        Returns:
            Callable[[~.GetAnywhereCacheRequest],
                    ~.AnywhereCache]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_anywhere_cache" not in self._stubs:
            self._stubs["get_anywhere_cache"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetAnywhereCache",
                request_serializer=storage_control.GetAnywhereCacheRequest.serialize,
                response_deserializer=storage_control.AnywhereCache.deserialize,
            )
        return self._stubs["get_anywhere_cache"]

    @property
    def list_anywhere_caches(
        self,
    ) -> Callable[
        [storage_control.ListAnywhereCachesRequest],
        storage_control.ListAnywhereCachesResponse,
    ]:
        r"""Return a callable for the list anywhere caches method over gRPC.

        Lists Anywhere Cache instances for a given bucket.

        Returns:
            Callable[[~.ListAnywhereCachesRequest],
                    ~.ListAnywhereCachesResponse]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "list_anywhere_caches" not in self._stubs:
            self._stubs["list_anywhere_caches"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/ListAnywhereCaches",
                request_serializer=storage_control.ListAnywhereCachesRequest.serialize,
                response_deserializer=storage_control.ListAnywhereCachesResponse.deserialize,
            )
        return self._stubs["list_anywhere_caches"]

    @property
    def get_project_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetProjectIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        r"""Return a callable for the get project intelligence
        config method over gRPC.

        Returns the Project scoped singleton
        IntelligenceConfig resource.

        Returns:
            Callable[[~.GetProjectIntelligenceConfigRequest],
                    ~.IntelligenceConfig]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_project_intelligence_config" not in self._stubs:
            self._stubs[
                "get_project_intelligence_config"
            ] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetProjectIntelligenceConfig",
                request_serializer=storage_control.GetProjectIntelligenceConfigRequest.serialize,
                response_deserializer=storage_control.IntelligenceConfig.deserialize,
            )
        return self._stubs["get_project_intelligence_config"]

    @property
    def update_project_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateProjectIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        r"""Return a callable for the update project intelligence
        config method over gRPC.

        Updates the Project scoped singleton
        IntelligenceConfig resource.

        Returns:
            Callable[[~.UpdateProjectIntelligenceConfigRequest],
                    ~.IntelligenceConfig]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "update_project_intelligence_config" not in self._stubs:
            self._stubs[
                "update_project_intelligence_config"
            ] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/UpdateProjectIntelligenceConfig",
                request_serializer=storage_control.UpdateProjectIntelligenceConfigRequest.serialize,
                response_deserializer=storage_control.IntelligenceConfig.deserialize,
            )
        return self._stubs["update_project_intelligence_config"]

    @property
    def get_folder_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetFolderIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        r"""Return a callable for the get folder intelligence config method over gRPC.

        Returns the Folder scoped singleton
        IntelligenceConfig resource.

        Returns:
            Callable[[~.GetFolderIntelligenceConfigRequest],
                    ~.IntelligenceConfig]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_folder_intelligence_config" not in self._stubs:
            self._stubs[
                "get_folder_intelligence_config"
            ] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetFolderIntelligenceConfig",
                request_serializer=storage_control.GetFolderIntelligenceConfigRequest.serialize,
                response_deserializer=storage_control.IntelligenceConfig.deserialize,
            )
        return self._stubs["get_folder_intelligence_config"]

    @property
    def update_folder_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateFolderIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        r"""Return a callable for the update folder intelligence
        config method over gRPC.

        Updates the Folder scoped singleton
        IntelligenceConfig resource.

        Returns:
            Callable[[~.UpdateFolderIntelligenceConfigRequest],
                    ~.IntelligenceConfig]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "update_folder_intelligence_config" not in self._stubs:
            self._stubs[
                "update_folder_intelligence_config"
            ] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/UpdateFolderIntelligenceConfig",
                request_serializer=storage_control.UpdateFolderIntelligenceConfigRequest.serialize,
                response_deserializer=storage_control.IntelligenceConfig.deserialize,
            )
        return self._stubs["update_folder_intelligence_config"]

    @property
    def get_organization_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetOrganizationIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        r"""Return a callable for the get organization intelligence
        config method over gRPC.

        Returns the Organization scoped singleton
        IntelligenceConfig resource.

        Returns:
            Callable[[~.GetOrganizationIntelligenceConfigRequest],
                    ~.IntelligenceConfig]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_organization_intelligence_config" not in self._stubs:
            self._stubs[
                "get_organization_intelligence_config"
            ] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetOrganizationIntelligenceConfig",
                request_serializer=storage_control.GetOrganizationIntelligenceConfigRequest.serialize,
                response_deserializer=storage_control.IntelligenceConfig.deserialize,
            )
        return self._stubs["get_organization_intelligence_config"]

    @property
    def update_organization_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateOrganizationIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        r"""Return a callable for the update organization
        intelligence config method over gRPC.

        Updates the Organization scoped singleton
        IntelligenceConfig resource.

        Returns:
            Callable[[~.UpdateOrganizationIntelligenceConfigRequest],
                    ~.IntelligenceConfig]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "update_organization_intelligence_config" not in self._stubs:
            self._stubs[
                "update_organization_intelligence_config"
            ] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/UpdateOrganizationIntelligenceConfig",
                request_serializer=storage_control.UpdateOrganizationIntelligenceConfigRequest.serialize,
                response_deserializer=storage_control.IntelligenceConfig.deserialize,
            )
        return self._stubs["update_organization_intelligence_config"]

    @property
    def get_iam_policy(
        self,
    ) -> Callable[[iam_policy_pb2.GetIamPolicyRequest], policy_pb2.Policy]:
        r"""Return a callable for the get iam policy method over gRPC.

        Gets the IAM policy for a specified bucket. The ``resource``
        field in the request should be ``projects/_/buckets/{bucket}``
        for a bucket, or
        ``projects/_/buckets/{bucket}/managedFolders/{managedFolder}``
        for a managed folder.

        Returns:
            Callable[[~.GetIamPolicyRequest],
                    ~.Policy]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "get_iam_policy" not in self._stubs:
            self._stubs["get_iam_policy"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/GetIamPolicy",
                request_serializer=iam_policy_pb2.GetIamPolicyRequest.SerializeToString,
                response_deserializer=policy_pb2.Policy.FromString,
            )
        return self._stubs["get_iam_policy"]

    @property
    def set_iam_policy(
        self,
    ) -> Callable[[iam_policy_pb2.SetIamPolicyRequest], policy_pb2.Policy]:
        r"""Return a callable for the set iam policy method over gRPC.

        Updates an IAM policy for the specified bucket. The ``resource``
        field in the request should be ``projects/_/buckets/{bucket}``
        for a bucket, or
        ``projects/_/buckets/{bucket}/managedFolders/{managedFolder}``
        for a managed folder.

        Returns:
            Callable[[~.SetIamPolicyRequest],
                    ~.Policy]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "set_iam_policy" not in self._stubs:
            self._stubs["set_iam_policy"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/SetIamPolicy",
                request_serializer=iam_policy_pb2.SetIamPolicyRequest.SerializeToString,
                response_deserializer=policy_pb2.Policy.FromString,
            )
        return self._stubs["set_iam_policy"]

    @property
    def test_iam_permissions(
        self,
    ) -> Callable[
        [iam_policy_pb2.TestIamPermissionsRequest],
        iam_policy_pb2.TestIamPermissionsResponse,
    ]:
        r"""Return a callable for the test iam permissions method over gRPC.

        Tests a set of permissions on the given bucket, object, or
        managed folder to see which, if any, are held by the caller. The
        ``resource`` field in the request should be
        ``projects/_/buckets/{bucket}`` for a bucket,
        ``projects/_/buckets/{bucket}/objects/{object}`` for an object,
        or
        ``projects/_/buckets/{bucket}/managedFolders/{managedFolder}``
        for a managed folder.

        Returns:
            Callable[[~.TestIamPermissionsRequest],
                    ~.TestIamPermissionsResponse]:
                A function that, when called, will call the underlying RPC
                on the server.
        """
        # Generate a "stub function" on-the-fly which will actually make
        # the request.
        # gRPC handles serialization and deserialization, so we just need
        # to pass in the functions for each.
        if "test_iam_permissions" not in self._stubs:
            self._stubs["test_iam_permissions"] = self._logged_channel.unary_unary(
                "/google.storage.control.v2.StorageControl/TestIamPermissions",
                request_serializer=iam_policy_pb2.TestIamPermissionsRequest.SerializeToString,
                response_deserializer=iam_policy_pb2.TestIamPermissionsResponse.FromString,
            )
        return self._stubs["test_iam_permissions"]

    def close(self):
        self._logged_channel.close()

    @property
    def kind(self) -> str:
        return "grpc"


__all__ = ("StorageControlGrpcTransport",)
