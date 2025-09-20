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

"""An async client for interacting with Google Cloud Storage using the gRPC API."""

from google.cloud import _storage_v2 as storage_v2


class AsyncGrpcClient:
    """An asynchronous client for interacting with Google Cloud Storage using the gRPC API.

    :type credentials: :class:`~google.auth.credentials.Credentials`
    :param credentials: (Optional) The OAuth2 Credentials to use for this
                        client. If not passed, falls back to the default
                        inferred from the environment.

    :type client_info: :class:`~google.api_core.client_info.ClientInfo`
    :param client_info:
        The client info used to send a user-agent string along with API
        requests. If ``None``, then default info will be used.

    :type client_options: :class:`~google.api_core.client_options.ClientOptions` or :class:`dict`
    :param client_options: (Optional) Client options used to set user options
        on the client.

    :type attempt_direct_path: bool
    :param attempt_direct_path:
        (Optional) Whether to attempt to use DirectPath for gRPC connections.
        Defaults to ``True``.
    """

    def __init__(
        self,
        credentials=None,
        client_info=None,
        client_options=None,
        *,
        attempt_direct_path=True,
    ):
        self._grpc_client = self._create_async_grpc_client(
            credentials=credentials,
            client_info=client_info,
            client_options=client_options,
            attempt_direct_path=attempt_direct_path,
        )

    def _create_async_grpc_client(
        self,
        credentials=None,
        client_info=None,
        client_options=None,
        attempt_direct_path=True,
    ):
        transport_cls = storage_v2.StorageAsyncClient.get_transport_class(
            "grpc_asyncio"
        )
        channel = transport_cls.create_channel(attempt_direct_path=attempt_direct_path)
        transport = transport_cls(credentials=credentials, channel=channel)

        return storage_v2.StorageAsyncClient(
            credentials=credentials,
            transport=transport,
            client_info=client_info,
            client_options=client_options,
        )

    @property
    def grpc_client(self):
        """The underlying gRPC client.

        This property gives users direct access to the `_storage_v2.StorageAsyncClient`
         instance. This can be useful for accessing
        newly added or experimental RPCs that are not yet exposed through
        the high-level GrpcClient.
        Returns:
            google.cloud._storage_v2.StorageAsyncClient: The configured GAPIC client.
        """
        return self._grpc_client
