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

"""A client for interacting with Google Cloud Storage using the gRPC API."""

from google.cloud.client import ClientWithProject
from google.cloud import _storage_v2 as storage_v2

_marker = object()


class GrpcClient(ClientWithProject):
    """A client for interacting with Google Cloud Storage using the gRPC API.

    :type project: str or None
    :param project: The project which the client acts on behalf of. If not
                    passed, falls back to the default inferred from the
                    environment.

    :type credentials: :class:`~google.auth.credentials.Credentials`
    :param credentials: (Optional) The OAuth2 Credentials to use for this
                        client. If not passed, falls back to the default
                        inferred from the environment.

    :type client_info: :class:`~google.api_core.client_info.ClientInfo`
    :param client_info:
        The client info used to send a user-agent string along with API
        requests. If ``None``, then default info will be used. Generally,
        you only need to set this if you're developing your own library
        or partner tool.

    :type client_options: :class:`~google.api_core.client_options.ClientOptions` or :class:`dict`
    :param client_options: (Optional) Client options used to set user options
        on the client. A non-default universe domain or API endpoint should be
        set through client_options.

    :type api_key: string
    :param api_key:
        (Optional) An API key. Mutually exclusive with any other credentials.
        This parameter is an alias for setting `client_options.api_key` and
        will supersede any API key set in the `client_options` parameter.

    :type attempt_direct_path: bool
    :param attempt_direct_path:
        (Optional) Whether to attempt to use DirectPath for gRPC connections.
        This provides a direct, unproxied connection to GCS for lower latency
        and higher throughput, and is highly recommended when running on Google
        Cloud infrastructure. Defaults to ``True``.
    """

    def __init__(
        self,
        project=_marker,
        credentials=None,
        client_info=None,
        client_options=None,
        *,
        api_key=None,
        attempt_direct_path=True,
    ):
        super(GrpcClient, self).__init__(project=project, credentials=credentials)

        if isinstance(client_options, dict):
            if api_key:
                client_options["api_key"] = api_key
        elif client_options is None:
            client_options = {} if not api_key else {"api_key": api_key}
        elif api_key:
            client_options.api_key = api_key

        self._grpc_client = self._create_gapic_client(
            credentials=credentials,
            client_info=client_info,
            client_options=client_options,
            attempt_direct_path=attempt_direct_path,
        )

    def _create_gapic_client(
        self,
        credentials=None,
        client_info=None,
        client_options=None,
        attempt_direct_path=True,
    ):
        """Creates and configures the low-level GAPIC `storage_v2` client."""
        transport_cls = storage_v2.StorageClient.get_transport_class("grpc")

        channel = transport_cls.create_channel(attempt_direct_path=attempt_direct_path)

        transport = transport_cls(credentials=credentials, channel=channel)

        return storage_v2.StorageClient(
            credentials=credentials,
            transport=transport,
            client_info=client_info,
            client_options=client_options,
        )

    @property
    def grpc_client(self):
        """The underlying gRPC client.

        This property gives users direct access to the `storage_v2.StorageClient`
         instance. This can be useful for accessing
        newly added or experimental RPCs that are not yet exposed through
        the high-level GrpcClient.

        Returns:
            google.cloud.storage_v2.StorageClient: The configured GAPIC client.
        """
        return self._grpc_client
