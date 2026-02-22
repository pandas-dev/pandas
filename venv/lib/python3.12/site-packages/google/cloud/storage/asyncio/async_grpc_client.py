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
from google.cloud._storage_v2.services.storage.transports.base import (
    DEFAULT_CLIENT_INFO,
)
from google.cloud.storage import __version__


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

    :type client_options: :class:`~google.api_core.client_options.ClientOptions`
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
        if client_info is None:
            client_info = DEFAULT_CLIENT_INFO
        client_info.client_library_version = __version__
        if client_info.user_agent is None:
            client_info.user_agent = ""
        agent_version = f"gcloud-python/{__version__}"
        if agent_version not in client_info.user_agent:
            client_info.user_agent += f" {agent_version} "

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

        primary_user_agent = client_info.to_user_agent()

        channel = transport_cls.create_channel(
            attempt_direct_path=attempt_direct_path,
            credentials=credentials,
            options=(("grpc.primary_user_agent", primary_user_agent),),
        )
        transport = transport_cls(channel=channel)

        return storage_v2.StorageAsyncClient(
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

    async def delete_object(
        self,
        bucket_name,
        object_name,
        generation=None,
        if_generation_match=None,
        if_generation_not_match=None,
        if_metageneration_match=None,
        if_metageneration_not_match=None,
        **kwargs,
    ):
        """Deletes an object and its metadata.

        :type bucket_name: str
        :param bucket_name: The name of the bucket in which the object resides.

        :type object_name: str
        :param object_name: The name of the object to delete.

        :type generation: int
        :param generation:
            (Optional) If present, permanently deletes a specific generation
            of an object.

        :type if_generation_match: int
        :param if_generation_match: (Optional)

        :type if_generation_not_match: int
        :param if_generation_not_match: (Optional)

        :type if_metageneration_match: int
        :param if_metageneration_match: (Optional)

        :type if_metageneration_not_match: int
        :param if_metageneration_not_match: (Optional)


        """
        # The gRPC API requires the bucket name to be in the format "projects/_/buckets/bucket_name"
        bucket_path = f"projects/_/buckets/{bucket_name}"
        request = storage_v2.DeleteObjectRequest(
            bucket=bucket_path,
            object=object_name,
            generation=generation,
            if_generation_match=if_generation_match,
            if_generation_not_match=if_generation_not_match,
            if_metageneration_match=if_metageneration_match,
            if_metageneration_not_match=if_metageneration_not_match,
            **kwargs,
        )
        await self._grpc_client.delete_object(request=request)

    async def get_object(
        self,
        bucket_name,
        object_name,
        generation=None,
        if_generation_match=None,
        if_generation_not_match=None,
        if_metageneration_match=None,
        if_metageneration_not_match=None,
        soft_deleted=None,
        **kwargs,
    ):
        """Retrieves an object's metadata.

        In the gRPC API, this is performed by the GetObject RPC, which
        returns the object resource (metadata) without the object's data.

        :type bucket_name: str
        :param bucket_name: The name of the bucket in which the object resides.

        :type object_name: str
        :param object_name: The name of the object.

        :type generation: int
        :param generation:
            (Optional) If present, selects a specific generation of an object.

        :type if_generation_match: int
        :param if_generation_match: (Optional) Precondition for object generation match.

        :type if_generation_not_match: int
        :param if_generation_not_match: (Optional) Precondition for object generation mismatch.

        :type if_metageneration_match: int
        :param if_metageneration_match: (Optional) Precondition for metageneration match.

        :type if_metageneration_not_match: int
        :param if_metageneration_not_match: (Optional) Precondition for metageneration mismatch.

        :type soft_deleted: bool
        :param soft_deleted:
            (Optional) If True, return the soft-deleted version of this object.

        :rtype: :class:`google.cloud._storage_v2.types.Object`
        :returns: The object metadata resource.
        """
        bucket_path = f"projects/_/buckets/{bucket_name}"

        request = storage_v2.GetObjectRequest(
            bucket=bucket_path,
            object=object_name,
            generation=generation,
            if_generation_match=if_generation_match,
            if_generation_not_match=if_generation_not_match,
            if_metageneration_match=if_metageneration_match,
            if_metageneration_not_match=if_metageneration_not_match,
            soft_deleted=soft_deleted or False,
            **kwargs,
        )

        # Calls the underlying GAPIC StorageAsyncClient.get_object method
        return await self._grpc_client.get_object(request=request)
