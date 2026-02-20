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

from typing import List, Optional, Tuple
import grpc
from google.cloud import _storage_v2
from google.cloud.storage.asyncio import _utils
from google.cloud.storage.asyncio.async_grpc_client import AsyncGrpcClient
from google.cloud.storage.asyncio.async_abstract_object_stream import (
    _AsyncAbstractObjectStream,
)
from google.api_core.bidi_async import AsyncBidiRpc


class _AsyncWriteObjectStream(_AsyncAbstractObjectStream):
    """Class representing a gRPC bidi-stream for writing data from a GCS
      ``Appendable Object``.

    This class provides a unix socket-like interface to a GCS ``Object``, with
    methods like ``open``, ``close``, ``send``, and ``recv``.

    :type client: :class:`~google.cloud.storage.asyncio.async_grpc_client.AsyncGrpcClient.grpc_client`
    :param client: async grpc client to use for making API requests.

    :type bucket_name: str
    :param bucket_name: The name of the GCS ``bucket`` containing the object.

    :type object_name: str
    :param object_name: The name of the GCS ``Appendable Object`` to be write.

    :type generation_number: int
    :param generation_number: (Optional) If present, creates writer for that
        specific revision of that object. Use this to append data to an
        existing Appendable Object.

        Setting to ``0`` makes the `writer.open()` succeed only if
        object doesn't exist in the bucket (useful for not accidentally
        overwriting existing objects).

        Warning: If `None`, a new object is created. If an object with the
        same name already exists, it will be overwritten the moment
        `writer.open()` is called.

    :type write_handle: _storage_v2.BidiWriteHandle
    :param write_handle: (Optional) An existing handle for writing the object.
                        If provided, opening the bidi-gRPC connection will be faster.
    """

    def __init__(
        self,
        client: AsyncGrpcClient.grpc_client,
        bucket_name: str,
        object_name: str,
        generation_number: Optional[int] = None,  # None means new object
        write_handle: Optional[_storage_v2.BidiWriteHandle] = None,
        routing_token: Optional[str] = None,
    ) -> None:
        if client is None:
            raise ValueError("client must be provided")
        if bucket_name is None:
            raise ValueError("bucket_name must be provided")
        if object_name is None:
            raise ValueError("object_name must be provided")

        super().__init__(
            bucket_name=bucket_name,
            object_name=object_name,
            generation_number=generation_number,
        )
        self.client: AsyncGrpcClient.grpc_client = client
        self.write_handle: Optional[_storage_v2.BidiWriteHandle] = write_handle
        self.routing_token: Optional[str] = routing_token

        self._full_bucket_name = f"projects/_/buckets/{self.bucket_name}"

        self.rpc = self.client._client._transport._wrapped_methods[
            self.client._client._transport.bidi_write_object
        ]

        self.metadata = (("x-goog-request-params", f"bucket={self._full_bucket_name}"),)
        self.socket_like_rpc: Optional[AsyncBidiRpc] = None
        self._is_stream_open: bool = False
        self.first_bidi_write_req = None
        self.persisted_size = 0
        self.object_resource: Optional[_storage_v2.Object] = None

    async def open(self, metadata: Optional[List[Tuple[str, str]]] = None) -> None:
        """
        Opens the bidi-gRPC connection to write to the object.

        This method sends an initial request to start the stream and receives
        the first response containing metadata and a write handle.

        :rtype: None
        :raises ValueError: If the stream is already open.
        :raises google.api_core.exceptions.FailedPrecondition:
            if `generation_number` is 0 and object already exists.
        """
        if self._is_stream_open:
            raise ValueError("Stream is already open")

        # Create a new object or overwrite existing one if generation_number
        # is None. This makes it consistent with GCS JSON API behavior.
        # Created object type would be Appendable Object.
        # if `generation_number` == 0 new object will be created only if there
        # isn't any existing object.
        if self.generation_number is None or self.generation_number == 0:
            self.first_bidi_write_req = _storage_v2.BidiWriteObjectRequest(
                write_object_spec=_storage_v2.WriteObjectSpec(
                    resource=_storage_v2.Object(
                        name=self.object_name, bucket=self._full_bucket_name
                    ),
                    appendable=True,
                    if_generation_match=self.generation_number,
                ),
            )
        else:
            self.first_bidi_write_req = _storage_v2.BidiWriteObjectRequest(
                append_object_spec=_storage_v2.AppendObjectSpec(
                    bucket=self._full_bucket_name,
                    object=self.object_name,
                    generation=self.generation_number,
                    write_handle=self.write_handle if self.write_handle else None,
                    routing_token=self.routing_token if self.routing_token else None,
                ),
            )

        request_param_values = [f"bucket={self._full_bucket_name}"]
        final_metadata = []
        if metadata:
            for key, value in metadata:
                if key == "x-goog-request-params":
                    request_param_values.append(value)
                else:
                    final_metadata.append((key, value))

        final_metadata.append(("x-goog-request-params", ",".join(request_param_values)))

        self.socket_like_rpc = AsyncBidiRpc(
            self.rpc,
            initial_request=self.first_bidi_write_req,
            metadata=final_metadata,
        )

        await self.socket_like_rpc.open()  # this is actually 1 send
        response = await self.socket_like_rpc.recv()
        self._is_stream_open = True

        if response.persisted_size:
            self.persisted_size = response.persisted_size

        if response.resource:
            if not response.resource.size:
                # Appending to a 0 byte appendable object.
                self.persisted_size = 0
            else:
                self.persisted_size = response.resource.size

            self.generation_number = response.resource.generation

        if response.write_handle:
            self.write_handle = response.write_handle

    async def close(self) -> None:
        """Closes the bidi-gRPC connection."""
        if not self._is_stream_open:
            raise ValueError("Stream is not open")
        await self.requests_done()
        await self.socket_like_rpc.close()
        self._is_stream_open = False

    async def requests_done(self):
        """Signals that all requests have been sent."""
        await self.socket_like_rpc.send(None)

        # The server may send a final "EOF" response immediately, or it may
        # first send an intermediate response followed by the EOF response depending on whether the object was finalized or not.
        first_resp = await self.socket_like_rpc.recv()
        _utils.update_write_handle_if_exists(self, first_resp)

        if first_resp != grpc.aio.EOF:
            self.persisted_size = first_resp.persisted_size
            second_resp = await self.socket_like_rpc.recv()
            assert second_resp == grpc.aio.EOF

    async def send(
        self, bidi_write_object_request: _storage_v2.BidiWriteObjectRequest
    ) -> None:
        """Sends a request message on the stream.

        Args:
            bidi_write_object_request (:class:`~google.cloud._storage_v2.types.BidiReadObjectRequest`):
                The request message to send. This is typically used to specify
                the read offset and limit.
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open")
        await self.socket_like_rpc.send(bidi_write_object_request)

    async def recv(self) -> _storage_v2.BidiWriteObjectResponse:
        """Receives a response from the stream.

        This method waits for the next message from the server, which could
        contain object data or metadata.

        Returns:
            :class:`~google.cloud._storage_v2.types.BidiWriteObjectResponse`:
                The response message from the server.
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open")
        response = await self.socket_like_rpc.recv()
        # Update write_handle if present in response
        if response:
            if response.write_handle:
                self.write_handle = response.write_handle
            if response.persisted_size is not None:
                self.persisted_size = response.persisted_size
            if response.resource and response.resource.size:
                self.persisted_size = response.resource.size
        return response

    @property
    def is_stream_open(self) -> bool:
        return self._is_stream_open
