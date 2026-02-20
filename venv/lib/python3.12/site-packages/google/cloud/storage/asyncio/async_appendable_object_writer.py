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

from io import BufferedReader
import io
import logging
from typing import List, Optional, Tuple, Union

from google.api_core import exceptions
from google.api_core.retry_async import AsyncRetry
from google.rpc import status_pb2
from google.cloud._storage_v2.types import BidiWriteObjectRedirectedError
from google.cloud._storage_v2.types.storage import BidiWriteObjectRequest


from . import _utils
from google.cloud import _storage_v2
from google.cloud.storage.asyncio.async_grpc_client import (
    AsyncGrpcClient,
)
from google.cloud.storage.asyncio.async_write_object_stream import (
    _AsyncWriteObjectStream,
)
from google.cloud.storage.asyncio.retry.bidi_stream_retry_manager import (
    _BidiStreamRetryManager,
)
from google.cloud.storage.asyncio.retry.writes_resumption_strategy import (
    _WriteResumptionStrategy,
    _WriteState,
)
from google.cloud.storage.asyncio.retry._helpers import (
    _extract_bidi_writes_redirect_proto,
)


_MAX_CHUNK_SIZE_BYTES = 2 * 1024 * 1024  # 2 MiB
_DEFAULT_FLUSH_INTERVAL_BYTES = 16 * 1024 * 1024  # 16 MiB
_BIDI_WRITE_REDIRECTED_TYPE_URL = (
    "type.googleapis.com/google.storage.v2.BidiWriteObjectRedirectedError"
)
logger = logging.getLogger(__name__)


def _is_write_retryable(exc):
    """Predicate to determine if a write operation should be retried."""

    if isinstance(
        exc,
        (
            exceptions.InternalServerError,
            exceptions.ServiceUnavailable,
            exceptions.DeadlineExceeded,
            exceptions.TooManyRequests,
            BidiWriteObjectRedirectedError,
        ),
    ):
        logger.warning(f"Retryable write exception encountered: {exc}")
        return True

    grpc_error = None
    if isinstance(exc, exceptions.Aborted) and exc.errors:
        grpc_error = exc.errors[0]
        if isinstance(grpc_error, BidiWriteObjectRedirectedError):
            return True

        trailers = grpc_error.trailing_metadata()
        if not trailers:
            return False

        status_details_bin = None
        for key, value in trailers:
            if key == "grpc-status-details-bin":
                status_details_bin = value
                break

        if status_details_bin:
            status_proto = status_pb2.Status()
            try:
                status_proto.ParseFromString(status_details_bin)
                for detail in status_proto.details:
                    if detail.type_url == _BIDI_WRITE_REDIRECTED_TYPE_URL:
                        return True
            except Exception:
                logger.error(
                    "Error unpacking redirect details from gRPC error. Exception: ",
                    {exc},
                )
                return False
    return False


class AsyncAppendableObjectWriter:
    """Class for appending data to a GCS Appendable Object asynchronously."""

    def __init__(
        self,
        client: AsyncGrpcClient,
        bucket_name: str,
        object_name: str,
        generation: Optional[int] = None,
        write_handle: Optional[_storage_v2.BidiWriteHandle] = None,
        writer_options: Optional[dict] = None,
    ):
        """
        Class for appending data to a GCS Appendable Object.

        Example usage:

        ```

        from google.cloud.storage.asyncio.async_grpc_client import AsyncGrpcClient
        from google.cloud.storage.asyncio.async_appendable_object_writer import AsyncAppendableObjectWriter
        import asyncio

        client = AsyncGrpcClient()
        bucket_name = "my-bucket"
        object_name = "my-appendable-object"

        # instantiate the writer
        writer = AsyncAppendableObjectWriter(client, bucket_name, object_name)
        # open the writer, (underlying gRPC bidi-stream will be opened)
        await writer.open()

        # append data, it can be called multiple times.
        await writer.append(b"hello world")
        await writer.append(b"some more data")

        # optionally flush data to persist.
        await writer.flush()

        # close the gRPC stream.
        # Please note closing the program will also close the stream,
        # however it's recommended to close the stream if no more data to append
        # to clean up gRPC connection (which means CPU/memory/network resources)
        await writer.close()
        ```

        :type client: :class:`~google.cloud.storage.asyncio.async_grpc_client.AsyncGrpcClient`
        :param client: async grpc client to use for making API requests.

        :type bucket_name: str
        :param bucket_name: The name of the GCS bucket containing the object.

        :type object_name: str
        :param object_name: The name of the GCS Appendable Object to be written.

        :type generation: Optional[int]
        :param generation: (Optional) If present, creates writer for that
            specific revision of that object. Use this to append data to an
            existing Appendable Object.

            Setting to ``0`` makes the `writer.open()` succeed only if
            object doesn't exist in the bucket (useful for not accidentally
            overwriting existing objects).

            Warning: If `None`, a new object is created. If an object with the
            same name already exists, it will be overwritten the moment
            `writer.open()` is called.

        :type write_handle: _storage_v2.BidiWriteHandle
        :param write_handle: (Optional) An handle for writing the object.
            If provided, opening the bidi-gRPC connection will be faster.

        :type writer_options: dict
        :param writer_options: (Optional) A dictionary of writer options.
            Supported options:
            - "FLUSH_INTERVAL_BYTES": int
                The number of bytes to append before "persisting" data in GCS
                servers. Default is `_DEFAULT_FLUSH_INTERVAL_BYTES`.
                Must be a multiple of `_MAX_CHUNK_SIZE_BYTES`.
        """
        _utils.raise_if_no_fast_crc32c()
        self.client = client
        self.bucket_name = bucket_name
        self.object_name = object_name
        self.write_handle = write_handle
        self.generation = generation

        self.write_obj_stream: Optional[_AsyncWriteObjectStream] = None
        self._is_stream_open: bool = False
        # `offset` is the latest size of the object without staleless.
        self.offset: Optional[int] = None
        # `persisted_size` is the total_bytes persisted in the GCS server.
        # Please note: `offset` and `persisted_size` are same when the stream is
        # opened.
        self.persisted_size: Optional[int] = None
        if writer_options is None:
            writer_options = {}
        self.flush_interval = writer_options.get(
            "FLUSH_INTERVAL_BYTES", _DEFAULT_FLUSH_INTERVAL_BYTES
        )
        if self.flush_interval < _MAX_CHUNK_SIZE_BYTES:
            raise exceptions.OutOfRange(
                f"flush_interval must be >= {_MAX_CHUNK_SIZE_BYTES} , but provided {self.flush_interval}"
            )
        if self.flush_interval % _MAX_CHUNK_SIZE_BYTES != 0:
            raise exceptions.OutOfRange(
                f"flush_interval must be a multiple of {_MAX_CHUNK_SIZE_BYTES}, but provided {self.flush_interval}"
            )
        self.bytes_appended_since_last_flush = 0
        self._routing_token: Optional[str] = None
        self.object_resource: Optional[_storage_v2.Object] = None

    async def state_lookup(self) -> int:
        """Returns the persisted_size

        :rtype: int
        :returns: persisted size.

        :raises ValueError: If the stream is not open (i.e., `open()` has not
            been called).
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open. Call open() before state_lookup().")

        await self.write_obj_stream.send(
            _storage_v2.BidiWriteObjectRequest(
                state_lookup=True,
            )
        )
        response = await self.write_obj_stream.recv()
        self.persisted_size = response.persisted_size
        return self.persisted_size

    def _on_open_error(self, exc):
        """Extracts routing token and write handle on redirect error during open."""
        redirect_proto = _extract_bidi_writes_redirect_proto(exc)
        if redirect_proto:
            if redirect_proto.routing_token:
                self._routing_token = redirect_proto.routing_token
            if redirect_proto.write_handle:
                self.write_handle = redirect_proto.write_handle
            if redirect_proto.generation:
                self.generation = redirect_proto.generation

    async def open(
        self,
        retry_policy: Optional[AsyncRetry] = None,
        metadata: Optional[List[Tuple[str, str]]] = None,
    ) -> None:
        """Opens the underlying bidi-gRPC stream.

        :raises ValueError: If the stream is already open.

        """
        if self._is_stream_open:
            raise ValueError("Underlying bidi-gRPC stream is already open")

        if retry_policy is None:
            retry_policy = AsyncRetry(
                predicate=_is_write_retryable, on_error=self._on_open_error
            )
        else:
            original_on_error = retry_policy._on_error

            def combined_on_error(exc):
                self._on_open_error(exc)
                if original_on_error:
                    original_on_error(exc)

            retry_policy = AsyncRetry(
                predicate=_is_write_retryable,
                initial=retry_policy._initial,
                maximum=retry_policy._maximum,
                multiplier=retry_policy._multiplier,
                deadline=retry_policy._deadline,
                on_error=combined_on_error,
            )

        async def _do_open():
            current_metadata = list(metadata) if metadata else []

            # Cleanup stream from previous failed attempt, if any.
            if self.write_obj_stream:
                if self.write_obj_stream.is_stream_open:
                    try:
                        await self.write_obj_stream.close()
                    except Exception as e:
                        logger.warning(
                            "Error closing previous write stream during open retry. Got exception: ",
                            {e},
                        )
                self.write_obj_stream = None
                self._is_stream_open = False

            self.write_obj_stream = _AsyncWriteObjectStream(
                client=self.client.grpc_client,
                bucket_name=self.bucket_name,
                object_name=self.object_name,
                generation_number=self.generation,
                write_handle=self.write_handle,
                routing_token=self._routing_token,
            )

            if self._routing_token:
                current_metadata.append(
                    ("x-goog-request-params", f"routing_token={self._routing_token}")
                )

            await self.write_obj_stream.open(
                metadata=current_metadata if metadata else None
            )

            if self.write_obj_stream.generation_number:
                self.generation = self.write_obj_stream.generation_number
            if self.write_obj_stream.write_handle:
                self.write_handle = self.write_obj_stream.write_handle
            if self.write_obj_stream.persisted_size is not None:
                self.persisted_size = self.write_obj_stream.persisted_size

            self._is_stream_open = True
            self._routing_token = None

        await retry_policy(_do_open)()

    async def append(
        self,
        data: bytes,
        retry_policy: Optional[AsyncRetry] = None,
        metadata: Optional[List[Tuple[str, str]]] = None,
    ) -> None:
        """Appends data to the Appendable object with automatic retries.

        calling `self.append` will append bytes at the end of the current size
        ie. `self.offset` bytes relative to the begining of the object.

        This method sends the provided `data` to the GCS server in chunks.
        and persists data in GCS at every `_DEFAULT_FLUSH_INTERVAL_BYTES` bytes
        or at the last chunk whichever is earlier. Persisting is done by setting
        `flush=True` on request.

        :type data: bytes
        :param data: The bytes to append to the object.

        :type retry_policy: :class:`~google.api_core.retry_async.AsyncRetry`
        :param retry_policy: (Optional) The retry policy to use for the operation.

        :type metadata: List[Tuple[str, str]]
        :param metadata: (Optional) The metadata to be sent with the request.

        :raises ValueError: If the stream is not open.
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open. Call open() before append().")
        if not data:
            logger.debug("No data provided to append; returning without action.")
            return

        if retry_policy is None:
            retry_policy = AsyncRetry(predicate=_is_write_retryable)

        strategy = _WriteResumptionStrategy()
        buffer = io.BytesIO(data)
        attempt_count = 0

        def send_and_recv_generator(
            requests: List[BidiWriteObjectRequest],
            state: dict[str, _WriteState],
            metadata: Optional[List[Tuple[str, str]]] = None,
        ):
            async def generator():
                nonlocal attempt_count
                nonlocal requests
                attempt_count += 1
                resp = None
                write_state = state["write_state"]
                # If this is a retry or redirect, we must re-open the stream
                if attempt_count > 1 or write_state.routing_token:
                    logger.info(
                        f"Re-opening the stream with attempt_count: {attempt_count}"
                    )
                    if self.write_obj_stream and self.write_obj_stream.is_stream_open:
                        await self.write_obj_stream.close()

                    current_metadata = list(metadata) if metadata else []
                    if write_state.routing_token:
                        current_metadata.append(
                            (
                                "x-goog-request-params",
                                f"routing_token={write_state.routing_token}",
                            )
                        )
                        self._routing_token = write_state.routing_token

                    self._is_stream_open = False
                    await self.open(metadata=current_metadata)

                    write_state.persisted_size = self.persisted_size
                    write_state.write_handle = self.write_handle
                    write_state.routing_token = None

                    write_state.user_buffer.seek(write_state.persisted_size)
                    write_state.bytes_sent = write_state.persisted_size
                    write_state.bytes_since_last_flush = 0

                    requests = strategy.generate_requests(state)

                num_requests = len(requests)
                for i, chunk_req in enumerate(requests):
                    if i == num_requests - 1:
                        chunk_req.state_lookup = True
                        chunk_req.flush = True
                    await self.write_obj_stream.send(chunk_req)

                resp = await self.write_obj_stream.recv()
                if resp:
                    if resp.persisted_size is not None:
                        self.persisted_size = resp.persisted_size
                        state["write_state"].persisted_size = resp.persisted_size
                        self.offset = self.persisted_size
                    if resp.write_handle:
                        self.write_handle = resp.write_handle
                        state["write_state"].write_handle = resp.write_handle
                    self.bytes_appended_since_last_flush = 0

                yield resp

            return generator()

        # State initialization
        write_state = _WriteState(_MAX_CHUNK_SIZE_BYTES, buffer, self.flush_interval)
        write_state.write_handle = self.write_handle
        write_state.persisted_size = self.persisted_size
        write_state.bytes_sent = self.persisted_size
        write_state.bytes_since_last_flush = self.bytes_appended_since_last_flush

        retry_manager = _BidiStreamRetryManager(
            _WriteResumptionStrategy(),
            lambda r, s: send_and_recv_generator(r, s, metadata),
        )
        await retry_manager.execute({"write_state": write_state}, retry_policy)

        # Sync local markers
        self.write_obj_stream.persisted_size = write_state.persisted_size
        self.write_obj_stream.write_handle = write_state.write_handle
        self.bytes_appended_since_last_flush = write_state.bytes_since_last_flush
        self.persisted_size = write_state.persisted_size
        self.offset = write_state.persisted_size

    async def simple_flush(self) -> None:
        """Flushes the data to the server.
        Please note: Unlike `flush` it does not do `state_lookup`

        :rtype: None

        :raises ValueError: If the stream is not open (i.e., `open()` has not
            been called).
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open. Call open() before simple_flush().")

        await self.write_obj_stream.send(
            _storage_v2.BidiWriteObjectRequest(
                flush=True,
            )
        )
        self.bytes_appended_since_last_flush = 0

    async def flush(self) -> int:
        """Flushes the data to the server.

        :rtype: int
        :returns: The persisted size after flush.

        :raises ValueError: If the stream is not open (i.e., `open()` has not
            been called).
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open. Call open() before flush().")

        await self.write_obj_stream.send(
            _storage_v2.BidiWriteObjectRequest(
                flush=True,
                state_lookup=True,
            )
        )
        response = await self.write_obj_stream.recv()
        self.persisted_size = response.persisted_size
        self.offset = self.persisted_size
        self.bytes_appended_since_last_flush = 0
        return self.persisted_size

    async def close(self, finalize_on_close=False) -> Union[int, _storage_v2.Object]:
        """Closes the underlying bidi-gRPC stream.

        :type finalize_on_close: bool
        :param finalize_on_close: Finalizes the Appendable Object. No more data
          can be appended.

        rtype: Union[int, _storage_v2.Object]
        returns: Updated `self.persisted_size` by default after closing the
            bidi-gRPC stream. However, if `finalize_on_close=True` is passed,
            returns the finalized object resource.

        :raises ValueError: If the stream is not open (i.e., `open()` has not
            been called).

        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open. Call open() before close().")

        if finalize_on_close:
            return await self.finalize()

        await self.write_obj_stream.close()

        self._is_stream_open = False
        self.offset = None
        return self.persisted_size

    async def finalize(self) -> _storage_v2.Object:
        """Finalizes the Appendable Object.

        Note: Once finalized no more data can be appended.
        This method is different from `close`. if `.close()` is called data may
        still be appended to object at a later point in time by opening with
        generation number.
        (i.e. `open(..., generation=<object_generation_number>)`.
        However if `.finalize()` is called no more data can be appended to the
        object.

        rtype: google.cloud.storage_v2.types.Object
        returns: The finalized object resource.

        :raises ValueError: If the stream is not open (i.e., `open()` has not
            been called).
        """
        if not self._is_stream_open:
            raise ValueError("Stream is not open. Call open() before finalize().")

        await self.write_obj_stream.send(
            _storage_v2.BidiWriteObjectRequest(finish_write=True)
        )
        response = await self.write_obj_stream.recv()
        self.object_resource = response.resource
        self.persisted_size = self.object_resource.size
        await self.write_obj_stream.close()

        self._is_stream_open = False
        self.offset = None
        return self.object_resource

    @property
    def is_stream_open(self) -> bool:
        return self._is_stream_open

    # helper methods.
    async def append_from_string(self, data: str):
        """
        str data will be encoded to bytes using utf-8 encoding calling

        self.append(data.encode("utf-8"))
        """
        raise NotImplementedError("append_from_string is not implemented yet.")

    async def append_from_stream(self, stream_obj):
        """
        At a time read a chunk of data (16MiB) from `stream_obj`
        and call self.append(chunk)
        """
        raise NotImplementedError("append_from_stream is not implemented yet.")

    async def append_from_file(
        self, file_obj: BufferedReader, block_size: int = _DEFAULT_FLUSH_INTERVAL_BYTES
    ):
        """
        Appends data to an Appendable Object using file_handle which is opened
        for reading in binary mode.

        :type file_obj: file
        :param file_obj: A file handle opened in binary mode for reading.

        """
        while block := file_obj.read(block_size):
            await self.append(block)
