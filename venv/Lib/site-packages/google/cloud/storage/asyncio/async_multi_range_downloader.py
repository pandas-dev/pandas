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

from __future__ import annotations

import asyncio
import logging
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple

from google.api_core import exceptions
from google.api_core.retry_async import AsyncRetry
from google.rpc import status_pb2

from google.cloud import _storage_v2
from google.cloud.storage._helpers import generate_random_56_bit_integer
from google.cloud.storage.asyncio._stream_multiplexer import (
    _StreamEnd,
    _StreamError,
    _StreamMultiplexer,
)
from google.cloud.storage.asyncio.async_grpc_client import (
    AsyncGrpcClient,
)
from google.cloud.storage.asyncio.async_read_object_stream import (
    _AsyncReadObjectStream,
)
from google.cloud.storage.asyncio.retry._helpers import _handle_redirect
from google.cloud.storage.asyncio.retry.bidi_stream_retry_manager import (
    _BidiStreamRetryManager,
)
from google.cloud.storage.asyncio.retry.reads_resumption_strategy import (
    _DownloadState,
    _ReadResumptionStrategy,
)
from google.cloud.storage.exceptions import DataCorruption

from ._utils import raise_if_no_fast_crc32c

_MAX_READ_RANGES_PER_BIDI_READ_REQUEST = 100
_BIDI_READ_REDIRECTED_TYPE_URL = (
    "type.googleapis.com/google.storage.v2.BidiReadObjectRedirectedError"
)

logger = logging.getLogger(__name__)


def _is_read_retryable(exc):
    """Predicate to determine if a read operation should be retried."""
    if isinstance(
        exc,
        (
            exceptions.InternalServerError,
            exceptions.ServiceUnavailable,
            exceptions.DeadlineExceeded,
            exceptions.TooManyRequests,
        ),
    ):
        return True

    if not isinstance(exc, exceptions.Aborted) or not exc.errors:
        return False

    try:
        grpc_error = exc.errors[0]
        trailers = grpc_error.trailing_metadata()
        if not trailers:
            return False

        status_details_bin = next(
            (v for k, v in trailers if k == "grpc-status-details-bin"), None
        )

        if not status_details_bin:
            return False

        status_proto = status_pb2.Status()
        status_proto.ParseFromString(status_details_bin)
        return any(
            detail.type_url == _BIDI_READ_REDIRECTED_TYPE_URL
            for detail in status_proto.details
        )
    except Exception as e:
        logger.error(f"Error parsing status_details_bin: {e}")
        return False


class AsyncMultiRangeDownloader:
    """Provides an interface for downloading multiple ranges of a GCS ``Object``
    concurrently.

    Example usage:

    .. code-block:: python

        client = AsyncGrpcClient()
        mrd = await AsyncMultiRangeDownloader.create_mrd(
            client, bucket_name="chandrasiri-rs", object_name="test_open9"
        )
        my_buff1 = open('my_fav_file.txt', 'wb')
        my_buff2 = BytesIO()
        my_buff3 = BytesIO()
        my_buff4 = any_object_which_provides_BytesIO_like_interface()
        await mrd.download_ranges(
            [
                # (start_byte, bytes_to_read, writeable_buffer)
                (0, 100, my_buff1),
                (100, 20, my_buff2),
                (200, 123, my_buff3),
                (300, 789, my_buff4),
            ]
        )

        # verify data in buffers...
        assert my_buff2.getbuffer().nbytes == 20


    """

    @classmethod
    async def create_mrd(
        cls,
        client: AsyncGrpcClient,
        bucket_name: str,
        object_name: str,
        generation: Optional[int] = None,
        read_handle: Optional[_storage_v2.BidiReadHandle] = None,
        retry_policy: Optional[AsyncRetry] = None,
        metadata: Optional[List[Tuple[str, str]]] = None,
        **kwargs,
    ) -> AsyncMultiRangeDownloader:
        """Initializes a MultiRangeDownloader and opens the underlying bidi-gRPC
        object for reading.

        :type client: :class:`~google.cloud.storage.asyncio.async_grpc_client.AsyncGrpcClient`
        :param client: The asynchronous client to use for making API requests.

        :type bucket_name: str
        :param bucket_name: The name of the bucket containing the object.

        :type object_name: str
        :param object_name: The name of the object to be read.

        :type generation: int
        :param generation: (Optional) If present, selects a specific
                                  revision of this object.

        :type read_handle: _storage_v2.BidiReadHandle
        :param read_handle: (Optional) An existing handle for reading the object.
                            If provided, opening the bidi-gRPC connection will be faster.

        :type retry_policy: :class:`~google.api_core.retry_async.AsyncRetry`
        :param retry_policy: (Optional) The retry policy to use for the ``open`` operation.

        :type metadata: List[Tuple[str, str]]
        :param metadata: (Optional) The metadata to be sent with the ``open`` request.

        :rtype: :class:`~google.cloud.storage.asyncio.async_multi_range_downloader.AsyncMultiRangeDownloader`
        :returns: An initialized AsyncMultiRangeDownloader instance for reading.
        """
        mrd = cls(
            client,
            bucket_name,
            object_name,
            generation=generation,
            read_handle=read_handle,
            **kwargs,
        )
        await mrd.open(retry_policy=retry_policy, metadata=metadata)
        return mrd

    def __init__(
        self,
        client: AsyncGrpcClient,
        bucket_name: str,
        object_name: str,
        generation: Optional[int] = None,
        read_handle: Optional[_storage_v2.BidiReadHandle] = None,
        **kwargs,
    ) -> None:
        """Constructor for AsyncMultiRangeDownloader, clients are not adviced to
         use it directly. Instead it's adviced to use the classmethod `create_mrd`.

        :type client: :class:`~google.cloud.storage.asyncio.async_grpc_client.AsyncGrpcClient`
        :param client: The asynchronous client to use for making API requests.

        :type bucket_name: str
        :param bucket_name: The name of the bucket containing the object.

        :type object_name: str
        :param object_name: The name of the object to be read.

        :type generation: int
        :param generation: (Optional) If present, selects a specific revision of
                                  this object.

        :type read_handle: _storage_v2.BidiReadHandle
        :param read_handle: (Optional) An existing read handle.
        """
        if "generation_number" in kwargs:
            if generation is not None:
                raise TypeError(
                    "Cannot set both 'generation' and 'generation_number'. "
                    "Use 'generation' for new code."
                )
            logger.warning(
                "'generation_number' is deprecated and will be removed in a future "
                "major release. Please use 'generation' instead."
            )
            generation = kwargs.pop("generation_number")

        self.client = client
        self.bucket_name = bucket_name
        self.object_name = object_name
        self.generation = generation
        self.read_handle: Optional[_storage_v2.BidiReadHandle] = read_handle
        self.read_obj_str: Optional[_AsyncReadObjectStream] = None
        self._is_stream_open: bool = False
        self._routing_token: Optional[str] = None
        self._multiplexer: Optional[_StreamMultiplexer] = None
        self.persisted_size: Optional[int] = None  # updated after opening the stream
        self._open_retries: int = 0
        self.is_finalized: bool = False
        self.full_obj_server_crc32c: Optional[int] = None

    async def __aenter__(self):
        """Opens the underlying bidi-gRPC connection to read from the object."""
        await self.open()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Closes the underlying bidi-gRPC connection."""
        if self.is_stream_open:
            await self.close()

    def _on_open_error(self, exc):
        """Extracts routing token and read handle on redirect error during open."""
        logger.warning(f"Error occurred while opening MRD: {exc}")
        routing_token, read_handle = _handle_redirect(exc)
        if routing_token:
            self._routing_token = routing_token
        if read_handle:
            self.read_handle = read_handle

    async def open(
        self,
        retry_policy: Optional[AsyncRetry] = None,
        metadata: Optional[List[Tuple[str, str]]] = None,
    ) -> None:
        """Opens the bidi-gRPC connection to read from the object."""
        if self._is_stream_open:
            raise ValueError("Underlying bidi-gRPC stream is already open")

        if retry_policy is None:

            def on_error_wrapper(exc):
                self._open_retries += 1
                self._on_open_error(exc)

            retry_policy = AsyncRetry(
                predicate=_is_read_retryable, on_error=on_error_wrapper
            )
        else:
            original_on_error = retry_policy._on_error

            def combined_on_error(exc):
                self._open_retries += 1
                self._on_open_error(exc)
                if original_on_error:
                    original_on_error(exc)

            retry_policy = AsyncRetry(
                predicate=_is_read_retryable,
                initial=retry_policy._initial,
                maximum=retry_policy._maximum,
                multiplier=retry_policy._multiplier,
                deadline=retry_policy._deadline,
                on_error=combined_on_error,
            )

        async def _do_open():
            current_metadata = list(metadata) if metadata else []

            # Cleanup stream from previous failed attempt, if any.
            if self.read_obj_str:
                if self.read_obj_str.is_stream_open:
                    try:
                        await self.read_obj_str.close()
                    except exceptions.GoogleAPICallError as e:
                        logger.warning(
                            f"Failed to close existing stream during resumption: {e}"
                        )
                self.read_obj_str = None
                self._is_stream_open = False

            self.read_obj_str = _AsyncReadObjectStream(
                client=self.client.grpc_client,
                bucket_name=self.bucket_name,
                object_name=self.object_name,
                generation_number=self.generation,
                read_handle=self.read_handle,
            )

            if self._routing_token:
                current_metadata.append(
                    ("x-goog-request-params", f"routing_token={self._routing_token}")
                )
                self._routing_token = None

            await self.read_obj_str.open(
                metadata=current_metadata if current_metadata else None
            )

            if self.read_obj_str.generation_number:
                self.generation = self.read_obj_str.generation_number
            if self.read_obj_str.read_handle:
                self.read_handle = self.read_obj_str.read_handle
            if self.read_obj_str.persisted_size is not None:
                self.persisted_size = self.read_obj_str.persisted_size
            self.is_finalized = self.read_obj_str.is_finalized
            self.full_obj_server_crc32c = self.read_obj_str.full_obj_server_crc32c

            self._is_stream_open = True

        await retry_policy(_do_open)()
        self._multiplexer = _StreamMultiplexer(self.read_obj_str)

    def _create_stream_factory(self, state, metadata):
        """Create a factory that opens a new stream with current routing state."""

        async def factory():
            current_handle = state.get("read_handle")
            current_token = state.get("routing_token")

            stream = _AsyncReadObjectStream(
                client=self.client.grpc_client,
                bucket_name=self.bucket_name,
                object_name=self.object_name,
                generation_number=self.generation,
                read_handle=current_handle,
            )

            current_metadata = list(metadata) if metadata else []
            if current_token:
                current_metadata.append(
                    (
                        "x-goog-request-params",
                        f"routing_token={current_token}",
                    )
                )

            await stream.open(metadata=current_metadata if current_metadata else None)

            if stream.generation_number:
                self.generation = stream.generation_number
            if stream.read_handle:
                self.read_handle = stream.read_handle
            self.is_finalized = stream.is_finalized
            self.full_obj_server_crc32c = stream.full_obj_server_crc32c

            self.read_obj_str = stream
            self._is_stream_open = True

            return stream

        return factory

    async def download_ranges(
        self,
        read_ranges: List[Tuple[int, int, BytesIO]],
        lock: asyncio.Lock = None,
        retry_policy: Optional[AsyncRetry] = None,
        metadata: Optional[List[Tuple[str, str]]] = None,
        enable_checksum: bool = True,
    ) -> None:
        """Downloads multiple byte ranges from the object into the buffers
        provided by user with automatic retries.

        :type read_ranges: List[Tuple[int, int, "BytesIO"]]
        :param read_ranges: A list of tuples, where each tuple represents a
            combination of byte_range and writeable buffer in format -
            (`start_byte`, `bytes_to_read`, `writeable_buffer`). Buffer has
            to be provided by the user, and user has to make sure appropriate
            memory is available in the application to avoid out-of-memory crash.

            Special cases:
            if the value of `bytes_to_read` is 0, it'll be interpreted as
            download all contents until the end of the file from `start_byte`.
            Examples:
                * (0, 0, buffer) : downloads 0 to end , i.e. entire object.
                * (100, 0, buffer) : downloads from 100 to end.

        :type lock: asyncio.Lock
        :param lock: (Deprecated) This parameter is deprecated and has no effect.

        :type retry_policy: :class:`~google.api_core.retry_async.AsyncRetry`
        :param retry_policy: (Optional) The retry policy to use for the operation.

        :type metadata: List[Tuple[str, str]]
        :param metadata: (Optional) The metadata to be sent with the request.

        :type enable_checksum: bool
        :param enable_checksum: (Optional) If True, checksums are verified for downloaded data. Defaults to True.

        :raises ValueError: if the underlying bidi-GRPC stream is not open.
        :raises ValueError: if the length of read_ranges is more than 1000.
        :raises DataCorruption: if a checksum mismatch is detected while reading data.

        """

        if len(read_ranges) > 1000:
            raise ValueError(
                "Invalid input - length of read_ranges cannot be more than 1000"
            )

        if enable_checksum:
            raise_if_no_fast_crc32c()

        if not self._is_stream_open:
            raise ValueError("Underlying bidi-gRPC stream is not open")

        if retry_policy is None:
            retry_policy = AsyncRetry(predicate=_is_read_retryable)

        # Initialize Global State for Retry Strategy
        download_states = {}
        for read_range in read_ranges:
            read_id = generate_random_56_bit_integer()
            # Unpack tuple into self-documenting variable names to improve readability.
            offset, length, user_buffer = read_range

            # Heuristic to detect full object reads:
            # - Implicit full object read: start offset is 0 and length is 0 (read all).
            # - Explicit full object read: start offset is 0 and length matches the exact persisted size.
            is_full_object_read = (offset == 0 and length == 0) or (
                self.persisted_size is not None
                and offset == 0
                and length == self.persisted_size
            )
            download_states[read_id] = _DownloadState(
                initial_offset=offset,
                initial_length=length,
                user_buffer=user_buffer,
                is_full_object_read=is_full_object_read,
            )

        initial_state = {
            "download_states": download_states,
            "read_handle": self.read_handle,
            "routing_token": None,
            "enable_checksum": enable_checksum,
            "full_obj_server_crc32c": self.full_obj_server_crc32c,
        }

        read_ids = set(download_states.keys())
        queue = self._multiplexer.register(read_ids)

        try:
            attempt_count = 0
            last_broken_generation = None

            def send_and_recv_via_multiplexer(
                requests: List[_storage_v2.ReadRange],
                state: Dict[str, Any],
            ):
                async def generator():
                    nonlocal attempt_count, last_broken_generation
                    attempt_count += 1

                    if attempt_count > 1:
                        logger.info(
                            f"Resuming download (attempt {attempt_count}) for {len(requests)} ranges."
                        )

                    # Reopen stream if needed
                    should_reopen = (
                        attempt_count > 1 and last_broken_generation is not None
                    ) or (attempt_count == 1 and metadata is not None)
                    if should_reopen:
                        broken_gen = (
                            last_broken_generation
                            if attempt_count > 1
                            else self._multiplexer.stream_generation
                        )
                        stream_factory = self._create_stream_factory(state, metadata)
                        await self._multiplexer.reopen_stream(
                            broken_gen, stream_factory
                        )

                    stream_generation = self._multiplexer.stream_generation

                    # Send Requests
                    pending_read_ids = {r.read_id for r in requests}
                    for i in range(
                        0, len(requests), _MAX_READ_RANGES_PER_BIDI_READ_REQUEST
                    ):
                        batch = requests[i : i + _MAX_READ_RANGES_PER_BIDI_READ_REQUEST]
                        try:
                            await self._multiplexer.send(
                                _storage_v2.BidiReadObjectRequest(read_ranges=batch)
                            )
                        except Exception:
                            last_broken_generation = stream_generation
                            raise

                    # Receive Responses
                    while pending_read_ids:
                        item = await queue.get()

                        if isinstance(item, _StreamEnd):
                            if pending_read_ids:
                                last_broken_generation = stream_generation
                                raise exceptions.ServiceUnavailable(
                                    "Stream ended with pending read_ids"
                                )
                            break

                        if isinstance(item, _StreamError):
                            if item.generation < stream_generation:
                                continue  # stale error, skip
                            last_broken_generation = item.generation
                            raise item.exception

                        # Track completion
                        if item.object_data_ranges:
                            for data_range in item.object_data_ranges:
                                if data_range.range_end:
                                    pending_read_ids.discard(
                                        data_range.read_range.read_id
                                    )
                        yield item

                return generator()

            strategy = _ReadResumptionStrategy()
            retry_manager = _BidiStreamRetryManager(
                strategy, send_and_recv_via_multiplexer
            )

            try:
                await retry_manager.execute(initial_state, retry_policy)
            except DataCorruption:
                if self.is_stream_open:
                    await self.close()
                raise

            if initial_state.get("read_handle"):
                self.read_handle = initial_state["read_handle"]
        finally:
            if self._multiplexer is not None:
                self._multiplexer.unregister(read_ids)

    async def close(self):
        """
        Closes the underlying bidi-gRPC connection.
        """
        if not self._is_stream_open:
            raise ValueError("Underlying bidi-gRPC stream is not open")

        if self._multiplexer:
            await self._multiplexer.close()
            self._multiplexer = None

        if self.read_obj_str:
            try:
                await self.read_obj_str.close()
            except (asyncio.CancelledError, exceptions.GoogleAPICallError):
                pass
        self.read_obj_str = None
        self._is_stream_open = False

    @property
    def is_stream_open(self) -> bool:
        return self._is_stream_open
