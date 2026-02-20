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

from typing import Any, Dict, IO, List, Optional, Union

import google_crc32c
from google.cloud._storage_v2.types import storage as storage_type
from google.cloud._storage_v2.types.storage import BidiWriteObjectRedirectedError
from google.cloud.storage.asyncio.retry.base_strategy import (
    _BaseResumptionStrategy,
)
from google.cloud.storage.asyncio.retry._helpers import (
    _extract_bidi_writes_redirect_proto,
)


class _WriteState:
    """A helper class to track the state of a single upload operation.

    :type chunk_size: int
    :param chunk_size: The size of chunks to write to the server.

    :type user_buffer: IO[bytes]
    :param user_buffer: The data source.

    :type flush_interval: int
    :param flush_interval: The flush interval at which the data is flushed.
    """

    def __init__(
        self,
        chunk_size: int,
        user_buffer: IO[bytes],
        flush_interval: int,
    ):
        self.chunk_size = chunk_size
        self.user_buffer = user_buffer
        self.persisted_size: int = 0
        self.bytes_sent: int = 0
        self.bytes_since_last_flush: int = 0
        self.flush_interval: int = flush_interval
        self.write_handle: Union[bytes, storage_type.BidiWriteHandle, None] = None
        self.routing_token: Optional[str] = None
        self.is_finalized: bool = False


class _WriteResumptionStrategy(_BaseResumptionStrategy):
    """The concrete resumption strategy for bidi writes."""

    def generate_requests(
        self, state: Dict[str, Any]
    ) -> List[storage_type.BidiWriteObjectRequest]:
        """Generates BidiWriteObjectRequests to resume or continue the upload.

        This method is not applicable for `open` methods.
        """
        write_state: _WriteState = state["write_state"]

        requests = []
        # The buffer should already be seeked to the correct position (persisted_size)
        # by the `recover_state_on_failure` method before this is called.
        while not write_state.is_finalized:
            chunk = write_state.user_buffer.read(write_state.chunk_size)

            # End of File detection
            if not chunk:
                break

            checksummed_data = storage_type.ChecksummedData(content=chunk)
            checksum = google_crc32c.Checksum(chunk)
            checksummed_data.crc32c = int.from_bytes(checksum.digest(), "big")

            request = storage_type.BidiWriteObjectRequest(
                write_offset=write_state.bytes_sent,
                checksummed_data=checksummed_data,
            )
            chunk_len = len(chunk)
            write_state.bytes_sent += chunk_len
            write_state.bytes_since_last_flush += chunk_len

            if write_state.bytes_since_last_flush >= write_state.flush_interval:
                request.flush = True
                # reset counter after marking flush
                write_state.bytes_since_last_flush = 0

            requests.append(request)
        return requests

    def update_state_from_response(
        self, response: storage_type.BidiWriteObjectResponse, state: Dict[str, Any]
    ) -> None:
        """Processes a server response and updates the write state."""
        write_state: _WriteState = state["write_state"]
        if response is None:
            return
        if response.persisted_size:
            write_state.persisted_size = response.persisted_size

        if response.write_handle:
            write_state.write_handle = response.write_handle

        if response.resource:
            write_state.persisted_size = response.resource.size
            if response.resource.finalize_time:
                write_state.is_finalized = True

    async def recover_state_on_failure(
        self, error: Exception, state: Dict[str, Any]
    ) -> None:
        """
        Handles errors, specifically BidiWriteObjectRedirectedError, and rewinds state.

        This method rewinds the user buffer and internal byte tracking to the
        last confirmed 'persisted_size' from the server.
        """
        write_state: _WriteState = state["write_state"]

        redirect_proto = None

        if isinstance(error, BidiWriteObjectRedirectedError):
            redirect_proto = error
        else:
            redirect_proto = _extract_bidi_writes_redirect_proto(error)

        # Extract routing token and potentially a new write handle for redirection.
        if redirect_proto:
            if redirect_proto.routing_token:
                write_state.routing_token = redirect_proto.routing_token
            if redirect_proto.write_handle:
                write_state.write_handle = redirect_proto.write_handle

        # We must assume any data sent beyond 'persisted_size' was lost.
        # Reset the user buffer to the last known good byte confirmed by the server.
        write_state.user_buffer.seek(write_state.persisted_size)
        write_state.bytes_sent = write_state.persisted_size
        write_state.bytes_since_last_flush = 0
