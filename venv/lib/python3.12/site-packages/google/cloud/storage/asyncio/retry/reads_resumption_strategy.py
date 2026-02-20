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

from typing import Any, Dict, List, IO
import logging

from google_crc32c import Checksum
from google.cloud import _storage_v2 as storage_v2
from google.cloud.storage.exceptions import DataCorruption
from google.cloud.storage.asyncio.retry._helpers import (
    _handle_redirect,
)
from google.cloud.storage.asyncio.retry.base_strategy import (
    _BaseResumptionStrategy,
)


_BIDI_READ_REDIRECTED_TYPE_URL = (
    "type.googleapis.com/google.storage.v2.BidiReadObjectRedirectedError"
)
logger = logging.getLogger(__name__)


class _DownloadState:
    """A helper class to track the state of a single range download."""

    def __init__(
        self, initial_offset: int, initial_length: int, user_buffer: IO[bytes]
    ):
        self.initial_offset = initial_offset
        self.initial_length = initial_length
        self.user_buffer = user_buffer
        self.bytes_written = 0
        self.next_expected_offset = initial_offset
        self.is_complete = False


class _ReadResumptionStrategy(_BaseResumptionStrategy):
    """The concrete resumption strategy for bidi reads."""

    def generate_requests(self, state: Dict[str, Any]) -> List[storage_v2.ReadRange]:
        """Generates new ReadRange requests for all incomplete downloads.

        :type state: dict
        :param state: A dictionary mapping a read_id to its corresponding
                  _DownloadState object.
        """
        pending_requests = []
        download_states: Dict[int, _DownloadState] = state["download_states"]

        for read_id, read_state in download_states.items():
            if not read_state.is_complete:
                new_offset = read_state.initial_offset + read_state.bytes_written

                # Calculate remaining length. If initial_length is 0 (read to end),
                # it stays 0. Otherwise, subtract bytes_written.
                new_length = 0
                if read_state.initial_length > 0:
                    new_length = read_state.initial_length - read_state.bytes_written

                new_request = storage_v2.ReadRange(
                    read_offset=new_offset,
                    read_length=new_length,
                    read_id=read_id,
                )
                pending_requests.append(new_request)
        return pending_requests

    def update_state_from_response(
        self, response: storage_v2.BidiReadObjectResponse, state: Dict[str, Any]
    ) -> None:
        """Processes a server response, performs integrity checks, and updates state."""

        # Capture read_handle if provided.
        if response.read_handle:
            state["read_handle"] = response.read_handle

        download_states = state["download_states"]

        for object_data_range in response.object_data_ranges:
            # Ignore empty ranges or ranges for IDs not in our state
            # (e.g., from a previously cancelled request on the same stream).
            if not object_data_range.read_range:
                logger.warning(
                    "Received response with missing read_range field; ignoring."
                )
                continue

            read_id = object_data_range.read_range.read_id
            if read_id not in download_states:
                logger.warning(
                    f"Received data for unknown or stale read_id {read_id}; ignoring."
                )
                continue

            read_state = download_states[read_id]

            # Offset Verification
            chunk_offset = object_data_range.read_range.read_offset
            if chunk_offset != read_state.next_expected_offset:
                raise DataCorruption(
                    response,
                    f"Offset mismatch for read_id {read_id}. "
                    f"Expected {read_state.next_expected_offset}, got {chunk_offset}",
                )

            # Checksum Verification
            # We must validate data before updating state or writing to buffer.
            data = object_data_range.checksummed_data.content
            server_checksum = object_data_range.checksummed_data.crc32c

            if server_checksum is not None:
                client_checksum = int.from_bytes(Checksum(data).digest(), "big")
                if server_checksum != client_checksum:
                    raise DataCorruption(
                        response,
                        f"Checksum mismatch for read_id {read_id}. "
                        f"Server sent {server_checksum}, client calculated {client_checksum}.",
                    )

            # Update State & Write Data
            chunk_size = len(data)
            read_state.user_buffer.write(data)
            read_state.bytes_written += chunk_size
            read_state.next_expected_offset += chunk_size

            # Final Byte Count Verification
            if object_data_range.range_end:
                read_state.is_complete = True
                if (
                    read_state.initial_length != 0
                    and read_state.bytes_written > read_state.initial_length
                ):
                    raise DataCorruption(
                        response,
                        f"Byte count mismatch for read_id {read_id}. "
                        f"Expected {read_state.initial_length}, got {read_state.bytes_written}",
                    )

    async def recover_state_on_failure(self, error: Exception, state: Any) -> None:
        """Handles BidiReadObjectRedirectedError for reads."""
        routing_token, read_handle = _handle_redirect(error)
        if routing_token:
            state["routing_token"] = routing_token
        if read_handle:
            state["read_handle"] = read_handle
