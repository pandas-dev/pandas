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

import logging
from typing import Tuple, Optional

from google.api_core import exceptions
from google.cloud._storage_v2.types import (
    BidiReadObjectRedirectedError,
    BidiWriteObjectRedirectedError,
)
from google.rpc import status_pb2

_BIDI_READ_REDIRECTED_TYPE_URL = (
    "type.googleapis.com/google.storage.v2.BidiReadObjectRedirectedError"
)
_BIDI_WRITE_REDIRECTED_TYPE_URL = (
    "type.googleapis.com/google.storage.v2.BidiWriteObjectRedirectedError"
)
logger = logging.getLogger(__name__)


def _handle_redirect(
    exc: Exception,
) -> Tuple[Optional[str], Optional[bytes]]:
    """
    Extracts routing token and read handle from a gRPC error.

    :type exc: Exception
    :param exc: The exception to parse.

    :rtype: Tuple[Optional[str], Optional[bytes]]
    :returns: A tuple of (routing_token, read_handle).
    """
    routing_token = None
    read_handle = None

    grpc_error = None
    if isinstance(exc, exceptions.Aborted) and exc.errors:
        grpc_error = exc.errors[0]

    if grpc_error:
        if isinstance(grpc_error, BidiReadObjectRedirectedError):
            routing_token = grpc_error.routing_token
            if grpc_error.read_handle:
                read_handle = grpc_error.read_handle
            return routing_token, read_handle

        if hasattr(grpc_error, "trailing_metadata"):
            trailers = grpc_error.trailing_metadata()
            if not trailers:
                return None, None

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
                        if detail.type_url == _BIDI_READ_REDIRECTED_TYPE_URL:
                            redirect_proto = BidiReadObjectRedirectedError.deserialize(
                                detail.value
                            )
                            if redirect_proto.routing_token:
                                routing_token = redirect_proto.routing_token
                            if redirect_proto.read_handle:
                                read_handle = redirect_proto.read_handle
                            break
                except Exception as e:
                    logger.error(f"Error unpacking redirect: {e}")

    return routing_token, read_handle


def _extract_bidi_writes_redirect_proto(exc: Exception):
    grpc_error = None
    if isinstance(exc, exceptions.Aborted) and exc.errors:
        grpc_error = exc.errors[0]

    if grpc_error:
        if isinstance(grpc_error, BidiWriteObjectRedirectedError):
            return grpc_error

        if hasattr(grpc_error, "trailing_metadata"):
            trailers = grpc_error.trailing_metadata()
            if not trailers:
                return

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
                            redirect_proto = BidiWriteObjectRedirectedError.deserialize(
                                detail.value
                            )
                            return redirect_proto
                except Exception:
                    logger.error("Error unpacking redirect details from gRPC error.")
                    pass
