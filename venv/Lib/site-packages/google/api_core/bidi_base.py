# Copyright 2025, Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# You may obtain a copy of the License at
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Base class for bi-directional streaming RPC helpers."""


class BidiRpcBase:
    """A base class for consuming a bi-directional streaming RPC.

    This maps gRPC's built-in interface which uses a request iterator and a
    response iterator into a socket-like :func:`send` and :func:`recv`. This
    is a more useful pattern for long-running or asymmetric streams (streams
    where there is not a direct correlation between the requests and
    responses).

    This does *not* retry the stream on errors.

    Args:
        start_rpc (Union[grpc.StreamStreamMultiCallable,
                    grpc.aio.StreamStreamMultiCallable]): The gRPC method used
                    to start the RPC.
        initial_request (Union[protobuf.Message,
                Callable[[], protobuf.Message]]): The initial request to
            yield. This is useful if an initial request is needed to start the
            stream.
        metadata (Sequence[Tuple(str, str)]): RPC metadata to include in
            the request.
    """

    def __init__(self, start_rpc, initial_request=None, metadata=None):
        self._start_rpc = start_rpc
        self._initial_request = initial_request
        self._rpc_metadata = metadata
        self._request_queue = self._create_queue()
        self._request_generator = None
        self._callbacks = []
        self.call = None

    def _create_queue(self):
        """Create a queue for requests."""
        raise NotImplementedError("`_create_queue` is not implemented.")

    def add_done_callback(self, callback):
        """Adds a callback that will be called when the RPC terminates.

        This occurs when the RPC errors or is successfully terminated.

        Args:
            callback (Union[Callable[[grpc.Future], None], Callable[[Any], None]]):
                The callback to execute after gRPC call completed (success or
                failure).

                For sync streaming gRPC: Callable[[grpc.Future], None]

                For async streaming gRPC: Callable[[Any], None]
        """
        self._callbacks.append(callback)

    def _on_call_done(self, future):
        # This occurs when the RPC errors or is successfully terminated.
        # Note that grpc's "future" here can also be a grpc.RpcError.
        # See note in https://github.com/grpc/grpc/issues/10885#issuecomment-302651331
        # that `grpc.RpcError` is also `grpc.Call`.
        # for asynchronous gRPC call it would be `grpc.aio.AioRpcError`

        # Note: sync callbacks can be limiting for async code, because you can't
        # await anything in a sync callback.
        for callback in self._callbacks:
            callback(future)

    @property
    def is_active(self):
        """True if the gRPC call is not done yet."""
        raise NotImplementedError("`is_active` is not implemented.")

    @property
    def pending_requests(self):
        """Estimate of the number of queued requests."""
        return self._request_queue.qsize()
