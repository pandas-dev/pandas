# Copyright 2025, Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Asynchronous bi-directional streaming RPC helpers."""

import asyncio
import logging
from typing import Callable, Optional, Union

from grpc import aio

from google.api_core import exceptions
from google.api_core.bidi_base import BidiRpcBase

from google.protobuf.message import Message as ProtobufMessage


_LOGGER = logging.getLogger(__name__)


class _AsyncRequestQueueGenerator:
    """_AsyncRequestQueueGenerator is a helper class for sending asynchronous
      requests to a gRPC stream from a Queue.

    This generator takes asynchronous requests off a given `asyncio.Queue` and
    yields them to gRPC.

    It's useful when you have an indeterminate, indefinite, or otherwise
    open-ended set of requests to send through a request-streaming (or
    bidirectional) RPC.

    Example::

        requests = _AsyncRequestQueueGenerator(q)
        call = await stub.StreamingRequest(requests)
        requests.call = call

        async for response in call:
            print(response)
            await q.put(...)

    Args:
        queue (asyncio.Queue): The request queue.
        initial_request (Union[ProtobufMessage,
                Callable[[], ProtobufMessage]]): The initial request to
            yield. This is done independently of the request queue to allow for
            easily restarting streams that require some initial configuration
            request.
    """

    def __init__(
        self,
        queue: asyncio.Queue,
        initial_request: Optional[
            Union[ProtobufMessage, Callable[[], ProtobufMessage]]
        ] = None,
    ) -> None:
        self._queue = queue
        self._initial_request = initial_request
        self.call: Optional[aio.Call] = None

    def _is_active(self) -> bool:
        """Returns true if the call is not set or not completed."""
        # Note: there is a possibility that this starts *before* the call
        # property is set. So we have to check if self.call is set before
        # seeing if it's active. We need to return True if self.call is None.
        # See https://github.com/googleapis/python-api-core/issues/560.
        return self.call is None or not self.call.done()

    async def __aiter__(self):
        # The reason this is necessary is because it lets the user have
        # control on when they would want to send requests proto messages
        # instead of sending all of them initially.
        #
        # This is achieved via asynchronous queue (asyncio.Queue),
        # gRPC awaits until there's a message in the queue.
        #
        # Finally, it allows for retrying without swapping queues because if
        # it does pull an item off the queue when the RPC is inactive, it'll
        # immediately put it back and then exit. This is necessary because
        # yielding the item in this case will cause gRPC to discard it. In
        # practice, this means that the order of messages is not guaranteed.
        # If preserving order is necessary it would be easy to use a priority
        # queue.
        if self._initial_request is not None:
            if callable(self._initial_request):
                yield self._initial_request()
            else:
                yield self._initial_request

        while True:
            item = await self._queue.get()

            # The consumer explicitly sent "None", indicating that the request
            # should end.
            if item is None:
                _LOGGER.debug("Cleanly exiting request generator.")
                return

            if not self._is_active():
                # We have an item, but the call is closed. We should put the
                # item back on the queue so that the next call can consume it.
                await self._queue.put(item)
                _LOGGER.debug(
                    "Inactive call, replacing item on queue and exiting "
                    "request generator."
                )
                return

            yield item


class AsyncBidiRpc(BidiRpcBase):
    """A helper for consuming a async bi-directional streaming RPC.

    This maps gRPC's built-in interface which uses a request iterator and a
    response iterator into a socket-like :func:`send` and :func:`recv`. This
    is a more useful pattern for long-running or asymmetric streams (streams
    where there is not a direct correlation between the requests and
    responses).

    Example::

        initial_request = example_pb2.StreamingRpcRequest(
            setting='example')
        rpc = AsyncBidiRpc(
            stub.StreamingRpc,
            initial_request=initial_request,
            metadata=[('name', 'value')]
        )

        await rpc.open()

        while rpc.is_active:
            print(await rpc.recv())
            await rpc.send(example_pb2.StreamingRpcRequest(
                data='example'))

        await rpc.close()

    This does *not* retry the stream on errors.

    Args:
        start_rpc (grpc.aio.StreamStreamMultiCallable): The gRPC method used to
            start the RPC.
        initial_request (Union[ProtobufMessage,
                Callable[[], ProtobufMessage]]): The initial request to
            yield. This is useful if an initial request is needed to start the
            stream.
        metadata (Sequence[Tuple(str, str)]): RPC metadata to include in
            the request.
    """

    def _create_queue(self) -> asyncio.Queue:
        """Create a queue for requests."""
        return asyncio.Queue()

    async def open(self) -> None:
        """Opens the stream."""
        if self.is_active:
            raise ValueError("Cannot open an already open stream.")

        request_generator = _AsyncRequestQueueGenerator(
            self._request_queue, initial_request=self._initial_request
        )
        try:
            call = await self._start_rpc(request_generator, metadata=self._rpc_metadata)
        except exceptions.GoogleAPICallError as exc:
            # The original `grpc.aio.AioRpcError` (which is usually also a
            # `grpc.aio.Call`) is available from the ``response`` property on
            # the mapped exception.
            self._on_call_done(exc.response)
            raise

        request_generator.call = call

        # TODO: api_core should expose the future interface for wrapped
        # callables as well.
        if hasattr(call, "_wrapped"):  # pragma: NO COVER
            call._wrapped.add_done_callback(self._on_call_done)
        else:
            call.add_done_callback(self._on_call_done)

        self._request_generator = request_generator
        self.call = call

    async def close(self) -> None:
        """Closes the stream."""
        if self.call is not None:
            self.call.cancel()

        # Put None in request queue to signal termination.
        await self._request_queue.put(None)
        self._request_generator = None
        self._initial_request = None
        self._callbacks = []
        # Don't set self.call to None. Keep it around so that send/recv can
        # raise the error.

    async def send(self, request: ProtobufMessage) -> None:
        """Queue a message to be sent on the stream.

        If the underlying RPC has been closed, this will raise.

        Args:
            request (ProtobufMessage): The request to send.
        """
        if self.call is None:
            raise ValueError("Cannot send on an RPC stream that has never been opened.")

        if not self.call.done():
            await self._request_queue.put(request)
        else:
            # calling read should cause the call to raise.
            await self.call.read()

    async def recv(self) -> ProtobufMessage:
        """Wait for a message to be returned from the stream.

        If the underlying RPC has been closed, this will raise.

        Returns:
            ProtobufMessage: The received message.
        """
        if self.call is None:
            raise ValueError("Cannot recv on an RPC stream that has never been opened.")

        return await self.call.read()

    @property
    def is_active(self) -> bool:
        """Whether the stream is currently open and active."""
        return self.call is not None and not self.call.done()
