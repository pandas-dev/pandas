# Copyright 2026 Google LLC
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
from typing import Awaitable, Callable, Dict, Optional, Set

import grpc

from google.cloud import _storage_v2
from google.cloud.storage.asyncio.async_read_object_stream import (
    _AsyncReadObjectStream,
)

logger = logging.getLogger(__name__)

_DEFAULT_QUEUE_MAX_SIZE = 100
_DEFAULT_PUT_TIMEOUT_SECONDS = 20.0


class _StreamError:
    """Wraps an error with the stream generation that produced it."""

    def __init__(self, exception: Exception, generation: int):
        self.exception = exception
        self.generation = generation


class _StreamEnd:
    """Signals the stream closed normally."""

    pass


class _StreamMultiplexer:
    """Multiplexes concurrent download tasks over a single bidi-gRPC stream.

    Routes responses from a background recv loop to per-task asyncio.Queues
    keyed by read_id. Coordinates stream reopening via generation-gated
    locking.

    A slow consumer on one task will slow down the entire shared connection
    due to bounded queue backpressure propagating through gRPC flow control.
    """

    def __init__(
        self,
        stream: _AsyncReadObjectStream,
        queue_max_size: int = _DEFAULT_QUEUE_MAX_SIZE,
    ):
        self._stream = stream
        self._stream_generation: int = 0
        self._queues: Dict[int, asyncio.Queue] = {}
        self._reopen_lock = asyncio.Lock()
        self._recv_task: Optional[asyncio.Task] = None
        self._queue_max_size = queue_max_size

    @property
    def stream_generation(self) -> int:
        return self._stream_generation

    def register(self, read_ids: Set[int]) -> asyncio.Queue:
        """Register read_ids for a task and return its response queue."""
        queue = asyncio.Queue(maxsize=self._queue_max_size)
        for read_id in read_ids:
            self._queues[read_id] = queue
        return queue

    def unregister(self, read_ids: Set[int]) -> None:
        """Remove read_ids from routing."""
        for read_id in read_ids:
            self._queues.pop(read_id, None)

    def _get_unique_queues(self) -> Set[asyncio.Queue]:
        return set(self._queues.values())

    async def _put_with_timeout(self, queue: asyncio.Queue, item) -> None:
        """Slow-path put: wait up to _DEFAULT_PUT_TIMEOUT_SECONDS, else drop.

        Callers should attempt ``queue.put_nowait(item)`` first and only call
        this when it raises :class:`asyncio.QueueFull`.
        """
        try:
            await asyncio.wait_for(
                queue.put(item), timeout=_DEFAULT_PUT_TIMEOUT_SECONDS
            )
        except asyncio.TimeoutError:
            if queue not in self._get_unique_queues():
                logger.debug("Dropped item for unregistered queue.")
            else:
                logger.warning(
                    "Queue full for too long. Dropping item to prevent multiplexer hang."
                )

    async def _put_to_queues(self, queues, item) -> None:
        """Deliver ``item`` to each queue.

        Fast path: ``put_nowait`` for queues with capacity (no Task, no
        timer handle, no coroutine yield). Slow path: ``_put_with_timeout``
        only for queues that were full, and a single direct ``await`` when
        exactly one queue needs the slow path (skips ``asyncio.gather``).
        """
        slow_queues = None
        for q in queues:
            try:
                q.put_nowait(item)
            except asyncio.QueueFull:
                if slow_queues is None:
                    slow_queues = [q]
                else:
                    slow_queues.append(q)
        if slow_queues is None:
            return
        if len(slow_queues) == 1:
            await self._put_with_timeout(slow_queues[0], item)
        else:
            await asyncio.gather(
                *(self._put_with_timeout(q, item) for q in slow_queues)
            )

    def _ensure_recv_loop(self) -> None:
        if self._recv_task is None or self._recv_task.done():
            self._recv_task = asyncio.create_task(self._recv_loop())

    def _stop_recv_loop(self) -> None:
        if self._recv_task and not self._recv_task.done():
            self._recv_task.cancel()

    def _put_error_nowait(self, queue: asyncio.Queue, error: _StreamError) -> None:
        while True:
            try:
                queue.put_nowait(error)
                break
            except asyncio.QueueFull:
                try:
                    queue.get_nowait()
                except asyncio.QueueEmpty:
                    pass

    async def _recv_loop(self) -> None:
        try:
            while True:
                response = await self._stream.recv()
                if response == grpc.aio.EOF:
                    await self._put_to_queues(self._get_unique_queues(), _StreamEnd())
                    return

                if response.object_data_ranges:
                    queues_to_notify: Set[asyncio.Queue] = set()
                    for data_range in response.object_data_ranges:
                        read_id = data_range.read_range.read_id
                        queue = self._queues.get(read_id)
                        if queue:
                            queues_to_notify.add(queue)
                        else:
                            logger.warning(
                                f"Received data for unregistered read_id: {read_id}"
                            )
                    await self._put_to_queues(queues_to_notify, response)
                else:
                    await self._put_to_queues(self._get_unique_queues(), response)
        except asyncio.CancelledError:
            raise
        except Exception as e:
            logger.warning(f"Stream multiplexer recv loop failed: {e}", exc_info=True)
            error = _StreamError(e, self._stream_generation)
            for queue in self._get_unique_queues():
                self._put_error_nowait(queue, error)

    async def send(self, request: _storage_v2.BidiReadObjectRequest) -> int:
        self._ensure_recv_loop()
        await self._stream.send(request)
        return self._stream_generation

    async def reopen_stream(
        self,
        broken_generation: int,
        stream_factory: Callable[[], Awaitable[_AsyncReadObjectStream]],
    ) -> None:
        async with self._reopen_lock:
            if self._stream_generation != broken_generation:
                return
            self._stop_recv_loop()
            if self._recv_task:
                try:
                    await self._recv_task
                except (asyncio.CancelledError, Exception):
                    pass
            error = _StreamError(Exception("Stream reopening"), self._stream_generation)
            for queue in self._get_unique_queues():
                self._put_error_nowait(queue, error)
            try:
                await self._stream.close()
            except Exception:
                pass
            self._stream = await stream_factory()
            self._stream_generation += 1
            self._ensure_recv_loop()

    async def close(self) -> None:
        self._stop_recv_loop()
        if self._recv_task:
            try:
                await self._recv_task
            except (asyncio.CancelledError, Exception):
                pass
        error = _StreamError(Exception("Multiplexer closed"), self._stream_generation)
        for queue in self._get_unique_queues():
            self._put_error_nowait(queue, error)
