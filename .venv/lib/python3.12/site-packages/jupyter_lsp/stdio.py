""" Language Server stdio-mode readers

Parts of this code are derived from:

> https://github.com/palantir/python-jsonrpc-server/blob/0.2.0/pyls_jsonrpc/streams.py#L83   # noqa
> https://github.com/palantir/python-jsonrpc-server/blob/45ed1931e4b2e5100cc61b3992c16d6f68af2e80/pyls_jsonrpc/streams.py  # noqa
> > MIT License   https://github.com/palantir/python-jsonrpc-server/blob/0.2.0/LICENSE
> > Copyright 2018 Palantir Technologies, Inc.
"""

# pylint: disable=broad-except
import asyncio
import io
import os
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional, Text

from tornado.concurrent import run_on_executor
from tornado.gen import convert_yielded
from tornado.httputil import HTTPHeaders
from tornado.ioloop import IOLoop
from tornado.queues import Queue
from traitlets import Float, Instance, default
from traitlets.config import LoggingConfigurable

from .non_blocking import make_non_blocking


class LspStdIoBase(LoggingConfigurable):
    """Non-blocking, queued base for communicating with stdio Language Servers"""

    executor = None

    stream = Instance(  # type:ignore[assignment]
        io.RawIOBase, help="the stream to read/write"
    )  # type: io.RawIOBase
    queue = Instance(Queue, help="queue to get/put")

    def __repr__(self):  # pragma: no cover
        return "<{}(parent={})>".format(self.__class__.__name__, self.parent)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.log.debug("%s initialized", self)
        self.executor = ThreadPoolExecutor(max_workers=1)

    def close(self):
        self.stream.close()
        self.log.debug("%s closed", self)


class LspStdIoReader(LspStdIoBase):
    """Language Server stdio Reader

    Because non-blocking (but still synchronous) IO is used, rudimentary
    exponential backoff is used.
    """

    max_wait = Float(help="maximum time to wait on idle stream").tag(config=True)
    min_wait = Float(0.05, help="minimum time to wait on idle stream").tag(config=True)
    next_wait = Float(0.05, help="next time to wait on idle stream").tag(config=True)

    @default("max_wait")
    def _default_max_wait(self):
        return 0.1 if os.name == "nt" else self.min_wait * 2

    async def sleep(self):
        """Simple exponential backoff for sleeping"""
        if self.stream.closed:  # pragma: no cover
            return
        self.next_wait = min(self.next_wait * 2, self.max_wait)
        try:
            await asyncio.sleep(self.next_wait)
        except Exception:  # pragma: no cover
            pass

    def wake(self):
        """Reset the wait time"""
        self.wait = self.min_wait

    async def read(self) -> None:
        """Read from a Language Server until it is closed"""
        make_non_blocking(self.stream)

        while not self.stream.closed:
            message = None
            try:
                message = await self.read_one()

                if not message:
                    await self.sleep()
                    continue
                else:
                    self.wake()

                IOLoop.current().add_callback(self.queue.put_nowait, message)
            except Exception as e:  # pragma: no cover
                self.log.exception(
                    "%s couldn't enqueue message: %s (%s)", self, message, e
                )
                await self.sleep()

    async def _read_content(
        self, length: int, max_parts=1000, max_empties=200
    ) -> Optional[bytes]:
        """Read the full length of the message unless exceeding max_parts or
           max_empties empty reads occur.

        See https://github.com/jupyter-lsp/jupyterlab-lsp/issues/450

        Crucial docs or read():
            "If the argument is positive, and the underlying raw
             stream is not interactive, multiple raw reads may be issued
             to satisfy the byte count (unless EOF is reached first)"

        Args:
           - length: the content length
           - max_parts: prevent absurdly long messages (1000 parts is several MBs):
             1 part is usually sufficient but not enough for some long
             messages 2 or 3 parts are often needed.
        """
        raw = None
        raw_parts: List[bytes] = []
        received_size = 0
        while received_size < length and len(raw_parts) < max_parts and max_empties > 0:
            part = None
            try:
                part = self.stream.read(length - received_size)
            except OSError:  # pragma: no cover
                pass
            if part is None:
                max_empties -= 1
                await self.sleep()
                continue
            received_size += len(part)
            raw_parts.append(part)

        if raw_parts:
            raw = b"".join(raw_parts)
            if len(raw) != length:  # pragma: no cover
                self.log.warning(
                    f"Readout and content-length mismatch: {len(raw)} vs {length};"
                    f"remaining empties: {max_empties}; remaining parts: {max_parts}"
                )

        return raw

    async def read_one(self) -> Text:
        """Read a single message"""
        message = ""
        headers = HTTPHeaders()

        line = await convert_yielded(self._readline())

        if line:
            while line and line.strip():
                headers.parse_line(line)
                line = await convert_yielded(self._readline())

            content_length = int(headers.get("content-length", "0"))

            if content_length:
                raw = await self._read_content(length=content_length)
                if raw is not None:
                    message = raw.decode("utf-8").strip()
                else:  # pragma: no cover
                    self.log.warning(
                        "%s failed to read message of length %s",
                        self,
                        content_length,
                    )

        return message

    @run_on_executor
    def _readline(self) -> Text:
        """Read a line (or immediately return None)"""
        try:
            return self.stream.readline().decode("utf-8").strip()
        except OSError:  # pragma: no cover
            return ""


class LspStdIoWriter(LspStdIoBase):
    """Language Server stdio Writer"""

    async def write(self) -> None:
        """Write to a Language Server until it closes"""
        while not self.stream.closed:
            message = await self.queue.get()
            try:
                body = message.encode("utf-8")
                response = "Content-Length: {}\r\n\r\n{}".format(len(body), message)
                await convert_yielded(self._write_one(response.encode("utf-8")))
            except Exception:  # pragma: no cover
                self.log.exception("%s couldn't write message: %s", self, response)
            finally:
                self.queue.task_done()

    @run_on_executor
    def _write_one(self, message) -> None:
        self.stream.write(message)
        self.stream.flush()
