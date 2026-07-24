import asyncio
from typing import TYPE_CHECKING, Any, cast

from .client_exceptions import ClientConnectionResetError
from .helpers import set_exception
from .tcp_helpers import tcp_nodelay

if TYPE_CHECKING:
    from .http_parser import HttpParser

# Raised by transport.pause_reading()/resume_reading() when the transport
# does not support flow control; safe to ignore.
# NOTE: Catch these with a plain try/except/pass, never contextlib.suppress():
# pause/resume run on the hot read path and suppress() is ~6x slower than
# try/except here (it builds a context manager and unpacks this tuple per call).
PAUSE_RESUME_READING_ERRORS = (AttributeError, NotImplementedError, RuntimeError)


class BaseProtocol(asyncio.Protocol):
    __slots__ = (
        "_loop",
        "_paused",
        "_parser",
        "_drain_waiter",
        "_connection_lost",
        "_reading_paused",
        "_upgraded",
        "transport",
    )

    def __init__(
        self, loop: asyncio.AbstractEventLoop, parser: "HttpParser[Any] | None" = None
    ) -> None:
        self._loop: asyncio.AbstractEventLoop = loop
        self._paused = False
        self._drain_waiter: asyncio.Future[None] | None = None
        self._reading_paused = False
        self._parser = parser
        self._upgraded = False

        self.transport: asyncio.Transport | None = None

    @property
    def connected(self) -> bool:
        """Return True if the connection is open."""
        return self.transport is not None

    @property
    def writing_paused(self) -> bool:
        return self._paused

    def pause_writing(self) -> None:
        assert not self._paused
        self._paused = True

    def resume_writing(self) -> None:
        assert self._paused
        self._paused = False

        waiter = self._drain_waiter
        if waiter is not None:
            self._drain_waiter = None
            if not waiter.done():
                waiter.set_result(None)

    def pause_reading(self) -> None:
        self._reading_paused = True
        # Parser shouldn't be paused on websockets.
        if not self._upgraded:
            assert self._parser is not None
            self._parser.pause_reading()
        if self.transport is not None:
            try:
                self.transport.pause_reading()
            except PAUSE_RESUME_READING_ERRORS:
                # Transport lacks flow control; nothing to pause. Intentionally
                # ignored (see PAUSE_RESUME_READING_ERRORS; do not use suppress).
                pass

    def _reading_paused_for_msg_queue(self) -> bool:
        """Keep the transport paused for protocol-specific reasons (overridden)."""
        return False

    def resume_reading(self, resume_parser: bool = True) -> None:
        self._reading_paused = False

        # This will resume parsing any unprocessed data from the last pause.
        if not self._upgraded and resume_parser:
            self.data_received(b"")

        # Reading may have been paused again in the above call if there was a lot of
        # compressed data still pending.
        if (
            not self._reading_paused
            and not self._reading_paused_for_msg_queue()
            and self.transport is not None
        ):
            try:
                self.transport.resume_reading()
            except PAUSE_RESUME_READING_ERRORS:
                # Transport lacks flow control; nothing to resume. Intentionally
                # ignored (see PAUSE_RESUME_READING_ERRORS; do not use suppress).
                pass
            self._reading_paused = False

    def connection_made(self, transport: asyncio.BaseTransport) -> None:
        tr = cast(asyncio.Transport, transport)
        tcp_nodelay(tr, True)
        self.transport = tr

    def connection_lost(self, exc: BaseException | None) -> None:
        # Wake up the writer if currently paused.
        self.transport = None
        if not self._paused:
            return
        waiter = self._drain_waiter
        if waiter is None:
            return
        self._drain_waiter = None
        if waiter.done():
            return
        if exc is None:
            waiter.set_result(None)
        else:
            set_exception(
                waiter,
                ConnectionError("Connection lost"),
                exc,
            )

    async def _drain_helper(self) -> None:
        if self.transport is None:
            raise ClientConnectionResetError("Connection lost")
        if not self._paused:
            return
        waiter = self._drain_waiter
        if waiter is None:
            waiter = self._loop.create_future()
            self._drain_waiter = waiter
        await waiter
