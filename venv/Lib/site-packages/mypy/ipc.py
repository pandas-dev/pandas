"""Cross-platform abstractions for inter-process communication

On Unix, this uses AF_UNIX sockets.
On Windows, this uses NamedPipes.
"""

from __future__ import annotations

import json
import os
import shutil
import struct
import sys
import tempfile
from abc import abstractmethod
from collections.abc import Callable, Sequence
from select import select
from types import TracebackType
from typing import Final
from typing_extensions import Self

from librt.base64 import urlsafe_b64encode
from librt.internal import ReadBuffer, WriteBuffer

if sys.platform == "win32":
    # This may be private, but it is needed for IPC on Windows, and is basically stable
    import _winapi
    import ctypes

    _IPCHandle = int

    kernel32 = ctypes.windll.kernel32
    DisconnectNamedPipe: Callable[[_IPCHandle], int] = kernel32.DisconnectNamedPipe
    FlushFileBuffers: Callable[[_IPCHandle], int] = kernel32.FlushFileBuffers
else:
    import socket

    _IPCHandle = socket.socket

# Size of the message packed as !L, i.e. 4 bytes in network order (big-endian).
HEADER_SIZE: Final = 4

# This is Linux default socket buffer size (for 64 bit), so we will not
# introduce an additional obstacle when exchanging a large IPC message.
MAX_READ: Final = 212992


# TODO: we should make sure consistent exceptions are raised on different platforms.
# Currently we raise either IPCException or OSError for equivalent conditions.
class IPCException(Exception):
    """Exception for IPC issues."""


class IPCBase:
    """Base class for communication between the dmypy client and server.

    This contains logic shared between the client and server, such as reading
    and writing.
    We want to be able to send multiple "messages" over a single connection and
    to be able to separate the messages. We do this by prefixing each message
    with its size in a fixed format.
    """

    connection: _IPCHandle

    def __init__(self, name: str, timeout: float | None) -> None:
        self.name = name
        self.timeout = timeout
        self.message_size: int | None = None
        self.buffer = bytearray()

    def frame_from_buffer(self) -> bytes | None:
        """Return a full frame from the bytes we have in the buffer."""
        size = len(self.buffer)
        if size < HEADER_SIZE:
            return None
        if self.message_size is None:
            self.message_size = struct.unpack("!L", self.buffer[:HEADER_SIZE])[0]
        if size < self.message_size + HEADER_SIZE:
            return None
        # We have a full frame, avoid extra copy in case we get a large frame.
        bdata = memoryview(self.buffer)[HEADER_SIZE : HEADER_SIZE + self.message_size]
        self.buffer = self.buffer[HEADER_SIZE + self.message_size :]
        self.message_size = None
        return bytes(bdata)

    def read(self, size: int = MAX_READ) -> str:
        return self.read_bytes(size).decode("utf-8")

    def read_bytes(self, size: int = MAX_READ) -> bytes:
        """Read bytes from an IPC connection until we have a full frame."""
        if sys.platform == "win32":
            while True:
                # Check if we already have a message in the buffer before
                # receiving any more data from the socket.
                bdata = self.frame_from_buffer()
                if bdata is not None:
                    break

                # Receive more data into the buffer.
                ov, err = _winapi.ReadFile(self.connection, size, overlapped=True)
                try:
                    if err == _winapi.ERROR_IO_PENDING:
                        timeout = int(self.timeout * 1000) if self.timeout else _winapi.INFINITE
                        res = _winapi.WaitForSingleObject(ov.event, timeout)
                        if res != _winapi.WAIT_OBJECT_0:
                            raise IPCException(f"Bad result from I/O wait: {res}")
                except BaseException:
                    ov.cancel()
                    raise
                _, err = ov.GetOverlappedResult(True)
                more = ov.getbuffer()
                if more:
                    self.buffer.extend(more)
                    bdata = self.frame_from_buffer()
                    if bdata is not None:
                        break
                if err == 0:
                    # we are done!
                    break
                elif err == _winapi.ERROR_MORE_DATA:
                    # read again
                    continue
                elif err == _winapi.ERROR_OPERATION_ABORTED:
                    raise IPCException("ReadFile operation aborted.")
        else:
            while True:
                # Check if we already have a message in the buffer before
                # receiving any more data from the socket.
                bdata = self.frame_from_buffer()
                if bdata is not None:
                    break

                # Receive more data into the buffer.
                more = self.connection.recv(size)
                if not more:
                    # Connection closed
                    break
                self.buffer.extend(more)

        if not bdata:
            # Socket was empty, and we didn't get any frame.
            # This should only happen if the socket was closed.
            return b""
        return bdata

    def write(self, data: str) -> None:
        self.write_bytes(data.encode("utf-8"))

    def write_bytes(self, data: bytes) -> None:
        """Write to an IPC connection."""

        # Frame the data by adding fixed size header.
        encoded_data = struct.pack("!L", len(data)) + data

        if sys.platform == "win32":
            try:
                ov, err = _winapi.WriteFile(self.connection, encoded_data, overlapped=True)
                try:
                    if err == _winapi.ERROR_IO_PENDING:
                        timeout = int(self.timeout * 1000) if self.timeout else _winapi.INFINITE
                        res = _winapi.WaitForSingleObject(ov.event, timeout)
                        if res != _winapi.WAIT_OBJECT_0:
                            raise IPCException(f"Bad result from I/O wait: {res}")
                    elif err != 0:
                        raise IPCException(f"Failed writing to pipe with error: {err}")
                except BaseException:
                    ov.cancel()
                    raise
                bytes_written, err = ov.GetOverlappedResult(True)
                assert err == 0, err
                assert bytes_written == len(encoded_data)
            except OSError as e:
                raise IPCException(f"Failed to write with error: {e.winerror}") from e
        else:
            self.connection.sendall(encoded_data)

    def close(self) -> None:
        if sys.platform == "win32":
            if self.connection != _winapi.NULL:
                _winapi.CloseHandle(self.connection)
        else:
            self.connection.close()


class IPCClient(IPCBase):
    """The client side of an IPC connection."""

    def __init__(self, name: str, timeout: float | None) -> None:
        super().__init__(name, timeout)
        if sys.platform == "win32":
            timeout = int(self.timeout * 1000) if self.timeout else _winapi.NMPWAIT_WAIT_FOREVER
            try:
                _winapi.WaitNamedPipe(self.name, timeout)
            except FileNotFoundError as e:
                raise IPCException(f"The NamedPipe at {self.name} was not found.") from e
            except OSError as e:
                if e.winerror == _winapi.ERROR_SEM_TIMEOUT:
                    raise IPCException("Timed out waiting for connection.") from e
                else:
                    raise
            try:
                self.connection = _winapi.CreateFile(
                    self.name,
                    _winapi.GENERIC_READ | _winapi.GENERIC_WRITE,
                    0,
                    _winapi.NULL,
                    _winapi.OPEN_EXISTING,
                    _winapi.FILE_FLAG_OVERLAPPED,
                    _winapi.NULL,
                )
            except OSError as e:
                if e.winerror == _winapi.ERROR_PIPE_BUSY:
                    raise IPCException("The connection is busy.") from e
                else:
                    raise
            _winapi.SetNamedPipeHandleState(
                self.connection, _winapi.PIPE_READMODE_MESSAGE, None, None
            )
        else:
            self.connection = socket.socket(socket.AF_UNIX)
            # This is already default on Linux, we set same buffer size
            # for macOS vs Linux consistency to simplify reasoning.
            self.connection.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, MAX_READ)
            self.connection.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, MAX_READ)
            self.connection.settimeout(timeout)
            self.connection.connect(name)

    def __enter__(self) -> IPCClient:
        return self

    def __exit__(
        self,
        exc_ty: type[BaseException] | None = None,
        exc_val: BaseException | None = None,
        exc_tb: TracebackType | None = None,
    ) -> None:
        self.close()


class IPCServer(IPCBase):
    BUFFER_SIZE: Final = 2**16

    def __init__(self, name: str, timeout: float | None = None) -> None:
        if sys.platform == "win32":
            name = r"\\.\pipe\{}-{}.pipe".format(name, urlsafe_b64encode(os.urandom(6)).decode())
        else:
            name = f"{name}.sock"
        super().__init__(name, timeout)
        if sys.platform == "win32":
            self.connection = _winapi.CreateNamedPipe(
                self.name,
                _winapi.PIPE_ACCESS_DUPLEX
                | _winapi.FILE_FLAG_FIRST_PIPE_INSTANCE
                | _winapi.FILE_FLAG_OVERLAPPED,
                _winapi.PIPE_READMODE_MESSAGE
                | _winapi.PIPE_TYPE_MESSAGE
                | _winapi.PIPE_WAIT
                | 0x8,  # PIPE_REJECT_REMOTE_CLIENTS
                1,  # one instance
                self.BUFFER_SIZE,
                self.BUFFER_SIZE,
                _winapi.NMPWAIT_WAIT_FOREVER,
                0,  # Use default security descriptor
            )
            if self.connection == -1:  # INVALID_HANDLE_VALUE
                err = _winapi.GetLastError()
                raise IPCException(f"Invalid handle to pipe: {err}")
        else:
            self.sock_directory = tempfile.mkdtemp()
            sockfile = os.path.join(self.sock_directory, self.name)
            self.sock = socket.socket(socket.AF_UNIX)
            self.sock.bind(sockfile)
            self.sock.listen(1)
            if timeout is not None:
                self.sock.settimeout(timeout)

    def __enter__(self) -> IPCServer:
        if sys.platform == "win32":
            # NOTE: It is theoretically possible that this will hang forever if the
            # client never connects, though this can be "solved" by killing the server
            try:
                ov = _winapi.ConnectNamedPipe(self.connection, overlapped=True)
            except OSError as e:
                # Don't raise if the client already exists, or the client already connected
                if e.winerror not in (_winapi.ERROR_PIPE_CONNECTED, _winapi.ERROR_NO_DATA):
                    raise
            else:
                try:
                    timeout = int(self.timeout * 1000) if self.timeout else _winapi.INFINITE
                    res = _winapi.WaitForSingleObject(ov.event, timeout)
                    assert res == _winapi.WAIT_OBJECT_0
                except BaseException:
                    ov.cancel()
                    _winapi.CloseHandle(self.connection)
                    raise
                _, err = ov.GetOverlappedResult(True)
                assert err == 0
        else:
            try:
                self.connection, _ = self.sock.accept()
                # This is already default on Linux, we set same buffer size
                # for macOS vs Linux consistency to simplify reasoning.
                self.connection.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, MAX_READ)
                self.connection.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, MAX_READ)
            except TimeoutError as e:
                raise IPCException("The socket timed out") from e
        return self

    def __exit__(
        self,
        exc_ty: type[BaseException] | None = None,
        exc_val: BaseException | None = None,
        exc_tb: TracebackType | None = None,
    ) -> None:
        if sys.platform == "win32":
            try:
                # Wait for the client to finish reading the last write before disconnecting
                if not FlushFileBuffers(self.connection):
                    raise IPCException(
                        "Failed to flush NamedPipe buffer, maybe the client hung up?"
                    )
            finally:
                DisconnectNamedPipe(self.connection)
        else:
            self.close()

    def cleanup(self) -> None:
        if sys.platform == "win32":
            self.close()
        else:
            shutil.rmtree(self.sock_directory)

    @property
    def connection_name(self) -> str:
        if sys.platform == "win32":
            return self.name
        elif sys.platform == "gnu0":
            # GNU/Hurd returns empty string from getsockname()
            # for AF_UNIX sockets
            return os.path.join(self.sock_directory, self.name)
        else:
            name = self.sock.getsockname()
            assert isinstance(name, str)
            return name


class BadStatus(Exception):
    """Exception raised when there is something wrong with the status file.

    For example:
    - No status file found
    - Status file malformed
    - Process whose pid is in the status file does not exist
    """


def read_status(status_file: str) -> dict[str, object]:
    """Read status file.

    Raise BadStatus if the status file doesn't exist or contains
    invalid JSON or the JSON is not a dict.
    """
    if not os.path.isfile(status_file):
        raise BadStatus("No status file found")
    with open(status_file) as f:
        try:
            data = json.load(f)
        except Exception as e:
            raise BadStatus(f"Malformed status file: {str(e)}") from e
    if not isinstance(data, dict):
        raise BadStatus(f"Invalid status file (not a dict): {data}")
    return data


def ready_to_read(conns: Sequence[IPCBase], timeout: float | None = None) -> list[int]:
    """Wait until some connections are readable.

    Return index of each readable connection in the original list.
    """
    unread_messages = [i for i, conn in enumerate(conns) if conn.buffer]
    if unread_messages:
        # If we already have unread messages in the buffer, return those first.
        return unread_messages
    if sys.platform == "win32":
        # Windows doesn't support select() on named pipes. Instead, start an overlapped
        # ReadFile on each pipe (which internally creates an event via CreateEventW),
        # then WaitForMultipleObjects on those events for efficient OS-level waiting.
        # Any data consumed by the probe reads is stored into each connection's buffer
        # so the subsequent read_bytes() call will find it via frame_from_buffer().
        WAIT_FAILED = 0xFFFFFFFF
        pending: list[tuple[int, _winapi.Overlapped]] = []
        events: list[int] = []
        ready: list[int] = []

        for i, conn in enumerate(conns):
            try:
                ov, err = _winapi.ReadFile(conn.connection, 1, overlapped=True)
            except OSError:
                # Broken/closed pipe. Mimic Linux behavior here, caller will get
                # the exception when trying to read from this socket.
                ready.append(i)
                continue
            if err == _winapi.ERROR_IO_PENDING:
                events.append(ov.event)
                pending.append((i, ov))
            else:
                # Data was immediately available (err == 0 or ERROR_MORE_DATA)
                _, err = ov.GetOverlappedResult(True)
                data = ov.getbuffer()
                if data:
                    conn.buffer.extend(data)
                ready.append(i)

        # Wait only if nothing is immediately ready and there are pending operations
        if not ready and events:
            timeout_ms = int(timeout * 1000) if timeout is not None else _winapi.INFINITE
            res = _winapi.WaitForMultipleObjects(events, False, timeout_ms)
            if res == WAIT_FAILED:
                for _, ov in pending:
                    ov.cancel()
                raise IPCException(f"Failed to wait for connections: {_winapi.GetLastError()}")

        # Cancel all pending operations. CancelIoEx is asynchronous, so an
        # operation may have completed before the cancel took effect. We then
        # wait for all operations to finalize and check each result: completed
        # reads get their data saved and are marked ready; cancelled ones are
        # simply skipped. This avoids a race between checking if an operation
        # is signaled and cancelling it.
        for _, ov in pending:
            ov.cancel()
        for i, ov in pending:
            try:
                _, err = ov.GetOverlappedResult(True)
            except OSError as e:
                err = e.winerror
                # Cancellation is expected here; broken/disconnected pipes should be
                # surfaced as readable so the follow-up receive observes EOF/closure.
                if err not in (
                    _winapi.ERROR_OPERATION_ABORTED,
                    _winapi.ERROR_BROKEN_PIPE,
                    _winapi.ERROR_NETNAME_DELETED,
                ):
                    # Anything else is a real IPC failure, not part of the probe race.
                    raise
            if err == _winapi.ERROR_OPERATION_ABORTED:
                # Operation was successfully cancelled -- no data consumed.
                continue
            if err in (0, _winapi.ERROR_MORE_DATA):
                data = ov.getbuffer()
                if data:
                    conns[i].buffer.extend(data)
            ready.append(i)

        return ready

    else:
        connections = [conn.connection for conn in conns]
        ready, _, _ = select(connections, [], [], timeout)
        return [connections.index(r) for r in ready]


def send(connection: IPCBase, data: IPCMessage) -> None:
    """Send data to a connection encoded and framed.

    The data must be a non-abstract IPCMessage. We assume that a single send call is a
    single frame to be sent.
    """
    buf = WriteBuffer()
    data.write(buf)
    connection.write_bytes(buf.getvalue())


def receive(connection: IPCBase) -> ReadBuffer:
    """Receive single encoded IPCMessage frame from a connection.

    Raise OSError if the data received is not valid.
    """
    bdata = connection.read_bytes()
    if not bdata:
        raise OSError("No data received")
    return ReadBuffer(bdata)


class IPCMessage:
    @classmethod
    @abstractmethod
    def read(cls, buf: ReadBuffer) -> Self:
        raise NotImplementedError

    @abstractmethod
    def write(self, buf: WriteBuffer) -> None:
        raise NotImplementedError
