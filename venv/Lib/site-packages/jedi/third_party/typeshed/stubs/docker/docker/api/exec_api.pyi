from _io import _BufferedReaderStream
from _typeshed import Incomplete
from socket import SocketIO
from typing import Literal, overload

from docker.transport.sshconn import SSHSocket
from docker.types.daemon import CancellableStream

class ExecApiMixin:
    def exec_create(
        self,
        container,
        cmd,
        stdout: bool = True,
        stderr: bool = True,
        stdin: bool = False,
        tty: bool = False,
        privileged: bool = False,
        user: str = "",
        environment: dict[str, str] | list[str] | None = None,
        workdir: str | None = None,
        detach_keys: str | None = None,
    ) -> dict[str, Incomplete]: ...
    def exec_inspect(self, exec_id: str) -> dict[str, Incomplete]: ...
    def exec_resize(self, exec_id: str, height: int | None = None, width: int | None = None) -> None: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[True],
        tty: bool = False,
        stream: bool = False,
        socket: bool = False,
        demux: bool = False,
    ) -> bytes: ...
    @overload
    def exec_start(
        self, exec_id: str, detach: Literal[False], tty: bool, stream: bool, socket: Literal[True], demux: bool = False
    ) -> SocketIO | _BufferedReaderStream | SSHSocket: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False] = False,
        tty: bool = False,
        stream: bool = False,
        *,
        socket: Literal[True],
        demux: bool = False,
    ) -> SocketIO | _BufferedReaderStream | SSHSocket: ...
    @overload
    def exec_start(
        self, exec_id: str, detach: Literal[False], tty: bool, stream: Literal[True], socket: Literal[False], demux: Literal[True]
    ) -> CancellableStream[tuple[bytes | None, bytes | None]]: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False] = False,
        tty: bool = False,
        socket: Literal[False] = False,
        *,
        stream: Literal[True],
        demux: Literal[True],
    ) -> CancellableStream[tuple[bytes | None, bytes | None]]: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False],
        tty: bool,
        stream: Literal[True],
        socket: Literal[False],
        demux: Literal[False],
    ) -> CancellableStream[bytes]: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False] = False,
        tty: bool = False,
        *,
        stream: Literal[True],
        socket: Literal[False] = False,
        demux: Literal[False] = False,
    ) -> CancellableStream[bytes]: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False],
        tty: bool,
        stream: Literal[False],
        socket: Literal[False],
        demux: Literal[True],
    ) -> tuple[bytes | None, bytes | None]: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False] = False,
        tty: bool = False,
        stream: Literal[False] = False,
        socket: Literal[False] = False,
        *,
        demux: Literal[True],
    ) -> tuple[bytes | None, bytes | None]: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: Literal[False] = False,
        tty: bool = False,
        stream: Literal[False] = False,
        socket: Literal[False] = False,
        demux: Literal[False] = False,
    ) -> bytes: ...
    @overload
    def exec_start(
        self,
        exec_id: str,
        detach: bool = False,
        tty: bool = False,
        stream: bool = False,
        socket: bool = False,
        demux: bool = False,
    ) -> (
        str
        | SocketIO
        | _BufferedReaderStream
        | SSHSocket
        | CancellableStream[bytes]
        | CancellableStream[tuple[bytes | None, bytes | None]]
        | tuple[bytes | None, bytes | None]
        | bytes
    ): ...
