import subprocess
import sys
from _typeshed import AnyPath
from asyncio import events, protocols, streams, transports
from typing import IO, Any, Callable, Optional, Tuple, Union
from typing_extensions import Literal

if sys.version_info >= (3, 8):
    _ExecArg = AnyPath
else:
    _ExecArg = Union[str, bytes]

PIPE: int
STDOUT: int
DEVNULL: int

class SubprocessStreamProtocol(streams.FlowControlMixin, protocols.SubprocessProtocol):
    stdin: Optional[streams.StreamWriter]
    stdout: Optional[streams.StreamReader]
    stderr: Optional[streams.StreamReader]
    def __init__(self, limit: int, loop: events.AbstractEventLoop) -> None: ...
    def connection_made(self, transport: transports.BaseTransport) -> None: ...
    def pipe_data_received(self, fd: int, data: Union[bytes, str]) -> None: ...
    def pipe_connection_lost(self, fd: int, exc: Optional[Exception]) -> None: ...
    def process_exited(self) -> None: ...

class Process:
    stdin: Optional[streams.StreamWriter]
    stdout: Optional[streams.StreamReader]
    stderr: Optional[streams.StreamReader]
    pid: int
    def __init__(
        self, transport: transports.BaseTransport, protocol: protocols.BaseProtocol, loop: events.AbstractEventLoop
    ) -> None: ...
    @property
    def returncode(self) -> Optional[int]: ...
    async def wait(self) -> int: ...
    def send_signal(self, signal: int) -> None: ...
    def terminate(self) -> None: ...
    def kill(self) -> None: ...
    async def communicate(self, input: Optional[bytes] = ...) -> Tuple[bytes, bytes]: ...

if sys.version_info >= (3, 10):
    async def create_subprocess_shell(
        cmd: Union[str, bytes],
        stdin: Union[int, IO[Any], None] = ...,
        stdout: Union[int, IO[Any], None] = ...,
        stderr: Union[int, IO[Any], None] = ...,
        limit: int = ...,
        *,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        text: Literal[False, None] = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        executable: Optional[AnyPath] = ...,
        preexec_fn: Optional[Callable[[], Any]] = ...,
        close_fds: bool = ...,
        cwd: Optional[AnyPath] = ...,
        env: Optional[subprocess._ENV] = ...,
        startupinfo: Optional[Any] = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...
    async def create_subprocess_exec(
        program: _ExecArg,
        *args: _ExecArg,
        stdin: Union[int, IO[Any], None] = ...,
        stdout: Union[int, IO[Any], None] = ...,
        stderr: Union[int, IO[Any], None] = ...,
        limit: int = ...,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        text: Optional[bool] = ...,
        executable: Optional[AnyPath] = ...,
        preexec_fn: Optional[Callable[[], Any]] = ...,
        close_fds: bool = ...,
        cwd: Optional[AnyPath] = ...,
        env: Optional[subprocess._ENV] = ...,
        startupinfo: Optional[Any] = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...

else:
    async def create_subprocess_shell(
        cmd: Union[str, bytes],
        stdin: Union[int, IO[Any], None] = ...,
        stdout: Union[int, IO[Any], None] = ...,
        stderr: Union[int, IO[Any], None] = ...,
        loop: Optional[events.AbstractEventLoop] = ...,
        limit: int = ...,
        *,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        text: Literal[False, None] = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        executable: Optional[AnyPath] = ...,
        preexec_fn: Optional[Callable[[], Any]] = ...,
        close_fds: bool = ...,
        cwd: Optional[AnyPath] = ...,
        env: Optional[subprocess._ENV] = ...,
        startupinfo: Optional[Any] = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...
    async def create_subprocess_exec(
        program: _ExecArg,
        *args: _ExecArg,
        stdin: Union[int, IO[Any], None] = ...,
        stdout: Union[int, IO[Any], None] = ...,
        stderr: Union[int, IO[Any], None] = ...,
        loop: Optional[events.AbstractEventLoop] = ...,
        limit: int = ...,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        text: Optional[bool] = ...,
        executable: Optional[AnyPath] = ...,
        preexec_fn: Optional[Callable[[], Any]] = ...,
        close_fds: bool = ...,
        cwd: Optional[AnyPath] = ...,
        env: Optional[subprocess._ENV] = ...,
        startupinfo: Optional[Any] = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...
