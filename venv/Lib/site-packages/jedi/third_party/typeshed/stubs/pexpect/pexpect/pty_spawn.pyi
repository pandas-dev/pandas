from _typeshed import FileDescriptorOrPath
from collections.abc import Callable, Mapping
from typing import AnyStr

from .spawnbase import SpawnBase, _Logfile

PY3: bool

class spawn(SpawnBase[AnyStr]):
    use_native_pty_fork: bool
    STDIN_FILENO: int
    STDOUT_FILENO: int
    STDERR_FILENO: int
    str_last_chars: int
    cwd: FileDescriptorOrPath | None
    env: Mapping[str, str] | None
    echo: bool
    ignore_sighup: bool
    command: str
    args: list[str]
    name: str
    use_poll: bool
    def __init__(
        self,
        command: str,
        args: list[str] = [],
        timeout: float | None = 30,
        maxread: int = 2000,
        searchwindowsize: int | None = None,
        logfile: _Logfile | None = None,
        cwd: FileDescriptorOrPath | None = None,
        env: Mapping[str, str] | None = None,
        ignore_sighup: bool = False,
        echo: bool = True,
        preexec_fn: Callable[[], None] | None = None,
        encoding: str | None = None,
        codec_errors: str = "strict",
        dimensions: tuple[int, int] | None = None,
        use_poll: bool = False,
    ) -> None: ...
    child_fd: int
    closed: bool
    def close(self, force: bool = True) -> None: ...
    def isatty(self) -> bool: ...
    def waitnoecho(self, timeout: float | None = -1) -> None: ...
    def getecho(self) -> bool: ...
    def setecho(self, state: bool) -> None: ...
    def read_nonblocking(self, size: int = 1, timeout: float | None = -1) -> AnyStr: ...
    def write(self, s: str | bytes) -> None: ...
    def writelines(self, sequence: list[str | bytes]) -> None: ...
    def send(self, s: str | bytes) -> int: ...
    def sendline(self, s: str | bytes = "") -> int: ...
    def sendcontrol(self, char: str) -> int: ...
    def sendeof(self) -> None: ...
    def sendintr(self) -> None: ...
    @property
    def flag_eof(self) -> bool: ...
    @flag_eof.setter
    def flag_eof(self, value: bool) -> None: ...
    def eof(self) -> bool: ...
    def terminate(self, force: bool = False) -> bool: ...
    status: int | None
    exitstatus: int | None
    signalstatus: int | None
    terminated: bool
    def wait(self) -> int: ...
    def isalive(self) -> bool: ...
    def kill(self, sig: int) -> None: ...
    def getwinsize(self) -> tuple[int, int]: ...
    def setwinsize(self, rows, cols) -> None: ...
    def interact(
        self,
        escape_character="\x1d",
        input_filter: Callable[[AnyStr], AnyStr] | None = None,
        output_filter: Callable[[AnyStr], AnyStr] | None = None,
    ) -> None: ...

def spawnu(
    command: str,
    args: list[str] = [],
    timeout: float | None = 30,
    maxread: int = 2000,
    searchwindowsize: int | None = None,
    logfile: _Logfile | None = None,
    cwd: FileDescriptorOrPath | None = None,
    env: Mapping[str, str] | None = None,
    ignore_sighup: bool = False,
    echo: bool = True,
    preexec_fn: Callable[[], None] | None = None,
    encoding: str | None = "utf-8",
    codec_errors: str = "strict",
    dimensions: tuple[int, int] | None = None,
    use_poll: bool = False,
) -> spawn[str]: ...
