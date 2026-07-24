from __future__ import annotations

__lazy_modules__ = {"atexit", "contextlib", "heapq", "itertools"}

import atexit
import contextlib
import enum
import heapq
import itertools
import sys
import time
import typing
from queue import Empty as QueueEmpty
from queue import Queue
from threading import Lock, Thread
from typing import Any

from plumbum.lib import IS_WIN32

if typing.TYPE_CHECKING:
    import subprocess
    from collections.abc import Callable, Container, Generator
    from typing import IO, Literal

    from plumbum.machines.base import PopenWithAddons


# ===================================================================================================
# utility functions
# ===================================================================================================
def _check_process(
    proc: PopenWithAddons[Any],
    retcode: int | Container[int] | None,
    timeout: float | None,
    stdout: str | bytes,
    stderr: str | bytes,
) -> tuple[int | None, str | bytes, str | bytes]:
    proc.verify(retcode, timeout, stdout, stderr)
    return proc.returncode, stdout, stderr


def _get_piped_streams(proc: PopenWithAddons[Any] | None) -> list[tuple[int, IO[Any]]]:
    """Get a list of all valid standard streams for proc that were opened with PIPE option.

    If proc was started from a Pipeline command, this function assumes it will have a
    "srcproc" member pointing to the previous command in the pipeline. That link will
    be used to traverse all started processes started from the pipeline, the list will
    include stdout/stderr streams opened as PIPE for all commands in the pipeline.
    If that was not the case, some processes could write to pipes no one reads from
    which would result in process stalling after the pipe's buffer is filled.

    Streams that were closed (because they were redirected to the input of a subsequent command)
    are not included in the result
    """
    streams = []

    def add_stream(type_: int, stream: IO[Any] | None) -> None:
        if stream is None or stream.closed:
            return
        streams.append((type_, stream))

    while proc:
        add_stream(1, proc.stderr)
        add_stream(0, proc.stdout)
        proc = getattr(proc, "srcproc", None)

    return streams


def _iter_lines_posix(
    proc: PopenWithAddons[Any],
    decode: Callable[[bytes], str],
    linesize: int,
    line_timeout: float | None = None,
) -> Generator[tuple[int, str], None, None]:
    from selectors import EVENT_READ, DefaultSelector

    streams = _get_piped_streams(proc)

    # Python 3.4+ implementation
    def selector() -> Generator[tuple[int, str], None, None]:
        sel = DefaultSelector()
        for stream_type, stream in streams:
            sel.register(stream, EVENT_READ, stream_type)
        while True:
            poll_timeout = line_timeout if line_timeout is not None else 0.1
            ready = sel.select(poll_timeout)
            if not ready and line_timeout:
                raise ProcessLineTimedOut(
                    "popen line timeout expired",
                    getattr(proc, "argv", None),
                    getattr(proc, "machine", None),
                )
            if not ready and proc.poll() is not None:
                return
            for key, _mask in ready:
                # We pass the stream to the selector, so we get a stream out
                line = key.fileobj.readline(linesize)  # type: ignore[union-attr]
                if not line:
                    with contextlib.suppress(Exception):
                        sel.unregister(key.fileobj)
                    if not sel.get_map():
                        return
                    continue
                yield key.data, decode(line)

    for ret in selector():
        yield ret
        if proc.poll() is not None:
            break
    for stream_type, stream in streams:
        for line in stream:
            yield stream_type, decode(line)


def _iter_lines_win32(
    proc: PopenWithAddons[Any],
    decode: Callable[[bytes], str],
    linesize: int,
    line_timeout: float | None = None,
) -> Generator[tuple[int, str], None, None]:
    class Piper(Thread):
        __slots__ = ("empty", "fd", "pipe")

        def __init__(self, fd: int, pipe: IO[bytes]) -> None:
            super().__init__(name=f"PlumbumPiper{fd}Thread")
            self.pipe = pipe
            self.fd = fd
            self.empty = False
            self.daemon = True
            super().start()

        def read_from_pipe(self) -> bytes:
            return self.pipe.readline(linesize)

        def run(self) -> None:
            for line in iter(self.read_from_pipe, b""):
                queue.put((self.fd, decode(line)))
            # self.pipe.close()

    if line_timeout is None:
        line_timeout = float("inf")

    queue = Queue[tuple[int, str]]()
    assert proc.stdout is not None
    assert proc.stderr is not None
    pipers = [Piper(0, proc.stdout), Piper(1, proc.stderr)]
    last_line_ts = time.time()
    empty = True
    while True:
        try:
            yield queue.get_nowait()
            last_line_ts = time.time()
            empty = False
        except QueueEmpty:
            empty = True
        if time.time() - last_line_ts > line_timeout:
            raise ProcessLineTimedOut(
                "popen line timeout expired",
                getattr(proc, "argv", None),
                getattr(proc, "machine", None),
            )
        if proc.poll() is not None:
            break
        if empty:
            time.sleep(0.1)

    for piper in pipers:
        piper.join()

    with contextlib.suppress(QueueEmpty):
        while True:
            yield queue.get_nowait()


_iter_lines = _iter_lines_win32 if IS_WIN32 else _iter_lines_posix


# ===================================================================================================
# Exceptions
# ===================================================================================================
class ProcessExecutionError(OSError):
    """Represents the failure of a process. When the exit code of a terminated process does not
    match the expected result, this exception is raised by :func:`run_proc
    <plumbum.commands.processes.run_proc>`. It contains the process' return code, stdout, and stderr, as
    well as the command line used to create the process (``argv``)
    """

    def __init__(
        self,
        argv: list[str],
        retcode: int | str | None,
        stdout: str | bytes,
        stderr: str | bytes,
        message: str | None = None,
        *,
        host: str | None = None,
    ):
        # we can't use 'super' here since OSError only keeps the first 2 args,
        # which leads to failing in loading this object from a pickle.dumps.
        # pylint: disable-next=non-parent-init-called
        Exception.__init__(self, argv, retcode, stdout, stderr)

        self.message = message
        self.host = host
        self.argv = argv
        self.retcode = retcode
        if isinstance(stdout, bytes):
            stdout = ascii(stdout)
        if isinstance(stderr, bytes):
            stderr = ascii(stderr)
        self.stdout = stdout
        self.stderr = stderr

    def __str__(self) -> str:
        # avoid an import cycle
        from plumbum.commands.base import shquote_list

        stdout = "\n              | ".join(str(self.stdout).splitlines())
        stderr = "\n              | ".join(str(self.stderr).splitlines())
        cmd = " ".join(shquote_list(self.argv))
        lines = []
        if self.message:
            lines = [self.message, "\nReturn code:  | ", str(self.retcode)]
        else:
            lines = ["Unexpected exit code: ", str(self.retcode)]
        cmd = "\n              | ".join(cmd.splitlines())
        lines += ["\nCommand line: | ", cmd]
        if self.host:
            lines += ["\nHost:         | ", self.host]
        if stdout:
            lines += ["\nStdout:       | ", stdout]
        if stderr:
            lines += ["\nStderr:       | ", stderr]
        return "".join(lines)


class ProcessTimedOut(Exception):
    """Raises by :func:`run_proc <plumbum.commands.processes.run_proc>` when a ``timeout`` has been
    specified and it has elapsed before the process terminated"""

    def __init__(self, msg: str, argv: Any):
        Exception.__init__(self, msg, argv)
        self.argv = argv


class ProcessLineTimedOut(Exception):
    """Raises by :func:`iter_lines <plumbum.commands.processes.iter_lines>` when a ``line_timeout`` has been
    specified and it has elapsed before the process yielded another line"""

    def __init__(self, msg: str, argv: list[str] | None, machine: str | None):
        Exception.__init__(self, msg, argv, machine)
        self.argv = argv
        self.machine = machine


class CommandNotFound(AttributeError):
    """Raised by :func:`local.which <plumbum.machines.local.LocalMachine.which>` and
    :func:`BaseRemoteMachine.which <plumbum.machines.remote.BaseRemoteMachine.which>` when a
    command was not found in the system's ``PATH``"""

    def __init__(self, program: object, path: object):
        super().__init__(program, path)
        self.program = program
        self.path = path


# ===================================================================================================
# Timeout thread
# ===================================================================================================
class MinHeap:
    """Deadline-ordered heap of (deadline, proc) pairs.

    An internal monotonic counter breaks deadline ties so heapq never falls
    back to comparing the (uncomparable) Popen objects, which would raise
    TypeError.
    """

    __slots__ = ("_counter", "_items")

    def __init__(self) -> None:
        self._items: list[tuple[float, int, subprocess.Popen[str]]] = []
        self._counter = itertools.count()

    def __len__(self) -> int:
        return len(self._items)

    def push(self, deadline: float, proc: subprocess.Popen[str]) -> None:
        heapq.heappush(self._items, (deadline, next(self._counter), proc))

    def pop(self) -> None:
        heapq.heappop(self._items)

    def peek(self) -> tuple[float, subprocess.Popen[str]]:
        deadline, _, proc = self._items[0]
        return deadline, proc


_timeout_queue = Queue[tuple[Any, float]]()
_shutting_down = False
_timeout_thread_lock = Lock()
_shutdown_registered = False
bgthd: Thread | None = None


def _timeout_thread_func() -> None:
    waiting = MinHeap()
    try:
        while not _shutting_down:
            if waiting:
                ttk, _ = waiting.peek()
                timeout = max(0, ttk - time.time())
            else:
                timeout = None
            with contextlib.suppress(QueueEmpty):
                proc, time_to_kill = _timeout_queue.get(timeout=timeout)
                if proc is SystemExit:
                    # terminate
                    return
                waiting.push(time_to_kill, proc)
            now = time.time()
            while waiting:
                ttk, proc = waiting.peek()
                if ttk > now:
                    break
                waiting.pop()
                with contextlib.suppress(OSError):
                    if proc.poll() is None:
                        proc.kill()
                        proc._timed_out = True  # type: ignore[attr-defined]
    except Exception:
        if _shutting_down:
            # to prevent all sorts of exceptions during interpreter shutdown
            pass
        else:
            raise


def _ensure_timeout_thread_started() -> None:
    global bgthd, _shutdown_registered  # noqa: PLW0603
    if bgthd is not None and bgthd.is_alive():
        return

    with _timeout_thread_lock:
        if bgthd is not None and bgthd.is_alive():
            return
        bgthd = Thread(target=_timeout_thread_func, name="PlumbumTimeoutThread")
        bgthd.daemon = True
        bgthd.start()
        if not _shutdown_registered:
            atexit.register(_shutdown_bg_threads)
            _shutdown_registered = True


def _register_proc_timeout(proc: PopenWithAddons[Any], timeout: float | None) -> None:
    if timeout is not None:
        _ensure_timeout_thread_started()
        _timeout_queue.put((proc, time.time() + timeout))


def _shutdown_bg_threads() -> None:
    global _shutting_down  # noqa: PLW0603
    _shutting_down = True
    # _timeout_queue could be deleted by this point
    if _timeout_queue and bgthd and bgthd.is_alive():  # type: ignore[truthy-bool]
        _timeout_queue.put((SystemExit, 0))
        # grace period
        bgthd.join(0.1)


def _close_streams(*streams: IO[Any] | None) -> None:
    """Close the given streams, ignoring ``None`` entries and close errors."""
    for stream in streams:
        if stream is not None:
            with contextlib.suppress(Exception):
                stream.close()


def _terminate_and_reap(proc: PopenWithAddons[Any], grace: float = 1.0) -> None:
    """Terminate *proc*, give it *grace* seconds to exit, then kill and reap it.

    Polls instead of ``wait(timeout=...)`` because not every popen-like object
    (remote/session procs) supports a wait timeout.
    """
    if proc.poll() is not None:
        return
    with contextlib.suppress(Exception):
        proc.terminate()
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            return
        time.sleep(0.05)
    with contextlib.suppress(Exception):
        proc.kill()
    with contextlib.suppress(Exception):
        proc.wait()


# ===================================================================================================
# run_proc
# ===================================================================================================
def run_proc(
    proc: PopenWithAddons[Any],
    retcode: int | None | Container[int],
    timeout: float | None = None,
) -> tuple[int | None, str, str]:
    """Waits for the given process to terminate, with the expected exit code

    :param proc: a running Popen-like object, with all the expected methods.

    :param retcode: the expected return (exit) code of the process. It defaults to 0 (the
                    convention for success). If ``None``, the return code is ignored.
                    It may also be a tuple (or any object that supports ``__contains__``)
                    of expected return codes.

    :param timeout: the number of seconds (a ``float``) to allow the process to run, before
                    forcefully terminating it. If ``None``, not timeout is imposed; otherwise
                    the process is expected to terminate within that timeout value, or it will
                    be killed and :class:`ProcessTimedOut <plumbum.commands.processes.ProcessTimedOut>`
                    will be raised

    :returns: A tuple of (return code, stdout, stderr)
    """
    _register_proc_timeout(proc, timeout)
    stdout: bytes | str
    stderr: bytes | str
    try:
        stdout, stderr = proc.communicate()
        proc._end_time = time.time()  # type: ignore[attr-defined]
        if not stdout:
            stdout = b""
        if not stderr:
            stderr = b""
        if custom_encoding := getattr(proc, "custom_encoding", None):
            assert isinstance(stdout, bytes)
            assert isinstance(stderr, bytes)
            stdout = stdout.decode(custom_encoding, "ignore")
            stderr = stderr.decode(custom_encoding, "ignore")

        return _check_process(proc, retcode, timeout, stdout, stderr)  # type: ignore[return-value]
    finally:
        if getattr(proc, "close_streams_after_communicate", True):
            _close_streams(proc.stdin, proc.stdout, proc.stderr)


# ===================================================================================================
# iter_lines
# ===================================================================================================


class Mode(enum.Enum):
    BY_POSITION = enum.auto()
    BY_TYPE = enum.auto()


BY_POSITION = Mode.BY_POSITION
BY_TYPE = Mode.BY_TYPE
DEFAULT_ITER_LINES_MODE = BY_POSITION

DEFAULT_BUFFER_SIZE = sys.maxsize


@typing.overload
def iter_lines(
    proc: PopenWithAddons[Any],
    retcode: int = ...,
    timeout: float | None = ...,
    linesize: int = ...,
    line_timeout: float | None = ...,
    buffer_size: int | None = ...,
    *,
    mode: Literal[Mode.BY_POSITION] | None = ...,
    _iter_lines: Callable[
        [PopenWithAddons[Any], Callable[[bytes], str], int, float | None],
        Generator[tuple[int, str], None, None],
    ] = _iter_lines,
) -> Generator[tuple[int, str], None, None]: ...


@typing.overload
def iter_lines(
    proc: PopenWithAddons[Any],
    retcode: int = ...,
    timeout: float | None = ...,
    linesize: int = ...,
    line_timeout: float | None = ...,
    buffer_size: int | None = ...,
    *,
    mode: Literal[Mode.BY_TYPE],
    _iter_lines: Callable[
        [PopenWithAddons[Any], Callable[[bytes], str], int, float | None],
        Generator[tuple[int, str], None, None],
    ] = _iter_lines,
) -> Generator[tuple[None, str] | tuple[str, None], None, None]: ...


def iter_lines(
    proc: PopenWithAddons[Any],
    retcode: int = 0,
    timeout: float | None = None,
    linesize: int = -1,
    line_timeout: float | None = None,
    buffer_size: int | None = None,
    *,
    mode: Mode | None = None,
    _iter_lines: Callable[
        [PopenWithAddons[Any], Callable[[bytes], str], int, float | None],
        Generator[tuple[int, str], None, None],
    ] = _iter_lines,
) -> Generator[tuple[int, str] | tuple[None, str] | tuple[str, None], None, None]:
    """Runs the given process (equivalent to run_proc()) and yields a tuples of (out, err) line pairs.
    If the exit code of the process does not match the expected one, :class:`ProcessExecutionError
    <plumbum.commands.processes.ProcessExecutionError>` is raised.

    :param retcode: The expected return code of this process (defaults to 0).
                    In order to disable exit-code validation, pass ``None``. It may also
                    be a tuple (or any iterable) of expected exit codes.

    :param timeout: The maximal amount of time (in seconds) to allow the process to run.
                    ``None`` means no timeout is imposed; otherwise, if the process hasn't
                    terminated after that many seconds, the process will be forcefully
                    terminated an exception will be raised

    :param linesize: Maximum number of characters to read from stdout/stderr at each iteration.
                    ``-1`` (default) reads until a b'\\n' is encountered.

    :param line_timeout: The maximal amount of time (in seconds) to allow between consecutive lines in either stream.
                    Raise an :class:`ProcessLineTimedOut <plumbum.commands.processes.ProcessLineTimedOut>` if the timeout has
                    been reached. ``None`` means no timeout is imposed.

    :param buffer_size: Maximum number of lines to keep in the stdout/stderr buffers, in case of a ProcessExecutionError.
                    Default is ``None``, which defaults to DEFAULT_BUFFER_SIZE (which is infinite by default).
                    ``0`` will disable buffering completely.

    :param mode: Controls what the generator yields. Defaults to DEFAULT_ITER_LINES_MODE (which is BY_POSITION by default)
                - BY_POSITION (default): yields ``(out, err)`` line tuples, where either item may be ``None``
                - BY_TYPE: yields ``(fd, line)`` tuples, where ``fd`` is 1 (stdout) or 2 (stderr)

    :returns: An iterator of (out, err) line tuples.
    """
    if mode is None:
        mode = DEFAULT_ITER_LINES_MODE

    if buffer_size is None:
        buffer_size = DEFAULT_BUFFER_SIZE

    assert mode in Mode

    encoding = getattr(proc, "custom_encoding", None) or "utf-8"
    decode = lambda s: s.decode(encoding, errors="replace").rstrip()  # noqa: E731

    _register_proc_timeout(proc, timeout)

    buffers: list[list[tuple[str | None, str | None] | str]] = [[], []]
    completed = False
    timed_out = False

    try:
        for t, line in _iter_lines(proc, decode, linesize, line_timeout):
            # verify that the proc hasn't timed out yet
            proc.verify(timeout=timeout, retcode=None, stdout=None, stderr=None)

            buffer = buffers[t]
            if buffer_size > 0:
                buffer.append(line)
                if buffer_size < sys.maxsize:
                    del buffer[:-buffer_size]

            if mode is BY_POSITION:
                if t == 0:
                    yield (line, None)
                else:
                    yield (None, line)

            elif mode is BY_TYPE:
                yield (t + 1), line  # 1=stdout, 2=stderr

        completed = True
    except ProcessLineTimedOut:
        timed_out = True
        raise
    finally:
        process_timed_out = timed_out or getattr(proc, "_timed_out", False)

        proc_running = proc.poll() is None

        if proc_running and process_timed_out:
            with contextlib.suppress(Exception):
                proc.kill()
            with contextlib.suppress(Exception):
                proc.wait()
        elif completed:
            # Generator consumed all lines; ensure process is reaped
            with contextlib.suppress(Exception):
                proc.wait()
        elif not proc_running:
            # Process already finished even though generator did not complete;
            # reap it to avoid zombies
            with contextlib.suppress(Exception):
                proc.wait()

        # Recompute running state after possible kill/wait above
        proc_running = proc.poll() is None

        # Only close streams once the process is known to have finished
        if not proc_running:
            _close_streams(proc.stdin, proc.stdout, proc.stderr)
    if completed:
        # this will take care of checking return code and timeouts
        _check_process(proc, retcode, timeout, *("\n".join(s) + "\n" for s in buffers))  # type: ignore[arg-type]


__all__ = [
    "BY_POSITION",
    "BY_TYPE",
    "DEFAULT_BUFFER_SIZE",
    "DEFAULT_ITER_LINES_MODE",
    "CommandNotFound",
    "Mode",
    "ProcessExecutionError",
    "ProcessLineTimedOut",
    "ProcessTimedOut",
    "iter_lines",
    "run_proc",
]


def __dir__() -> list[str]:
    return list(__all__)
