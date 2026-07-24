from __future__ import annotations

__lazy_modules__ = {"contextlib", "plumbum.commands", "random", "threading"}

import contextlib
import logging
import random
import threading
import time
import typing
from typing import IO, Any

from plumbum.commands import BaseCommand, run_proc
from plumbum.commands.processes import ProcessExecutionError
from plumbum.machines.base import PopenAddons

if typing.TYPE_CHECKING:
    from collections.abc import Container

    from plumbum._compat.typing import Self
    from plumbum.machines.base import PopenWithAddons
    from plumbum.machines.paramiko_machine import ParamikoPopen


class ShellSessionError(Exception):
    """Raises when something goes wrong when calling
    :func:`ShellSession.popen <plumbum.machines.session.ShellSession.popen>`"""


class SSHCommsError(ProcessExecutionError, EOFError):
    """Raises when the communication channel can't be created on the
    remote host or it times out."""


class SSHCommsChannel2Error(SSHCommsError):
    """Raises when channel 2 (stderr) is not available"""


class IncorrectLogin(SSHCommsError):
    """Raises when incorrect login credentials are provided"""


class HostPublicKeyUnknown(SSHCommsError):
    """Raises when the host public key isn't known"""


shell_logger = logging.getLogger("plumbum.shell")


# ===================================================================================================
# Shell Session Popen
# ===================================================================================================
class MarkedPipe:
    """A pipe-like object from which you can read lines; the pipe will return report EOF (the
    empty string) when a special marker is detected"""

    __slots__ = ["__weakref__", "marker", "pipe"]

    def __init__(self, pipe: IO[bytes], marker: str) -> None:
        self.pipe: IO[bytes] | None = pipe
        self.marker = bytes(marker, "ascii")

    def close(self) -> None:
        """'Closes' the marked pipe; following calls to ``readline`` will return "" """
        # consume everything
        while self.readline():
            pass
        self.pipe = None

    def readline(self) -> bytes:
        """Reads the next line from the pipe; returns "" when the special marker is reached.
        Raises ``EOFError`` if the underlying pipe has closed"""
        if self.pipe is None:
            return b""
        line = self.pipe.readline()
        if not line:
            raise EOFError()
        if line.strip() == self.marker:
            self.pipe = None
            return b""
        return line


class SessionPopen(PopenAddons):
    """A shell-session-based ``Popen``-like object (has the following attributes: ``stdin``,
    ``stdout``, ``stderr``, ``returncode``)"""

    returncode: int | None

    def __init__(
        self,
        proc: SessionPopen | PopenWithAddons[Any],
        argv: bytes,
        isatty: bool,
        stdin: IO[bytes],
        stdout: MarkedPipe,
        stderr: MarkedPipe,
        encoding: str,
        *,
        host: str | None,
    ) -> None:
        self.host = host
        self.proc: PopenWithAddons[Any] = proc  # type: ignore[assignment]
        self.argv = argv
        self.isatty = isatty
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        self.custom_encoding = encoding
        self.close_streams_after_communicate = False
        self.returncode = None
        self._done = False

    def poll(self) -> int | None:
        """Returns the process' exit code or ``None`` if it's still running"""
        return self.returncode if self._done else None

    def wait(self) -> int:
        """Waits for the process to terminate and returns its exit code"""
        self.communicate()
        assert self.returncode is not None
        return self.returncode

    def communicate(
        self,
        input: bytes | bytearray | None = None,  # pylint: disable=redefined-builtin
    ) -> tuple[bytes, bytes]:
        """Consumes the process' stdout and stderr until the it terminates.

        :param input: An optional bytes/buffer object to send to the process over stdin
        :returns: A tuple of (stdout, stderr)
        """
        stdout: list[bytes] = []
        stderr: list[bytes] = []
        sources = [("1", stdout, self.stdout)]
        if not self.isatty:
            # in tty mode, stdout and stderr are unified
            sources.append(("2", stderr, self.stderr))
        i = 0
        while sources:
            if input:
                chunk = input[:1000]
                self.stdin.write(chunk)
                self.stdin.flush()
                input = input[1000:]
            i = (i + 1) % len(sources)
            name, coll, pipe = sources[i]
            try:
                line = pipe.readline()
                shell_logger.debug("%s> %r", name, line)
            except EOFError as err:
                shell_logger.debug("%s> Nothing returned.", name)

                self.proc.poll()
                returncode = self.proc.returncode
                assert self.custom_encoding is not None
                stdout_txt = b"".join(stdout).decode(self.custom_encoding, "ignore")
                stderr_txt = b"".join(stderr).decode(self.custom_encoding, "ignore")
                argv = self.argv.decode(self.custom_encoding, "ignore").split(";")[:1]

                if returncode == 5:
                    raise IncorrectLogin(
                        argv,
                        returncode,
                        stdout_txt,
                        stderr_txt,
                        message="Incorrect username or password provided",
                        host=self.host,
                    ) from None
                if returncode == 6:
                    raise HostPublicKeyUnknown(
                        argv,
                        returncode,
                        stdout_txt,
                        stderr_txt,
                        message="The authenticity of the host can't be established",
                        host=self.host,
                    ) from None
                if returncode != 0:
                    raise SSHCommsError(
                        argv,
                        returncode,
                        stdout_txt,
                        stderr_txt,
                        message="SSH communication failed",
                        host=self.host,
                    ) from None
                if name == "2":
                    raise SSHCommsChannel2Error(
                        argv,
                        returncode,
                        stdout_txt,
                        stderr_txt,
                        message="No stderr result detected. Does the remote have Bash as the default shell?",
                        host=self.host,
                    ) from None

                raise SSHCommsError(
                    argv,
                    returncode,
                    stdout_txt,
                    stderr_txt,
                    message="No communication channel detected. Does the remote exist?",
                    host=self.host,
                ) from err
            if not line:
                del sources[i]
            else:
                coll.append(line)
        if self.isatty:
            stdout.pop(0)  # discard first line of prompt
        try:
            self.returncode = int(stdout.pop(-1))
        except (IndexError, ValueError):
            self.returncode = "Unknown"  # type: ignore[assignment]
        self._done = True
        stdout_bytes = b"".join(stdout)
        # the command trailer prints the return code (stdout) and the marker
        # (stderr) on their own lines, so it always appends a newline to the real
        # output; strip it to return the output verbatim (e.g. `printf abc` yields
        # b"abc", not b"abc\n"). Under a tty the injected newline is translated to
        # CRLF, so strip that instead of eating a real "\r".
        newline = b"\r\n" if self.isatty else b"\n"
        stdout_bytes = stdout_bytes.removesuffix(newline)
        stderr_bytes = b"".join(stderr).removesuffix(b"\n")
        return stdout_bytes, stderr_bytes


class ShellSession:
    """An abstraction layer over *shell sessions*. A shell session is the execution of an
    interactive shell (``/bin/sh`` or something compatible), over which you may run commands
    (sent over stdin). The output of is then read from stdout and stderr. Shell sessions are
    less "robust" than executing a process on its own, and they are susseptible to all sorts
    of malformatted-strings attacks, and there is little benefit from using them locally.
    However, they can greatly speed up remote connections, and are required for the implementation
    of :class:`SshMachine <plumbum.machines.ssh_machine.SshMachine>`, as they allow us to send multiple
    commands over a single SSH connection (setting up separate SSH connections incurs a high
    overhead). Try to avoid using shell sessions, unless you know what you're doing.

    Instances of this class may be used as *context-managers*.

    :param proc: The underlying shell process (with open stdin, stdout and stderr)
    :param encoding: The encoding to use for the shell session. If ``"auto"``, the underlying
                     process' encoding is used.
    :param isatty: If true, assume the shell has a TTY and that stdout and stderr are unified
    :param connect_timeout: The timeout to connect to the shell, after which, if no prompt
                            is seen, the shell process is killed
    """

    __slots__ = (
        "_current",
        "_lock",
        "_startup_result",
        "custom_encoding",
        "host",
        "isatty",
        "proc",
    )

    def __init__(
        self,
        proc: PopenWithAddons[Any] | SessionPopen | ParamikoPopen,
        encoding: str = "auto",
        isatty: bool = False,
        connect_timeout: float | None = 5,
        *,
        host: str | None = None,
    ) -> None:
        self.host = host
        self.proc: PopenWithAddons[Any] | None = proc  # type: ignore[assignment]
        self.custom_encoding: str | None = (
            proc.custom_encoding if encoding == "auto" else encoding
        )
        self.isatty = isatty
        self._lock = threading.RLock()
        self._current: SessionPopen | None = None
        self._startup_result: tuple[int | None, str, str] | None = None
        if connect_timeout:

            def closer() -> None:
                shell_logger.error(
                    "Connection to %s timed out (%d sec)", proc, connect_timeout
                )
                self.close()

            timer = threading.Timer(connect_timeout, closer)
            timer.start()
        try:
            self._startup_result = self.run("")
        finally:
            if connect_timeout:
                timer.cancel()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, t: object, v: object, tb: object) -> None:
        self.close()

    def __del__(self) -> None:
        with contextlib.suppress(Exception):
            self.close()

    def alive(self) -> bool:
        """Returns ``True`` if the underlying shell process is alive, ``False`` otherwise"""
        return bool(self.proc and self.proc.poll() is None)

    def close(self) -> None:
        """Closes (terminates) the shell session"""
        if not self.alive():
            return
        assert self.proc is not None
        with contextlib.suppress(ValueError, OSError):
            self.proc.stdin.write(b"\nexit\n\n\nexit\n\n")  # type: ignore[union-attr]
            self.proc.stdin.flush()  # type: ignore[union-attr]
            time.sleep(0.05)
        for p in (self.proc.stdin, self.proc.stdout, self.proc.stderr):
            with contextlib.suppress(Exception):
                p.close()  # type: ignore[union-attr]
        with contextlib.suppress(OSError):
            self.proc.kill()
        self.proc = None

    def popen(self, cmd: str | BaseCommand) -> PopenWithAddons[Any]:
        """Runs the given command in the shell, adding some decoration around it. Only a single
        command can be executed at any given time.

        :param cmd: The command (string or :class:`Command <plumbum.commands.base.BaseCommand>` object)
                    to run
        :returns: A :class:`SessionPopen <plumbum.machines.session.SessionPopen>` instance
        """
        if self.proc is None:
            raise ShellSessionError("Shell session has already been closed")
        if self._current and not self._current._done:
            raise ShellSessionError("Each shell may start only one process at a time")

        full_cmd = " ".join(cmd.formulate(1)) if isinstance(cmd, BaseCommand) else cmd
        marker = f"--.END{time.time() * random.random()}.--"
        if full_cmd.strip():
            full_cmd += " ; "
        else:
            full_cmd = "true ; "
        full_cmd += f"printf '\\n%s\\n' \"$?\" ; echo '{marker}'"
        if not self.isatty:
            # the leading newline keeps the marker off the last line of stderr
            # output, just like the return code on stdout; otherwise unterminated
            # stderr output would hide the marker and hang communicate()
            full_cmd += f" ; printf '\\n%s\\n' '{marker}' 1>&2"
        if self.custom_encoding:
            full_cmd_bytes = full_cmd.encode(self.custom_encoding)
        else:
            # TODO: I don't think bytes are supported above
            full_cmd_bytes = full_cmd  # type: ignore[assignment]
        shell_logger.debug("Running %r", full_cmd)
        assert self.proc.stdin is not None
        self.proc.stdin.write(full_cmd_bytes + b"\n")
        self.proc.stdin.flush()
        assert self.proc.stdout is not None
        assert self.proc.stderr is not None
        assert self.custom_encoding is not None
        _current = SessionPopen(
            self.proc,
            full_cmd_bytes,
            self.isatty,
            self.proc.stdin,
            MarkedPipe(self.proc.stdout, marker),
            MarkedPipe(self.proc.stderr, marker),
            self.custom_encoding,
            host=self.host,
        )
        self._current = _current
        return _current  # type: ignore[return-value]

    def run(
        self, cmd: str | BaseCommand, retcode: int | None | Container[int] = 0
    ) -> tuple[int | None, str, str]:
        """Runs the given command

        :param cmd: The command (string or :class:`Command <plumbum.commands.base.BaseCommand>` object)
                    to run
        :param retcode: The expected return code (0 by default). Set to ``None`` in order to
                        ignore erroneous return codes
        :returns: A tuple of (return code, stdout, stderr)
        """
        with self._lock:
            return run_proc(self.popen(cmd), retcode)


__all__ = [
    "HostPublicKeyUnknown",
    "IncorrectLogin",
    "MarkedPipe",
    "SSHCommsChannel2Error",
    "SSHCommsError",
    "SessionPopen",
    "ShellSession",
    "ShellSessionError",
]


def __dir__() -> list[str]:
    return list(__all__)
