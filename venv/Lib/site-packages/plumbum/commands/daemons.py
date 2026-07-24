from __future__ import annotations

__lazy_modules__ = {"contextlib", "pickle", "plumbum.commands.processes", "subprocess"}

import contextlib
import os
import pickle
import subprocess
import sys
import time
import typing

from plumbum.commands.processes import ProcessExecutionError

if typing.TYPE_CHECKING:
    from plumbum.commands.base import BaseCommand
    from plumbum.machines.base import PopenWithAddons


class _fake_lock:
    """Needed to allow normal os.exit() to work without error"""

    __slots__ = ()

    @staticmethod
    def acquire(_: bool) -> bool:
        return True

    @staticmethod
    def release() -> None:
        pass


def posix_daemonize(
    command: BaseCommand,
    cwd: str,
    stdout: str | None = None,
    stderr: str | None = None,
    append: bool = True,
) -> PopenWithAddons[str]:
    if stdout is None:
        stdout = os.devnull
    if stderr is None:
        stderr = stdout

    argv = command.formulate()
    payload = pickle.dumps(
        {
            "command": command,
            "cwd": cwd,
            "stdout": stdout,
            "stderr": stderr,
            "append": append,
        },
        protocol=pickle.HIGHEST_PROTOCOL,
    )

    rfd, wfd = os.pipe()
    try:
        try:
            launcher = subprocess.Popen(
                [sys.executable, "-m", "plumbum.commands._daemon_launcher", str(wfd)],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                close_fds=True,
                pass_fds=(wfd,),
            )
        finally:
            os.close(wfd)

        launch_stdout, launch_stderr = launcher.communicate(payload)
        output: str | bytes = os.read(rfd, 16384)
        assert isinstance(output, bytes)
    finally:
        os.close(rfd)
    with contextlib.suppress(UnicodeError):
        output = output.decode("utf8")

    if launcher.returncode == 0 and output.isdigit():
        secondpid = int(output)
    else:
        launch_stdout_text = launch_stdout.decode("utf8", "ignore")
        launch_stderr_text = launch_stderr.decode("utf8", "ignore")
        launcher_output = "\n".join(
            part
            for part in (str(output), launch_stdout_text, launch_stderr_text)
            if part
        )
        raise ProcessExecutionError(
            argv,
            launcher.returncode,
            "",
            launcher_output,
        )
    proc: subprocess.Popen[bytes] = subprocess.Popen.__new__(subprocess.Popen)  # type: ignore[arg-type]
    proc._child_created = False  # type: ignore[attr-defined]
    proc.returncode = None
    proc.stdout = None
    proc.stdin = None
    proc.stderr = None
    proc.pid = secondpid
    proc.universal_newlines = False
    proc._input = None  # type: ignore[attr-defined]
    proc._waitpid_lock = _fake_lock()  # type: ignore[attr-defined]
    proc._communication_started = False  # type: ignore[attr-defined]
    proc.args = argv
    proc.argv = argv  # type: ignore[attr-defined]

    def poll(self: subprocess.Popen[bytes] = proc) -> int | None:
        if self.returncode is None:
            try:
                os.kill(self.pid, 0)
            except ProcessLookupError:
                # process does not exist (ESRCH)
                self.returncode = 0
            except PermissionError:
                # the process exists but we may not signal it (EPERM); treat it
                # as still running
                pass
        return self.returncode

    def wait(self: subprocess.Popen[bytes] = proc, timeout: float | None = None) -> int:
        deadline = None if timeout is None else time.monotonic() + timeout
        while self.poll() is None:
            if deadline is None:
                time.sleep(0.5)
            else:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise subprocess.TimeoutExpired(self.args, timeout)  # type: ignore[arg-type]
                time.sleep(min(0.5, remaining))
        return self.returncode

    proc.poll = poll  # type: ignore[method-assign]
    proc.wait = wait  # type: ignore[method-assign, assignment]
    return proc  # type: ignore[return-value]


def win32_daemonize(
    command: BaseCommand,
    cwd: str,
    stdout: str | None = None,
    stderr: str | None = None,
    append: bool = True,
) -> PopenWithAddons[str]:
    if stdout is None:
        stdout = os.devnull
    if stderr is None:
        stderr = stdout
    stdin_file = open(os.devnull, encoding="utf-8")
    stdout_file = open(stdout, "a" if append else "w", encoding="utf-8")
    stderr_file = open(stderr, "a" if append else "w", encoding="utf-8")
    try:
        return command.popen(
            cwd=cwd,
            stdin=stdin_file.fileno(),
            stdout=stdout_file.fileno(),
            stderr=stderr_file.fileno(),
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP  # type: ignore[attr-defined]
            | subprocess.DETACHED_PROCESS,  # type: ignore[attr-defined]
        )
    finally:
        for stream in (stdin_file, stdout_file, stderr_file):
            with contextlib.suppress(Exception):
                stream.close()


__all__ = [
    "posix_daemonize",
    "win32_daemonize",
]


def __dir__() -> list[str]:
    return list(__all__)
