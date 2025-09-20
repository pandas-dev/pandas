"""JupyterLab Server process handler"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import atexit
import logging
import os
import re
import signal
import subprocess
import sys
import threading
import time
import weakref
from logging import Logger
from shutil import which as _which
from typing import Any

from tornado import gen

try:
    import pty
except ImportError:
    pty = None  # type:ignore[assignment]

if sys.platform == "win32":
    list2cmdline = subprocess.list2cmdline
else:

    def list2cmdline(cmd_list: list[str]) -> str:
        """Shim for list2cmdline on posix."""
        import shlex

        return " ".join(map(shlex.quote, cmd_list))


def which(command: str, env: dict[str, str] | None = None) -> str:
    """Get the full path to a command.

    Parameters
    ----------
    command: str
        The command name or path.
    env: dict, optional
        The environment variables, defaults to `os.environ`.
    """
    env = env or os.environ  # type:ignore[assignment]
    path = env.get("PATH") or os.defpath  # type:ignore[union-attr]
    command_with_path = _which(command, path=path)

    # Allow nodejs as an alias to node.
    if command == "node" and not command_with_path:
        command = "nodejs"
        command_with_path = _which("nodejs", path=path)

    if not command_with_path:
        if command in ["nodejs", "node", "npm"]:
            msg = "Please install Node.js and npm before continuing installation. You may be able to install Node.js from your package manager, from conda, or directly from the Node.js website (https://nodejs.org)."
            raise ValueError(msg)
        raise ValueError("The command was not found or was not " + "executable: %s." % command)
    return os.path.abspath(command_with_path)


class Process:
    """A wrapper for a child process."""

    _procs: weakref.WeakSet = weakref.WeakSet()
    _pool = None

    def __init__(
        self,
        cmd: list[str],
        logger: Logger | None = None,
        cwd: str | None = None,
        kill_event: threading.Event | None = None,
        env: dict[str, str] | None = None,
        quiet: bool = False,
    ) -> None:
        """Start a subprocess that can be run asynchronously.

        Parameters
        ----------
        cmd: list
            The command to run.
        logger: :class:`~logger.Logger`, optional
            The logger instance.
        cwd: string, optional
            The cwd of the process.
        env: dict, optional
            The environment for the process.
        kill_event: :class:`~threading.Event`, optional
            An event used to kill the process operation.
        quiet: bool, optional
            Whether to suppress output.
        """
        if not isinstance(cmd, (list, tuple)):
            msg = "Command must be given as a list"  # type:ignore[unreachable]
            raise ValueError(msg)

        if kill_event and kill_event.is_set():
            msg = "Process aborted"
            raise ValueError(msg)

        self.logger = logger or self.get_log()
        self._last_line = ""
        if not quiet:
            self.logger.info("> %s", list2cmdline(cmd))
        self.cmd = cmd

        kwargs = {}
        if quiet:
            kwargs["stdout"] = subprocess.DEVNULL

        self.proc = self._create_process(cwd=cwd, env=env, **kwargs)
        self._kill_event = kill_event or threading.Event()

        Process._procs.add(self)

    def terminate(self) -> int:
        """Terminate the process and return the exit code."""
        proc = self.proc

        # Kill the process.
        if proc.poll() is None:
            os.kill(proc.pid, signal.SIGTERM)

        # Wait for the process to close.
        try:
            proc.wait(timeout=2.0)
        except subprocess.TimeoutExpired:
            if os.name == "nt":  # noqa: SIM108
                sig = signal.SIGBREAK  # type:ignore[attr-defined]
            else:
                sig = signal.SIGKILL

            if proc.poll() is None:
                os.kill(proc.pid, sig)

        finally:
            if self in Process._procs:
                Process._procs.remove(self)

        return proc.wait()

    def wait(self) -> int:
        """Wait for the process to finish.

        Returns
        -------
        The process exit code.
        """
        proc = self.proc
        kill_event = self._kill_event
        while proc.poll() is None:
            if kill_event.is_set():
                self.terminate()
                msg = "Process was aborted"
                raise ValueError(msg)
            time.sleep(1.0)
        return self.terminate()

    @gen.coroutine
    def wait_async(self) -> Any:
        """Asynchronously wait for the process to finish."""
        proc = self.proc
        kill_event = self._kill_event
        while proc.poll() is None:
            if kill_event.is_set():
                self.terminate()
                msg = "Process was aborted"
                raise ValueError(msg)
            yield gen.sleep(1.0)

        raise gen.Return(self.terminate())

    def _create_process(self, **kwargs: Any) -> subprocess.Popen[str]:
        """Create the process."""
        cmd = list(self.cmd)
        kwargs.setdefault("stderr", subprocess.STDOUT)

        cmd[0] = which(cmd[0], kwargs.get("env"))

        if os.name == "nt":
            kwargs["shell"] = True

        return subprocess.Popen(cmd, **kwargs)  # noqa: S603

    @classmethod
    def _cleanup(cls: type[Process]) -> None:
        """Clean up the started subprocesses at exit."""
        for proc in list(cls._procs):
            proc.terminate()

    def get_log(self) -> Logger:
        """Get our logger."""
        if hasattr(self, "logger") and self.logger is not None:
            return self.logger
        # fallback logger
        self.logger = logging.getLogger("jupyterlab")
        self.logger.setLevel(logging.INFO)
        return self.logger


class WatchHelper(Process):
    """A process helper for a watch process."""

    def __init__(
        self,
        cmd: list[str],
        startup_regex: str,
        logger: Logger | None = None,
        cwd: str | None = None,
        kill_event: threading.Event | None = None,
        env: dict[str, str] | None = None,
    ) -> None:
        """Initialize the process helper.

        Parameters
        ----------
        cmd: list
            The command to run.
        startup_regex: string
            The regex to wait for at startup.
        logger: :class:`~logger.Logger`, optional
            The logger instance.
        cwd: string, optional
            The cwd of the process.
        env: dict, optional
            The environment for the process.
        kill_event: callable, optional
            A function to call to check if we should abort.
        """
        super().__init__(cmd, logger=logger, cwd=cwd, kill_event=kill_event, env=env)

        if pty is None:
            self._stdout = self.proc.stdout  # type:ignore[unreachable]

        while 1:
            line = self._stdout.readline().decode("utf-8")  # type:ignore[has-type]
            if not line:
                msg = "Process ended improperly"
                raise RuntimeError(msg)
            print(line.rstrip())
            if re.match(startup_regex, line):
                break

        self._read_thread = threading.Thread(target=self._read_incoming, daemon=True)
        self._read_thread.start()

    def terminate(self) -> int:
        """Terminate the process."""
        proc = self.proc

        if proc.poll() is None:
            if os.name != "nt":
                # Kill the process group if we started a new session.
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            else:
                os.kill(proc.pid, signal.SIGTERM)

        # Wait for the process to close.
        try:
            proc.wait()
        finally:
            if self in Process._procs:
                Process._procs.remove(self)

        return proc.returncode

    def _read_incoming(self) -> None:
        """Run in a thread to read stdout and print"""
        fileno = self._stdout.fileno()  # type:ignore[has-type]
        while 1:
            try:
                buf = os.read(fileno, 1024)
            except OSError as e:
                self.logger.debug("Read incoming error %s", e)
                return

            if not buf:
                return

            print(buf.decode("utf-8"), end="")

    def _create_process(self, **kwargs: Any) -> subprocess.Popen[str]:
        """Create the watcher helper process."""
        kwargs["bufsize"] = 0

        if pty is not None:
            master, slave = pty.openpty()
            kwargs["stderr"] = kwargs["stdout"] = slave
            kwargs["start_new_session"] = True
            self._stdout = os.fdopen(master, "rb")  # type:ignore[has-type]
        else:
            kwargs["stdout"] = subprocess.PIPE  # type:ignore[unreachable]

            if os.name == "nt":
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                kwargs["startupinfo"] = startupinfo
                kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
                kwargs["shell"] = True

        return super()._create_process(**kwargs)


# Register the cleanup handler.
atexit.register(Process._cleanup)
