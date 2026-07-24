from __future__ import annotations

__lazy_modules__ = {
    "contextlib",
    "platform",
    "plumbum.commands.daemons",
    "plumbum.commands.processes",
    "plumbum.machines._windows",
    "plumbum.machines.session",
    "plumbum.path",
    "plumbum.path.remote",
    "re",
    "subprocess",
    "tempfile",
}

import contextlib
import logging
import os
import platform
import re
import subprocess
import sys
import threading
import time
import typing
from contextlib import AbstractContextManager, contextmanager
from subprocess import PIPE, Popen
from tempfile import mkdtemp
from typing import Any, ClassVar

from plumbum.commands import CommandNotFound, ConcreteCommand
from plumbum.commands.daemons import posix_daemonize, win32_daemonize
from plumbum.commands.processes import iter_lines
from plumbum.lib import IS_WIN32, ProcInfo, StaticProperty
from plumbum.machines.base import BaseMachine, PopenAddons, PopenWithAddons
from plumbum.machines.env import BaseEnv
from plumbum.machines.session import ShellSession
from plumbum.path.local import LocalPath, LocalWorkdir
from plumbum.path.remote import RemotePath

if typing.TYPE_CHECKING:
    from collections.abc import Generator, Iterator, Sequence
    from types import TracebackType

    from plumbum.commands.async_ import AsyncLocalCommand
    from plumbum.commands.base import BaseCommand


class PlumbumLocalPopen(PopenAddons):
    iter_lines = iter_lines

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._proc = Popen(*args, **kwargs)  # pylint: disable=consider-using-with

    def __iter__(self) -> Iterator[tuple[int, str]]:
        return self.iter_lines()

    def __enter__(self) -> Popen[bytes]:
        return self._proc.__enter__()

    def __exit__(
        self,
        t: type[BaseException] | None,
        v: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        return self._proc.__exit__(t, v, tb)

    def __getattr__(self, name: str) -> Any:
        try:
            proc = object.__getattribute__(self, "_proc")
        except AttributeError:
            raise AttributeError(name) from None
        return getattr(proc, name)

    def __del__(self) -> None:
        with contextlib.suppress(AttributeError):
            proc = object.__getattribute__(self, "_proc")
            with contextlib.suppress(Exception):
                proc.poll()
            if proc.returncode is None:
                with contextlib.suppress(subprocess.TimeoutExpired, Exception):
                    proc.wait(timeout=0.2)
            for stream in (proc.stdin, proc.stdout, proc.stderr):
                if stream is not None:
                    with contextlib.suppress(Exception):
                        stream.close()


if IS_WIN32:
    from plumbum.machines._windows import IMAGE_SUBSYSTEM_WINDOWS_CUI, get_pe_subsystem

logger = logging.getLogger("plumbum.local")

# ``expand``/``expanduser`` below temporarily swap out the global ``os.environ``
# while expanding. Without serialization, two threads interleaving their
# save/restore can permanently leave ``os.environ`` pointing at a throwaway
# dict. This matters now that the async layer runs sync commands in executor
# threads. The lock keeps the swap/restore atomic with respect to each other.
_environ_lock = threading.RLock()


# ===================================================================================================
# Environment
# ===================================================================================================
class LocalEnv(BaseEnv[LocalPath]):
    """The local machine's environment; exposes a dict-like interface"""

    __slots__ = ()
    CASE_SENSITIVE = not IS_WIN32

    def __init__(self) -> None:
        # os.environ already takes care of upper'ing on windows
        super().__init__(LocalPath, os.path.pathsep, _curr=os.environ.copy())
        if IS_WIN32 and "HOME" not in self and self.home is not None:
            self["HOME"] = self.home

    def expand(self, expr: str) -> str:
        """Expands any environment variables and home shortcuts found in ``expr``
        (like ``os.path.expanduser`` combined with ``os.path.expandvars``)

        :param expr: An expression containing environment variables (as ``$FOO``) or
                     home shortcuts (as ``~/.bashrc``)

        :returns: The expanded string"""
        with _environ_lock:
            prev = os.environ
            os.environ = self.getdict()  # type: ignore[assignment] # noqa: B003
            try:
                output = os.path.expanduser(os.path.expandvars(expr))
            finally:
                os.environ = prev  # noqa: B003
        return output

    def expanduser(self, expr: str) -> str:
        """Expand home shortcuts (e.g., ``~/foo/bar`` or ``~john/foo/bar``)

        :param expr: An expression containing home shortcuts

        :returns: The expanded string"""
        with _environ_lock:
            prev = os.environ
            os.environ = self.getdict()  # type: ignore[assignment] # noqa: B003
            try:
                output = os.path.expanduser(expr)
            finally:
                os.environ = prev  # noqa: B003
        return output


# ===================================================================================================
# Local Commands
# ===================================================================================================
class LocalCommand(ConcreteCommand):
    __slots__ = ()
    QUOTE_LEVEL = 2

    def __init__(self, executable: LocalPath | str, encoding: str = "auto") -> None:
        ConcreteCommand.__init__(
            self, executable, local.custom_encoding if encoding == "auto" else encoding
        )

    @property
    def machine(self) -> LocalMachine:
        return local

    def popen(
        self,
        args: Sequence[Any] = (),
        cwd: str | LocalPath | None = None,
        env: dict[str, str] | BaseEnv[LocalPath] | None = None,
        **kwargs: Any,
    ) -> PopenWithAddons[str]:
        if isinstance(args, str):
            args = (args,)
        return self.machine._popen(
            self.executable,
            self.formulate(0, args),
            cwd=self.cwd if cwd is None else cwd,
            env=self.env if env is None else env,
            **kwargs,
        )


# ===================================================================================================
# Local Machine
# ===================================================================================================


class LocalMachine(BaseMachine):
    """The *local machine* (a singleton object). It serves as an entry point to everything
    related to the local machine, such as working directory and environment manipulation,
    command creation, etc.

    Attributes:

    * ``cwd`` - the local working directory
    * ``env`` - the local environment
    * ``custom_encoding`` - the local machine's default encoding (``sys.getfilesystemencoding()``)
    """

    __slots__ = ("_as_user_stack", "_start_time")

    cwd = StaticProperty(LocalWorkdir)
    env = LocalEnv()

    custom_encoding = sys.getfilesystemencoding()
    uname = platform.uname()[0]
    _program_cache: ClassVar[dict[tuple[str, str], LocalPath]] = {}

    def __init__(self) -> None:
        self._as_user_stack: list[Any] = []

    @classmethod
    def clear_program_cache(cls) -> None:
        cls._program_cache.clear()

    if IS_WIN32:
        _EXTENSIONS: list[str] = [
            "",
            *env.get("PATHEXT", ":.exe:.bat").lower().split(os.path.pathsep),
        ]

        @classmethod
        def _which(cls, progname: str) -> LocalPath | None:
            progname = progname.lower()
            for p in cls.env.path:
                for ext in cls._EXTENSIONS:
                    fn = p / (progname + ext)
                    if fn.access("x") and not fn.is_dir():
                        return fn
            return None

    else:

        @classmethod
        def _which(cls, progname: str) -> LocalPath | None:
            for p in cls.env.path:
                assert isinstance(p, LocalPath)
                fn = p / progname
                if fn.access("x") and not fn.is_dir():
                    return fn
            return None

    @classmethod
    def which(cls, progname: str) -> LocalPath:
        """Looks up a program in the ``PATH``. If the program is not found, raises
        :class:`CommandNotFound <plumbum.commands.processes.CommandNotFound>`

        :param progname: The program's name. Note that if underscores (``_``) are present
                         in the name, and the exact name is not found, they will be replaced
                         in turn by hyphens (``-``) then periods (``.``), and the name will
                         be looked up again for each alternative

        :returns: A :class:`LocalPath <plumbum.path.local.LocalPath>`
        """

        key = (progname, cls.env.get("PATH", ""))

        with contextlib.suppress(KeyError):
            return cls._program_cache[key]

        alternatives = [progname]
        if "_" in progname:
            alternatives += [progname.replace("_", "-"), progname.replace("_", ".")]
        for pn in alternatives:
            path = cls._which(pn)
            if path:
                cls._program_cache[key] = path
                return path
        raise CommandNotFound(progname, list(cls.env.path))

    def path(self, *parts: str) -> LocalPath:
        """A factory for :class:`LocalPaths <plumbum.path.local.LocalPath>`.
        Usage: ``p = local.path("/usr", "lib", "python2.7")``
        """
        parts2 = [str(self.cwd)]
        for p in parts:
            if isinstance(p, RemotePath):
                raise TypeError(f"Cannot construct LocalPath from {p!r}")
            parts2.append(self.env.expanduser(str(p)))
        return LocalPath(os.path.join(*parts2))

    def __getitem__(self, cmd: str | LocalPath) -> LocalCommand:
        """Returns a `Command` object representing the given program. ``cmd`` can be a string or
        a :class:`LocalPath <plumbum.path.local.LocalPath>`; if it is a path, a command
        representing this path will be returned; otherwise, the program name will be looked up
        in the system's ``PATH`` (using ``which``). Usage::

            ls = local["ls"]
        """

        if isinstance(cmd, LocalPath):
            return LocalCommand(cmd)

        if not isinstance(cmd, RemotePath):
            # handle "path-like" (pathlib.Path) objects
            cmd = os.fspath(cmd)
            if "/" in cmd or "\\" in cmd:
                # assume path
                return LocalCommand(local.path(cmd))
            # search for command
            return LocalCommand(self.which(cmd))

        raise TypeError(f"cmd must not be a RemotePath: {cmd!r}")

    def _popen(
        self,
        executable: LocalPath,
        argv: list[str],
        stdin: Any = PIPE,
        stdout: Any = PIPE,
        stderr: Any = PIPE,
        cwd: str | LocalPath | None = None,
        env: dict[str, str] | BaseEnv[LocalPath] | None = None,
        new_session: bool = False,
        **kwargs: Any,
    ) -> PlumbumLocalPopen:
        if new_session:
            kwargs["start_new_session"] = True

        if IS_WIN32 and "startupinfo" not in kwargs and stdin not in (sys.stdin, None):
            # pylint: disable-next=used-before-assignment
            subsystem = get_pe_subsystem(str(executable))

            # pylint: disable-next=used-before-assignment
            if subsystem == IMAGE_SUBSYSTEM_WINDOWS_CUI:
                # don't open a new console
                sui = subprocess.STARTUPINFO()  # type: ignore[attr-defined]
                kwargs["startupinfo"] = sui
                if hasattr(subprocess, "_subprocess"):
                    sui.dwFlags |= subprocess._subprocess.STARTF_USESHOWWINDOW
                    sui.wShowWindow = subprocess._subprocess.SW_HIDE
                else:
                    sui.dwFlags |= subprocess.STARTF_USESHOWWINDOW  # type: ignore[attr-defined]
                    sui.wShowWindow = subprocess.SW_HIDE  # type: ignore[attr-defined]

        if cwd is None:
            cwd = self.cwd

        envs = [self.env, env]
        env = {}
        for _env in envs:
            if not _env:
                continue
            if isinstance(_env, BaseEnv):
                _env = _env.getdict()
            env.update(_env)

        if self._as_user_stack:
            argv, executable = self._as_user_stack[-1](argv)

        logger.debug("Running %r", argv)
        proc = PlumbumLocalPopen(
            argv,
            executable=str(executable),
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            cwd=str(cwd),
            env=env,
            **kwargs,
        )  # bufsize = 4096
        proc._start_time = time.time()  # type: ignore[attr-defined]
        proc.custom_encoding = self.custom_encoding
        proc.argv = argv  # type: ignore[attr-defined]
        return proc

    def daemonic_popen(
        self,
        command: BaseCommand,
        cwd: str = "/",
        stdout: str | None = None,
        stderr: str | None = None,
        append: bool = True,
    ) -> PopenWithAddons[str]:
        """
        On POSIX systems:

        Run ``command`` as a UNIX daemon: fork a child process to setpid, redirect std handles to /dev/null,
        umask, close all fds, chdir to ``cwd``, then fork and exec ``command``. Returns a ``Popen`` process that
        can be used to poll/wait for the executed command (but keep in mind that you cannot access std handles)

        On Windows:

        Run ``command`` as a "Windows daemon": detach from controlling console and create a new process group.
        This means that the command will not receive console events and would survive its parent's termination.
        Returns a ``Popen`` object.

        .. note:: this does not run ``command`` as a system service, only detaches it from its parent.

        .. versionadded:: 1.3
        """
        if IS_WIN32:
            return win32_daemonize(command, cwd, stdout, stderr, append)

        return posix_daemonize(command, cwd, stdout, stderr, append)

    if IS_WIN32:

        def list_processes(self) -> Generator[ProcInfo, None, None]:  # pylint: disable=no-self-use
            """
            Returns information about all running processes (on Windows: using ``tasklist``)

            .. versionadded:: 1.3
            """
            import csv

            tasklist = local["tasklist"]
            output = tasklist("/V", "/FO", "CSV")
            lines = output.splitlines()
            rows = csv.reader(lines)
            try:
                header = next(rows)
            except StopIteration:
                raise RuntimeError("tasklist must at least have header") from None
            imgidx = header.index("Image Name")
            pididx = header.index("PID")
            statidx = header.index("Status")
            useridx = header.index("User Name")
            for row in rows:
                yield ProcInfo(
                    int(row[pididx]), row[useridx], row[statidx], row[imgidx]
                )

    else:

        def list_processes(self) -> Generator[ProcInfo, None, None]:
            """
            Returns information about all running processes (on POSIX systems: using ``ps``)

            .. versionadded:: 1.3
            """
            ps = self["ps"]
            lines = ps("-e", "-o", "pid,uid,stat,args").splitlines()
            lines.pop(0)  # header
            for line in lines:
                parts = line.strip().split()
                yield ProcInfo(
                    int(parts[0]), int(parts[1]), parts[2], " ".join(parts[3:])
                )

    def pgrep(self, pattern: str) -> Generator[ProcInfo, None, None]:
        """
        Process grep: return information about all processes whose command-line args match the given regex pattern
        """
        pat = re.compile(pattern)
        for procinfo in self.list_processes():
            if pat.search(procinfo.args):
                yield procinfo

    def session(self, new_session: bool = False) -> ShellSession:
        """Creates a new :class:`ShellSession <plumbum.machines.session.ShellSession>` object; this
        invokes ``/bin/sh`` and executes commands on it over stdin/stdout/stderr"""
        return ShellSession(self["sh"].popen(new_session=new_session))

    @contextmanager
    def tempdir(self) -> Generator[LocalPath, None, None]:
        """A context manager that creates a temporary directory, which is removed when the context
        exits"""
        new_dir = self.path(mkdtemp())
        try:
            yield new_dir
        finally:
            new_dir.delete()

    @contextmanager
    def as_user(self, username: str | None = None) -> Generator[None, None, None]:
        """Run nested commands as the given user. For example::

            head = local["head"]
            head("-n1", "/dev/sda1")    # this will fail...
            with local.as_user():
                head("-n1", "/dev/sda1")

        :param username: The user to run commands as. If not given, root (or Administrator) is assumed
        """
        if IS_WIN32:
            if username is None:
                username = "Administrator"
            self._as_user_stack.append(
                lambda argv: (
                    [
                        "runas",
                        "/savecred",
                        f"/user:{username}",
                        '"' + " ".join(str(a) for a in argv) + '"',
                    ],
                    self.which("runas"),
                )
            )
        elif username is None:
            self._as_user_stack.append(
                lambda argv: (["sudo", *list(argv)], self.which("sudo"))
            )
        else:
            self._as_user_stack.append(
                lambda argv: (
                    ["sudo", "-u", username, *list(argv)],
                    self.which("sudo"),
                )
            )
        try:
            yield
        finally:
            self._as_user_stack.pop(-1)

    def as_root(self) -> AbstractContextManager[None]:
        """A shorthand for :func:`as_user("root") <plumbum.machines.local.LocalMachine.as_user>`"""
        return self.as_user()

    @property
    def python(self) -> LocalCommand:
        """A command that represents the current python interpreter (``sys.executable``)."""
        return LocalCommand(sys.executable, self.custom_encoding)


local = LocalMachine()
"""The *local machine* (a singleton object). It serves as an entry point to everything
related to the local machine, such as working directory and environment manipulation,
command creation, etc.

Attributes:

* ``cwd`` - the local working directory
* ``env`` - the local environment
* ``custom_encoding`` - the local machine's default encoding (``sys.getfilesystemencoding()``)
"""


class AsyncLocalMachine:
    """Async version of LocalMachine.

    This class provides async access to local commands and utilities.
    It delegates to the sync LocalMachine for command lookup and wraps
    the results in AsyncLocalCommand.

    Example::

        from plumbum.machines.local import async_local

        # Command lookup delegates to sync local machine
        ls = async_local["ls"]

        # Execution is async
        result = await ls("-la")

    .. versionadded:: 2.0
    """

    def __getitem__(self, cmd: str | LocalPath) -> AsyncLocalCommand:
        """Get an async command by name or path.

        This delegates to local[cmd] to get the sync command, then wraps it.

        Args:
            cmd: Command name (will be looked up in PATH) or LocalPath

        Returns:
            AsyncLocalCommand instance

        Raises:
            CommandNotFound: If command is not found in PATH
        """
        from plumbum.commands.async_ import AsyncLocalCommand

        # Delegate to sync local machine for command lookup
        sync_cmd = local[cmd]
        return AsyncLocalCommand(sync_cmd)

    def __contains__(self, cmd: str) -> bool:
        """Check if a command exists in PATH."""
        # Delegate to sync local machine
        return cmd in local

    @property
    def cwd(self) -> LocalPath:
        """Current working directory."""
        return local.cwd  # type: ignore[no-any-return]

    @property
    def env(self) -> LocalEnv:
        """Environment variables."""
        return local.env

    @staticmethod
    def path(*parts: str) -> LocalPath:
        """Create a LocalPath from parts."""
        return local.path(*parts)


async_local = AsyncLocalMachine()
"""Async version of the local machine singleton.

Use this to access async commands::

    from plumbum import async_local
    # or
    from plumbum.machines.local import async_local

    async def main():
        result = await async_local["ls"]("-la")
        print(result)

.. versionadded:: 2.0
"""

__all__ = [
    "AsyncLocalMachine",
    "LocalCommand",
    "LocalEnv",
    "LocalMachine",
    "PlumbumLocalPopen",
    "async_local",
    "local",
]


def __dir__() -> list[str]:
    return list(__all__)
