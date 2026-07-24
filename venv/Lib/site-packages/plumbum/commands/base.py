from __future__ import annotations

__lazy_modules__ = {
    "contextlib",
    "functools",
    "plumbum.commands.modifiers",
    "plumbum.commands.processes",
    "shlex",
    "tempfile",
    "types",
}

import contextlib
import functools
import shlex
import subprocess
import typing
from subprocess import PIPE, Popen
from tempfile import TemporaryFile
from types import MethodType
from typing import ClassVar

import plumbum.commands.modifiers
from plumbum.commands.processes import (
    ProcessTimedOut,
    _close_streams,
    _terminate_and_reap,
    iter_lines,
    run_proc,
)

if typing.TYPE_CHECKING:
    from collections.abc import Container, Generator, Sequence
    from typing import Any

    from plumbum._compat.typing import Self
    from plumbum.machines.base import BaseMachine, PopenWithAddons

__all__ = (
    "ERROUT",
    "AppendingStdoutRedirection",
    "BaseCommand",
    "BaseRedirection",
    "BoundCommand",
    "BoundEnvCommand",
    "ConcreteCommand",
    "Pipeline",
    "RedirectionError",
    "StderrRedirection",
    "StdinDataRedirection",
    "StdinRedirection",
    "StdoutRedirection",
    "iter_lines",
    "run_proc",
    "shquote",
    "shquote_list",
)


def __dir__() -> list[str]:
    return list(__all__)


class RedirectionError(Exception):
    """Raised when an attempt is made to redirect an process' standard handle,
    which was already redirected to/from a file"""


# ===================================================================================================
# Utilities
# ===================================================================================================
def shquote(text: Any) -> str:
    """Quotes the given text with shell escaping (assumes as syntax similar to ``sh``)"""
    text = str(text)
    return shlex.quote(text)


def shquote_list(seq: Sequence[Any]) -> list[str]:
    return [shquote(item) for item in seq]


# ===================================================================================================
# Commands
# ===================================================================================================
class BaseCommand:
    """Base of all command objects"""

    __slots__ = ("__weakref__", "custom_encoding", "cwd", "env")

    custom_encoding: str | None
    cwd: str | None
    env: dict[str, str] | None

    def __str__(self) -> str:
        return " ".join(self.formulate(1))

    def __or__(self, other: BaseCommand) -> Pipeline:
        """Creates a pipe with the other command"""
        return Pipeline(self, other)

    def __gt__(self, file: Any) -> StdoutRedirection:
        """Redirects the process' stdout to the given file"""
        return StdoutRedirection(self, file)

    def __rshift__(self, file: Any) -> AppendingStdoutRedirection:
        """Redirects the process' stdout to the given file (appending)"""
        return AppendingStdoutRedirection(self, file)

    def __ge__(self, file: Any) -> StderrRedirection:
        """Redirects the process' stderr to the given file"""
        return StderrRedirection(self, file)

    def __lt__(self, file: Any) -> StdinRedirection:
        """Redirects the given file into the process' stdin"""
        return StdinRedirection(self, file)

    def __lshift__(self, data: Any) -> StdinDataRedirection:
        """Redirects the given data into the process' stdin"""
        return StdinDataRedirection(self, data)

    def __getitem__(self, args: Any) -> BoundCommand:
        """Creates a bound-command with the given arguments. Shortcut for
        bound_command."""
        if not isinstance(args, (tuple, list)):
            args = [
                args,
            ]
        return self.bound_command(*args)

    @typing.overload
    def bound_command(self) -> Self: ...

    @typing.overload
    def bound_command(self, *args: Any) -> BoundCommand: ...

    def bound_command(self, *args: Any) -> Self | BoundCommand:
        """Creates a bound-command with the given arguments"""
        if not args:
            return self

        if isinstance(self, BoundCommand):
            return BoundCommand(self.cmd, self.args + list(args))

        return BoundCommand(self, args)

    def __call__(self, *args: Any, **kwargs: Any) -> str:
        """A shortcut for `run(args)`, returning only the process' stdout"""
        return self.run(args, **kwargs)[1]

    def _get_encoding(self) -> str | None:
        raise NotImplementedError()

    def with_env(self, **env: str) -> BaseCommand:
        """Returns a BoundEnvCommand with the given environment variables"""
        if not env:
            return self
        return BoundEnvCommand(self, env=env)

    def with_cwd(self, path: Any) -> BaseCommand:
        """
        Returns a BoundEnvCommand with the specified working directory.
        This overrides a cwd specified in a wrapping `machine.cwd()` context manager.
        """
        if not path:
            return self
        return BoundEnvCommand(self, cwd=path)

    setenv = with_env

    @property
    def machine(self) -> BaseMachine:
        raise NotImplementedError()

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        """Formulates the command into a command-line, i.e., a list of shell-quoted strings
        that can be executed by ``Popen`` or shells.

        :param level: The nesting level of the formulation; it dictates how much shell-quoting
                      (if any) should be performed

        :param args: The arguments passed to this command (a tuple)

        :returns: A list of strings
        """
        raise NotImplementedError()

    def popen(self, args: Sequence[Any] = (), **kwargs: Any) -> PopenWithAddons[str]:
        """Spawns the given command, returning a ``Popen``-like object.

        .. note::

           When processes run in the **background** (either via ``popen`` or
           :class:`& BG <plumbum.commands.modifiers.BG>`), their stdout/stderr pipes might fill up,
           causing them to hang. If you know a process produces output, be sure to consume it
           every once in a while, using a monitoring thread/reactor in the background.
           For more info, see `#48 <https://github.com/tomerfiliba/plumbum/issues/48>`_

        :param args: Any arguments to be passed to the process (a tuple)

        :param kwargs: Any keyword-arguments to be passed to the ``Popen`` constructor

        :returns: A ``Popen``-like object
        """
        raise NotImplementedError()

    def nohup(
        self,
        cwd: str = ".",
        stdout: str = "nohup.out",
        stderr: str | None = None,
        append: bool = True,
    ) -> PopenWithAddons[str]:
        """Runs a command detached."""
        return self.machine.daemonic_popen(self, cwd, stdout, stderr, append)

    @contextlib.contextmanager
    def bgrun(
        self, args: Sequence[Any] = (), **kwargs: Any
    ) -> Generator[PopenWithAddons[str], None, None]:
        """Runs the given command as a context manager, allowing you to create a
        `pipeline <https://en.wikipedia.org/wiki/Pipeline_(computing)>`_ (not in the UNIX sense)
        of programs, parallelizing their work. In other words, instead of running programs
        one after the other, you can start all of them at the same time and wait for them to
        finish. For a more thorough review, see
        `Lightweight Asynchronism <https://tomerfiliba.com/blog/Toying-with-Context-Managers/>`_.

        Example::

            from plumbum.cmd import mkfs

            with mkfs["-t", "ext3", "/dev/sda1"] as p1:
                with mkfs["-t", "ext3", "/dev/sdb1"] as p2:
                    pass

        .. note::

           When processes run in the **background** (either via ``popen`` or
           :data:`& BG <plumbum.commands.modifiers.BG>`), their stdout/stderr pipes might fill up,
           causing them to hang. If you know a process produces output, be sure to consume it
           every once in a while, using a monitoring thread/reactor in the background.
           For more info, see `#48 <https://github.com/tomerfiliba/plumbum/issues/48>`_

        For the arguments, see :func:`run <BaseCommand.run>`.

        :returns: A Popen object, augmented with a ``.run()`` method, which returns a tuple of
                  (return code, stdout, stderr)
        """
        retcode = kwargs.pop("retcode", 0)
        timeout = kwargs.pop("timeout", None)
        p = self.popen(args, **kwargs)
        was_run = [False]

        def cleanup() -> None:
            del p.run  # type: ignore[attr-defined]
            _close_streams(p.stdin, p.stdout, p.stderr)

        def runner() -> tuple[int | None, str | bytes, str | bytes] | None:
            if was_run[0]:
                return None  # already done
            was_run[0] = True
            try:
                return run_proc(p, retcode, timeout)
            finally:
                cleanup()

        p.run = runner  # type: ignore[attr-defined]
        try:
            yield p
        except BaseException:
            # The body raised (including KeyboardInterrupt): waiting for the
            # process to finish could delay the exception indefinitely, and
            # exit-code validation could replace it. Terminate and reap
            # without validation instead.
            if not was_run[0]:
                was_run[0] = True
                try:
                    _terminate_and_reap(p)
                finally:
                    cleanup()
            raise
        # Reap/cleanup on normal exit. ``runner`` is guarded by ``was_run`` so
        # an explicit ``p.run()`` in the body won't run it a second time.
        runner()

    def run(self, args: Sequence[Any] = (), **kwargs: Any) -> tuple[int, str, str]:
        """Runs the given command (equivalent to popen() followed by
        :func:`run_proc <plumbum.commands.processes.run_proc>`). If the exit code of the process does
        not match the expected one, :class:`ProcessExecutionError
        <plumbum.commands.processes.ProcessExecutionError>` is raised.

        :param args: Any arguments to be passed to the process (a tuple)

        :param retcode: The expected return code of this process (defaults to 0).
                        In order to disable exit-code validation, pass ``None``. It may also
                        be a tuple (or any iterable) of expected exit codes.

                        .. note:: this argument must be passed as a keyword argument.

        :param timeout: The maximal amount of time (in seconds) to allow the process to run.
                        ``None`` means no timeout is imposed; otherwise, if the process hasn't
                        terminated after that many seconds, the process will be forcefully
                        terminated an exception will be raised

                        .. note:: this argument must be passed as a keyword argument.

        :param kwargs: Any keyword-arguments to be passed to the ``Popen`` constructor

        :returns: A tuple of (return code, stdout, stderr)
        """
        with self.bgrun(args, **kwargs) as p:
            return p.run()  # type: ignore[attr-defined, no-any-return]

    def run_bg(self, **kwargs: Any) -> plumbum.commands.modifiers.Future:
        """
        Run this command in the background. Uses all arguments from the BG construct
        :py:class: `plumbum.commands.modifiers.BG`
        """
        return self & plumbum.commands.modifiers.BG(**kwargs)

    def run_fg(self, **kwargs: Any) -> None:
        """
        Run this command in the foreground. Uses all arguments from the FG construct
        :py:class: `plumbum.commands.modifiers.FG`
        """
        return self & plumbum.commands.modifiers.FG(**kwargs)

    def run_tee(self, **kwargs: Any) -> tuple[int, str, str]:
        """
        Run this command using the TEE construct. Inherits all arguments from TEE
        :py:class: `plumbum.commands.modifiers.TEE`
        """
        return self & plumbum.commands.modifiers.TEE(**kwargs)

    def run_tf(self, **kwargs: Any) -> bool:
        """
        Run this command using the TF construct. Inherits all arguments from TF
        :py:class: `plumbum.commands.modifiers.TF`
        """
        return self & plumbum.commands.modifiers.TF(**kwargs)

    def run_retcode(self, **kwargs: Any) -> int:
        """
        Run this command using the RETCODE construct. Inherits all arguments from RETCODE
        :py:class: `plumbum.commands.modifiers.RETCODE`
        """
        return self & plumbum.commands.modifiers.RETCODE(**kwargs)

    def run_nohup(self, **kwargs: Any) -> PopenWithAddons[str]:
        """
        Run this command using the NOHUP construct. Inherits all arguments from NOHUP
        :py:class: `plumbum.commands.modifiers.NOHUP`
        """
        return self & plumbum.commands.modifiers.NOHUP(**kwargs)


class BoundCommand(BaseCommand):
    __slots__ = ("args", "cmd")

    cmd: BaseCommand
    args: list[Any]

    def __init__(self, cmd: BaseCommand, args: Sequence[Any]) -> None:
        self.cmd = cmd
        self.args = list(args)

    def __repr__(self) -> str:
        return f"BoundCommand({self.cmd!r}, {self.args!r})"

    def _get_encoding(self) -> str | None:
        return self.cmd._get_encoding()

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        return self.cmd.formulate(level + 1, self.args + list(args))

    @property
    def machine(self) -> BaseMachine:
        return self.cmd.machine

    def popen(
        self, args: Sequence[Any] | str = (), **kwargs: Any
    ) -> PopenWithAddons[str]:
        if isinstance(args, str):
            args = [
                args,
            ]
        return self.cmd.popen(self.args + list(args), **kwargs)


class BoundEnvCommand(BaseCommand):
    __slots__ = ("cmd",)

    cmd: BaseCommand
    env: dict[str, str]
    cwd: str | None

    def __init__(
        self,
        cmd: BaseCommand,
        env: dict[str, str] | None = None,
        cwd: str | None = None,
    ) -> None:
        self.cmd = cmd
        self.env = env or {}
        self.cwd = cwd

    def __repr__(self) -> str:
        return f"BoundEnvCommand({self.cmd!r}, {self.env!r})"

    def _get_encoding(self) -> str | None:
        return self.cmd._get_encoding()

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        return self.cmd.formulate(level, args)

    @property
    def machine(self) -> BaseMachine:
        return self.cmd.machine

    def popen(
        self,
        args: Sequence[Any] = (),
        cwd: str | None = None,
        env: dict[str, str] | None = None,
        **kwargs: Any,
    ) -> PopenWithAddons[str]:
        env = env or {}
        return self.cmd.popen(
            args,
            cwd=self.cwd if cwd is None else cwd,
            env={**self.env, **env},
            **kwargs,
        )


class Pipeline(BaseCommand):
    __slots__ = ("dstcmd", "srccmd")

    srccmd: BaseCommand
    dstcmd: BaseCommand

    def __init__(self, srccmd: BaseCommand, dstcmd: BaseCommand) -> None:
        self.srccmd = srccmd
        self.dstcmd = dstcmd

    def __repr__(self) -> str:
        return f"Pipeline({self.srccmd!r}, {self.dstcmd!r})"

    def _get_encoding(self) -> str | None:
        return self.srccmd._get_encoding() or self.dstcmd._get_encoding()

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        # Call-time args are bound to the *source* command, matching ``popen``
        # (which passes them to ``self.srccmd.popen``); e.g. ``(a | b)("x")``
        # runs ``a x | b``.
        return [
            *self.srccmd.formulate(level + 1, args),
            "|",
            *self.dstcmd.formulate(level + 1),
        ]

    @property
    def machine(self) -> BaseMachine:
        return self.srccmd.machine

    def popen(self, args: Sequence[Any] = (), **kwargs: Any) -> PopenWithAddons[str]:
        src_kwargs = kwargs.copy()
        src_kwargs["stdout"] = PIPE
        if "stdin" in kwargs:
            src_kwargs["stdin"] = kwargs["stdin"]

        srcproc = self.srccmd.popen(args, **src_kwargs)
        kwargs["stdin"] = srcproc.stdout
        dstproc = self.dstcmd.popen(**kwargs)
        # allow p1 to receive a SIGPIPE if p2 exits
        srcproc.stdout.close()  # type: ignore[union-attr]
        if srcproc.stdin and src_kwargs.get("stdin") != PIPE:
            srcproc.stdin.close()
        dstproc.srcproc = srcproc  # type: ignore[attr-defined]

        # monkey-patch .wait() to wait on srcproc as well (it's expected to die when dstproc dies)
        dstproc_wait = dstproc.wait

        @functools.wraps(Popen.wait)
        def wait2(*args: Any, **kwargs: Any) -> int:
            rc_dst = dstproc_wait(*args, **kwargs)
            rc_src = srcproc.wait(*args, **kwargs)
            dstproc.returncode = rc_dst or rc_src
            # The source's stdout was already closed (redirected into dstproc's
            # stdin) above; its stderr/stdin pipes, however, are left open and
            # would leak. Now that both stages have been reaped, close them.
            _close_streams(srcproc.stderr, srcproc.stdin)
            return dstproc.returncode

        dstproc._proc.wait = wait2  # type: ignore[attr-defined]

        dstproc_verify = dstproc.verify

        def verify(
            proc: Any,
            retcode: int | Container[int] | None,
            timeout: float | None,
            stdout: str,
            stderr: str,
        ) -> None:
            # Check if any process in the pipeline has timed out
            # If so, raise ProcessTimedOut instead of ProcessExecutionError
            if getattr(proc.srcproc, "_timed_out", False) or getattr(
                proc, "_timed_out", False
            ):
                raise ProcessTimedOut(
                    f"Process did not terminate within {timeout} seconds",
                    getattr(proc, "argv", None),
                )

            # TODO: right now it's impossible to specify different expected
            # return codes for different stages of the pipeline
            try:
                or_retcode: list[int] | None = [0, *list(retcode)]  # type: ignore[call-overload]
            except TypeError:
                # no-retcode-verification acts "greedily"
                or_retcode = None if retcode is None else [0, retcode]  # type: ignore[list-item]
            proc.srcproc.verify(or_retcode, timeout, stdout, stderr)
            dstproc_verify(retcode, timeout, stdout, stderr)

        dstproc.verify = MethodType(verify, dstproc)  # type: ignore[method-assign]

        dstproc.stdin = srcproc.stdin
        return dstproc


class BaseRedirection(BaseCommand):
    __slots__ = ("cmd", "file")

    # These must be defined by subclasses
    SYM: ClassVar[str]  # pylint: disable=declare-non-slot
    KWARG: ClassVar[str]  # pylint: disable=declare-non-slot
    MODE: ClassVar[str]  # pylint: disable=declare-non-slot

    cmd: BaseCommand
    file: Any

    def __init__(self, cmd: BaseCommand, file: Any) -> None:
        self.cmd = cmd
        self.file = file

    def _get_encoding(self) -> str | None:
        return self.cmd._get_encoding()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.cmd!r}, {self.file!r})"

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        return [
            *self.cmd.formulate(level + 1, args),
            self.SYM,
            shquote(getattr(self.file, "name", self.file)),
        ]

    @property
    def machine(self) -> BaseMachine:
        return self.cmd.machine

    def popen(self, args: Sequence[Any] = (), **kwargs: Any) -> PopenWithAddons[str]:
        from plumbum.path.local import LocalPath
        from plumbum.path.remote import RemotePath

        if self.KWARG in kwargs and kwargs[self.KWARG] not in (PIPE, None):
            raise RedirectionError(f"{self.KWARG} is already redirected")
        if isinstance(self.file, RemotePath):
            raise TypeError("Cannot redirect to/from remote paths")
        if isinstance(self.file, (str, LocalPath)):
            f = kwargs[self.KWARG] = open(str(self.file), self.MODE, encoding="utf-8")
        else:
            kwargs[self.KWARG] = self.file
            f = None
        try:
            return self.cmd.popen(args, **kwargs)
        finally:
            if f:
                f.close()


class StdinRedirection(BaseRedirection):
    __slots__ = ()
    SYM = "<"
    KWARG = "stdin"
    MODE = "r"


class StdoutRedirection(BaseRedirection):
    __slots__ = ()
    SYM = ">"
    KWARG = "stdout"
    MODE = "w"


class AppendingStdoutRedirection(BaseRedirection):
    __slots__ = ()
    SYM = ">>"
    KWARG = "stdout"
    MODE = "a"


class StderrRedirection(BaseRedirection):
    __slots__ = ()
    SYM = "2>"
    KWARG = "stderr"
    MODE = "w"


class _ERROUT(int):
    def __repr__(self) -> str:
        return "ERROUT"

    def __str__(self) -> str:
        return "&1"


ERROUT = _ERROUT(subprocess.STDOUT)


class StdinDataRedirection(BaseCommand):
    __slots__ = ("cmd", "data")

    cmd: BaseCommand
    data: bytes | str

    def __init__(self, cmd: BaseCommand, data: bytes | str) -> None:
        self.cmd = cmd
        self.data = data

    def _get_encoding(self) -> str | None:
        return self.cmd._get_encoding()

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        return [
            f"echo {shquote(self.data)}",
            "|",
            *self.cmd.formulate(level + 1, args),
        ]

    @property
    def machine(self) -> BaseMachine:
        return self.cmd.machine

    def popen(self, args: Sequence[Any] = (), **kwargs: Any) -> PopenWithAddons[str]:
        if kwargs.get("stdin") not in (PIPE, None):
            raise RedirectionError("stdin is already redirected")
        data = self.data
        encoding = self._get_encoding()
        if isinstance(data, str) and encoding is not None:
            data = data.encode(encoding)
        f = TemporaryFile()
        f.write(data)  # type: ignore[arg-type]
        f.seek(0)
        kwargs["stdin"] = f
        try:
            return self.cmd.popen(args, **kwargs)
        finally:
            f.close()


class ConcreteCommand(BaseCommand):
    __slots__ = ("executable",)

    # These must be defined by subclasses
    QUOTE_LEVEL: ClassVar[int]  # pylint: disable=declare-non-slot

    executable: Any

    def __init__(self, executable: Any, encoding: str | None) -> None:
        self.executable = executable
        self.custom_encoding = encoding
        self.cwd = None
        self.env = None

    def __str__(self) -> str:
        return str(self.executable)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.executable})"

    def _get_encoding(self) -> str | None:
        return self.custom_encoding

    def formulate(self, level: int = 0, args: Sequence[Any] = ()) -> list[str]:
        argv = [str(self.executable)]
        for a in args:
            if a is None:
                continue
            if isinstance(a, BaseCommand):
                if level >= self.QUOTE_LEVEL:
                    argv.extend(shquote_list(a.formulate(level + 1)))
                else:
                    argv.extend(a.formulate(level + 1))
            elif isinstance(a, (list, tuple)):
                argv.extend(
                    shquote(b) if level >= self.QUOTE_LEVEL else str(b) for b in a
                )
            else:
                argv.append(shquote(a) if level >= self.QUOTE_LEVEL else str(a))
        return argv

    @property
    def machine(self) -> BaseMachine:
        raise NotImplementedError()

    def popen(self, args: Sequence[Any] = (), **kwargs: Any) -> PopenWithAddons[str]:
        raise NotImplementedError()
