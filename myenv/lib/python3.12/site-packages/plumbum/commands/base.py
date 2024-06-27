import contextlib
import functools
import shlex
import subprocess
from subprocess import PIPE, Popen
from tempfile import TemporaryFile
from types import MethodType
from typing import ClassVar

import plumbum.commands.modifiers
from plumbum.commands.processes import iter_lines, run_proc

__all__ = (
    "iter_lines",
    "run_proc",
    "shquote",
    "shquote_list",
    "RedirectionError",
    "BaseCommand",
    "Pipeline",
    "BaseRedirection",
    "BoundCommand",
    "BoundEnvCommand",
    "ConcreteCommand",
    "ERROUT",
    "StdinRedirection",
    "StdoutRedirection",
    "StderrRedirection",
    "AppendingStdoutRedirection",
    "StdinDataRedirection",
)


class RedirectionError(Exception):
    """Raised when an attempt is made to redirect an process' standard handle,
    which was already redirected to/from a file"""


# ===================================================================================================
# Utilities
# ===================================================================================================
# modified from the stdlib pipes module for windows
_safechars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@%_-+=:,./"
_funnychars = '"`$\\'


def shquote(text):
    """Quotes the given text with shell escaping (assumes as syntax similar to ``sh``)"""
    text = str(text)
    return shlex.quote(text)


def shquote_list(seq):
    return [shquote(item) for item in seq]


# ===================================================================================================
# Commands
# ===================================================================================================
class BaseCommand:
    """Base of all command objects"""

    __slots__ = ("cwd", "env", "custom_encoding", "__weakref__")

    def __str__(self):
        return " ".join(self.formulate())

    def __or__(self, other):
        """Creates a pipe with the other command"""
        return Pipeline(self, other)

    def __gt__(self, file):
        """Redirects the process' stdout to the given file"""
        return StdoutRedirection(self, file)

    def __rshift__(self, file):
        """Redirects the process' stdout to the given file (appending)"""
        return AppendingStdoutRedirection(self, file)

    def __ge__(self, file):
        """Redirects the process' stderr to the given file"""
        return StderrRedirection(self, file)

    def __lt__(self, file):
        """Redirects the given file into the process' stdin"""
        return StdinRedirection(self, file)

    def __lshift__(self, data):
        """Redirects the given data into the process' stdin"""
        return StdinDataRedirection(self, data)

    def __getitem__(self, args):
        """Creates a bound-command with the given arguments. Shortcut for
        bound_command."""
        if not isinstance(args, (tuple, list)):
            args = [
                args,
            ]
        return self.bound_command(*args)

    def bound_command(self, *args):
        """Creates a bound-command with the given arguments"""
        if not args:
            return self

        if isinstance(self, BoundCommand):
            return BoundCommand(self.cmd, self.args + list(args))

        return BoundCommand(self, args)

    def __call__(self, *args, **kwargs):
        """A shortcut for `run(args)`, returning only the process' stdout"""
        return self.run(args, **kwargs)[1]

    def _get_encoding(self):
        raise NotImplementedError()

    def with_env(self, **env):
        """Returns a BoundEnvCommand with the given environment variables"""
        if not env:
            return self
        return BoundEnvCommand(self, env=env)

    def with_cwd(self, path):
        """
        Returns a BoundEnvCommand with the specified working directory.
        This overrides a cwd specified in a wrapping `machine.cwd()` context manager.
        """
        if not path:
            return self
        return BoundEnvCommand(self, cwd=path)

    setenv = with_env

    @property
    def machine(self):
        raise NotImplementedError()

    def formulate(self, level=0, args=()):
        """Formulates the command into a command-line, i.e., a list of shell-quoted strings
        that can be executed by ``Popen`` or shells.

        :param level: The nesting level of the formulation; it dictates how much shell-quoting
                      (if any) should be performed

        :param args: The arguments passed to this command (a tuple)

        :returns: A list of strings
        """
        raise NotImplementedError()

    def popen(self, args=(), **kwargs):
        """Spawns the given command, returning a ``Popen``-like object.

        .. note::

           When processes run in the **background** (either via ``popen`` or
           :class:`& BG <plumbum.commands.BG>`), their stdout/stderr pipes might fill up,
           causing them to hang. If you know a process produces output, be sure to consume it
           every once in a while, using a monitoring thread/reactor in the background.
           For more info, see `#48 <https://github.com/tomerfiliba/plumbum/issues/48>`_

        :param args: Any arguments to be passed to the process (a tuple)

        :param kwargs: Any keyword-arguments to be passed to the ``Popen`` constructor

        :returns: A ``Popen``-like object
        """
        raise NotImplementedError()

    def nohup(self, cwd=".", stdout="nohup.out", stderr=None, append=True):
        """Runs a command detached."""
        return self.machine.daemonic_popen(self, cwd, stdout, stderr, append)

    @contextlib.contextmanager
    def bgrun(self, args=(), **kwargs):
        """Runs the given command as a context manager, allowing you to create a
        `pipeline <http://en.wikipedia.org/wiki/Pipeline_(computing)>`_ (not in the UNIX sense)
        of programs, parallelizing their work. In other words, instead of running programs
        one after the other, you can start all of them at the same time and wait for them to
        finish. For a more thorough review, see
        `Lightweight Asynchronism <http://tomerfiliba.com/blog/Toying-with-Context-Managers/>`_.

        Example::

            from plumbum.cmd import mkfs

            with mkfs["-t", "ext3", "/dev/sda1"] as p1:
                with mkfs["-t", "ext3", "/dev/sdb1"] as p2:
                    pass

        .. note::

           When processes run in the **background** (either via ``popen`` or
           :class:`& BG <plumbum.commands.BG>`), their stdout/stderr pipes might fill up,
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

        def runner():
            if was_run[0]:
                return None  # already done
            was_run[0] = True
            try:
                return run_proc(p, retcode, timeout)
            finally:
                del p.run  # to break cyclic reference p -> cell -> p
                for f in (p.stdin, p.stdout, p.stderr):
                    with contextlib.suppress(Exception):
                        f.close()

        p.run = runner
        yield p
        runner()

    def run(self, args=(), **kwargs):
        """Runs the given command (equivalent to popen() followed by
        :func:`run_proc <plumbum.commands.run_proc>`). If the exit code of the process does
        not match the expected one, :class:`ProcessExecutionError
        <plumbum.commands.ProcessExecutionError>` is raised.

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
            return p.run()

    def _use_modifier(self, modifier, args):
        """
        Applies a modifier to the current object (e.g. FG, NOHUP)
        :param modifier: The modifier class to apply (e.g. FG)
        :param args: A dictionary of arguments to pass to this modifier
        :return:
        """
        modifier_instance = modifier(**args)
        return self & modifier_instance

    def run_bg(self, **kwargs):
        """
        Run this command in the background. Uses all arguments from the BG construct
        :py:class: `plumbum.commands.modifiers.BG`
        """
        return self._use_modifier(plumbum.commands.modifiers.BG, kwargs)

    def run_fg(self, **kwargs):
        """
        Run this command in the foreground. Uses all arguments from the FG construct
        :py:class: `plumbum.commands.modifiers.FG`
        """
        return self._use_modifier(plumbum.commands.modifiers.FG, kwargs)

    def run_tee(self, **kwargs):
        """
        Run this command using the TEE construct. Inherits all arguments from TEE
        :py:class: `plumbum.commands.modifiers.TEE`
        """
        return self._use_modifier(plumbum.commands.modifiers.TEE, kwargs)

    def run_tf(self, **kwargs):
        """
        Run this command using the TF construct. Inherits all arguments from TF
        :py:class: `plumbum.commands.modifiers.TF`
        """
        return self._use_modifier(plumbum.commands.modifiers.TF, kwargs)

    def run_retcode(self, **kwargs):
        """
        Run this command using the RETCODE construct. Inherits all arguments from RETCODE
        :py:class: `plumbum.commands.modifiers.RETCODE`
        """
        return self._use_modifier(plumbum.commands.modifiers.RETCODE, kwargs)

    def run_nohup(self, **kwargs):
        """
        Run this command using the NOHUP construct. Inherits all arguments from NOHUP
        :py:class: `plumbum.commands.modifiers.NOHUP`
        """
        return self._use_modifier(plumbum.commands.modifiers.NOHUP, kwargs)


class BoundCommand(BaseCommand):
    __slots__ = ("cmd", "args")

    def __init__(self, cmd, args):
        self.cmd = cmd
        self.args = list(args)

    def __repr__(self):
        return f"BoundCommand({self.cmd!r}, {self.args!r})"

    def _get_encoding(self):
        return self.cmd._get_encoding()

    def formulate(self, level=0, args=()):
        return self.cmd.formulate(level + 1, self.args + list(args))

    @property
    def machine(self):
        return self.cmd.machine

    def popen(self, args=(), **kwargs):
        if isinstance(args, str):
            args = [
                args,
            ]
        return self.cmd.popen(self.args + list(args), **kwargs)


class BoundEnvCommand(BaseCommand):
    __slots__ = ("cmd",)

    def __init__(self, cmd, env=None, cwd=None):
        self.cmd = cmd
        self.env = env or {}
        self.cwd = cwd

    def __repr__(self):
        return f"BoundEnvCommand({self.cmd!r}, {self.env!r})"

    def _get_encoding(self):
        return self.cmd._get_encoding()

    def formulate(self, level=0, args=()):
        return self.cmd.formulate(level, args)

    @property
    def machine(self):
        return self.cmd.machine

    def popen(self, args=(), cwd=None, env=None, **kwargs):
        env = env or {}
        return self.cmd.popen(
            args,
            cwd=self.cwd if cwd is None else cwd,
            env=dict(self.env, **env),
            **kwargs,
        )


class Pipeline(BaseCommand):
    __slots__ = ("srccmd", "dstcmd")

    def __init__(self, srccmd, dstcmd):
        self.srccmd = srccmd
        self.dstcmd = dstcmd

    def __repr__(self):
        return f"Pipeline({self.srccmd!r}, {self.dstcmd!r})"

    def _get_encoding(self):
        return self.srccmd._get_encoding() or self.dstcmd._get_encoding()

    def formulate(self, level=0, args=()):
        return [
            *self.srccmd.formulate(level + 1),
            "|",
            *self.dstcmd.formulate(level + 1, args),
        ]

    @property
    def machine(self):
        return self.srccmd.machine

    def popen(self, args=(), **kwargs):
        src_kwargs = kwargs.copy()
        src_kwargs["stdout"] = PIPE
        if "stdin" in kwargs:
            src_kwargs["stdin"] = kwargs["stdin"]

        srcproc = self.srccmd.popen(args, **src_kwargs)
        kwargs["stdin"] = srcproc.stdout
        dstproc = self.dstcmd.popen(**kwargs)
        # allow p1 to receive a SIGPIPE if p2 exits
        srcproc.stdout.close()
        if srcproc.stdin and src_kwargs.get("stdin") != PIPE:
            srcproc.stdin.close()
        dstproc.srcproc = srcproc

        # monkey-patch .wait() to wait on srcproc as well (it's expected to die when dstproc dies)
        dstproc_wait = dstproc.wait

        @functools.wraps(Popen.wait)
        def wait2(*args, **kwargs):
            rc_dst = dstproc_wait(*args, **kwargs)
            rc_src = srcproc.wait(*args, **kwargs)
            dstproc.returncode = rc_dst or rc_src
            return dstproc.returncode

        dstproc._proc.wait = wait2

        dstproc_verify = dstproc.verify

        def verify(proc, retcode, timeout, stdout, stderr):
            # TODO: right now it's impossible to specify different expected
            # return codes for different stages of the pipeline
            try:
                or_retcode = [0, *list(retcode)]
            except TypeError:
                # no-retcode-verification acts "greedily"
                or_retcode = None if retcode is None else [0, retcode]
            proc.srcproc.verify(or_retcode, timeout, stdout, stderr)
            dstproc_verify(retcode, timeout, stdout, stderr)

        dstproc.verify = MethodType(verify, dstproc)

        dstproc.stdin = srcproc.stdin
        return dstproc


class BaseRedirection(BaseCommand):
    __slots__ = ("cmd", "file")
    SYM: ClassVar[str]
    KWARG: ClassVar[str]
    MODE: ClassVar[str]

    def __init__(self, cmd, file):
        self.cmd = cmd
        self.file = file

    def _get_encoding(self):
        return self.cmd._get_encoding()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.cmd!r}, {self.file!r})"

    def formulate(self, level=0, args=()):
        return [
            *self.cmd.formulate(level + 1, args),
            self.SYM,
            shquote(getattr(self.file, "name", self.file)),
        ]

    @property
    def machine(self):
        return self.cmd.machine

    def popen(self, args=(), **kwargs):
        from plumbum.machines.local import LocalPath
        from plumbum.machines.remote import RemotePath

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
    def __repr__(self):
        return "ERROUT"

    def __str__(self):
        return "&1"


ERROUT = _ERROUT(subprocess.STDOUT)


class StdinDataRedirection(BaseCommand):
    __slots__ = ("cmd", "data")
    CHUNK_SIZE = 16000

    def __init__(self, cmd, data):
        self.cmd = cmd
        self.data = data

    def _get_encoding(self):
        return self.cmd._get_encoding()

    def formulate(self, level=0, args=()):
        return [
            f"echo {shquote(self.data)}",
            "|",
            *self.cmd.formulate(level + 1, args),
        ]

    @property
    def machine(self):
        return self.cmd.machine

    def popen(self, args=(), **kwargs):
        if kwargs.get("stdin") not in (PIPE, None):
            raise RedirectionError("stdin is already redirected")
        data = self.data
        if isinstance(data, str) and self._get_encoding() is not None:
            data = data.encode(self._get_encoding())
        f = TemporaryFile()
        while data:
            chunk = data[: self.CHUNK_SIZE]
            f.write(chunk)
            data = data[self.CHUNK_SIZE :]
        f.seek(0)
        kwargs["stdin"] = f
        # try:
        return self.cmd.popen(args, **kwargs)
        # finally:
        #    f.close()


class ConcreteCommand(BaseCommand):
    QUOTE_LEVEL: ClassVar[int]
    __slots__ = ("executable",)

    def __init__(self, executable, encoding):
        self.executable = executable
        self.custom_encoding = encoding
        self.cwd = None
        self.env = None

    def __str__(self):
        return str(self.executable)

    def __repr__(self):
        return f"{type(self).__name__}({self.executable})"

    def _get_encoding(self):
        return self.custom_encoding

    def formulate(self, level=0, args=()):
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
        # if self.custom_encoding:
        #    argv = [a.encode(self.custom_encoding) for a in argv if isinstance(a, str)]
        return argv

    @property
    def machine(self):
        raise NotImplementedError()

    def popen(self, args=(), **kwargs):
        raise NotImplementedError()
