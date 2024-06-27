# Copyright (C) 2008, 2009 Michael Trier (mtrier@gmail.com) and contributors
#
# This module is part of GitPython and is released under the
# 3-Clause BSD License: https://opensource.org/license/bsd-3-clause/

from __future__ import annotations

__all__ = ["GitMeta", "Git"]

import contextlib
import io
import itertools
import logging
import os
import re
import signal
import subprocess
from subprocess import DEVNULL, PIPE, Popen
import sys
from textwrap import dedent
import threading
import warnings

from git.compat import defenc, force_bytes, safe_decode
from git.exc import (
    CommandError,
    GitCommandError,
    GitCommandNotFound,
    UnsafeOptionError,
    UnsafeProtocolError,
)
from git.util import (
    cygpath,
    expand_path,
    is_cygwin_git,
    patch_env,
    remove_password_if_present,
    stream_copy,
)

# typing ---------------------------------------------------------------------------

from typing import (
    Any,
    AnyStr,
    BinaryIO,
    Callable,
    Dict,
    IO,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    TYPE_CHECKING,
    TextIO,
    Tuple,
    Union,
    cast,
    overload,
)

from git.types import Literal, PathLike, TBD

if TYPE_CHECKING:
    from git.diff import DiffIndex
    from git.repo.base import Repo

# ---------------------------------------------------------------------------------

execute_kwargs = {
    "istream",
    "with_extended_output",
    "with_exceptions",
    "as_process",
    "output_stream",
    "stdout_as_string",
    "kill_after_timeout",
    "with_stdout",
    "universal_newlines",
    "shell",
    "env",
    "max_chunk_size",
    "strip_newline_in_stdout",
}

_logger = logging.getLogger(__name__)


# ==============================================================================
## @name Utilities
# ------------------------------------------------------------------------------
# Documentation
## @{


def handle_process_output(
    process: "Git.AutoInterrupt" | Popen,
    stdout_handler: Union[
        None,
        Callable[[AnyStr], None],
        Callable[[List[AnyStr]], None],
        Callable[[bytes, "Repo", "DiffIndex"], None],
    ],
    stderr_handler: Union[None, Callable[[AnyStr], None], Callable[[List[AnyStr]], None]],
    finalizer: Union[None, Callable[[Union[Popen, "Git.AutoInterrupt"]], None]] = None,
    decode_streams: bool = True,
    kill_after_timeout: Union[None, float] = None,
) -> None:
    R"""Register for notifications to learn that process output is ready to read, and
    dispatch lines to the respective line handlers.

    This function returns once the finalizer returns.

    :param process:
        :class:`subprocess.Popen` instance.

    :param stdout_handler:
        f(stdout_line_string), or ``None``.

    :param stderr_handler:
        f(stderr_line_string), or ``None``.

    :param finalizer:
        f(proc) - wait for proc to finish.

    :param decode_streams:
        Assume stdout/stderr streams are binary and decode them before pushing their
        contents to handlers.

        This defaults to ``True``. Set it to ``False`` if:

        - ``universal_newlines == True``, as then streams are in text mode, or
        - decoding must happen later, such as for :class:`~git.diff.Diff`\s.

    :param kill_after_timeout:
        :class:`float` or ``None``, Default = ``None``

        To specify a timeout in seconds for the git command, after which the process
        should be killed.
    """

    # Use 2 "pump" threads and wait for both to finish.
    def pump_stream(
        cmdline: List[str],
        name: str,
        stream: Union[BinaryIO, TextIO],
        is_decode: bool,
        handler: Union[None, Callable[[Union[bytes, str]], None]],
    ) -> None:
        try:
            for line in stream:
                if handler:
                    if is_decode:
                        assert isinstance(line, bytes)
                        line_str = line.decode(defenc)
                        handler(line_str)
                    else:
                        handler(line)

        except Exception as ex:
            _logger.error(f"Pumping {name!r} of cmd({remove_password_if_present(cmdline)}) failed due to: {ex!r}")
            if "I/O operation on closed file" not in str(ex):
                # Only reraise if the error was not due to the stream closing.
                raise CommandError([f"<{name}-pump>"] + remove_password_if_present(cmdline), ex) from ex
        finally:
            stream.close()

    if hasattr(process, "proc"):
        process = cast("Git.AutoInterrupt", process)
        cmdline: str | Tuple[str, ...] | List[str] = getattr(process.proc, "args", "")
        p_stdout = process.proc.stdout if process.proc else None
        p_stderr = process.proc.stderr if process.proc else None
    else:
        process = cast(Popen, process)  # type: ignore[redundant-cast]
        cmdline = getattr(process, "args", "")
        p_stdout = process.stdout
        p_stderr = process.stderr

    if not isinstance(cmdline, (tuple, list)):
        cmdline = cmdline.split()

    pumps: List[Tuple[str, IO, Callable[..., None] | None]] = []
    if p_stdout:
        pumps.append(("stdout", p_stdout, stdout_handler))
    if p_stderr:
        pumps.append(("stderr", p_stderr, stderr_handler))

    threads: List[threading.Thread] = []

    for name, stream, handler in pumps:
        t = threading.Thread(target=pump_stream, args=(cmdline, name, stream, decode_streams, handler))
        t.daemon = True
        t.start()
        threads.append(t)

    # FIXME: Why join? Will block if stdin needs feeding...
    for t in threads:
        t.join(timeout=kill_after_timeout)
        if t.is_alive():
            if isinstance(process, Git.AutoInterrupt):
                process._terminate()
            else:  # Don't want to deal with the other case.
                raise RuntimeError(
                    "Thread join() timed out in cmd.handle_process_output()."
                    f" kill_after_timeout={kill_after_timeout} seconds"
                )
            if stderr_handler:
                error_str: Union[str, bytes] = (
                    "error: process killed because it timed out." f" kill_after_timeout={kill_after_timeout} seconds"
                )
                if not decode_streams and isinstance(p_stderr, BinaryIO):
                    # Assume stderr_handler needs binary input.
                    error_str = cast(str, error_str)
                    error_str = error_str.encode()
                # We ignore typing on the next line because mypy does not like the way
                # we inferred that stderr takes str or bytes.
                stderr_handler(error_str)  # type: ignore[arg-type]

    if finalizer:
        finalizer(process)


safer_popen: Callable[..., Popen]

if sys.platform == "win32":

    def _safer_popen_windows(
        command: Union[str, Sequence[Any]],
        *,
        shell: bool = False,
        env: Optional[Mapping[str, str]] = None,
        **kwargs: Any,
    ) -> Popen:
        """Call :class:`subprocess.Popen` on Windows but don't include a CWD in the
        search.

        This avoids an untrusted search path condition where a file like ``git.exe`` in
        a malicious repository would be run when GitPython operates on the repository.
        The process using GitPython may have an untrusted repository's working tree as
        its current working directory. Some operations may temporarily change to that
        directory before running a subprocess. In addition, while by default GitPython
        does not run external commands with a shell, it can be made to do so, in which
        case the CWD of the subprocess, which GitPython usually sets to a repository
        working tree, can itself be searched automatically by the shell. This wrapper
        covers all those cases.

        :note:
            This currently works by setting the
            :envvar:`NoDefaultCurrentDirectoryInExePath` environment variable during
            subprocess creation. It also takes care of passing Windows-specific process
            creation flags, but that is unrelated to path search.

        :note:
            The current implementation contains a race condition on :attr:`os.environ`.
            GitPython isn't thread-safe, but a program using it on one thread should
            ideally be able to mutate :attr:`os.environ` on another, without
            unpredictable results. See comments in:
            https://github.com/gitpython-developers/GitPython/pull/1650
        """
        # CREATE_NEW_PROCESS_GROUP is needed for some ways of killing it afterwards.
        # https://docs.python.org/3/library/subprocess.html#subprocess.Popen.send_signal
        # https://docs.python.org/3/library/subprocess.html#subprocess.CREATE_NEW_PROCESS_GROUP
        creationflags = subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP

        # When using a shell, the shell is the direct subprocess, so the variable must
        # be set in its environment, to affect its search behavior.
        if shell:
            # The original may be immutable, or the caller may reuse it. Mutate a copy.
            env = {} if env is None else dict(env)
            env["NoDefaultCurrentDirectoryInExePath"] = "1"  # The "1" can be an value.

        # When not using a shell, the current process does the search in a
        # CreateProcessW API call, so the variable must be set in our environment. With
        # a shell, that's unnecessary if https://github.com/python/cpython/issues/101283
        # is patched. In Python versions where it is unpatched, and in the rare case the
        # ComSpec environment variable is unset, the search for the shell itself is
        # unsafe. Setting NoDefaultCurrentDirectoryInExePath in all cases, as done here,
        # is simpler and protects against that. (As above, the "1" can be any value.)
        with patch_env("NoDefaultCurrentDirectoryInExePath", "1"):
            return Popen(
                command,
                shell=shell,
                env=env,
                creationflags=creationflags,
                **kwargs,
            )

    safer_popen = _safer_popen_windows
else:
    safer_popen = Popen


def dashify(string: str) -> str:
    return string.replace("_", "-")


def slots_to_dict(self: "Git", exclude: Sequence[str] = ()) -> Dict[str, Any]:
    return {s: getattr(self, s) for s in self.__slots__ if s not in exclude}


def dict_to_slots_and__excluded_are_none(self: object, d: Mapping[str, Any], excluded: Sequence[str] = ()) -> None:
    for k, v in d.items():
        setattr(self, k, v)
    for k in excluded:
        setattr(self, k, None)


## -- End Utilities -- @}

_USE_SHELL_DEFAULT_MESSAGE = (
    "Git.USE_SHELL is deprecated, because only its default value of False is safe. "
    "It will be removed in a future release."
)

_USE_SHELL_DANGER_MESSAGE = (
    "Setting Git.USE_SHELL to True is unsafe and insecure, as the effect of special "
    "shell syntax cannot usually be accounted for. This can result in a command "
    "injection vulnerability and arbitrary code execution. Git.USE_SHELL is deprecated "
    "and will be removed in a future release."
)


def _warn_use_shell(extra_danger: bool) -> None:
    warnings.warn(
        _USE_SHELL_DANGER_MESSAGE if extra_danger else _USE_SHELL_DEFAULT_MESSAGE,
        DeprecationWarning,
        stacklevel=3,
    )


class _GitMeta(type):
    """Metaclass for :class:`Git`.

    This helps issue :class:`DeprecationWarning` if :attr:`Git.USE_SHELL` is used.
    """

    def __getattribute(cls, name: str) -> Any:
        if name == "USE_SHELL":
            _warn_use_shell(False)
        return super().__getattribute__(name)

    def __setattr(cls, name: str, value: Any) -> Any:
        if name == "USE_SHELL":
            _warn_use_shell(value)
        super().__setattr__(name, value)

    if not TYPE_CHECKING:
        # To preserve static checking for undefined/misspelled attributes while letting
        # the methods' bodies be type-checked, these are defined as non-special methods,
        # then bound to special names out of view of static type checkers. (The original
        # names invoke name mangling (leading "__") to avoid confusion in other scopes.)
        __getattribute__ = __getattribute
        __setattr__ = __setattr


GitMeta = _GitMeta
"""Alias of :class:`Git`'s metaclass, whether it is :class:`type` or a custom metaclass.

Whether the :class:`Git` class has the default :class:`type` as its metaclass or uses a
custom metaclass is not documented and may change at any time. This statically checkable
metaclass alias is equivalent at runtime to ``type(Git)``. This should almost never be
used. Code that benefits from it is likely to be remain brittle even if it is used.

In view of the :class:`Git` class's intended use and :class:`Git` objects' dynamic
callable attributes representing git subcommands, it rarely makes sense to inherit from
:class:`Git` at all. Using :class:`Git` in multiple inheritance can be especially tricky
to do correctly. Attempting uses of :class:`Git` where its metaclass is relevant, such
as when a sibling class has an unrelated metaclass and a shared lower bound metaclass
might have to be introduced to solve a metaclass conflict, is not recommended.

:note:
    The correct static type of the :class:`Git` class itself, and any subclasses, is
    ``Type[Git]``. (This can be written as ``type[Git]`` in Python 3.9 later.)

    :class:`GitMeta` should never be used in any annotation where ``Type[Git]`` is
    intended or otherwise possible to use. This alias is truly only for very rare and
    inherently precarious situations where it is necessary to deal with the metaclass
    explicitly.
"""


class Git(metaclass=_GitMeta):
    """The Git class manages communication with the Git binary.

    It provides a convenient interface to calling the Git binary, such as in::

     g = Git( git_dir )
     g.init()                   # calls 'git init' program
     rval = g.ls_files()        # calls 'git ls-files' program

    Debugging:

    * Set the :envvar:`GIT_PYTHON_TRACE` environment variable to print each invocation
      of the command to stdout.
    * Set its value to ``full`` to see details about the returned values.
    """

    __slots__ = (
        "_working_dir",
        "cat_file_all",
        "cat_file_header",
        "_version_info",
        "_version_info_token",
        "_git_options",
        "_persistent_git_options",
        "_environment",
    )

    _excluded_ = (
        "cat_file_all",
        "cat_file_header",
        "_version_info",
        "_version_info_token",
    )

    re_unsafe_protocol = re.compile(r"(.+)::.+")

    def __getstate__(self) -> Dict[str, Any]:
        return slots_to_dict(self, exclude=self._excluded_)

    def __setstate__(self, d: Dict[str, Any]) -> None:
        dict_to_slots_and__excluded_are_none(self, d, excluded=self._excluded_)

    # CONFIGURATION

    git_exec_name = "git"
    """Default git command that should work on Linux, Windows, and other systems."""

    GIT_PYTHON_TRACE = os.environ.get("GIT_PYTHON_TRACE", False)
    """Enables debugging of GitPython's git commands."""

    USE_SHELL: bool = False
    """Deprecated. If set to ``True``, a shell will be used when executing git commands.

    Code that uses ``USE_SHELL = True`` or that passes ``shell=True`` to any GitPython
    functions should be updated to use the default value of ``False`` instead. ``True``
    is unsafe unless the effect of syntax treated specially by the shell is fully
    considered and accounted for, which is not possible under most circumstances. As
    detailed below, it is also no longer needed, even where it had been in the past.

    It is in many if not most cases a command injection vulnerability for an application
    to set :attr:`USE_SHELL` to ``True``. Any attacker who can cause a specially crafted
    fragment of text to make its way into any part of any argument to any git command
    (including paths, branch names, etc.) can cause the shell to read and write
    arbitrary files and execute arbitrary commands. Innocent input may also accidentally
    contain special shell syntax, leading to inadvertent malfunctions.

    In addition, how a value of ``True`` interacts with some aspects of GitPython's
    operation is not precisely specified and may change without warning, even before
    GitPython 4.0.0 when :attr:`USE_SHELL` may be removed. This includes:

    * Whether or how GitPython automatically customizes the shell environment.

    * Whether, outside of Windows (where :class:`subprocess.Popen` supports lists of
      separate arguments even when ``shell=True``), this can be used with any GitPython
      functionality other than direct calls to the :meth:`execute` method.

    * Whether any GitPython feature that runs git commands ever attempts to partially
      sanitize data a shell may treat specially. Currently this is not done.

    Prior to GitPython 2.0.8, this had a narrow purpose in suppressing console windows
    in graphical Windows applications. In 2.0.8 and higher, it provides no benefit, as
    GitPython solves that problem more robustly and safely by using the
    ``CREATE_NO_WINDOW`` process creation flag on Windows.

    Because Windows path search differs subtly based on whether a shell is used, in rare
    cases changing this from ``True`` to ``False`` may keep an unusual git "executable",
    such as a batch file, from being found. To fix this, set the command name or full
    path in the :envvar:`GIT_PYTHON_GIT_EXECUTABLE` environment variable or pass the
    full path to :func:`git.refresh` (or invoke the script using a ``.exe`` shim).

    Further reading:

    * :meth:`Git.execute` (on the ``shell`` parameter).
    * https://github.com/gitpython-developers/GitPython/commit/0d9390866f9ce42870d3116094cd49e0019a970a
    * https://learn.microsoft.com/en-us/windows/win32/procthread/process-creation-flags
    * https://github.com/python/cpython/issues/91558#issuecomment-1100942950
    * https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-createprocessw
    """

    _git_exec_env_var = "GIT_PYTHON_GIT_EXECUTABLE"
    _refresh_env_var = "GIT_PYTHON_REFRESH"

    GIT_PYTHON_GIT_EXECUTABLE = None
    """Provide the full path to the git executable. Otherwise it assumes git is in the
    executable search path.

    :note:
        The git executable is actually found during the refresh step in the top level
        ``__init__``. It can also be changed by explicitly calling :func:`git.refresh`.
    """

    _refresh_token = object()  # Since None would match an initial _version_info_token.

    @classmethod
    def refresh(cls, path: Union[None, PathLike] = None) -> bool:
        """Update information about the git executable :class:`Git` objects will use.

        Called by the :func:`git.refresh` function in the top level ``__init__``.

        :param path:
            Optional path to the git executable. If not absolute, it is resolved
            immediately, relative to the current directory. (See note below.)

        :note:
            The top-level :func:`git.refresh` should be preferred because it calls this
            method and may also update other state accordingly.

        :note:
            There are three different ways to specify the command that refreshing causes
            to be used for git:

            1. Pass no `path` argument and do not set the
               :envvar:`GIT_PYTHON_GIT_EXECUTABLE` environment variable. The command
               name ``git`` is used. It is looked up in a path search by the system, in
               each command run (roughly similar to how git is found when running
               ``git`` commands manually). This is usually the desired behavior.

            2. Pass no `path` argument but set the :envvar:`GIT_PYTHON_GIT_EXECUTABLE`
               environment variable. The command given as the value of that variable is
               used. This may be a simple command or an arbitrary path. It is looked up
               in each command run. Setting :envvar:`GIT_PYTHON_GIT_EXECUTABLE` to
               ``git`` has the same effect as not setting it.

            3. Pass a `path` argument. This path, if not absolute, is immediately
               resolved, relative to the current directory. This resolution occurs at
               the time of the refresh. When git commands are run, they are run using
               that previously resolved path. If a `path` argument is passed, the
               :envvar:`GIT_PYTHON_GIT_EXECUTABLE` environment variable is not
               consulted.

        :note:
            Refreshing always sets the :attr:`Git.GIT_PYTHON_GIT_EXECUTABLE` class
            attribute, which can be read on the :class:`Git` class or any of its
            instances to check what command is used to run git. This attribute should
            not be confused with the related :envvar:`GIT_PYTHON_GIT_EXECUTABLE`
            environment variable. The class attribute is set no matter how refreshing is
            performed.
        """
        # Discern which path to refresh with.
        if path is not None:
            new_git = os.path.expanduser(path)
            new_git = os.path.abspath(new_git)
        else:
            new_git = os.environ.get(cls._git_exec_env_var, cls.git_exec_name)

        # Keep track of the old and new git executable path.
        old_git = cls.GIT_PYTHON_GIT_EXECUTABLE
        old_refresh_token = cls._refresh_token
        cls.GIT_PYTHON_GIT_EXECUTABLE = new_git
        cls._refresh_token = object()

        # Test if the new git executable path is valid. A GitCommandNotFound error is
        # raised by us. A PermissionError is raised if the git executable cannot be
        # executed for whatever reason.
        has_git = False
        try:
            cls().version()
            has_git = True
        except (GitCommandNotFound, PermissionError):
            pass

        # Warn or raise exception if test failed.
        if not has_git:
            err = (
                dedent(
                    """\
                Bad git executable.
                The git executable must be specified in one of the following ways:
                    - be included in your $PATH
                    - be set via $%s
                    - explicitly set via git.refresh(<full-path-to-git-executable>)
                """
                )
                % cls._git_exec_env_var
            )

            # Revert to whatever the old_git was.
            cls.GIT_PYTHON_GIT_EXECUTABLE = old_git
            cls._refresh_token = old_refresh_token

            if old_git is None:
                # On the first refresh (when GIT_PYTHON_GIT_EXECUTABLE is None) we only
                # are quiet, warn, or error depending on the GIT_PYTHON_REFRESH value.

                # Determine what the user wants to happen during the initial refresh. We
                # expect GIT_PYTHON_REFRESH to either be unset or be one of the
                # following values:
                #
                #   0|q|quiet|s|silence|silent|n|none
                #   1|w|warn|warning|l|log
                #   2|r|raise|e|error|exception

                mode = os.environ.get(cls._refresh_env_var, "raise").lower()

                quiet = ["quiet", "q", "silence", "s", "silent", "none", "n", "0"]
                warn = ["warn", "w", "warning", "log", "l", "1"]
                error = ["error", "e", "exception", "raise", "r", "2"]

                if mode in quiet:
                    pass
                elif mode in warn or mode in error:
                    err = dedent(
                        """\
                        %s
                        All git commands will error until this is rectified.

                        This initial message can be silenced or aggravated in the future by setting the
                        $%s environment variable. Use one of the following values:
                            - %s: for no message or exception
                            - %s: for a warning message (logging level CRITICAL, displayed by default)
                            - %s: for a raised exception

                        Example:
                            export %s=%s
                        """
                    ) % (
                        err,
                        cls._refresh_env_var,
                        "|".join(quiet),
                        "|".join(warn),
                        "|".join(error),
                        cls._refresh_env_var,
                        quiet[0],
                    )

                    if mode in warn:
                        _logger.critical(err)
                    else:
                        raise ImportError(err)
                else:
                    err = dedent(
                        """\
                        %s environment variable has been set but it has been set with an invalid value.

                        Use only the following values:
                            - %s: for no message or exception
                            - %s: for a warning message (logging level CRITICAL, displayed by default)
                            - %s: for a raised exception
                        """
                    ) % (
                        cls._refresh_env_var,
                        "|".join(quiet),
                        "|".join(warn),
                        "|".join(error),
                    )
                    raise ImportError(err)

                # We get here if this was the initial refresh and the refresh mode was
                # not error. Go ahead and set the GIT_PYTHON_GIT_EXECUTABLE such that we
                # discern the difference between the first refresh at import time
                # and subsequent calls to git.refresh or this refresh method.
                cls.GIT_PYTHON_GIT_EXECUTABLE = cls.git_exec_name
            else:
                # After the first refresh (when GIT_PYTHON_GIT_EXECUTABLE is no longer
                # None) we raise an exception.
                raise GitCommandNotFound(new_git, err)

        return has_git

    @classmethod
    def is_cygwin(cls) -> bool:
        return is_cygwin_git(cls.GIT_PYTHON_GIT_EXECUTABLE)

    @overload
    @classmethod
    def polish_url(cls, url: str, is_cygwin: Literal[False] = ...) -> str: ...

    @overload
    @classmethod
    def polish_url(cls, url: str, is_cygwin: Union[None, bool] = None) -> str: ...

    @classmethod
    def polish_url(cls, url: str, is_cygwin: Union[None, bool] = None) -> PathLike:
        """Remove any backslashes from URLs to be written in config files.

        Windows might create config files containing paths with backslashes, but git
        stops liking them as it will escape the backslashes. Hence we undo the escaping
        just to be sure.
        """
        if is_cygwin is None:
            is_cygwin = cls.is_cygwin()

        if is_cygwin:
            url = cygpath(url)
        else:
            url = os.path.expandvars(url)
            if url.startswith("~"):
                url = os.path.expanduser(url)
            url = url.replace("\\\\", "\\").replace("\\", "/")
        return url

    @classmethod
    def check_unsafe_protocols(cls, url: str) -> None:
        """Check for unsafe protocols.

        Apart from the usual protocols (http, git, ssh), Git allows "remote helpers"
        that have the form ``<transport>::<address>``. One of these helpers (``ext::``)
        can be used to invoke any arbitrary command.

        See:

        - https://git-scm.com/docs/gitremote-helpers
        - https://git-scm.com/docs/git-remote-ext
        """
        match = cls.re_unsafe_protocol.match(url)
        if match:
            protocol = match.group(1)
            raise UnsafeProtocolError(
                f"The `{protocol}::` protocol looks suspicious, use `allow_unsafe_protocols=True` to allow it."
            )

    @classmethod
    def check_unsafe_options(cls, options: List[str], unsafe_options: List[str]) -> None:
        """Check for unsafe options.

        Some options that are passed to ``git <command>`` can be used to execute
        arbitrary commands. These are blocked by default.
        """
        # Options can be of the form `foo`, `--foo bar`, or `--foo=bar`, so we need to
        # check if they start with "--foo" or if they are equal to "foo".
        bare_unsafe_options = [option.lstrip("-") for option in unsafe_options]
        for option in options:
            for unsafe_option, bare_option in zip(unsafe_options, bare_unsafe_options):
                if option.startswith(unsafe_option) or option == bare_option:
                    raise UnsafeOptionError(
                        f"{unsafe_option} is not allowed, use `allow_unsafe_options=True` to allow it."
                    )

    class AutoInterrupt:
        """Process wrapper that terminates the wrapped process on finalization.

        This kills/interrupts the stored process instance once this instance goes out of
        scope. It is used to prevent processes piling up in case iterators stop reading.

        All attributes are wired through to the contained process object.

        The wait method is overridden to perform automatic status code checking and
        possibly raise.
        """

        __slots__ = ("proc", "args", "status")

        # If this is non-zero it will override any status code during _terminate, used
        # to prevent race conditions in testing.
        _status_code_if_terminate: int = 0

        def __init__(self, proc: Union[None, subprocess.Popen], args: Any) -> None:
            self.proc = proc
            self.args = args
            self.status: Union[int, None] = None

        def _terminate(self) -> None:
            """Terminate the underlying process."""
            if self.proc is None:
                return

            proc = self.proc
            self.proc = None
            if proc.stdin:
                proc.stdin.close()
            if proc.stdout:
                proc.stdout.close()
            if proc.stderr:
                proc.stderr.close()
            # Did the process finish already so we have a return code?
            try:
                if proc.poll() is not None:
                    self.status = self._status_code_if_terminate or proc.poll()
                    return
            except OSError as ex:
                _logger.info("Ignored error after process had died: %r", ex)

            # It can be that nothing really exists anymore...
            if os is None or getattr(os, "kill", None) is None:
                return

            # Try to kill it.
            try:
                proc.terminate()
                status = proc.wait()  # Ensure the process goes away.

                self.status = self._status_code_if_terminate or status
            except OSError as ex:
                _logger.info("Ignored error after process had died: %r", ex)
            # END exception handling

        def __del__(self) -> None:
            self._terminate()

        def __getattr__(self, attr: str) -> Any:
            return getattr(self.proc, attr)

        # TODO: Bad choice to mimic `proc.wait()` but with different args.
        def wait(self, stderr: Union[None, str, bytes] = b"") -> int:
            """Wait for the process and return its status code.

            :param stderr:
                Previously read value of stderr, in case stderr is already closed.

            :warn:
                May deadlock if output or error pipes are used and not handled
                separately.

            :raise git.exc.GitCommandError:
                If the return status is not 0.
            """
            if stderr is None:
                stderr_b = b""
            stderr_b = force_bytes(data=stderr, encoding="utf-8")
            status: Union[int, None]
            if self.proc is not None:
                status = self.proc.wait()
                p_stderr = self.proc.stderr
            else:  # Assume the underlying proc was killed earlier or never existed.
                status = self.status
                p_stderr = None

            def read_all_from_possibly_closed_stream(stream: Union[IO[bytes], None]) -> bytes:
                if stream:
                    try:
                        return stderr_b + force_bytes(stream.read())
                    except (OSError, ValueError):
                        return stderr_b or b""
                else:
                    return stderr_b or b""

            # END status handling

            if status != 0:
                errstr = read_all_from_possibly_closed_stream(p_stderr)
                _logger.debug("AutoInterrupt wait stderr: %r" % (errstr,))
                raise GitCommandError(remove_password_if_present(self.args), status, errstr)
            return status

    # END auto interrupt

    class CatFileContentStream:
        """Object representing a sized read-only stream returning the contents of
        an object.

        This behaves like a stream, but counts the data read and simulates an empty
        stream once our sized content region is empty.

        If not all data are read to the end of the object's lifetime, we read the
        rest to ensure the underlying stream continues to work.
        """

        __slots__ = ("_stream", "_nbr", "_size")

        def __init__(self, size: int, stream: IO[bytes]) -> None:
            self._stream = stream
            self._size = size
            self._nbr = 0  # Number of bytes read.

            # Special case: If the object is empty, has null bytes, get the final
            # newline right away.
            if size == 0:
                stream.read(1)
            # END handle empty streams

        def read(self, size: int = -1) -> bytes:
            bytes_left = self._size - self._nbr
            if bytes_left == 0:
                return b""
            if size > -1:
                # Ensure we don't try to read past our limit.
                size = min(bytes_left, size)
            else:
                # They try to read all, make sure it's not more than what remains.
                size = bytes_left
            # END check early depletion
            data = self._stream.read(size)
            self._nbr += len(data)

            # Check for depletion, read our final byte to make the stream usable by
            # others.
            if self._size - self._nbr == 0:
                self._stream.read(1)  # final newline
            # END finish reading
            return data

        def readline(self, size: int = -1) -> bytes:
            if self._nbr == self._size:
                return b""

            # Clamp size to lowest allowed value.
            bytes_left = self._size - self._nbr
            if size > -1:
                size = min(bytes_left, size)
            else:
                size = bytes_left
            # END handle size

            data = self._stream.readline(size)
            self._nbr += len(data)

            # Handle final byte.
            if self._size - self._nbr == 0:
                self._stream.read(1)
            # END finish reading

            return data

        def readlines(self, size: int = -1) -> List[bytes]:
            if self._nbr == self._size:
                return []

            # Leave all additional logic to our readline method, we just check the size.
            out = []
            nbr = 0
            while True:
                line = self.readline()
                if not line:
                    break
                out.append(line)
                if size > -1:
                    nbr += len(line)
                    if nbr > size:
                        break
                # END handle size constraint
            # END readline loop
            return out

        # skipcq: PYL-E0301
        def __iter__(self) -> "Git.CatFileContentStream":
            return self

        def __next__(self) -> bytes:
            line = self.readline()
            if not line:
                raise StopIteration

            return line

        next = __next__

        def __del__(self) -> None:
            bytes_left = self._size - self._nbr
            if bytes_left:
                # Read and discard - seeking is impossible within a stream.
                # This includes any terminating newline.
                self._stream.read(bytes_left + 1)
            # END handle incomplete read

    def __init__(self, working_dir: Union[None, PathLike] = None) -> None:
        """Initialize this instance with:

        :param working_dir:
            Git directory we should work in. If ``None``, we always work in the current
            directory as returned by :func:`os.getcwd`.
            This is meant to be the working tree directory if available, or the
            ``.git`` directory in case of bare repositories.
        """
        super().__init__()
        self._working_dir = expand_path(working_dir)
        self._git_options: Union[List[str], Tuple[str, ...]] = ()
        self._persistent_git_options: List[str] = []

        # Extra environment variables to pass to git commands
        self._environment: Dict[str, str] = {}

        # Cached version slots
        self._version_info: Union[Tuple[int, ...], None] = None
        self._version_info_token: object = None

        # Cached command slots
        self.cat_file_header: Union[None, TBD] = None
        self.cat_file_all: Union[None, TBD] = None

    def __getattribute__(self, name: str) -> Any:
        if name == "USE_SHELL":
            _warn_use_shell(False)
        return super().__getattribute__(name)

    def __getattr__(self, name: str) -> Any:
        """A convenience method as it allows to call the command as if it was an object.

        :return:
            Callable object that will execute call :meth:`_call_process` with your
            arguments.
        """
        if name.startswith("_"):
            return super().__getattribute__(name)
        return lambda *args, **kwargs: self._call_process(name, *args, **kwargs)

    def set_persistent_git_options(self, **kwargs: Any) -> None:
        """Specify command line options to the git executable for subsequent
        subcommand calls.

        :param kwargs:
            A dict of keyword arguments.
            These arguments are passed as in :meth:`_call_process`, but will be passed
            to the git command rather than the subcommand.
        """

        self._persistent_git_options = self.transform_kwargs(split_single_char_options=True, **kwargs)

    @property
    def working_dir(self) -> Union[None, PathLike]:
        """:return: Git directory we are working on"""
        return self._working_dir

    @property
    def version_info(self) -> Tuple[int, ...]:
        """
        :return: Tuple with integers representing the major, minor and additional
            version numbers as parsed from :manpage:`git-version(1)`. Up to four fields
            are used.

            This value is generated on demand and is cached.
        """
        # Refreshing is global, but version_info caching is per-instance.
        refresh_token = self._refresh_token  # Copy token in case of concurrent refresh.

        # Use the cached version if obtained after the most recent refresh.
        if self._version_info_token is refresh_token:
            assert self._version_info is not None, "Bug: corrupted token-check state"
            return self._version_info

        # Run "git version" and parse it.
        process_version = self._call_process("version")
        version_string = process_version.split(" ")[2]
        version_fields = version_string.split(".")[:4]
        leading_numeric_fields = itertools.takewhile(str.isdigit, version_fields)
        self._version_info = tuple(map(int, leading_numeric_fields))

        # This value will be considered valid until the next refresh.
        self._version_info_token = refresh_token
        return self._version_info

    @overload
    def execute(
        self,
        command: Union[str, Sequence[Any]],
        *,
        as_process: Literal[True],
    ) -> "AutoInterrupt": ...

    @overload
    def execute(
        self,
        command: Union[str, Sequence[Any]],
        *,
        as_process: Literal[False] = False,
        stdout_as_string: Literal[True],
    ) -> Union[str, Tuple[int, str, str]]: ...

    @overload
    def execute(
        self,
        command: Union[str, Sequence[Any]],
        *,
        as_process: Literal[False] = False,
        stdout_as_string: Literal[False] = False,
    ) -> Union[bytes, Tuple[int, bytes, str]]: ...

    @overload
    def execute(
        self,
        command: Union[str, Sequence[Any]],
        *,
        with_extended_output: Literal[False],
        as_process: Literal[False],
        stdout_as_string: Literal[True],
    ) -> str: ...

    @overload
    def execute(
        self,
        command: Union[str, Sequence[Any]],
        *,
        with_extended_output: Literal[False],
        as_process: Literal[False],
        stdout_as_string: Literal[False],
    ) -> bytes: ...

    def execute(
        self,
        command: Union[str, Sequence[Any]],
        istream: Union[None, BinaryIO] = None,
        with_extended_output: bool = False,
        with_exceptions: bool = True,
        as_process: bool = False,
        output_stream: Union[None, BinaryIO] = None,
        stdout_as_string: bool = True,
        kill_after_timeout: Union[None, float] = None,
        with_stdout: bool = True,
        universal_newlines: bool = False,
        shell: Union[None, bool] = None,
        env: Union[None, Mapping[str, str]] = None,
        max_chunk_size: int = io.DEFAULT_BUFFER_SIZE,
        strip_newline_in_stdout: bool = True,
        **subprocess_kwargs: Any,
    ) -> Union[str, bytes, Tuple[int, Union[str, bytes], str], AutoInterrupt]:
        R"""Handle executing the command, and consume and return the returned
        information (stdout).

        :param command:
            The command argument list to execute.
            It should be a sequence of program arguments, or a string. The
            program to execute is the first item in the args sequence or string.

        :param istream:
            Standard input filehandle passed to :class:`subprocess.Popen`.

        :param with_extended_output:
            Whether to return a (status, stdout, stderr) tuple.

        :param with_exceptions:
            Whether to raise an exception when git returns a non-zero status.

        :param as_process:
            Whether to return the created process instance directly from which
            streams can be read on demand. This will render `with_extended_output`
            and `with_exceptions` ineffective - the caller will have to deal with
            the details. It is important to note that the process will be placed
            into an :class:`AutoInterrupt` wrapper that will interrupt the process
            once it goes out of scope. If you use the command in iterators, you
            should pass the whole process instance instead of a single stream.

        :param output_stream:
            If set to a file-like object, data produced by the git command will be
            copied to the given stream instead of being returned as a string.
            This feature only has any effect if `as_process` is ``False``.

        :param stdout_as_string:
            If ``False``, the command's standard output will be bytes. Otherwise, it
            will be decoded into a string using the default encoding (usually UTF-8).
            The latter can fail, if the output contains binary data.

        :param kill_after_timeout:
            Specifies a timeout in seconds for the git command, after which the process
            should be killed. This will have no effect if `as_process` is set to
            ``True``. It is set to ``None`` by default and will let the process run
            until the timeout is explicitly specified. Uses of this feature should be
            carefully considered, due to the following limitations:

            1. This feature is not supported at all on Windows.
            2. Effectiveness may vary by operating system. ``ps --ppid`` is used to
               enumerate child processes, which is available on most GNU/Linux systems
               but not most others.
            3. Deeper descendants do not receive signals, though they may sometimes
               terminate as a consequence of their parent processes being killed.
            4. `kill_after_timeout` uses ``SIGKILL``, which can have negative side
               effects on a repository. For example, stale locks in case of
               :manpage:`git-gc(1)` could render the repository incapable of accepting
               changes until the lock is manually removed.

        :param with_stdout:
            If ``True``, default ``True``, we open stdout on the created process.

        :param universal_newlines:
            If ``True``, pipes will be opened as text, and lines are split at all known
            line endings.

        :param shell:
            Whether to invoke commands through a shell
            (see :class:`Popen(..., shell=True) <subprocess.Popen>`).
            If this is not ``None``, it overrides :attr:`USE_SHELL`.

            Passing ``shell=True`` to this or any other GitPython function should be
            avoided, as it is unsafe under most circumstances. This is because it is
            typically not feasible to fully consider and account for the effect of shell
            expansions, especially when passing ``shell=True`` to other methods that
            forward it to :meth:`Git.execute`. Passing ``shell=True`` is also no longer
            needed (nor useful) to work around any known operating system specific
            issues.

        :param env:
            A dictionary of environment variables to be passed to
            :class:`subprocess.Popen`.

        :param max_chunk_size:
            Maximum number of bytes in one chunk of data passed to the `output_stream`
            in one invocation of its ``write()`` method. If the given number is not
            positive then the default value is used.

        :param strip_newline_in_stdout:
            Whether to strip the trailing ``\n`` of the command stdout.

        :param subprocess_kwargs:
            Keyword arguments to be passed to :class:`subprocess.Popen`. Please note
            that some of the valid kwargs are already set by this method; the ones you
            specify may not be the same ones.

        :return:
            * str(output), if `extended_output` is ``False`` (Default)
            * tuple(int(status), str(stdout), str(stderr)),
              if `extended_output` is ``True``

            If `output_stream` is ``True``, the stdout value will be your output stream:

            * output_stream, if `extended_output` is ``False``
            * tuple(int(status), output_stream, str(stderr)),
              if `extended_output` is ``True``

            Note that git is executed with ``LC_MESSAGES="C"`` to ensure consistent
            output regardless of system language.

        :raise git.exc.GitCommandError:

        :note:
            If you add additional keyword arguments to the signature of this method, you
            must update the ``execute_kwargs`` variable housed in this module.
        """
        # Remove password for the command if present.
        redacted_command = remove_password_if_present(command)
        if self.GIT_PYTHON_TRACE and (self.GIT_PYTHON_TRACE != "full" or as_process):
            _logger.info(" ".join(redacted_command))

        # Allow the user to have the command executed in their working dir.
        try:
            cwd = self._working_dir or os.getcwd()  # type: Union[None, str]
            if not os.access(str(cwd), os.X_OK):
                cwd = None
        except FileNotFoundError:
            cwd = None

        # Start the process.
        inline_env = env
        env = os.environ.copy()
        # Attempt to force all output to plain ASCII English, which is what some parsing
        # code may expect.
        # According to https://askubuntu.com/a/311796, we are setting LANGUAGE as well
        # just to be sure.
        env["LANGUAGE"] = "C"
        env["LC_ALL"] = "C"
        env.update(self._environment)
        if inline_env is not None:
            env.update(inline_env)

        if sys.platform == "win32":
            if kill_after_timeout is not None:
                raise GitCommandError(
                    redacted_command,
                    '"kill_after_timeout" feature is not supported on Windows.',
                )
            cmd_not_found_exception = OSError
        else:
            cmd_not_found_exception = FileNotFoundError
        # END handle

        stdout_sink = PIPE if with_stdout else getattr(subprocess, "DEVNULL", None) or open(os.devnull, "wb")
        if shell is None:
            # Get the value of USE_SHELL with no deprecation warning. Do this without
            # warnings.catch_warnings, to avoid a race condition with application code
            # configuring warnings. The value could be looked up in type(self).__dict__
            # or Git.__dict__, but those can break under some circumstances. This works
            # the same as self.USE_SHELL in more situations; see Git.__getattribute__.
            shell = super().__getattribute__("USE_SHELL")
        _logger.debug(
            "Popen(%s, cwd=%s, stdin=%s, shell=%s, universal_newlines=%s)",
            redacted_command,
            cwd,
            "<valid stream>" if istream else "None",
            shell,
            universal_newlines,
        )
        try:
            proc = safer_popen(
                command,
                env=env,
                cwd=cwd,
                bufsize=-1,
                stdin=(istream or DEVNULL),
                stderr=PIPE,
                stdout=stdout_sink,
                shell=shell,
                universal_newlines=universal_newlines,
                **subprocess_kwargs,
            )
        except cmd_not_found_exception as err:
            raise GitCommandNotFound(redacted_command, err) from err
        else:
            # Replace with a typeguard for Popen[bytes]?
            proc.stdout = cast(BinaryIO, proc.stdout)
            proc.stderr = cast(BinaryIO, proc.stderr)

        if as_process:
            return self.AutoInterrupt(proc, command)

        if sys.platform != "win32" and kill_after_timeout is not None:
            # Help mypy figure out this is not None even when used inside communicate().
            timeout = kill_after_timeout

            def kill_process(pid: int) -> None:
                """Callback to kill a process.

                This callback implementation would be ineffective and unsafe on Windows.
                """
                p = Popen(["ps", "--ppid", str(pid)], stdout=PIPE)
                child_pids = []
                if p.stdout is not None:
                    for line in p.stdout:
                        if len(line.split()) > 0:
                            local_pid = (line.split())[0]
                            if local_pid.isdigit():
                                child_pids.append(int(local_pid))
                try:
                    os.kill(pid, signal.SIGKILL)
                    for child_pid in child_pids:
                        try:
                            os.kill(child_pid, signal.SIGKILL)
                        except OSError:
                            pass
                    # Tell the main routine that the process was killed.
                    kill_check.set()
                except OSError:
                    # It is possible that the process gets completed in the duration
                    # after timeout happens and before we try to kill the process.
                    pass
                return

            def communicate() -> Tuple[AnyStr, AnyStr]:
                watchdog.start()
                out, err = proc.communicate()
                watchdog.cancel()
                if kill_check.is_set():
                    err = 'Timeout: the command "%s" did not complete in %d ' "secs." % (
                        " ".join(redacted_command),
                        timeout,
                    )
                    if not universal_newlines:
                        err = err.encode(defenc)
                return out, err

            # END helpers

            kill_check = threading.Event()
            watchdog = threading.Timer(timeout, kill_process, args=(proc.pid,))
        else:
            communicate = proc.communicate

        # Wait for the process to return.
        status = 0
        stdout_value: Union[str, bytes] = b""
        stderr_value: Union[str, bytes] = b""
        newline = "\n" if universal_newlines else b"\n"
        try:
            if output_stream is None:
                stdout_value, stderr_value = communicate()
                # Strip trailing "\n".
                if stdout_value.endswith(newline) and strip_newline_in_stdout:  # type: ignore[arg-type]
                    stdout_value = stdout_value[:-1]
                if stderr_value.endswith(newline):  # type: ignore[arg-type]
                    stderr_value = stderr_value[:-1]

                status = proc.returncode
            else:
                max_chunk_size = max_chunk_size if max_chunk_size and max_chunk_size > 0 else io.DEFAULT_BUFFER_SIZE
                stream_copy(proc.stdout, output_stream, max_chunk_size)
                stdout_value = proc.stdout.read()
                stderr_value = proc.stderr.read()
                # Strip trailing "\n".
                if stderr_value.endswith(newline):  # type: ignore[arg-type]
                    stderr_value = stderr_value[:-1]
                status = proc.wait()
            # END stdout handling
        finally:
            proc.stdout.close()
            proc.stderr.close()

        if self.GIT_PYTHON_TRACE == "full":
            cmdstr = " ".join(redacted_command)

            def as_text(stdout_value: Union[bytes, str]) -> str:
                return not output_stream and safe_decode(stdout_value) or "<OUTPUT_STREAM>"

            # END as_text

            if stderr_value:
                _logger.info(
                    "%s -> %d; stdout: '%s'; stderr: '%s'",
                    cmdstr,
                    status,
                    as_text(stdout_value),
                    safe_decode(stderr_value),
                )
            elif stdout_value:
                _logger.info("%s -> %d; stdout: '%s'", cmdstr, status, as_text(stdout_value))
            else:
                _logger.info("%s -> %d", cmdstr, status)
        # END handle debug printing

        if with_exceptions and status != 0:
            raise GitCommandError(redacted_command, status, stderr_value, stdout_value)

        if isinstance(stdout_value, bytes) and stdout_as_string:  # Could also be output_stream.
            stdout_value = safe_decode(stdout_value)

        # Allow access to the command's status code.
        if with_extended_output:
            return (status, stdout_value, safe_decode(stderr_value))
        else:
            return stdout_value

    def environment(self) -> Dict[str, str]:
        return self._environment

    def update_environment(self, **kwargs: Any) -> Dict[str, Union[str, None]]:
        """Set environment variables for future git invocations. Return all changed
        values in a format that can be passed back into this function to revert the
        changes.

        Examples::

            old_env = self.update_environment(PWD='/tmp')
            self.update_environment(**old_env)

        :param kwargs:
            Environment variables to use for git processes.

        :return:
            Dict that maps environment variables to their old values
        """
        old_env = {}
        for key, value in kwargs.items():
            # Set value if it is None.
            if value is not None:
                old_env[key] = self._environment.get(key)
                self._environment[key] = value
            # Remove key from environment if its value is None.
            elif key in self._environment:
                old_env[key] = self._environment[key]
                del self._environment[key]
        return old_env

    @contextlib.contextmanager
    def custom_environment(self, **kwargs: Any) -> Iterator[None]:
        """A context manager around the above :meth:`update_environment` method to
        restore the environment back to its previous state after operation.

        Examples::

            with self.custom_environment(GIT_SSH='/bin/ssh_wrapper'):
                repo.remotes.origin.fetch()

        :param kwargs:
            See :meth:`update_environment`.
        """
        old_env = self.update_environment(**kwargs)
        try:
            yield
        finally:
            self.update_environment(**old_env)

    def transform_kwarg(self, name: str, value: Any, split_single_char_options: bool) -> List[str]:
        if len(name) == 1:
            if value is True:
                return ["-%s" % name]
            elif value not in (False, None):
                if split_single_char_options:
                    return ["-%s" % name, "%s" % value]
                else:
                    return ["-%s%s" % (name, value)]
        else:
            if value is True:
                return ["--%s" % dashify(name)]
            elif value is not False and value is not None:
                return ["--%s=%s" % (dashify(name), value)]
        return []

    def transform_kwargs(self, split_single_char_options: bool = True, **kwargs: Any) -> List[str]:
        """Transform Python-style kwargs into git command line options."""
        args = []
        for k, v in kwargs.items():
            if isinstance(v, (list, tuple)):
                for value in v:
                    args += self.transform_kwarg(k, value, split_single_char_options)
            else:
                args += self.transform_kwarg(k, v, split_single_char_options)
        return args

    @classmethod
    def _unpack_args(cls, arg_list: Sequence[str]) -> List[str]:
        outlist = []
        if isinstance(arg_list, (list, tuple)):
            for arg in arg_list:
                outlist.extend(cls._unpack_args(arg))
        else:
            outlist.append(str(arg_list))

        return outlist

    def __call__(self, **kwargs: Any) -> "Git":
        """Specify command line options to the git executable for a subcommand call.

        :param kwargs:
            A dict of keyword arguments.
            These arguments are passed as in :meth:`_call_process`, but will be passed
            to the git command rather than the subcommand.

        Examples::

            git(work_tree='/tmp').difftool()
        """
        self._git_options = self.transform_kwargs(split_single_char_options=True, **kwargs)
        return self

    @overload
    def _call_process(
        self, method: str, *args: None, **kwargs: None
    ) -> str: ...  # If no args were given, execute the call with all defaults.

    @overload
    def _call_process(
        self,
        method: str,
        istream: int,
        as_process: Literal[True],
        *args: Any,
        **kwargs: Any,
    ) -> "Git.AutoInterrupt": ...

    @overload
    def _call_process(
        self, method: str, *args: Any, **kwargs: Any
    ) -> Union[str, bytes, Tuple[int, Union[str, bytes], str], "Git.AutoInterrupt"]: ...

    def _call_process(
        self, method: str, *args: Any, **kwargs: Any
    ) -> Union[str, bytes, Tuple[int, Union[str, bytes], str], "Git.AutoInterrupt"]:
        """Run the given git command with the specified arguments and return the result
        as a string.

        :param method:
            The command. Contained ``_`` characters will be converted to hyphens, such
            as in ``ls_files`` to call ``ls-files``.

        :param args:
            The list of arguments. If ``None`` is included, it will be pruned.
            This allows your commands to call git more conveniently, as ``None`` is
            realized as non-existent.

        :param kwargs:
            Contains key-values for the following:

            - The :meth:`execute()` kwds, as listed in ``execute_kwargs``.
            - "Command options" to be converted by :meth:`transform_kwargs`.
            - The ``insert_kwargs_after`` key which its value must match one of
              ``*args``.

            It also contains any command options, to be appended after the matched arg.

        Examples::

            git.rev_list('master', max_count=10, header=True)

        turns into::

            git rev-list max-count 10 --header master

        :return:
            Same as :meth:`execute`. If no args are given, used :meth:`execute`'s
            default (especially ``as_process = False``, ``stdout_as_string = True``) and
            return :class:`str`.
        """
        # Handle optional arguments prior to calling transform_kwargs.
        # Otherwise these'll end up in args, which is bad.
        exec_kwargs = {k: v for k, v in kwargs.items() if k in execute_kwargs}
        opts_kwargs = {k: v for k, v in kwargs.items() if k not in execute_kwargs}

        insert_after_this_arg = opts_kwargs.pop("insert_kwargs_after", None)

        # Prepare the argument list.

        opt_args = self.transform_kwargs(**opts_kwargs)
        ext_args = self._unpack_args([a for a in args if a is not None])

        if insert_after_this_arg is None:
            args_list = opt_args + ext_args
        else:
            try:
                index = ext_args.index(insert_after_this_arg)
            except ValueError as err:
                raise ValueError(
                    "Couldn't find argument '%s' in args %s to insert cmd options after"
                    % (insert_after_this_arg, str(ext_args))
                ) from err
            # END handle error
            args_list = ext_args[: index + 1] + opt_args + ext_args[index + 1 :]
        # END handle opts_kwargs

        call = [self.GIT_PYTHON_GIT_EXECUTABLE]

        # Add persistent git options.
        call.extend(self._persistent_git_options)

        # Add the git options, then reset to empty to avoid side effects.
        call.extend(self._git_options)
        self._git_options = ()

        call.append(dashify(method))
        call.extend(args_list)

        return self.execute(call, **exec_kwargs)

    def _parse_object_header(self, header_line: str) -> Tuple[str, str, int]:
        """
        :param header_line:
            A line of the form::

                <hex_sha> type_string size_as_int

        :return:
            (hex_sha, type_string, size_as_int)

        :raise ValueError:
            If the header contains indication for an error due to incorrect input sha.
        """
        tokens = header_line.split()
        if len(tokens) != 3:
            if not tokens:
                err_msg = (
                    f"SHA is empty, possible dubious ownership in the repository "
                    f"""at {self._working_dir}.\n            If this is unintended run:\n\n         """
                    f"""             "git config --global --add safe.directory {self._working_dir}" """
                )
                raise ValueError(err_msg)
            else:
                raise ValueError("SHA %s could not be resolved, git returned: %r" % (tokens[0], header_line.strip()))
            # END handle actual return value
        # END error handling

        if len(tokens[0]) != 40:
            raise ValueError("Failed to parse header: %r" % header_line)
        return (tokens[0], tokens[1], int(tokens[2]))

    def _prepare_ref(self, ref: AnyStr) -> bytes:
        # Required for command to separate refs on stdin, as bytes.
        if isinstance(ref, bytes):
            # Assume 40 bytes hexsha - bin-to-ascii for some reason returns bytes, not text.
            refstr: str = ref.decode("ascii")
        elif not isinstance(ref, str):
            refstr = str(ref)  # Could be ref-object.
        else:
            refstr = ref

        if not refstr.endswith("\n"):
            refstr += "\n"
        return refstr.encode(defenc)

    def _get_persistent_cmd(self, attr_name: str, cmd_name: str, *args: Any, **kwargs: Any) -> "Git.AutoInterrupt":
        cur_val = getattr(self, attr_name)
        if cur_val is not None:
            return cur_val

        options = {"istream": PIPE, "as_process": True}
        options.update(kwargs)

        cmd = self._call_process(cmd_name, *args, **options)
        setattr(self, attr_name, cmd)
        cmd = cast("Git.AutoInterrupt", cmd)
        return cmd

    def __get_object_header(self, cmd: "Git.AutoInterrupt", ref: AnyStr) -> Tuple[str, str, int]:
        if cmd.stdin and cmd.stdout:
            cmd.stdin.write(self._prepare_ref(ref))
            cmd.stdin.flush()
            return self._parse_object_header(cmd.stdout.readline())
        else:
            raise ValueError("cmd stdin was empty")

    def get_object_header(self, ref: str) -> Tuple[str, str, int]:
        """Use this method to quickly examine the type and size of the object behind the
        given ref.

        :note:
            The method will only suffer from the costs of command invocation once and
            reuses the command in subsequent calls.

        :return:
            (hexsha, type_string, size_as_int)
        """
        cmd = self._get_persistent_cmd("cat_file_header", "cat_file", batch_check=True)
        return self.__get_object_header(cmd, ref)

    def get_object_data(self, ref: str) -> Tuple[str, str, int, bytes]:
        """Similar to :meth:`get_object_header`, but returns object data as well.

        :return:
            (hexsha, type_string, size_as_int, data_string)

        :note:
            Not threadsafe.
        """
        hexsha, typename, size, stream = self.stream_object_data(ref)
        data = stream.read(size)
        del stream
        return (hexsha, typename, size, data)

    def stream_object_data(self, ref: str) -> Tuple[str, str, int, "Git.CatFileContentStream"]:
        """Similar to :meth:`get_object_data`, but returns the data as a stream.

        :return:
            (hexsha, type_string, size_as_int, stream)

        :note:
            This method is not threadsafe. You need one independent :class:`Git`
            instance per thread to be safe!
        """
        cmd = self._get_persistent_cmd("cat_file_all", "cat_file", batch=True)
        hexsha, typename, size = self.__get_object_header(cmd, ref)
        cmd_stdout = cmd.stdout if cmd.stdout is not None else io.BytesIO()
        return (hexsha, typename, size, self.CatFileContentStream(size, cmd_stdout))

    def clear_cache(self) -> "Git":
        """Clear all kinds of internal caches to release resources.

        Currently persistent commands will be interrupted.

        :return:
            self
        """
        for cmd in (self.cat_file_all, self.cat_file_header):
            if cmd:
                cmd.__del__()

        self.cat_file_all = None
        self.cat_file_header = None
        return self
