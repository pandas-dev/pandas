from __future__ import annotations

__lazy_modules__ = {
    "contextlib",
    "plumbum.lib",
    "plumbum.path",
    "plumbum.path.local",
    "tempfile",
}

import contextlib
import re
import typing
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING

from plumbum.commands import CommandNotFound, ConcreteCommand, shquote
from plumbum.lib import ProcInfo
from plumbum.machines.base import BaseMachine, PopenWithAddons
from plumbum.machines.env import BaseEnv
from plumbum.path.local import LocalPath
from plumbum.path.remote import RemotePath, RemoteStatRes, RemoteWorkdir

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence

    from plumbum._compat.typing import Self
    from plumbum.commands.async_ import AsyncRemoteCommand
    from plumbum.machines.session import ShellSession


_ENV_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _check_env_name(name: str) -> str:
    """Validate a shell environment variable *name*.

    Variable names are interpolated unquoted into remote shell command lines,
    so a name containing shell metacharacters would allow command injection.
    Only POSIX identifiers (``[A-Za-z_][A-Za-z0-9_]*``) are permitted.

    :returns: the validated name (so it can be used inline)
    """
    if not _ENV_NAME_RE.match(name):
        raise ValueError(f"Invalid environment variable name: {name!r}")
    return name


class RemoteEnv(BaseEnv[RemotePath]):
    """The remote machine's environment; exposes a dict-like interface"""

    __slots__ = ["_orig", "remote"]

    def __init__(self, remote: BaseRemoteMachine) -> None:
        session = remote._session
        # GNU env has a -0 argument; use it if present. Otherwise,
        # fall back to calling printenv on each (possible) variable
        # from plain env.
        env0 = session.run("env -0; echo")
        if env0[0] == 0 and not env0[2].rstrip():
            _curr = dict(
                line.split("=", 1) for line in env0[1].split("\x00") if "=" in line
            )
        else:
            lines = session.run("env; echo")[1].splitlines()
            split = (line.split("=", 1) for line in lines)
            keys = (line[0] for line in split if len(line) > 1)
            runs = ((key, session.run(f'printenv "{key}"; echo')) for key in keys)
            _curr = {
                key: run[1].rstrip("\n")
                for (key, run) in runs
                if run[0] == 0 and run[1].rstrip("\n") and not run[2]
            }

        super().__init__(remote.path, ":", _curr=_curr)
        self.remote = remote
        self._orig = dict(self._curr)

    def __delitem__(self, name: str) -> None:
        _check_env_name(name)
        BaseEnv.__delitem__(self, name)
        self.remote._session.run(f"unset {name}")

    def __setitem__(self, name: str, value: str) -> None:
        _check_env_name(name)
        BaseEnv.__setitem__(self, name, value)
        self.remote._session.run(f"export {name}={shquote(value)}")

    def pop(self, name: str, *default: str) -> str | None:
        _check_env_name(name)
        value = BaseEnv.pop(self, name, *default)
        self.remote._session.run(f"unset {name}")
        return value

    def update(self, *args: typing.Any, **kwargs: typing.Any) -> None:
        BaseEnv.update(self, *args, **kwargs)
        self.remote._session.run(
            "export "
            + " ".join(
                f"{_check_env_name(k)}={shquote(v)}" for k, v in self.getdict().items()
            )
        )

    def expand(self, expr: str) -> str:
        """Expands any environment variables and home shortcuts found in ``expr``
        (like ``os.path.expanduser`` combined with ``os.path.expandvars``)

        :param expr: An expression containing environment variables (as ``$FOO``) or
                     home shortcuts (as ``~/.bashrc``)

        :returns: The expanded string"""
        return self.remote.expand(expr)

    def expanduser(self, expr: str) -> str:
        """Expand home shortcuts (e.g., ``~/foo/bar`` or ``~john/foo/bar``)

        :param expr: An expression containing home shortcuts

        :returns: The expanded string"""
        return self.remote.expanduser(expr)

    # def clear(self):
    #    BaseEnv.clear(self, *args, **kwargs)
    #    self.remote._session.run("export %s" % " ".join("%s=%s" % (k, v) for k, v in self.getdict()))

    def getdelta(self) -> dict[str, str]:
        """Returns the difference between the this environment and the original environment of
        the remote machine"""
        self._curr["PATH"] = self.path.join()

        delta = {}
        for k, v in self._curr.items():
            if k not in self._orig:
                delta[k] = str(v)
        for k, v in self._orig.items():
            if k not in self._curr:
                delta[k] = ""
            elif v != self._curr[k]:
                delta[k] = self._curr[k]

        return delta


class RemoteCommand(ConcreteCommand):
    __slots__ = ("remote",)
    QUOTE_LEVEL = 1

    def __init__(
        self, remote: BaseRemoteMachine, executable: RemotePath, encoding: str = "auto"
    ) -> None:
        self.remote = remote
        ConcreteCommand.__init__(
            self, executable, remote.custom_encoding if encoding == "auto" else encoding
        )

    @property
    def machine(self) -> BaseRemoteMachine:
        return self.remote

    def __repr__(self) -> str:
        return f"RemoteCommand({self.remote!r}, {self.executable!r})"

    def popen(
        self, args: Sequence[typing.Any] | str = (), **kwargs: typing.Any
    ) -> PopenWithAddons[str]:
        return self.remote.popen(self[args], **kwargs)  # type: ignore[arg-type]

    def nohup(
        self,
        cwd: str = ".",
        stdout: str = "nohup.out",
        stderr: str | None = None,
        append: bool = True,
    ) -> PopenWithAddons[str]:
        """Runs a command detached."""
        return self.machine.daemonic_popen(self, cwd, stdout, stderr, append)


class ClosedRemoteMachine(Exception):
    pass


class ClosedRemote:
    __slots__ = ["__weakref__", "_obj"]

    def __init__(self, obj: object) -> None:
        self._obj = obj

    def close(self) -> None:
        pass

    def __getattr__(
        self, name: str
    ) -> typing.NoReturn:  # pragma: no cover - always raises
        raise ClosedRemoteMachine(f"{self._obj!r} has been closed")


def _is_recursive_glob(pattern: str) -> bool:
    """Whether ``pattern`` uses ``**`` as a recursive wildcard.

    As in :mod:`glob`/:mod:`pathlib`, ``**`` is only special when it is a whole
    path segment; ``a**b`` is just two ``*`` wildcards within one segment.
    """
    return "**" in pattern.split("/")


def _segment_to_regex(segment: str) -> str:
    """Translate a single (slash-free) glob segment into a regex fragment.

    Supports ``*``, ``?`` and ``[...]`` character classes, all confined to a
    single path component (they never cross ``/``). A leading wildcard does not
    match a leading dot, mirroring :func:`glob.glob` (which skips dotfiles
    unless the pattern segment starts with a literal ``.``).
    """
    out = []
    # A wildcard at the start of a segment must not match a leading dot.
    if segment[:1] in ("*", "?", "["):
        out.append(r"(?!\.)")
    i, n = 0, len(segment)
    while i < n:
        c = segment[i]
        if c == "*":
            out.append("[^/]*")
            i += 1
        elif c == "?":
            out.append("[^/]")
            i += 1
        elif c == "[":
            j = i + 1
            if j < n and segment[j] == "!":
                j += 1
            if j < n and segment[j] == "]":
                j += 1
            while j < n and segment[j] != "]":
                j += 1
            if j >= n:  # no closing bracket -- treat "[" as a literal
                out.append(re.escape("["))
                i += 1
            else:
                stuff = segment[i + 1 : j].replace("\\", r"\\")
                if stuff[0] == "!":
                    stuff = "^" + stuff[1:]
                elif stuff[0] in ("^", "["):
                    stuff = "\\" + stuff
                out.append(f"[{stuff}]")
                i = j + 1
        else:
            out.append(re.escape(c))
            i += 1
    return "".join(out)


def _glob_to_regex(pattern: str) -> str:
    """Translate a glob pattern into an anchored regex with pathlib-like ``**``.

    ``**`` as a whole path segment matches any number of (non-hidden)
    directories; ``*``, ``?`` and ``[...]`` match within a single path segment.
    Dotfiles are not matched unless the relevant pattern segment starts with a
    literal ``.`` -- matching :func:`glob.glob`, which is the local backend.
    Used to match recursive globs in Python instead of relying on
    shell-specific recursion support.
    """
    segments = pattern.split("/")
    out = []
    need_sep = False  # whether a "/" must precede the next fragment
    for i, segment in enumerate(segments):
        if segment == "**":
            if need_sep:
                out.append("/")
            if i == len(segments) - 1:
                # trailing ``**``: one or more non-hidden path components
                out.append(r"(?!\.)[^/]+(?:/(?!\.)[^/]+)*")
                need_sep = False
            else:
                # ``**/``: zero or more non-hidden directories (slash included)
                out.append(r"(?:(?!\.)[^/]+/)*")
                need_sep = False
            continue
        if need_sep:
            out.append("/")
        out.append(_segment_to_regex(segment))
        need_sep = True
    return "(?s:" + "".join(out) + r")\Z"


class BaseRemoteMachine(BaseMachine):
    """Represents a *remote machine*; serves as an entry point to everything related to that
    remote machine, such as working directory and environment manipulation, command creation,
    etc.

    Attributes:

    * ``cwd`` - the remote working directory
    * ``env`` - the remote environment
    * ``custom_encoding`` - the remote machine's default encoding (assumed to be UTF8)
    * ``connect_timeout`` - the connection timeout


    There also is a _cwd attribute that exists if the cwd is not current (del if cwd is changed).
    """

    __slots__ = (
        "_cwd",
        "_program_cache",
        "_python",
        "_session",
        "connect_timeout",
        "env",
        "uname",
    )

    # allow inheritors to override the RemoteCommand class
    RemoteCommand = RemoteCommand

    @property
    def cwd(self) -> RemoteWorkdir:
        if not hasattr(self, "_cwd"):
            self._cwd = RemoteWorkdir(self)
        return self._cwd

    def __init__(
        self,
        encoding: str = "utf8",
        connect_timeout: float | None = 10,
        new_session: bool = False,
    ) -> None:
        self.custom_encoding = encoding
        self.connect_timeout = connect_timeout
        self._session: ShellSession | ClosedRemote = self.session(
            new_session=new_session
        )
        self.uname = self._get_uname()
        self.env = RemoteEnv(self)
        self._python: ConcreteCommand | None = None
        self._program_cache: dict[tuple[str, str], RemotePath] = {}

    def clear_program_cache(self) -> None:
        self._program_cache.clear()

    def _get_uname(self) -> str:
        rc, out, _ = self._session.run("uname", retcode=None)
        if rc == 0:
            return out.strip()

        rc, out, _ = self._session.run(
            "python3 -c 'import platform;print(platform.uname()[0])'", retcode=None
        )
        if rc == 0:
            return out.strip()

        # all POSIX systems should have uname. make an educated guess it's Windows
        return "Windows"

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} {self}>"

    def __enter__(self) -> Self:
        return self

    def __exit__(self, t: object, v: object, tb: object) -> None:
        self.close()

    def __del__(self) -> None:
        with contextlib.suppress(Exception):
            self.close()

    def close(self) -> None:
        """closes the connection to the remote machine; all paths and programs will
        become defunct"""
        self._session.close()
        self._session = ClosedRemote(self)

    def path(self, *parts: str | RemotePath | LocalPath) -> RemotePath:
        """A factory for :class:`RemotePaths <plumbum.path.remote.RemotePath>`.
        Usage: ``p = rem.path("/usr", "lib", "python2.7")``
        """
        parts2 = [str(self.cwd)]
        for p in parts:
            if isinstance(p, LocalPath):
                raise TypeError(f"Cannot construct RemotePath from {p!r}")
            parts2.append(self.expanduser(str(p)))
        return RemotePath(self, *parts2)

    def which(self, progname: str) -> RemotePath:
        """Looks up a program in the ``PATH``. If the program is not found, raises
        :class:`CommandNotFound <plumbum.commands.processes.CommandNotFound>`

        :param progname: The program's name. Note that if underscores (``_``) are present
                         in the name, and the exact name is not found, they will be replaced
                         in turn by hyphens (``-``) then periods (``.``), and the name will
                         be looked up again for each alternative

        :returns: A :class:`RemotePath <plumbum.path.remote.RemotePath>`
        """
        key = (progname, self.env.get("PATH", ""))

        with contextlib.suppress(KeyError):
            return self._program_cache[key]

        alternatives = [progname]
        if "_" in progname:
            alternatives += [progname.replace("_", "-"), progname.replace("_", ".")]
        for name in alternatives:
            for p in self.env.path:
                fn = p / name
                if fn.access("x") and not fn.is_dir():
                    self._program_cache[key] = fn
                    return fn

        raise CommandNotFound(progname, self.env.path)

    def __getitem__(self, cmd: str | RemotePath | LocalPath) -> ConcreteCommand:
        """Returns a `Command` object representing the given program. ``cmd`` can be a string or
        a :class:`RemotePath <plumbum.path.remote.RemotePath>`; if it is a path, a command
        representing this path will be returned; otherwise, the program name will be looked up in
        the system's ``PATH`` (using ``which``). Usage::

            r_ls = rem["ls"]
        """
        if isinstance(cmd, RemotePath):
            if cmd.remote is self:
                return self.RemoteCommand(self, cmd)

            raise TypeError(
                f"Given path does not belong to this remote machine: {cmd!r}"
            )

        if not isinstance(cmd, LocalPath):
            return self.RemoteCommand(
                self, self.path(cmd) if "/" in cmd or "\\" in cmd else self.which(cmd)
            )

        raise TypeError(f"cmd must not be a LocalPath: {cmd!r}")

    @property
    def python(self) -> ConcreteCommand:
        """A command that represents the default remote python interpreter"""
        if not self._python:
            self._python = self["python3"]
        return self._python

    def session(
        self, isatty: bool = False, *, new_session: bool = False
    ) -> ShellSession:
        """Creates a new :class:`ShellSession <plumbum.machines.session.ShellSession>` object; this invokes the user's
        shell on the remote machine and executes commands on it over stdin/stdout/stderr
        """
        raise NotImplementedError()

    def download(self, src: str | RemotePath, dst: str | LocalPath) -> None:
        """Downloads a remote file/directory (``src``) to a local destination (``dst``).
        ``src`` must be a string or a :class:`RemotePath <plumbum.path.remote.RemotePath>`
        pointing to this remote machine, and ``dst`` must be a string or a
        :class:`LocalPath <plumbum.path.local.LocalPath>`"""
        raise NotImplementedError()

    def upload(self, src: str | LocalPath, dst: str | RemotePath) -> None:
        """Uploads a local file/directory (``src``) to a remote destination (``dst``).
        ``src`` must be a string or a :class:`LocalPath <plumbum.path.local.LocalPath>`,
        and ``dst`` must be a string or a :class:`RemotePath <plumbum.path.remote.RemotePath>`
        pointing to this remote machine"""
        raise NotImplementedError()

    def popen(
        self, args: Sequence[typing.Any] | str, **kwargs: typing.Any
    ) -> PopenWithAddons[str]:
        """Spawns the given command on the remote machine, returning a ``Popen``-like object;
        do not use this method directly, unless you need "low-level" control on the remote
        process"""
        raise NotImplementedError()

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
            yield ProcInfo(int(parts[0]), int(parts[1]), parts[2], " ".join(parts[3:]))

    def pgrep(self, pattern: str) -> Generator[ProcInfo, None, None]:
        """
        Process grep: return information about all processes whose command-line args match the given regex pattern
        """
        pat = re.compile(pattern)
        for procinfo in self.list_processes():
            if pat.search(procinfo.args):
                yield procinfo

    @contextlib.contextmanager
    def tempdir(self) -> Generator[RemotePath, None, None]:
        """A context manager that creates a remote temporary directory, which is removed when
        the context exits"""
        _, out, _ = self._session.run(
            "mktemp -d 2>/dev/null || mktemp -d tmp.XXXXXXXXXX"
        )
        local_dir = self.path(out.strip())
        try:
            yield local_dir
        finally:
            local_dir.delete()

    #
    # Path implementation
    #
    def _path_listdir(self, fn: str) -> list[str]:
        files = self._session.run(f"ls -a {shquote(fn)}")[1].splitlines()
        files.remove(".")
        files.remove("..")
        return files

    def _path_glob(self, fn: str, pattern: str) -> list[str]:
        # Both branches enumerate the tree with POSIX ``find`` and match in
        # Python rather than relying on the remote shell's globbing. This is
        # portable across shells (``**`` only recurses in bash with
        # ``globstar``) and avoids the shell-glob quirks -- injection and
        # whitespace-mangling -- that bite directory names containing spaces or
        # glob/shell metacharacters.
        regex = re.compile(_glob_to_regex(pattern))
        # ``find`` exits non-zero on partial errors (e.g. an unreadable
        # subdirectory) while still printing the matches it did find, so match
        # whatever was printed rather than discarding it -- this mirrors
        # ``glob.glob``, which does not error on unreadable subdirs.
        if _is_recursive_glob(pattern):
            find_cmd = f"find {shquote(fn)}"
        else:
            # A non-recursive pattern spans a fixed number of path components,
            # so cap ``find``'s depth to avoid descending the whole tree.
            depth = pattern.count("/") + 1
            find_cmd = f"find {shquote(fn)} -maxdepth {depth}"
        _, out, _ = self._session.run(find_cmd, retcode=None)
        prefix = fn.rstrip("/") + "/"
        return sorted(
            line
            for line in out.splitlines()
            if line.startswith(prefix) and regex.match(line[len(prefix) :])
        )

    def _path_getuid(self, fn: str) -> list[str]:
        stat_cmd = (
            "stat -c '%u,%U' "
            if self.uname not in ("Darwin", "FreeBSD")
            else "stat -f '%u,%Su' "
        )
        return self._session.run(stat_cmd + shquote(fn))[1].strip().split(",")

    def _path_getgid(self, fn: str) -> list[str]:
        stat_cmd = (
            "stat -c '%g,%G' "
            if self.uname not in ("Darwin", "FreeBSD")
            else "stat -f '%g,%Sg' "
        )
        return self._session.run(stat_cmd + shquote(fn))[1].strip().split(",")

    def _path_stat(self, fn: str) -> RemoteStatRes | None:
        if self.uname not in ("Darwin", "FreeBSD"):
            stat_cmd = "stat -c '%F,%f,%i,%d,%h,%u,%g,%s,%X,%Y,%Z' "
        else:
            stat_cmd = "stat -f '%HT,%Xp,%i,%d,%l,%u,%g,%z,%a,%m,%c' "
        rc, out, _ = self._session.run(stat_cmd + shquote(fn), retcode=None)
        if rc != 0:
            return None
        statres = out.strip().split(",")
        text_mode = statres.pop(0).lower()
        res = RemoteStatRes(
            (int(statres[0], 16), *tuple(int(sr) for sr in statres[1:]))  # type: ignore[arg-type]
        )
        res.text_mode = text_mode
        return res

    def _path_lstat(self, fn: str) -> RemoteStatRes | None:
        # ``stat`` without -L/--dereference already reports the link itself.
        return self._path_stat(fn)

    def _path_delete(self, fn: str) -> None:
        self._session.run(f"rm -rf {shquote(fn)}")

    def _path_unlink(self, fn: str) -> None:
        # Non-recursive removal: refuses to descend into directories.
        self._session.run(f"rm -f {shquote(fn)}")

    def _path_move(self, src: str, dst: str) -> RemotePath:
        self._session.run(f"mv {shquote(src)} {shquote(dst)}")
        return RemotePath(self, dst)

    def _path_copy(self, src: str, dst: str) -> RemotePath:
        self._session.run(f"cp -r {shquote(src)} {shquote(dst)}")
        return RemotePath(self, dst)

    def _path_mkdir(
        self,
        fn: str,
        mode: int | None = None,  # noqa: ARG002
        minus_p: bool = True,
    ) -> None:
        p_str = "-p " if minus_p else ""
        cmd = f"mkdir {p_str}{shquote(fn)}"
        self._session.run(cmd)

    def _path_chmod(self, mode: int, fn: str) -> None:
        self._session.run(f"chmod {mode:o} {shquote(fn)}")

    def _path_touch(self, path: str) -> None:
        self._session.run(f"touch {shquote(path)}")

    def _path_chown(
        self,
        fn: str,
        owner: int | str | None,
        group: int | str | None,
        recursive: bool,
    ) -> None:
        args = ["chown"]
        if recursive:
            args.append("-R")
        if owner is not None and group is not None:
            args.append(f"{owner}:{group}")
        elif owner is not None:
            args.append(str(owner))
        elif group is not None:
            args.append(f":{group}")
        args.append(shquote(fn))
        self._session.run(" ".join(args))

    def _path_read(self, fn: str) -> bytes:
        data = self["cat"](fn)
        if self.custom_encoding and isinstance(data, str):
            return data.encode(self.custom_encoding)
        return typing.cast("bytes", data)

    def _path_write(self, fn: str, data: bytes | str) -> None:
        if self.custom_encoding and isinstance(data, str):
            data = data.encode(self.custom_encoding)
        assert isinstance(data, (bytes, bytearray))
        with NamedTemporaryFile() as f:
            f.write(data)
            f.flush()
            f.seek(0)
            self.upload(f.name, fn)

    def _path_link(self, src: str, dst: str, symlink: bool) -> None:
        symlink_str = "-s " if symlink else ""
        self._session.run(f"ln {symlink_str}{shquote(src)} {shquote(dst)}")

    def expand(self, expr: str) -> str:
        return self._session.run(f"echo {expr}")[1].strip()

    def expanduser(self, expr: str) -> str:
        # Only a leading ``~`` (i.e. ``~`` or ``~user`` as the first path
        # component) is meaningful, exactly like os.path.expanduser.
        if not expr.startswith("~"):
            return expr
        head, sep, tail = expr.partition("/")
        # Tilde expansion happens only on an *unquoted* leading ``~``/``~user``;
        # shquote would suppress it. A tilde-prefix can only name a login, so we
        # accept ``~`` or ``~user`` for a safe username and otherwise leave the
        # path untouched (matching os.path.expanduser, which returns the input
        # unchanged when the user is unknown/unexpandable).
        if head != "~" and not re.fullmatch(r"~[A-Za-z0-9_][A-Za-z0-9_.-]*", head):
            return expr
        expanded = self._session.run(f"echo {head}")[1].strip()
        # The remainder of the path is treated as a literal string; it is never
        # passed through the shell, so metacharacters cannot inject or expand.
        return expanded + sep + tail


class AsyncRemoteMachine:
    """Async version of BaseRemoteMachine.

    This class provides async access to remote commands via SSH.
    It wraps a sync RemoteMachine and provides async execution methods.

    Example::

        from plumbum.machines.ssh_machine import AsyncSshMachine

        async with AsyncSshMachine("host") as rem:
            ls = rem["ls"]
            result = await ls("-la")

    .. versionadded:: 2.0
    """

    __slots__ = ("_sync_machine",)

    def __init__(self, sync_machine: BaseRemoteMachine):
        """Initialize with a sync remote machine to wrap.

        Args:
            sync_machine: The sync remote machine to wrap
        """
        self._sync_machine = sync_machine

    def __getitem__(self, cmd: str | RemotePath | LocalPath) -> AsyncRemoteCommand:
        """Get an async remote command by name or path.

        This delegates to the sync machine for command lookup, then wraps it.

        Args:
            cmd: Command name (will be looked up in PATH) or RemotePath

        Returns:
            AsyncRemoteCommand instance

        Raises:
            CommandNotFound: If command is not found in PATH
        """
        from plumbum.commands.async_ import AsyncRemoteCommand

        sync_cmd = self._sync_machine[cmd]
        return AsyncRemoteCommand(sync_cmd)

    def __contains__(self, cmd: str) -> bool:
        """Check if a command exists in remote PATH."""
        return cmd in self._sync_machine

    @property
    def cwd(self) -> RemoteWorkdir:
        """Current working directory on remote machine."""
        return self._sync_machine.cwd

    @property
    def env(self) -> RemoteEnv:
        """Environment variables on remote machine."""
        return self._sync_machine.env

    def path(self, *parts: str | RemotePath | LocalPath) -> RemotePath:
        """Create a RemotePath from parts."""
        return self._sync_machine.path(*parts)

    def close(self) -> None:
        """Close the connection to the remote machine."""
        self._sync_machine.close()

    async def __aenter__(self) -> Self:
        """Async context manager entry."""
        return self

    async def __aexit__(self, t: object, v: object, tb: object) -> None:
        """Async context manager exit."""
        self.close()


__all__ = [
    "AsyncRemoteMachine",
    "BaseRemoteMachine",
    "ClosedRemote",
    "ClosedRemoteMachine",
    "RemoteCommand",
    "RemoteEnv",
]


def __dir__() -> list[str]:
    return list(__all__)
