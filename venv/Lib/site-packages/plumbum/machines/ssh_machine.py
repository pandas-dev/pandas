from __future__ import annotations

__lazy_modules__ = {
    "contextlib",
    "plumbum.commands",
    "plumbum.lib",
    "plumbum.machines.local",
    "plumbum.machines.session",
    "plumbum.path",
    "plumbum.path.local",
    "plumbum.path.remote",
    "re",
    "socket",
}

import re
import socket
from contextlib import closing
from typing import TYPE_CHECKING, Any

from plumbum.commands import BaseCommand, ProcessExecutionError, shquote
from plumbum.lib import IS_WIN32
from plumbum.machines.local import local
from plumbum.machines.remote import BaseRemoteMachine, RemoteEnv, _check_env_name
from plumbum.machines.session import ShellSession
from plumbum.path.local import LocalPath
from plumbum.path.remote import RemotePath, RemoteWorkdir

if TYPE_CHECKING:
    from collections.abc import Sequence

    from plumbum._compat.typing import Self
    from plumbum.commands.async_ import AsyncRemoteCommand
    from plumbum.machines.base import PopenWithAddons


def _get_free_port() -> int:
    """Attempts to find a free port."""
    s = socket.socket()
    with closing(s):
        s.bind(("localhost", 0))
        port: int = s.getsockname()[1]
        return port


class SshTunnel:
    """An object representing an SSH tunnel (created by
    :func:`SshMachine.tunnel <plumbum.machines.ssh_machine.SshMachine.tunnel>`)"""

    __slots__ = ["__weakref__", "_dport", "_lport", "_reverse", "_session"]

    def __init__(
        self,
        session: ShellSession,
        lport: int | str,
        dport: int | str,
        reverse: bool,
    ) -> None:
        self._session = session
        self._lport = lport
        self._dport = dport
        self._reverse = reverse
        if reverse and str(dport) == "0" and session._startup_result is not None:
            # Try to detect assigned remote port.
            regex = re.compile(
                r"^Allocated port (\d+) for remote forward to .+$", re.MULTILINE
            )
            match = regex.search(str(session._startup_result[2]))
            if match:
                self._dport = match.group(1)

    def __repr__(self) -> str:
        tunnel = self._session.proc if self._session.alive() else "(defunct)"
        return f"<SshTunnel {tunnel}>"

    def __enter__(self) -> Self:
        return self

    def __exit__(self, t: object, v: object, tb: object) -> None:
        self.close()

    def close(self) -> None:
        """Closes(terminates) the tunnel"""
        self._session.close()

    @property
    def lport(self) -> int | str:
        """Tunneled port or socket on the local machine."""
        return self._lport

    @property
    def dport(self) -> int | str:
        """Tunneled port or socket on the remote machine."""
        return self._dport

    @property
    def reverse(self) -> bool:
        """Represents if the tunnel is a reverse tunnel."""
        return self._reverse


class SshMachine(BaseRemoteMachine):
    """
    An implementation of :class:`remote machine <plumbum.machines.remote.BaseRemoteMachine>`
    over SSH. Invoking a remote command translates to invoking it over SSH ::

        with SshMachine("yourhostname") as rem:
            r_ls = rem["ls"]
            # r_ls is the remote `ls`
            # executing r_ls() translates to `ssh yourhostname ls`

    :param host: the host name to connect to (SSH server)

    :param user: the user to connect as (if ``None``, the default will be used)

    :param port: the server's port (if ``None``, the default will be used)

    :param keyfile: the path to the identity file (if ``None``, the default will be used)

    :param ssh_command: the ``ssh`` command to use; this has to be a ``Command`` object;
                        if ``None``, the default ssh client will be used.

    :param scp_command: the ``scp`` command to use; this has to be a ``Command`` object;
                        if ``None``, the default scp program will be used.

    :param ssh_opts: any additional options for ``ssh`` (a list of strings)

    :param scp_opts: any additional options for ``scp`` (a list of strings)

    :param password: the password to use; requires ``sshpass`` be installed. Cannot be used
                     in conjunction with ``ssh_command`` or ``scp_command`` (will be ignored).
                     NOTE: THIS IS A SECURITY RISK!

    :param encoding: the remote machine's encoding (defaults to UTF8)

    :param connect_timeout: specify a connection timeout (the time until shell prompt is seen).
                            The default is 10 seconds. Set to ``None`` to disable

    :param new_session: whether or not to start the background session as a new
                        session leader (setsid). This will prevent it from being killed on
                        Ctrl+C (SIGINT)
    """

    __slots__ = ("_fqhost", "_scp_command", "_scp_translate", "_ssh_command", "host")

    def __init__(
        self,
        host: str,
        user: str | None = None,
        port: int | None = None,
        keyfile: str | None = None,
        ssh_command: BaseCommand | None = None,
        scp_command: BaseCommand | None = None,
        ssh_opts: Sequence[str] = (),
        scp_opts: Sequence[str] = (),
        password: str | None = None,
        encoding: str = "utf8",
        connect_timeout: float | None = 10,
        new_session: bool = False,
    ) -> None:
        if ssh_command is None:
            if password is not None:
                ssh_command = local["sshpass"]["-p", password, "ssh"]
            else:
                ssh_command = local["ssh"]
        if scp_command is None:
            if password is not None:
                scp_command = local["sshpass"]["-p", password, "scp"]
            else:
                scp_command = local["scp"]

        scp_args = []
        ssh_args = []
        self.host = host
        if user:
            self._fqhost = f"{user}@{host}"
        else:
            self._fqhost = host
        if port:
            ssh_args.extend(["-p", str(port)])
            scp_args.extend(["-P", str(port)])
        if keyfile:
            ssh_args.extend(["-i", str(keyfile)])
            scp_args.extend(["-i", str(keyfile)])
        scp_args.append("-r")
        ssh_args.extend(ssh_opts)
        scp_args.extend(scp_opts)
        self._ssh_command = ssh_command[tuple(ssh_args)]
        self._scp_command = scp_command[tuple(scp_args)]
        self._scp_translate = IS_WIN32 and self._scp_uses_posix_paths(self._scp_command)
        BaseRemoteMachine.__init__(
            self,
            encoding=encoding,
            connect_timeout=connect_timeout,
            new_session=new_session,
        )

    def __str__(self) -> str:
        return f"ssh://{self._fqhost}"

    def popen(
        self,
        args: Sequence[str] | str,
        ssh_opts: Sequence[str] = (),
        env: dict[str, str] | None = None,
        cwd: str | LocalPath | None = None,
        **kwargs: Any,
    ) -> PopenWithAddons[str]:
        cmdline: list[str] = []
        cmdline.extend(ssh_opts)
        cmdline.append(self._fqhost)
        if args:
            envdelta = {}
            if hasattr(self, "env"):
                envdelta.update(self.env.getdelta())
            if env:
                envdelta.update(env)
            if cwd is None:
                cwd = getattr(self, "cwd", None)
            if cwd:
                cmdline.extend(["cd", shquote(str(cwd)), "&&"])
            if envdelta:
                cmdline.append("env")
                cmdline.extend(
                    f"{_check_env_name(k)}={shquote(v)}" for k, v in envdelta.items()
                )
            if isinstance(args, (tuple, list)):
                cmdline.extend(args)
            else:
                cmdline.append(args)  # type: ignore[arg-type]
        return self._ssh_command[tuple(cmdline)].popen(**kwargs)

    def daemonic_popen(
        self,
        command: BaseCommand,
        cwd: str = ".",
        stdout: str | None = None,
        stderr: str | None = None,
        append: bool = True,
    ) -> PopenWithAddons[str]:
        """
        Runs the given command using ``nohup`` and redirects std handles,
        allowing the command to run "detached" from its controlling TTY or parent.
        Does not return anything.

        .. versionadded:: 1.6.0

        """
        if stdout is None:
            stdout = "/dev/null"
        if stderr is None:
            stderr = "&1"

        args = [] if str(cwd) == "." else ["cd", shquote(str(cwd)), "&&"]
        args.append("nohup")
        args.extend(command.formulate())
        # ``&1`` (i.e. ``2>&1``) is a shell redirection operand, not a path, so
        # it must stay unquoted; everything else is a filename and is quoted.
        stderr_target = "&1" if stderr == "&1" else shquote(str(stderr))
        args.extend(
            [
                (">>" if append else ">") + shquote(str(stdout)),
                "2" + (">>" if (append and stderr != "&1") else ">") + stderr_target,
                "</dev/null",
            ]
        )
        proc = self.popen(args, ssh_opts=["-f"])
        rc = proc.wait()
        assert proc.stdin is not None
        assert proc.stdout is not None
        assert proc.stderr is not None
        try:
            if rc != 0:
                raise ProcessExecutionError(
                    args, rc, proc.stdout.read(), proc.stderr.read()
                )
        finally:
            proc.stdin.close()
            proc.stdout.close()
            proc.stderr.close()
        return proc

    def session(self, isatty: bool = False, new_session: bool = False) -> ShellSession:
        return ShellSession(
            self.popen(
                ["/bin/sh"], (["-tt"] if isatty else ["-T"]), new_session=new_session
            ),
            self.custom_encoding,
            isatty,
            self.connect_timeout,
            host=self.host,
        )

    def tunnel(
        self,
        lport: int | str,
        dport: int | str,
        lhost: str | None = "localhost",
        dhost: str | None = "localhost",
        connect_timeout: float | None = None,
        reverse: bool = False,
    ) -> SshTunnel:
        r"""Creates an SSH tunnel from the TCP port (``lport``) of the local machine
        (``lhost``, defaults to ``"localhost"``, but it can be any IP you can ``bind()``)
        to the remote TCP port (``dport``) of the destination machine (``dhost``, defaults
        to ``"localhost"``, which means *this remote machine*). This function also
        supports Unix sockets, in which case the local socket should be passed in as
        ``lport`` and the local bind address should be ``None``. The same can be done
        for a remote socket, by following the same pattern with ``dport`` and ``dhost``.
        The returned :class:`SshTunnel <plumbum.machines.ssh_machine.SshTunnel>` object can
        be used as a *context-manager*.

        The more conventional use case is the following::

            +---------+          +---------+
            | Your    |          | Remote  |
            | Machine |          | Machine |
            +----o----+          +---- ----+
                 |                    ^
                 |                    |
               lport                dport
                 |                    |
                 \______SSH TUNNEL____/
                        (secure)

        Here, you wish to communicate safely between port ``lport`` of your machine and
        port ``dport`` of the remote machine. Communication is tunneled over SSH, so the
        connection is authenticated and encrypted.

        The more general case is shown below (where ``dport != "localhost"``)::

            +---------+          +-------------+      +-------------+
            | Your    |          | Remote      |      | Destination |
            | Machine |          | Machine     |      | Machine     |
            +----o----+          +---- ----o---+      +---- --------+
                 |                    ^    |               ^
                 |                    |    |               |
            lhost:lport               |    |          dhost:dport
                 |                    |    |               |
                 \_____SSH TUNNEL_____/    \_____SOCKET____/
                        (secure)              (not secure)

        Usage::

            rem = SshMachine("megazord")

            with rem.tunnel(1234, "/var/lib/mysql/mysql.sock", dhost=None):
                sock = socket.socket()
                sock.connect(("localhost", 1234))
                # sock is now tunneled to the MySQL socket on megazord

        The ``connect_timeout`` is the time to wait for the tunnel's shell prompt;
        if not given, the machine-level ``connect_timeout`` is used.
        """
        formatted_lhost = "" if lhost is None else f"[{lhost}]:"
        formatted_dhost = "" if dhost is None else f"[{dhost}]:"
        if str(lport) == "0":
            lport = _get_free_port()
        ssh_opts = (
            [
                "-L",
                f"{formatted_lhost}{lport}:{formatted_dhost}{dport}",
            ]
            if not reverse
            else [
                "-R",
                f"{formatted_dhost}{dport}:{formatted_lhost}{lport}",
            ]
        )
        proc = self.popen((), ssh_opts=ssh_opts, new_session=True)
        if connect_timeout is None:
            connect_timeout = self.connect_timeout
        return SshTunnel(
            ShellSession(proc, self.custom_encoding, connect_timeout=connect_timeout),
            lport,
            dport,
            reverse,
        )

    @staticmethod
    def _translate_drive_letter(path: str | LocalPath) -> str:
        # replace c:\some\path with /c/some/path
        path = str(path)
        if ":" in path:
            return "/" + path.replace(":", "").replace("\\", "/")
        return path

    @staticmethod
    def _scp_uses_posix_paths(scp_command: BaseCommand) -> bool:
        """Whether ``scp_command`` expects POSIX-style ``/c/...`` local paths.

        Cygwin- and MSYS2-based builds of ``scp`` (such as the one bundled with
        Git for Windows) cannot make sense of native ``c:\\...`` paths and need
        drive letters rewritten as ``/c/...``. They are recognisable by the
        ``cygwin1.dll`` / ``msys-2.0.dll`` runtime that sits next to the
        executable. The native Windows OpenSSH ``scp`` accepts drive-letter
        paths directly, so for it (and on POSIX hosts) no translation is needed.
        """
        try:
            cmd: Any = scp_command
            while not hasattr(cmd, "executable"):
                cmd = getattr(cmd, "cmd", None)
                if cmd is None:
                    return False
            folder = local.path(cmd.executable).dirname
            return any(
                (folder / dll).exists() for dll in ("cygwin1.dll", "msys-2.0.dll")
            )
        except (AttributeError, OSError):
            return False

    def download(self, src: str | RemotePath, dst: str | LocalPath) -> None:
        if isinstance(src, LocalPath):
            raise TypeError(f"src of download cannot be {src!r}")
        if isinstance(src, RemotePath) and src.remote != self:
            raise TypeError(f"src {src!r} points to a different remote machine")
        if isinstance(dst, RemotePath):
            raise TypeError(f"dst of download cannot be {dst!r}")
        if self._scp_translate:
            # only the local path (dst) is parsed by the local scp binary
            dst = self._translate_drive_letter(dst)
        self._scp_command(f"{self._fqhost}:{shquote(src)}", dst)

    def upload(self, src: str | LocalPath, dst: str | RemotePath) -> None:
        if isinstance(src, RemotePath):
            raise TypeError(f"src of upload cannot be {src!r}")
        if isinstance(dst, LocalPath):
            raise TypeError(f"dst of upload cannot be {dst!r}")
        if isinstance(dst, RemotePath) and dst.remote != self:
            raise TypeError(f"dst {dst!r} points to a different remote machine")
        if self._scp_translate:
            # only the local path (src) is parsed by the local scp binary
            src = self._translate_drive_letter(src)
        self._scp_command(src, f"{self._fqhost}:{shquote(dst)}")


class PuttyMachine(SshMachine):
    """
    PuTTY-flavored SSH connection. The programs ``plink`` and ``pscp`` are expected to
    be in the path (or you may provide your own ``ssh_command`` and ``scp_command``)

    Arguments are the same as for :class:`plumbum.machines.ssh_machine.SshMachine`
    """

    __slots__ = ()

    def __init__(
        self,
        host: str,
        user: str | None = None,
        port: int | None = None,
        keyfile: str | None = None,
        ssh_command: BaseCommand | None = None,
        scp_command: BaseCommand | None = None,
        ssh_opts: Sequence[str] = (),
        scp_opts: Sequence[str] = (),
        encoding: str = "utf8",
        connect_timeout: float | None = 10,
        new_session: bool = False,
    ) -> None:
        if ssh_command is None:
            ssh_command = local["plink"]
        if scp_command is None:
            scp_command = local["pscp"]
        if not ssh_opts:
            ssh_opts = ["-ssh"]
        if user is None:
            user = local.env.user
        if port is not None:
            ssh_opts = [*ssh_opts, "-P", str(port)]
            scp_opts = [*scp_opts, "-P", str(port)]
            port = None
        SshMachine.__init__(
            self,
            host,
            user,
            port,
            keyfile=keyfile,
            ssh_command=ssh_command,
            scp_command=scp_command,
            ssh_opts=ssh_opts,
            scp_opts=scp_opts,
            encoding=encoding,
            connect_timeout=connect_timeout,
            new_session=new_session,
        )

    def __str__(self) -> str:
        return f"putty-ssh://{self._fqhost}"

    @staticmethod
    def _translate_drive_letter(path: str | LocalPath) -> str:
        # pscp takes care of windows paths automatically
        return str(path)

    def session(self, isatty: bool = False, new_session: bool = False) -> ShellSession:
        return ShellSession(
            self.popen((), (["-t"] if isatty else ["-T"]), new_session=new_session),
            self.custom_encoding,
            isatty,
            self.connect_timeout,
        )


class AsyncSshMachine:
    """Async version of SshMachine.

    This class provides async SSH command execution. It wraps a sync SshMachine
    and provides async execution methods.

    Example::

        from plumbum.machines.ssh_machine import AsyncSshMachine

        async with AsyncSshMachine("hostname") as rem:
            result = await rem["ls"]("-la")

            # Concurrent execution
            results = await asyncio.gather(
                rem["echo"]("task1"),
                rem["echo"]("task2"),
            )

    .. versionadded:: 2.0
    """

    __slots__ = ("_sync_machine",)

    def __init__(
        self,
        host: str,
        user: str | None = None,
        port: int | None = None,
        keyfile: str | None = None,
        ssh_command: BaseCommand | None = None,
        scp_command: BaseCommand | None = None,
        ssh_opts: Sequence[str] = (),
        scp_opts: Sequence[str] = (),
        password: str | None = None,
        encoding: str = "utf8",
        connect_timeout: float | None = 10,
        new_session: bool = False,
    ) -> None:
        """Initialize async SSH machine.

        Args:
            host: The host name to connect to (SSH server)
            user: The user to connect as (if None, the default will be used)
            port: The server's port (if None, the default will be used)
            keyfile: The path to the identity file (if None, the default will be used)
            ssh_command: The ssh command to use (if None, the default will be used)
            scp_command: The scp command to use (if None, the default will be used)
            ssh_opts: Any additional options for ssh (a list of strings)
            scp_opts: Any additional options for scp (a list of strings)
            password: The password to use (requires sshpass)
            encoding: The remote machine's encoding (defaults to UTF8)
            connect_timeout: Connection timeout (default 10 seconds)
            new_session: Whether to start as a new session leader
        """
        sync_machine = SshMachine(
            host=host,
            user=user,
            port=port,
            keyfile=keyfile,
            ssh_command=ssh_command,
            scp_command=scp_command,
            ssh_opts=ssh_opts,
            scp_opts=scp_opts,
            password=password,
            encoding=encoding,
            connect_timeout=connect_timeout,
            new_session=new_session,
        )
        self._sync_machine = sync_machine

    def __getitem__(self, cmd: str | RemotePath | LocalPath) -> AsyncRemoteCommand:
        """Get an async remote command by name or path."""
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
    "AsyncSshMachine",
    "PuttyMachine",
    "SshMachine",
    "SshTunnel",
]


def __dir__() -> list[str]:
    return list(__all__)
