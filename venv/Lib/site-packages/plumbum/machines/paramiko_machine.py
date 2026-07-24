from __future__ import annotations

__lazy_modules__ = {
    "contextlib",
    "plumbum.commands",
    "plumbum.commands.processes",
    "plumbum.machines.session",
    "plumbum.path",
    "plumbum.path.local",
    "plumbum.path.remote",
}

import contextlib
import errno
import logging
import os
import stat
import typing
from typing import IO, Any

from plumbum.commands.base import shquote
from plumbum.commands.processes import ProcessLineTimedOut, iter_lines
from plumbum.machines.base import PopenAddons
from plumbum.machines.remote import BaseRemoteMachine, _check_env_name
from plumbum.machines.session import ShellSession
from plumbum.path.local import LocalPath
from plumbum.path.remote import RemotePath, RemoteStatRes

if typing.TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterator, Sequence

    import paramiko
    from paramiko.channel import Channel

    from plumbum.commands.base import BaseCommand
    from plumbum.machines.base import PopenWithAddons
else:
    try:
        # Sigh... we need to gracefully-import paramiko for Sphinx builds, etc
        import paramiko
    except ImportError:

        class paramiko:
            __slots__ = ()

            def __bool__(self):
                return False

            def __getattr__(self, name):
                raise ImportError("No module named paramiko")

        paramiko = paramiko()

logger = logging.getLogger("plumbum.paramiko")


class ParamikoPopen(PopenAddons):
    def __init__(
        self,
        argv: Sequence[str],
        stdin: paramiko.channel.ChannelFile,
        stdout: paramiko.channel.ChannelFile,
        stderr: paramiko.channel.ChannelFile,
        encoding: str,
        stdin_file: IO[str] | None = None,
        stdout_file: IO[str] | None = None,
        stderr_file: IO[str] | None = None,
    ) -> None:
        self.argv = argv
        self.channel = stdout.channel
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        self.custom_encoding = encoding
        self.returncode: int | None = None
        self.pid: int | None = None
        self.stdin_file = stdin_file
        self.stdout_file = stdout_file
        self.stderr_file = stderr_file

    def poll(self) -> int | None:
        if self.returncode is None and self.channel.exit_status_ready():
            return self.wait()
        return self.returncode

    def wait(self) -> int | None:
        if self.returncode is None:
            self.channel.recv_exit_status()
            self.returncode = self.channel.exit_status
            self.close()
        return self.returncode

    def close(self) -> None:
        self.channel.shutdown_read()
        self.channel.shutdown_write()
        self.channel.close()

    @staticmethod
    def kill() -> None:
        # possible way to obtain pid:
        # "(cmd ; echo $?) & echo ?!"
        # and then client.exec_command("kill -9 %s" % (pid,))
        raise OSError("Cannot kill remote processes, we don't have their PIDs")

    terminate = kill

    def send_signal(self, sig: int) -> None:
        raise NotImplementedError()

    def communicate(self) -> tuple[bytes, bytes]:
        stdout: list[str] = []
        stderr: list[str] = []
        infile = self.stdin_file
        sources = [
            ("1", stdout, self.stdout, self.stdout_file),
            ("2", stderr, self.stderr, self.stderr_file),
        ]
        i = 0
        while sources:
            if infile:
                try:
                    line = infile.readline()
                except (ValueError, OSError):
                    line = None
                logger.debug("communicate: %r", line)
                if not line:
                    infile.close()
                    infile = None
                    self.stdin.close()
                else:
                    self.stdin.write(line)
                    self.stdin.flush()

            i = (i + 1) % len(sources)
            _name, coll, pipe, outfile = sources[i]
            line = pipe.readline()
            # logger.debug("%s> %r", name, line)
            if not line:
                del sources[i]
            elif outfile:
                outfile.write(line)
                outfile.flush()
            else:
                coll.append(line)
        self.wait()
        assert self.custom_encoding is not None
        stdout_bytes = "".join(stdout).encode(self.custom_encoding)
        stderr_bytes = "".join(stderr).encode(self.custom_encoding)
        return stdout_bytes, stderr_bytes

    def iter_lines(
        self, timeout: float | None = None, **kwargs: Any
    ) -> Iterator[tuple[int, str]]:
        if timeout is not None:
            raise NotImplementedError(
                "The 'timeout' parameter is not supported with ParamikoMachine"
            )
        self_view = typing.cast("PopenWithAddons[Any]", self)
        return iter_lines(self_view, _iter_lines=_iter_lines, **kwargs)

    __iter__ = iter_lines


class RemoteCommand(BaseRemoteMachine.RemoteCommand):  # type: ignore[valid-type, misc]
    __slots__ = ()

    def __or__(self, *_: Any) -> Any:
        return NotImplemented


class ParamikoMachine(BaseRemoteMachine):
    """
    An implementation of :class:`remote machine <plumbum.machines.remote.BaseRemoteMachine>`
    over Paramiko (a Python implementation of openSSH2 client/server). Invoking a remote command
    translates to invoking it over SSH ::

        with ParamikoMachine("yourhostname") as rem:
            r_ls = rem["ls"]
            # r_ls is the remote `ls`
            # executing r_ls() is equivalent to `ssh yourhostname ls`, only without
            # spawning a new ssh client

    :param host: the host name to connect to (SSH server)

    :param user: the user to connect as (if ``None``, the default will be used)

    :param port: the server's port (if ``None``, the default will be used)

    :param password: the user's password (if a password-based authentication is to be performed)
                     (if ``None``, key-based authentication will be used)

    :param keyfile: the path to the identity file (if ``None``, the default will be used)

    :param load_system_host_keys: whether or not to load the system's host keys (from ``/etc/ssh``
                                  and ``~/.ssh``). The default is ``True``, which means Paramiko
                                  behaves much like the ``ssh`` command-line client

    :param missing_host_policy: the value passed to the underlying ``set_missing_host_key_policy``
                                of the client. The default is ``None``, which means
                                ``set_missing_host_key_policy`` is not invoked and paramiko's
                                default behavior (reject) is employed

    :param encoding: the remote machine's encoding (defaults to UTF8)

    :param look_for_keys: set to False to disable searching for discoverable
                          private key files in ``~/.ssh``

    :param connect_timeout: timeout for TCP connection

    .. note:: If Paramiko 1.15 or above is installed, can use GSS_API authentication

    :param bool gss_auth: ``True`` if you want to use GSS-API authentication

    :param bool gss_kex: Perform GSS-API Key Exchange and user authentication

    :param bool gss_deleg_creds: Delegate GSS-API client credentials or not

    :param str gss_host: The targets name in the kerberos database. default: hostname

    :param bool get_pty: Execute remote commands with allocated pseudo-tty. default: False

    :param bool load_system_ssh_config: read system SSH config for ProxyCommand configuration. default: False
    """

    __slots__ = ("_client", "_fqhost", "_get_pty", "_keep_alive", "_sftp", "host")

    RemoteCommand = RemoteCommand

    def __gt__(self, *_: Any) -> Any:
        return NotImplemented

    def __rshift__(self, *_: Any) -> Any:
        return NotImplemented

    def __ge__(self, *_: Any) -> Any:
        return NotImplemented

    def __lt__(self, *_: Any) -> Any:
        return NotImplemented

    def __lshift__(self, *_: Any) -> Any:
        return NotImplemented

    def __init__(
        self,
        host: str,
        user: str | None = None,
        port: int | None = None,
        password: str | None = None,
        keyfile: str | None = None,
        load_system_host_keys: bool = True,
        missing_host_policy: Any | None = None,
        encoding: str = "utf8",
        look_for_keys: bool | None = None,
        connect_timeout: float | None = None,
        keep_alive: int = 0,
        gss_auth: bool = False,
        gss_kex: bool | None = None,
        gss_deleg_creds: bool | None = None,
        gss_host: str | None = None,
        get_pty: bool = False,
        load_system_ssh_config: bool = False,
    ) -> None:
        self.host = host
        kwargs: dict[str, Any] = {}
        if user:
            self._fqhost = f"{user}@{host}"
            kwargs["username"] = user
        else:
            self._fqhost = host
        self._client = paramiko.SSHClient()
        if load_system_host_keys:
            self._client.load_system_host_keys()
        if port is not None:
            kwargs["port"] = port
        if keyfile is not None:
            kwargs["key_filename"] = keyfile
        if password is not None:
            kwargs["password"] = password
        if missing_host_policy is not None:
            self._client.set_missing_host_key_policy(missing_host_policy)
        if look_for_keys is not None:
            kwargs["look_for_keys"] = look_for_keys
        if connect_timeout is not None:
            kwargs["timeout"] = connect_timeout
        if gss_auth:
            kwargs["gss_auth"] = gss_auth
            kwargs["gss_kex"] = gss_kex
            kwargs["gss_deleg_creds"] = gss_deleg_creds
            if not gss_host:
                gss_host = host
            kwargs["gss_host"] = gss_host
        if load_system_ssh_config:
            ssh_config = paramiko.SSHConfig()
            with open(os.path.expanduser("~/.ssh/config"), encoding="utf-8") as f:
                ssh_config.parse(f)
            with contextlib.suppress(KeyError):
                hostConfig = ssh_config.lookup(host)
                kwargs["sock"] = paramiko.ProxyCommand(hostConfig["proxycommand"])
        self._client.connect(host, **kwargs)
        self._keep_alive = keep_alive
        self._sftp: paramiko.SFTPClient | None = None
        self._get_pty = get_pty
        BaseRemoteMachine.__init__(self, encoding, connect_timeout)

    def __str__(self) -> str:
        return f"paramiko://{self._fqhost}"

    def close(self) -> None:
        BaseRemoteMachine.close(self)
        self._client.close()

    @property
    def sftp(self) -> paramiko.SFTPClient:
        """
        Returns an SFTP client on top of the current SSH connection; it can be used to manipulate
        files directly, much like an interactive FTP/SFTP session
        """
        if not self._sftp:
            self._sftp = self._client.open_sftp()
        return self._sftp

    def session(
        self,
        isatty: bool = False,
        term: str = "vt100",
        width: int = 80,
        height: int = 24,
        *,
        new_session: bool = False,  # noqa: ARG002
    ) -> ShellSession:
        # new_session is ignored for ParamikoMachine
        trans = self._client.get_transport()
        assert trans is not None
        trans.set_keepalive(self._keep_alive)
        chan = trans.open_session()
        if isatty:
            chan.get_pty(term, width, height)
            chan.set_combine_stderr(True)
        chan.invoke_shell()
        stdin = chan.makefile("wb", -1)
        stdout = chan.makefile("rb", -1)
        stderr = chan.makefile_stderr("rb", -1)
        proc = ParamikoPopen(["<shell>"], stdin, stdout, stderr, self.custom_encoding)
        return ShellSession(proc, self.custom_encoding, isatty)

    def popen(  # type: ignore[override]
        self,
        args: BaseCommand,
        stdin: IO[str] | None = None,
        stdout: IO[str] | None = None,
        stderr: IO[str] | None = None,
        new_session: bool = False,  # noqa: ARG002
        env: dict[str, str] | None = None,
        cwd: str | LocalPath | None = None,
    ) -> ParamikoPopen:
        # new_session is ignored for ParamikoMachine
        argv = []
        envdelta = self.env.getdelta()
        if env:
            envdelta.update(env)
        argv.extend(["cd", shquote(str(cwd or self.cwd)), "&&"])
        if envdelta:
            argv.append("env")
            argv.extend(
                f"{_check_env_name(k)}={shquote(v)}" for k, v in envdelta.items()
            )
        argv.extend(args.formulate())
        cmdline = " ".join(argv)
        logger.debug(cmdline)
        si, so, se = self._client.exec_command(cmdline, 1, get_pty=self._get_pty)
        return ParamikoPopen(
            argv,
            si,
            so,
            se,
            self.custom_encoding,
            stdin_file=stdin,
            stdout_file=stdout,
            stderr_file=stderr,
        )

    def download(self, src: str | RemotePath, dst: str | LocalPath) -> None:
        if isinstance(src, LocalPath):
            raise TypeError(f"src of download cannot be {src!r}")
        if isinstance(src, RemotePath) and src.remote != self:
            raise TypeError(f"src {src!r} points to a different remote machine")
        if isinstance(dst, RemotePath):
            raise TypeError(f"dst of download cannot be {dst!r}")
        return self._download(
            src if isinstance(src, RemotePath) else self.path(src),
            dst if isinstance(dst, LocalPath) else LocalPath(dst),
        )

    def _download(self, src: RemotePath, dst: LocalPath) -> None:
        if src.is_dir():
            if not dst.exists():
                self.sftp.mkdir(str(dst))
            for fn in src:
                self._download(fn, dst / fn.name)
        elif dst.is_dir():
            self.sftp.get(str(src), str(dst / src.name))
        else:
            self.sftp.get(str(src), str(dst))

    def upload(self, src: str | LocalPath, dst: str | RemotePath) -> None:
        if isinstance(src, RemotePath):
            raise TypeError(f"src of upload cannot be {src!r}")
        if isinstance(dst, LocalPath):
            raise TypeError(f"dst of upload cannot be {dst!r}")
        if isinstance(dst, RemotePath) and dst.remote != self:
            raise TypeError(f"dst {dst!r} points to a different remote machine")
        return self._upload(
            src if isinstance(src, LocalPath) else LocalPath(src),
            dst if isinstance(dst, RemotePath) else self.path(dst),
        )

    def _upload(self, src: LocalPath, dst: RemotePath) -> None:
        if src.is_dir():
            if not dst.exists():
                self.sftp.mkdir(str(dst))
            for fn in src:
                self._upload(fn, dst / fn.name)
        elif dst.is_dir():
            self.sftp.put(str(src), str(dst / src.name))
        else:
            self.sftp.put(str(src), str(dst))

    def connect_sock(
        self, dport: int, dhost: str = "localhost", ipv6: bool = False
    ) -> SocketCompatibleChannel:
        """Returns a Paramiko ``Channel``, connected to dhost:dport on the remote machine.
        The ``Channel`` behaves like a regular socket; you can ``send`` and ``recv`` on it
        and the data will pass encrypted over SSH. Usage::

            mach = ParamikoMachine("myhost")
            sock = mach.connect_sock(12345)
            data = sock.recv(100)
            sock.send("foobar")
            sock.close()
        """
        if ipv6 and dhost == "localhost":
            dhost = "::1"
        srcaddr = ("::1", 0, 0, 0) if ipv6 else ("127.0.0.1", 0)
        trans = self._client.get_transport()
        assert trans is not None
        trans.set_keepalive(self._keep_alive)
        # TODO: I think the two args need to match, so ipv6 might be broken here
        chan = trans.open_channel("direct-tcpip", (dhost, dport), srcaddr)
        return SocketCompatibleChannel(chan)

    #
    # Path implementation
    #
    def _path_listdir(self, fn: str) -> list[str]:
        return self.sftp.listdir(fn)

    def _path_read(self, fn: str) -> bytes:
        f = self.sftp.open(fn, "rb")
        data = f.read()
        f.close()
        return data

    def _path_write(self, fn: str, data: str | bytes) -> None:
        if self.custom_encoding and isinstance(data, str):
            data = data.encode(self.custom_encoding)
        f = self.sftp.open(fn, "wb")
        f.write(data)
        f.close()

    def _path_stat(self, fn: str) -> RemoteStatRes | None:
        return self._sftp_stat(self.sftp.stat, fn)

    def _path_lstat(self, fn: str) -> RemoteStatRes | None:
        return self._sftp_stat(self.sftp.lstat, fn)

    @staticmethod
    def _sftp_stat(
        stat_fn: Callable[[str], paramiko.SFTPAttributes], fn: str
    ) -> RemoteStatRes | None:
        try:
            st = stat_fn(fn)
        except OSError as e:
            if e.errno == errno.ENOENT:
                return None
            raise
        res = RemoteStatRes(
            (
                st.st_mode,
                0,
                0,
                0,
                st.st_uid,
                st.st_gid,
                st.st_size,
                st.st_atime,
                st.st_mtime,
                0,
            )  # type: ignore[arg-type]
        )

        assert st.st_mode is not None
        res.text_mode = ""
        if stat.S_ISLNK(st.st_mode):
            res.text_mode = "symbolic link"
        elif stat.S_ISDIR(st.st_mode):
            res.text_mode = "directory"
        elif stat.S_ISREG(st.st_mode):
            res.text_mode = "regular file"
        return res

    def daemonic_popen(
        self,
        command: BaseCommand,
        cwd: str = "/",
        stdout: str | None = None,
        stderr: str | None = None,
        append: bool = True,
    ) -> PopenWithAddons[str]:
        raise NotImplementedError("This is not implemented on ParamikoMachine!")


###################################################################################################
# Make paramiko.Channel adhere to the socket protocol, namely, send and recv should fail
# when the socket has been closed
###################################################################################################
class SocketCompatibleChannel:
    __slots__ = ("_chan",)

    def __init__(self, chan: Channel) -> None:
        self._chan = chan

    def __getattr__(self, name: str) -> Any:
        return getattr(self._chan, name)

    def send(self, s: bytes) -> int:
        if self._chan.closed:
            raise OSError(errno.EBADF, "Bad file descriptor")
        return self._chan.send(s)

    def recv(self, count: int) -> bytes:
        if self._chan.closed:
            raise OSError(errno.EBADF, "Bad file descriptor")
        return self._chan.recv(count)


###################################################################################################
# Custom iter_lines for paramiko.Channel
###################################################################################################
def _iter_lines(
    proc: PopenWithAddons[Any],
    decode: Callable[[bytes], str],  # noqa: ARG001
    linesize: int,
    line_timeout: float | None = None,
) -> Generator[tuple[int, str], None, None]:
    from selectors import EVENT_READ, DefaultSelector

    real_proc = typing.cast("ParamikoPopen", proc)

    # Python 3.4+ implementation
    def selector() -> Generator[None, None, None]:
        sel = DefaultSelector()
        sel.register(real_proc.stdout.channel, EVENT_READ)
        while True:
            ready = sel.select(line_timeout)
            if not ready and line_timeout:
                raise ProcessLineTimedOut(
                    "popen line timeout expired",
                    getattr(real_proc, "argv", None),
                    getattr(real_proc, "machine", None),
                )
            for _key, _mask in ready:
                yield

    for _ in selector():
        if real_proc.stdout.channel.recv_ready():
            yield 0, real_proc.stdout.readline(linesize)
        if real_proc.stdout.channel.recv_stderr_ready():
            yield 1, real_proc.stderr.readline(linesize)
        if real_proc.poll() is not None:
            break

    for line in real_proc.stdout:
        yield 0, line
    for line in real_proc.stderr:
        yield 1, line


__all__ = [
    "ParamikoMachine",
    "ParamikoPopen",
    "RemoteCommand",
    "SocketCompatibleChannel",
]


def __dir__() -> list[str]:
    return list(__all__)
