from __future__ import annotations

import re
import socket
import warnings
from contextlib import closing

from plumbum.commands import ProcessExecutionError, shquote
from plumbum.lib import IS_WIN32
from plumbum.machines.local import local
from plumbum.machines.remote import BaseRemoteMachine
from plumbum.machines.session import ShellSession
from plumbum.path.local import LocalPath
from plumbum.path.remote import RemotePath


def _get_free_port():
    """Attempts to find a free port."""
    s = socket.socket()
    with closing(s):
        s.bind(("localhost", 0))
        return s.getsockname()[1]


class SshTunnel:
    """An object representing an SSH tunnel (created by
    :func:`SshMachine.tunnel <plumbum.machines.remote.SshMachine.tunnel>`)"""

    __slots__ = ["_session", "_lport", "_dport", "_reverse", "__weakref__"]

    def __init__(self, session, lport, dport, reverse):
        self._session = session
        self._lport = lport
        self._dport = dport
        self._reverse = reverse
        if reverse and str(dport) == "0" and session._startup_result is not None:
            # Try to detect assigned remote port.
            regex = re.compile(
                r"^Allocated port (\d+) for remote forward to .+$", re.MULTILINE
            )
            match = regex.search(session._startup_result[2])
            if match:
                self._dport = match.group(1)

    def __repr__(self):
        tunnel = self._session.proc if self._session.alive() else "(defunct)"
        return f"<SshTunnel {tunnel}>"

    def __enter__(self):
        return self

    def __exit__(self, t, v, tb):
        self.close()

    def close(self):
        """Closes(terminates) the tunnel"""
        self._session.close()

    @property
    def lport(self):
        """Tunneled port or socket on the local machine."""
        return self._lport

    @property
    def dport(self):
        """Tunneled port or socket on the remote machine."""
        return self._dport

    @property
    def reverse(self):
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

    def __init__(
        self,
        host,
        user=None,
        port=None,
        keyfile=None,
        ssh_command=None,
        scp_command=None,
        ssh_opts=(),
        scp_opts=(),
        password=None,
        encoding="utf8",
        connect_timeout=10,
        new_session=False,
    ):
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
        BaseRemoteMachine.__init__(
            self,
            encoding=encoding,
            connect_timeout=connect_timeout,
            new_session=new_session,
        )

    def __str__(self):
        return f"ssh://{self._fqhost}"

    def popen(self, args, ssh_opts=(), env=None, cwd=None, **kwargs):
        cmdline = []
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
                cmdline.extend(["cd", str(cwd), "&&"])
            if envdelta:
                cmdline.append("env")
                cmdline.extend(f"{k}={shquote(v)}" for k, v in envdelta.items())
            if isinstance(args, (tuple, list)):
                cmdline.extend(args)
            else:
                cmdline.append(args)
        return self._ssh_command[tuple(cmdline)].popen(**kwargs)

    def nohup(self, command):
        """
        Runs the given command using ``nohup`` and redirects std handles,
        allowing the command to run "detached" from its controlling TTY or parent.
        Does not return anything. Depreciated (use command.nohup or daemonic_popen).
        """
        warnings.warn(
            "Use .nohup on the command or use daemonic_popen)",
            FutureWarning,
            stacklevel=2,
        )
        self.daemonic_popen(command, cwd=".", stdout=None, stderr=None, append=False)

    def daemonic_popen(self, command, cwd=".", stdout=None, stderr=None, append=True):
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

        args = [] if str(cwd) == "." else ["cd", str(cwd), "&&"]
        args.append("nohup")
        args.extend(command.formulate())
        args.extend(
            [
                (">>" if append else ">") + str(stdout),
                "2" + (">>" if (append and stderr != "&1") else ">") + str(stderr),
                "</dev/null",
            ]
        )
        proc = self.popen(args, ssh_opts=["-f"])
        rc = proc.wait()
        try:
            if rc != 0:
                raise ProcessExecutionError(
                    args, rc, proc.stdout.read(), proc.stderr.read()
                )
        finally:
            proc.stdin.close()
            proc.stdout.close()
            proc.stderr.close()

    def session(self, isatty=False, new_session=False):
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
        lport,
        dport,
        lhost="localhost",
        dhost="localhost",
        connect_timeout=5,  # noqa: ARG002
        reverse=False,
    ):
        r"""Creates an SSH tunnel from the TCP port (``lport``) of the local machine
        (``lhost``, defaults to ``"localhost"``, but it can be any IP you can ``bind()``)
        to the remote TCP port (``dport``) of the destination machine (``dhost``, defaults
        to ``"localhost"``, which means *this remote machine*). This function also
        supports Unix sockets, in which case the local socket should be passed in as
        ``lport`` and the local bind address should be ``None``. The same can be done
        for a remote socket, by following the same pattern with ``dport`` and ``dhost``.
        The returned :class:`SshTunnel <plumbum.machines.remote.SshTunnel>` object can
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
        return SshTunnel(
            ShellSession(
                proc, self.custom_encoding, connect_timeout=self.connect_timeout
            ),
            lport,
            dport,
            reverse,
        )

    @staticmethod
    def _translate_drive_letter(path):
        # replace c:\some\path with /c/some/path
        path = str(path)
        if ":" in path:
            return "/" + path.replace(":", "").replace("\\", "/")
        return path

    def download(self, src, dst):
        if isinstance(src, LocalPath):
            raise TypeError(f"src of download cannot be {src!r}")
        if isinstance(src, RemotePath) and src.remote != self:
            raise TypeError(f"src {src!r} points to a different remote machine")
        if isinstance(dst, RemotePath):
            raise TypeError(f"dst of download cannot be {dst!r}")
        if IS_WIN32:
            src = self._translate_drive_letter(src)
            dst = self._translate_drive_letter(dst)
        self._scp_command(f"{self._fqhost}:{shquote(src)}", dst)

    def upload(self, src, dst):
        if isinstance(src, RemotePath):
            raise TypeError(f"src of upload cannot be {src!r}")
        if isinstance(dst, LocalPath):
            raise TypeError(f"dst of upload cannot be {dst!r}")
        if isinstance(dst, RemotePath) and dst.remote != self:
            raise TypeError(f"dst {dst!r} points to a different remote machine")
        if IS_WIN32:
            src = self._translate_drive_letter(src)
            dst = self._translate_drive_letter(dst)
        self._scp_command(src, f"{self._fqhost}:{shquote(dst)}")


class PuttyMachine(SshMachine):
    """
    PuTTY-flavored SSH connection. The programs ``plink`` and ``pscp`` are expected to
    be in the path (or you may provide your own ``ssh_command`` and ``scp_command``)

    Arguments are the same as for :class:`plumbum.machines.remote.SshMachine`
    """

    def __init__(
        self,
        host,
        user=None,
        port=None,
        keyfile=None,
        ssh_command=None,
        scp_command=None,
        ssh_opts=(),
        scp_opts=(),
        encoding="utf8",
        connect_timeout=10,
        new_session=False,
    ):
        if ssh_command is None:
            ssh_command = local["plink"]
        if scp_command is None:
            scp_command = local["pscp"]
        if not ssh_opts:
            ssh_opts = ["-ssh"]
        if user is None:
            user = local.env.user
        if port is not None:
            ssh_opts.extend(["-P", str(port)])
            scp_opts = [*list(scp_opts), "-P", str(port)]
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

    def __str__(self):
        return f"putty-ssh://{self._fqhost}"

    def _translate_drive_letter(self, path):
        # pscp takes care of windows paths automatically
        return path

    def session(self, isatty=False, new_session=False):
        return ShellSession(
            self.popen((), (["-t"] if isatty else ["-T"]), new_session=new_session),
            self.custom_encoding,
            isatty,
            self.connect_timeout,
        )
