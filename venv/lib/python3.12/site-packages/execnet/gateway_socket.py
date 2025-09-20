from __future__ import annotations

import sys
from typing import cast

from execnet.gateway import Gateway
from execnet.gateway_base import ExecModel
from execnet.gateway_bootstrap import HostNotFound
from execnet.multi import Group
from execnet.xspec import XSpec


class SocketIO:
    remoteaddress: str

    def __init__(self, sock, execmodel: ExecModel) -> None:
        self.sock = sock
        self.execmodel = execmodel
        socket = execmodel.socket
        try:
            # IPTOS_LOWDELAY
            sock.setsockopt(socket.SOL_IP, socket.IP_TOS, 0x10)
            sock.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
        except (AttributeError, OSError):
            sys.stderr.write("WARNING: cannot set socketoption")

    def read(self, numbytes: int) -> bytes:
        "Read exactly 'bytes' bytes from the socket."
        buf = b""
        while len(buf) < numbytes:
            t = self.sock.recv(numbytes - len(buf))
            if not t:
                raise EOFError
            buf += t
        return buf

    def write(self, data: bytes) -> None:
        self.sock.sendall(data)

    def close_read(self) -> None:
        try:
            self.sock.shutdown(0)
        except self.execmodel.socket.error:
            pass

    def close_write(self) -> None:
        try:
            self.sock.shutdown(1)
        except self.execmodel.socket.error:
            pass

    def wait(self) -> None:
        pass

    def kill(self) -> None:
        pass


def start_via(
    gateway: Gateway, hostport: tuple[str, int] | None = None
) -> tuple[str, int]:
    """Instantiate a socketserver on the given gateway.

    Returns a host, port tuple.
    """
    if hostport is None:
        host, port = ("localhost", 0)
    else:
        host, port = hostport

    from execnet.script import socketserver

    # execute the above socketserverbootstrap on the other side
    channel = gateway.remote_exec(socketserver)
    channel.send((host, port))
    realhost, realport = cast("tuple[str, int]", channel.receive())
    # self._trace("new_remote received"
    #               "port=%r, hostname = %r" %(realport, hostname))
    if not realhost or realhost == "0.0.0.0":
        realhost = "localhost"
    return realhost, realport


def create_io(spec: XSpec, group: Group, execmodel: ExecModel) -> SocketIO:
    assert spec.socket is not None
    assert not spec.python, "socket: specifying python executables not yet supported"
    gateway_id = spec.installvia
    if gateway_id:
        host, port = start_via(group[gateway_id])
    else:
        host, port_str = spec.socket.split(":")
        port = int(port_str)

    socket = execmodel.socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    io = SocketIO(sock, execmodel)
    io.remoteaddress = "%s:%d" % (host, port)
    try:
        sock.connect((host, port))
    except execmodel.socket.gaierror as e:
        raise HostNotFound() from e
    return io
