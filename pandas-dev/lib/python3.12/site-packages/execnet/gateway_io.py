"""execnet IO initialization code.

Creates IO instances used for gateway IO.
"""

from __future__ import annotations

import shlex
import sys
from typing import TYPE_CHECKING
from typing import cast

if TYPE_CHECKING:
    from execnet.gateway_base import Channel
    from execnet.gateway_base import ExecModel
    from execnet.xspec import XSpec

try:
    from execnet.gateway_base import Message
    from execnet.gateway_base import Popen2IO
except ImportError:
    from __main__ import Message  # type: ignore[no-redef]
    from __main__ import Popen2IO  # type: ignore[no-redef]

from functools import partial


class Popen2IOMaster(Popen2IO):
    # Set externally, for some specs only.
    remoteaddress: str

    def __init__(self, args, execmodel: ExecModel) -> None:
        PIPE = execmodel.subprocess.PIPE
        self.popen = p = execmodel.subprocess.Popen(args, stdout=PIPE, stdin=PIPE)
        super().__init__(p.stdin, p.stdout, execmodel=execmodel)

    def wait(self) -> int | None:
        try:
            return self.popen.wait()  # type: ignore[no-any-return]
        except OSError:
            return None

    def kill(self) -> None:
        try:
            self.popen.kill()
        except OSError as e:
            sys.stderr.write("ERROR killing: %s\n" % e)
            sys.stderr.flush()


popen_bootstrapline = "import sys;exec(eval(sys.stdin.readline()))"


def shell_split_path(path: str) -> list[str]:
    """
    Use shell lexer to split the given path into a list of components,
    taking care to handle Windows' '\' correctly.
    """
    if sys.platform.startswith("win"):
        # replace \\ by / otherwise shlex will strip them out
        path = path.replace("\\", "/")
    return shlex.split(path)


def popen_args(spec: XSpec) -> list[str]:
    args = shell_split_path(spec.python) if spec.python else [sys.executable]
    args.append("-u")
    if spec.dont_write_bytecode:
        args.append("-B")
    args.extend(["-c", popen_bootstrapline])
    return args


def ssh_args(spec: XSpec) -> list[str]:
    # NOTE: If changing this, you need to sync those changes to vagrant_args
    # as well, or, take some time to further refactor the commonalities of
    # ssh_args and vagrant_args.
    remotepython = spec.python or "python"
    args = ["ssh", "-C"]
    if spec.ssh_config is not None:
        args.extend(["-F", str(spec.ssh_config)])

    assert spec.ssh is not None
    args.extend(spec.ssh.split())
    remotecmd = f'{remotepython} -c "{popen_bootstrapline}"'
    args.append(remotecmd)
    return args


def vagrant_ssh_args(spec: XSpec) -> list[str]:
    # This is the vagrant-wrapped version of SSH. Unfortunately the
    # command lines are incompatible to just channel through ssh_args
    # due to ordering/templating issues.
    # NOTE: This should be kept in sync with the ssh_args behaviour.
    # spec.vagrant is identical to spec.ssh in that they both carry
    # the remote host "address".
    assert spec.vagrant_ssh is not None
    remotepython = spec.python or "python"
    args = ["vagrant", "ssh", spec.vagrant_ssh, "--", "-C"]
    if spec.ssh_config is not None:
        args.extend(["-F", str(spec.ssh_config)])
    remotecmd = f'{remotepython} -c "{popen_bootstrapline}"'
    args.extend([remotecmd])
    return args


def create_io(spec: XSpec, execmodel: ExecModel) -> Popen2IOMaster:
    if spec.popen:
        args = popen_args(spec)
        return Popen2IOMaster(args, execmodel)
    if spec.ssh:
        args = ssh_args(spec)
        io = Popen2IOMaster(args, execmodel)
        io.remoteaddress = spec.ssh
        return io
    if spec.vagrant_ssh:
        args = vagrant_ssh_args(spec)
        io = Popen2IOMaster(args, execmodel)
        io.remoteaddress = spec.vagrant_ssh
        return io
    assert False


#
# Proxy Gateway handling code
#
# master: proxy initiator
# forwarder: forwards between master and sub
# sub: sub process that is proxied to the initiator

RIO_KILL = 1
RIO_WAIT = 2
RIO_REMOTEADDRESS = 3
RIO_CLOSE_WRITE = 4


class ProxyIO:
    """A Proxy IO object allows to instantiate a Gateway
    through another "via" gateway.

    A master:ProxyIO object provides an IO object effectively connected to the
    sub via the forwarder. To achieve this, master:ProxyIO interacts with
    forwarder:serve_proxy_io() which itself instantiates and interacts with the
    sub.
    """

    def __init__(self, proxy_channel: Channel, execmodel: ExecModel) -> None:
        # after exchanging the control channel we use proxy_channel
        # for messaging IO
        self.controlchan = proxy_channel.gateway.newchannel()
        proxy_channel.send(self.controlchan)
        self.iochan = proxy_channel
        self.iochan_file = self.iochan.makefile("r")
        self.execmodel = execmodel

    def read(self, nbytes: int) -> bytes:
        # TODO(typing): The IO protocol requires bytes here but ChannelFileRead
        # returns str.
        return self.iochan_file.read(nbytes)  # type: ignore[return-value]

    def write(self, data: bytes) -> None:
        self.iochan.send(data)

    def _controll(self, event: int) -> object:
        self.controlchan.send(event)
        return self.controlchan.receive()

    def close_write(self) -> None:
        self._controll(RIO_CLOSE_WRITE)

    def close_read(self) -> None:
        raise NotImplementedError()

    def kill(self) -> None:
        self._controll(RIO_KILL)

    def wait(self) -> int | None:
        response = self._controll(RIO_WAIT)
        assert response is None or isinstance(response, int)
        return response

    @property
    def remoteaddress(self) -> str:
        response = self._controll(RIO_REMOTEADDRESS)
        assert isinstance(response, str)
        return response

    def __repr__(self) -> str:
        return f"<RemoteIO via {self.iochan.gateway.id}>"


class PseudoSpec:
    def __init__(self, vars) -> None:
        self.__dict__.update(vars)

    def __getattr__(self, name: str) -> None:
        return None


def serve_proxy_io(proxy_channelX: Channel) -> None:
    execmodel = proxy_channelX.gateway.execmodel
    log = partial(
        proxy_channelX.gateway._trace, "serve_proxy_io:%s" % proxy_channelX.id
    )
    spec = cast("XSpec", PseudoSpec(proxy_channelX.receive()))
    # create sub IO object which we will proxy back to our proxy initiator
    sub_io = create_io(spec, execmodel)
    control_chan = cast("Channel", proxy_channelX.receive())
    log("got control chan", control_chan)

    # read data from master, forward it to the sub
    # XXX writing might block, thus blocking the receiver thread
    def forward_to_sub(data: bytes) -> None:
        log("forward data to sub, size %s" % len(data))
        sub_io.write(data)

    proxy_channelX.setcallback(forward_to_sub)

    def control(data: int) -> None:
        if data == RIO_WAIT:
            control_chan.send(sub_io.wait())
        elif data == RIO_KILL:
            sub_io.kill()
            control_chan.send(None)
        elif data == RIO_REMOTEADDRESS:
            control_chan.send(sub_io.remoteaddress)
        elif data == RIO_CLOSE_WRITE:
            sub_io.close_write()
            control_chan.send(None)

    control_chan.setcallback(control)

    # write data to the master coming from the sub
    forward_to_master_file = proxy_channelX.makefile("w")

    # read bootstrap byte from sub, send it on to master
    log("reading bootstrap byte from sub", spec.id)
    initial = sub_io.read(1)
    assert initial == b"1", initial
    log("forwarding bootstrap byte from sub", spec.id)
    forward_to_master_file.write(initial)

    # enter message forwarding loop
    while True:
        try:
            message = Message.from_io(sub_io)
        except EOFError:
            log("EOF from sub, terminating proxying loop", spec.id)
            break
        message.to_io(forward_to_master_file)
    # proxy_channelX will be closed from remote_exec's finalization code


if __name__ == "__channelexec__":
    serve_proxy_io(channel)  # type: ignore[name-defined] # noqa:F821
