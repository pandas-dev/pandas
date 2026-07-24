"""Code to initialize the remote side of a gateway once the IO is created."""

from __future__ import annotations

import inspect
import os

import execnet

from . import gateway_base
from .gateway_base import IO
from .xspec import XSpec

importdir = os.path.dirname(os.path.dirname(execnet.__file__))


class HostNotFound(Exception):
    pass


def bootstrap_import(io: IO, spec: XSpec) -> None:
    # Only insert the importdir into the path if we must.  This prevents
    # bugs where backports expect to be shadowed by the standard library on
    # newer versions of python but would instead shadow the standard library.
    sendexec(
        io,
        "import sys",
        "if %r not in sys.path:" % importdir,
        "    sys.path.insert(0, %r)" % importdir,
        "from execnet.gateway_base import serve, init_popen_io, get_execmodel",
        "sys.stdout.write('1')",
        "sys.stdout.flush()",
        "execmodel = get_execmodel(%r)" % spec.execmodel,
        "serve(init_popen_io(execmodel), id='%s-worker')" % spec.id,
    )
    s = io.read(1)
    assert s == b"1", repr(s)


def bootstrap_exec(io: IO, spec: XSpec) -> None:
    try:
        sendexec(
            io,
            inspect.getsource(gateway_base),
            "execmodel = get_execmodel(%r)" % spec.execmodel,
            "io = init_popen_io(execmodel)",
            "io.write('1'.encode('ascii'))",
            "serve(io, id='%s-worker')" % spec.id,
        )
        s = io.read(1)
        assert s == b"1"
    except EOFError:
        ret = io.wait()
        if ret == 255 and hasattr(io, "remoteaddress"):
            raise HostNotFound(io.remoteaddress) from None


def bootstrap_socket(io: IO, id) -> None:
    # XXX: switch to spec
    from execnet.gateway_socket import SocketIO

    sendexec(
        io,
        inspect.getsource(gateway_base),
        "import socket",
        inspect.getsource(SocketIO),
        "try: execmodel",
        "except NameError:",
        "   execmodel = get_execmodel('thread')",
        "io = SocketIO(clientsock, execmodel)",
        "io.write('1'.encode('ascii'))",
        "serve(io, id='%s-worker')" % id,
    )
    s = io.read(1)
    assert s == b"1"


def sendexec(io: IO, *sources: str) -> None:
    source = "\n".join(sources)
    io.write((repr(source) + "\n").encode("utf-8"))


def bootstrap(io: IO, spec: XSpec) -> execnet.Gateway:
    if spec.popen:
        if spec.via or spec.python:
            bootstrap_exec(io, spec)
        else:
            bootstrap_import(io, spec)
    elif spec.ssh or spec.vagrant_ssh:
        bootstrap_exec(io, spec)
    elif spec.socket:
        bootstrap_socket(io, spec)
    else:
        raise ValueError("unknown gateway type, can't bootstrap")
    gw = execnet.Gateway(io, spec)
    return gw
