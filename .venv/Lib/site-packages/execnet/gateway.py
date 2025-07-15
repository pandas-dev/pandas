"""Gateway code for initiating popen, socket and ssh connections.

(c) 2004-2013, Holger Krekel and others
"""

from __future__ import annotations

import inspect
import linecache
import textwrap
import types
from typing import TYPE_CHECKING
from typing import Any
from typing import Callable

from . import gateway_base
from .gateway_base import IO
from .gateway_base import Channel
from .gateway_base import Message
from .multi import Group
from .xspec import XSpec


class Gateway(gateway_base.BaseGateway):
    """Gateway to a local or remote Python Interpreter."""

    _group: Group

    def __init__(self, io: IO, spec: XSpec) -> None:
        """:private:"""
        super().__init__(io=io, id=spec.id, _startcount=1)
        self.spec = spec
        self._initreceive()

    @property
    def remoteaddress(self) -> str:
        # Only defined for remote IO types.
        return self._io.remoteaddress  # type: ignore[attr-defined,no-any-return]

    def __repr__(self) -> str:
        """A string representing gateway type and status."""
        try:
            r: str = self.hasreceiver() and "receive-live" or "not-receiving"
            i = str(len(self._channelfactory.channels()))
        except AttributeError:
            r = "uninitialized"
            i = "no"
        return f"<{self.__class__.__name__} id={self.id!r} {r}, {self.execmodel.backend} model, {i} active channels>"

    def exit(self) -> None:
        """Trigger gateway exit.

        Defer waiting for finishing of receiver-thread and subprocess activity
        to when group.terminate() is called.
        """
        self._trace("gateway.exit() called")
        if self not in self._group:
            self._trace("gateway already unregistered with group")
            return
        self._group._unregister(self)
        try:
            self._trace("--> sending GATEWAY_TERMINATE")
            self._send(Message.GATEWAY_TERMINATE)
            self._trace("--> io.close_write")
            self._io.close_write()
        except (ValueError, EOFError, OSError) as exc:
            self._trace("io-error: could not send termination sequence")
            self._trace(" exception: %r" % exc)

    def reconfigure(
        self, py2str_as_py3str: bool = True, py3str_as_py2str: bool = False
    ) -> None:
        """Set the string coercion for this gateway.

        The default is to try to convert py2 str as py3 str, but not to try and
        convert py3 str to py2 str.
        """
        self._strconfig = (py2str_as_py3str, py3str_as_py2str)
        data = gateway_base.dumps_internal(self._strconfig)
        self._send(Message.RECONFIGURE, data=data)

    def _rinfo(self, update: bool = False) -> RInfo:
        """Return some sys/env information from remote."""
        if update or not hasattr(self, "_cache_rinfo"):
            ch = self.remote_exec(rinfo_source)
            try:
                self._cache_rinfo = RInfo(ch.receive())
            finally:
                ch.waitclose()
        return self._cache_rinfo

    def hasreceiver(self) -> bool:
        """Whether gateway is able to receive data."""
        return self._receivepool.active_count() > 0

    def remote_status(self) -> RemoteStatus:
        """Obtain information about the remote execution status."""
        channel = self.newchannel()
        self._send(Message.STATUS, channel.id)
        statusdict = channel.receive()
        # the other side didn't actually instantiate a channel
        # so we just delete the internal id/channel mapping
        self._channelfactory._local_close(channel.id)
        return RemoteStatus(statusdict)

    def remote_exec(
        self,
        source: str | types.FunctionType | Callable[..., object] | types.ModuleType,
        **kwargs: object,
    ) -> Channel:
        """Return channel object and connect it to a remote
        execution thread where the given ``source`` executes.

        * ``source`` is a string: execute source string remotely
          with a ``channel`` put into the global namespace.
        * ``source`` is a pure function: serialize source and
          call function with ``**kwargs``, adding a
          ``channel`` object to the keyword arguments.
        * ``source`` is a pure module: execute source of module
          with a ``channel`` in its global namespace.

        In all cases the binding ``__name__='__channelexec__'``
        will be available in the global namespace of the remotely
        executing code.
        """
        call_name = None
        file_name = None
        if isinstance(source, types.ModuleType):
            file_name = inspect.getsourcefile(source)
            linecache.updatecache(file_name)  # type: ignore[arg-type]
            source = inspect.getsource(source)
        elif isinstance(source, types.FunctionType):
            call_name = source.__name__
            file_name = inspect.getsourcefile(source)
            source = _source_of_function(source)
        else:
            source = textwrap.dedent(str(source))

        if not call_name and kwargs:
            raise TypeError("can't pass kwargs to non-function remote_exec")

        channel = self.newchannel()
        self._send(
            Message.CHANNEL_EXEC,
            channel.id,
            gateway_base.dumps_internal((source, file_name, call_name, kwargs)),
        )
        return channel

    def remote_init_threads(self, num: int | None = None) -> None:
        """DEPRECATED.  Is currently a NO-OPERATION already."""
        print("WARNING: remote_init_threads() is a no-operation in execnet-1.2")


class RInfo:
    def __init__(self, kwargs) -> None:
        self.__dict__.update(kwargs)

    def __repr__(self) -> str:
        info = ", ".join(f"{k}={v}" for k, v in sorted(self.__dict__.items()))
        return "<RInfo %r>" % info

    if TYPE_CHECKING:

        def __getattr__(self, name: str) -> Any: ...


RemoteStatus = RInfo


def rinfo_source(channel) -> None:
    import os
    import sys

    channel.send(
        dict(
            executable=sys.executable,
            version_info=sys.version_info[:5],
            platform=sys.platform,
            cwd=os.getcwd(),
            pid=os.getpid(),
        )
    )


def _find_non_builtin_globals(source: str, codeobj: types.CodeType) -> list[str]:
    import ast
    import builtins

    vars = dict.fromkeys(codeobj.co_varnames)
    return [
        node.id
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Name)
        and node.id not in vars
        and node.id not in builtins.__dict__
    ]


def _source_of_function(function: types.FunctionType | Callable[..., object]) -> str:
    if function.__name__ == "<lambda>":
        raise ValueError("can't evaluate lambda functions'")
    # XXX: we dont check before remote instantiation
    #      if arguments are used properly
    try:
        sig = inspect.getfullargspec(function)
    except AttributeError:
        args = inspect.getargspec(function)[0]
    else:
        args = sig.args
    if not args or args[0] != "channel":
        raise ValueError("expected first function argument to be `channel`")

    closure = function.__closure__
    codeobj = function.__code__

    if closure is not None:
        raise ValueError("functions with closures can't be passed")

    try:
        source = inspect.getsource(function)
    except OSError as e:
        raise ValueError("can't find source file for %s" % function) from e

    source = textwrap.dedent(source)  # just for inner functions

    used_globals = _find_non_builtin_globals(source, codeobj)
    if used_globals:
        raise ValueError("the use of non-builtin globals isn't supported", used_globals)

    leading_ws = "\n" * (codeobj.co_firstlineno - 1)
    return leading_ws + source
