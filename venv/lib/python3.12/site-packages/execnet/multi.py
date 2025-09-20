"""
Managing Gateway Groups and interactions with multiple channels.

(c) 2008-2014, Holger Krekel and others
"""

from __future__ import annotations

import atexit
import types
from functools import partial
from threading import Lock
from typing import TYPE_CHECKING
from typing import Any
from typing import Callable
from typing import Iterable
from typing import Iterator
from typing import Literal
from typing import Sequence
from typing import overload

from . import gateway_bootstrap
from . import gateway_io
from .gateway_base import Channel
from .gateway_base import ExecModel
from .gateway_base import WorkerPool
from .gateway_base import get_execmodel
from .gateway_base import trace
from .xspec import XSpec

if TYPE_CHECKING:
    from .gateway import Gateway


NO_ENDMARKER_WANTED = object()


class Group:
    """Gateway Group."""

    defaultspec = "popen"

    def __init__(
        self, xspecs: Iterable[XSpec | str | None] = (), execmodel: str = "thread"
    ) -> None:
        """Initialize a group and make gateways as specified.

        execmodel can be one of the supported execution models.
        """
        self._gateways: list[Gateway] = []
        self._autoidcounter = 0
        self._autoidlock = Lock()
        self._gateways_to_join: list[Gateway] = []
        # we use the same execmodel for all of the Gateway objects
        # we spawn on our side.  Probably we should not allow different
        # execmodels between different groups but not clear.
        # Note that "other side" execmodels may differ and is typically
        # specified by the spec passed to makegateway.
        self.set_execmodel(execmodel)
        for xspec in xspecs:
            self.makegateway(xspec)
        atexit.register(self._cleanup_atexit)

    @property
    def execmodel(self) -> ExecModel:
        return self._execmodel

    @property
    def remote_execmodel(self) -> ExecModel:
        return self._remote_execmodel

    def set_execmodel(
        self, execmodel: str, remote_execmodel: str | None = None
    ) -> None:
        """Set the execution model for local and remote site.

        execmodel can be one of the supported execution models.
        It determines the execution model for any newly created gateway.
        If remote_execmodel is not specified it takes on the value of execmodel.

        NOTE: Execution models can only be set before any gateway is created.
        """
        if self._gateways:
            raise ValueError(
                "can not set execution models if " "gateways have been created already"
            )
        if remote_execmodel is None:
            remote_execmodel = execmodel
        self._execmodel = get_execmodel(execmodel)
        self._remote_execmodel = get_execmodel(remote_execmodel)

    def __repr__(self) -> str:
        idgateways = [gw.id for gw in self]
        return "<Group %r>" % idgateways

    def __getitem__(self, key: int | str | Gateway) -> Gateway:
        if isinstance(key, int):
            return self._gateways[key]
        for gw in self._gateways:
            if gw == key or gw.id == key:
                return gw
        raise KeyError(key)

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
            return True
        except KeyError:
            return False

    def __len__(self) -> int:
        return len(self._gateways)

    def __iter__(self) -> Iterator[Gateway]:
        return iter(list(self._gateways))

    def makegateway(self, spec: XSpec | str | None = None) -> Gateway:
        """Create and configure a gateway to a Python interpreter.

        The ``spec`` string encodes the target gateway type
        and configuration information. The general format is::

            key1=value1//key2=value2//...

        If you leave out the ``=value`` part a True value is assumed.
        Valid types: ``popen``, ``ssh=hostname``, ``socket=host:port``.
        Valid configuration::

            id=<string>     specifies the gateway id
            python=<path>   specifies which python interpreter to execute
            execmodel=model 'thread', 'main_thread_only', 'eventlet', 'gevent' execution model
            chdir=<path>    specifies to which directory to change
            nice=<path>     specifies process priority of new process
            env:NAME=value  specifies a remote environment variable setting.

        If no spec is given, self.defaultspec is used.
        """
        if not spec:
            spec = self.defaultspec
        if not isinstance(spec, XSpec):
            spec = XSpec(spec)
        self.allocate_id(spec)
        if spec.execmodel is None:
            spec.execmodel = self.remote_execmodel.backend
        if spec.via:
            assert not spec.socket
            master = self[spec.via]
            proxy_channel = master.remote_exec(gateway_io)
            proxy_channel.send(vars(spec))
            proxy_io_master = gateway_io.ProxyIO(proxy_channel, self.execmodel)
            gw = gateway_bootstrap.bootstrap(proxy_io_master, spec)
        elif spec.popen or spec.ssh or spec.vagrant_ssh:
            io = gateway_io.create_io(spec, execmodel=self.execmodel)
            gw = gateway_bootstrap.bootstrap(io, spec)
        elif spec.socket:
            from . import gateway_socket

            sio = gateway_socket.create_io(spec, self, execmodel=self.execmodel)
            gw = gateway_bootstrap.bootstrap(sio, spec)
        else:
            raise ValueError(f"no gateway type found for {spec._spec!r}")
        gw.spec = spec
        self._register(gw)
        if spec.chdir or spec.nice or spec.env:
            channel = gw.remote_exec(
                """
                import os
                path, nice, env = channel.receive()
                if path:
                    if not os.path.exists(path):
                        os.mkdir(path)
                    os.chdir(path)
                if nice and hasattr(os, 'nice'):
                    os.nice(nice)
                if env:
                    for name, value in env.items():
                        os.environ[name] = value
            """
            )
            nice = spec.nice and int(spec.nice) or 0
            channel.send((spec.chdir, nice, spec.env))
            channel.waitclose()
        return gw

    def allocate_id(self, spec: XSpec) -> None:
        """(re-entrant) allocate id for the given xspec object."""
        if spec.id is None:
            with self._autoidlock:
                id = "gw" + str(self._autoidcounter)
                self._autoidcounter += 1
                if id in self:
                    raise ValueError(f"already have gateway with id {id!r}")
                spec.id = id

    def _register(self, gateway: Gateway) -> None:
        assert not hasattr(gateway, "_group")
        assert gateway.id
        assert gateway.id not in self
        self._gateways.append(gateway)
        gateway._group = self

    def _unregister(self, gateway: Gateway) -> None:
        self._gateways.remove(gateway)
        self._gateways_to_join.append(gateway)

    def _cleanup_atexit(self) -> None:
        trace(f"=== atexit cleanup {self!r} ===")
        self.terminate(timeout=1.0)

    def terminate(self, timeout: float | None = None) -> None:
        """Trigger exit of member gateways and wait for termination
        of member gateways and associated subprocesses.

        After waiting timeout seconds try to to kill local sub processes of
        popen- and ssh-gateways.

        Timeout defaults to None meaning open-ended waiting and no kill
        attempts.
        """
        while self:
            vias: set[str] = set()
            for gw in self:
                if gw.spec.via:
                    vias.add(gw.spec.via)
            for gw in self:
                if gw.id not in vias:
                    gw.exit()

            def join_wait(gw: Gateway) -> None:
                gw.join()
                gw._io.wait()

            def kill(gw: Gateway) -> None:
                trace("Gateways did not come down after timeout: %r" % gw)
                gw._io.kill()

            safe_terminate(
                self.execmodel,
                timeout,
                [
                    (partial(join_wait, gw), partial(kill, gw))
                    for gw in self._gateways_to_join
                ],
            )
            self._gateways_to_join[:] = []

    def remote_exec(
        self,
        source: str | types.FunctionType | Callable[..., object] | types.ModuleType,
        **kwargs,
    ) -> MultiChannel:
        """remote_exec source on all member gateways and return
        a MultiChannel connecting to all sub processes."""
        channels = []
        for gw in self:
            channels.append(gw.remote_exec(source, **kwargs))
        return MultiChannel(channels)


class MultiChannel:
    def __init__(self, channels: Sequence[Channel]) -> None:
        self._channels = channels

    def __len__(self) -> int:
        return len(self._channels)

    def __iter__(self) -> Iterator[Channel]:
        return iter(self._channels)

    def __getitem__(self, key: int) -> Channel:
        return self._channels[key]

    def __contains__(self, chan: Channel) -> bool:
        return chan in self._channels

    def send_each(self, item: object) -> None:
        for ch in self._channels:
            ch.send(item)

    @overload
    def receive_each(self, withchannel: Literal[False] = ...) -> list[Any]:
        pass

    @overload
    def receive_each(self, withchannel: Literal[True]) -> list[tuple[Channel, Any]]:
        pass

    def receive_each(
        self, withchannel: bool = False
    ) -> list[tuple[Channel, Any]] | list[Any]:
        assert not hasattr(self, "_queue")
        l: list[object] = []
        for ch in self._channels:
            obj = ch.receive()
            if withchannel:
                l.append((ch, obj))
            else:
                l.append(obj)
        return l

    def make_receive_queue(self, endmarker: object = NO_ENDMARKER_WANTED):
        try:
            return self._queue  # type: ignore[has-type]
        except AttributeError:
            self._queue = None
            for ch in self._channels:
                if self._queue is None:
                    self._queue = ch.gateway.execmodel.queue.Queue()

                def putreceived(obj, channel: Channel = ch) -> None:
                    self._queue.put((channel, obj))  # type: ignore[union-attr]

                if endmarker is NO_ENDMARKER_WANTED:
                    ch.setcallback(putreceived)
                else:
                    ch.setcallback(putreceived, endmarker=endmarker)
            return self._queue

    def waitclose(self) -> None:
        first = None
        for ch in self._channels:
            try:
                ch.waitclose()
            except ch.RemoteError as exc:
                if first is None:
                    first = exc
        if first:
            raise first


def safe_terminate(
    execmodel: ExecModel, timeout: float | None, list_of_paired_functions
) -> None:
    workerpool = WorkerPool(execmodel)

    def termkill(termfunc, killfunc) -> None:
        termreply = workerpool.spawn(termfunc)
        try:
            termreply.get(timeout=timeout)
        except OSError:
            killfunc()

    replylist = []
    for termfunc, killfunc in list_of_paired_functions:
        reply = workerpool.spawn(termkill, termfunc, killfunc)
        replylist.append(reply)
    for reply in replylist:
        reply.get()
    workerpool.waitall(timeout=timeout)


default_group = Group()
makegateway = default_group.makegateway
set_execmodel = default_group.set_execmodel
