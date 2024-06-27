"""Base execnet gateway code send to the other side for bootstrapping.

:copyright: 2004-2015
:authors:
    - Holger Krekel
    - Armin Rigo
    - Benjamin Peterson
    - Ronny Pfannschmidt
    - many others
"""

from __future__ import annotations

import abc
import os
import struct
import sys
import traceback
import weakref
from _thread import interrupt_main
from io import BytesIO
from typing import Any
from typing import Callable
from typing import Iterator
from typing import Literal
from typing import MutableSet
from typing import Protocol
from typing import cast
from typing import overload


class WriteIO(Protocol):
    def write(self, data: bytes, /) -> None: ...


class ReadIO(Protocol):
    def read(self, numbytes: int, /) -> bytes: ...


class IO(Protocol):
    execmodel: ExecModel

    def read(self, numbytes: int, /) -> bytes: ...

    def write(self, data: bytes, /) -> None: ...

    def close_read(self) -> None: ...

    def close_write(self) -> None: ...

    def wait(self) -> int | None: ...

    def kill(self) -> None: ...


class Event(Protocol):
    """Protocol for types which look like threading.Event."""

    def is_set(self) -> bool: ...

    def set(self) -> None: ...

    def clear(self) -> None: ...

    def wait(self, timeout: float | None = None) -> bool: ...


class ExecModel(metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def backend(self) -> str:
        raise NotImplementedError()

    def __repr__(self) -> str:
        return "<ExecModel %r>" % self.backend

    @property
    @abc.abstractmethod
    def queue(self):
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def subprocess(self):
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def socket(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def start(self, func, args=()) -> None:
        raise NotImplementedError()

    @abc.abstractmethod
    def get_ident(self) -> int:
        raise NotImplementedError()

    @abc.abstractmethod
    def sleep(self, delay: float) -> None:
        raise NotImplementedError()

    @abc.abstractmethod
    def fdopen(self, fd, mode, bufsize=1, closefd=True):
        raise NotImplementedError()

    @abc.abstractmethod
    def Lock(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def RLock(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def Event(self) -> Event:
        raise NotImplementedError()


class ThreadExecModel(ExecModel):
    backend = "thread"

    @property
    def queue(self):
        import queue

        return queue

    @property
    def subprocess(self):
        import subprocess

        return subprocess

    @property
    def socket(self):
        import socket

        return socket

    def get_ident(self) -> int:
        import _thread

        return _thread.get_ident()

    def sleep(self, delay: float) -> None:
        import time

        time.sleep(delay)

    def start(self, func, args=()) -> None:
        import _thread

        _thread.start_new_thread(func, args)

    def fdopen(self, fd, mode, bufsize=1, closefd=True):
        import os

        return os.fdopen(fd, mode, bufsize, encoding="utf-8", closefd=closefd)

    def Lock(self):
        import threading

        return threading.RLock()

    def RLock(self):
        import threading

        return threading.RLock()

    def Event(self):
        import threading

        return threading.Event()


class MainThreadOnlyExecModel(ThreadExecModel):
    backend = "main_thread_only"


class EventletExecModel(ExecModel):
    backend = "eventlet"

    @property
    def queue(self):
        import eventlet

        return eventlet.queue

    @property
    def subprocess(self):
        import eventlet.green.subprocess

        return eventlet.green.subprocess

    @property
    def socket(self):
        import eventlet.green.socket

        return eventlet.green.socket

    def get_ident(self) -> int:
        import eventlet.green.thread

        return eventlet.green.thread.get_ident()  # type: ignore[no-any-return]

    def sleep(self, delay: float) -> None:
        import eventlet

        eventlet.sleep(delay)

    def start(self, func, args=()) -> None:
        import eventlet

        eventlet.spawn_n(func, *args)

    def fdopen(self, fd, mode, bufsize=1, closefd=True):
        import eventlet.green.os

        return eventlet.green.os.fdopen(fd, mode, bufsize, closefd=closefd)

    def Lock(self):
        import eventlet.green.threading

        return eventlet.green.threading.RLock()

    def RLock(self):
        import eventlet.green.threading

        return eventlet.green.threading.RLock()

    def Event(self):
        import eventlet.green.threading

        return eventlet.green.threading.Event()


class GeventExecModel(ExecModel):
    backend = "gevent"

    @property
    def queue(self):
        import gevent.queue

        return gevent.queue

    @property
    def subprocess(self):
        import gevent.subprocess

        return gevent.subprocess

    @property
    def socket(self):
        import gevent

        return gevent.socket

    def get_ident(self) -> int:
        import gevent.thread

        return gevent.thread.get_ident()  # type: ignore[no-any-return]

    def sleep(self, delay: float) -> None:
        import gevent

        gevent.sleep(delay)

    def start(self, func, args=()) -> None:
        import gevent

        gevent.spawn(func, *args)

    def fdopen(self, fd, mode, bufsize=1, closefd=True):
        # XXX
        import gevent.fileobject

        return gevent.fileobject.FileObjectThread(fd, mode, bufsize, closefd=closefd)

    def Lock(self):
        import gevent.lock

        return gevent.lock.RLock()

    def RLock(self):
        import gevent.lock

        return gevent.lock.RLock()

    def Event(self):
        import gevent.event

        return gevent.event.Event()


def get_execmodel(backend: str | ExecModel) -> ExecModel:
    if isinstance(backend, ExecModel):
        return backend
    if backend == "thread":
        return ThreadExecModel()
    elif backend == "main_thread_only":
        return MainThreadOnlyExecModel()
    elif backend == "eventlet":
        return EventletExecModel()
    elif backend == "gevent":
        return GeventExecModel()
    else:
        raise ValueError(f"unknown execmodel {backend!r}")


class Reply:
    """Provide access to the result of a function execution that got dispatched
    through WorkerPool.spawn()."""

    def __init__(self, task, threadmodel: ExecModel) -> None:
        self.task = task
        self._result_ready = threadmodel.Event()
        self.running = True

    def get(self, timeout: float | None = None):
        """get the result object from an asynchronous function execution.
        if the function execution raised an exception,
        then calling get() will reraise that exception
        including its traceback.
        """
        self.waitfinish(timeout)
        try:
            return self._result
        except AttributeError:
            raise self._exc from None

    def waitfinish(self, timeout: float | None = None) -> None:
        if not self._result_ready.wait(timeout):
            raise OSError(f"timeout waiting for {self.task!r}")

    def run(self) -> None:
        func, args, kwargs = self.task
        try:
            try:
                self._result = func(*args, **kwargs)
            except BaseException as exc:
                self._exc = exc
        finally:
            self._result_ready.set()
            self.running = False


class WorkerPool:
    """A WorkerPool allows to spawn function executions
    to threads, returning a reply object on which you
    can ask for the result (and get exceptions reraised).

    This implementation allows the main thread to integrate
    itself into performing function execution through
    calling integrate_as_primary_thread() which will return
    when the pool received a trigger_shutdown().

    By default allows unlimited number of spawns.
    """

    _primary_thread_task: Reply | None

    def __init__(self, execmodel: ExecModel, hasprimary: bool = False) -> None:
        self.execmodel = execmodel
        self._running_lock = self.execmodel.Lock()
        self._running: MutableSet[Reply] = set()
        self._shuttingdown = False
        self._waitall_events: list[Event] = []
        if hasprimary:
            if self.execmodel.backend not in ("thread", "main_thread_only"):
                raise ValueError("hasprimary=True requires thread model")
            self._primary_thread_task_ready: Event | None = self.execmodel.Event()
        else:
            self._primary_thread_task_ready = None

    def integrate_as_primary_thread(self) -> None:
        """Integrate the thread with which we are called as a primary
        thread for executing functions triggered with spawn()."""
        assert self.execmodel.backend in ("thread", "main_thread_only"), self.execmodel
        primary_thread_task_ready = self._primary_thread_task_ready
        assert primary_thread_task_ready is not None
        # interacts with code at REF1
        while 1:
            primary_thread_task_ready.wait()
            reply = self._primary_thread_task
            if reply is None:  # trigger_shutdown() woke us up
                break
            self._perform_spawn(reply)
            # we are concurrent with trigger_shutdown and spawn
            with self._running_lock:
                if self._shuttingdown:
                    break
                # Only clear if _try_send_to_primary_thread has not
                # yet set the next self._primary_thread_task reply
                # after waiting for this one to complete.
                if reply is self._primary_thread_task:
                    primary_thread_task_ready.clear()

    def trigger_shutdown(self) -> None:
        with self._running_lock:
            self._shuttingdown = True
            if self._primary_thread_task_ready is not None:
                self._primary_thread_task = None
                self._primary_thread_task_ready.set()

    def active_count(self) -> int:
        return len(self._running)

    def _perform_spawn(self, reply: Reply) -> None:
        reply.run()
        with self._running_lock:
            self._running.remove(reply)
            if not self._running:
                while self._waitall_events:
                    waitall_event = self._waitall_events.pop()
                    waitall_event.set()

    def _try_send_to_primary_thread(self, reply: Reply) -> bool:
        # REF1 in 'thread' model we give priority to running in main thread
        # note that we should be called with _running_lock hold
        primary_thread_task_ready = self._primary_thread_task_ready
        if primary_thread_task_ready is not None:
            if not primary_thread_task_ready.is_set():
                self._primary_thread_task = reply
                # wake up primary thread
                primary_thread_task_ready.set()
                return True
            elif (
                self.execmodel.backend == "main_thread_only"
                and self._primary_thread_task is not None
            ):
                self._primary_thread_task.waitfinish()
                self._primary_thread_task = reply
                # wake up primary thread (it's okay if this is already set
                # because we waited for the previous task to finish above
                # and integrate_as_primary_thread will not clear it when
                # it enters self._running_lock if it detects that a new
                # task is available)
                primary_thread_task_ready.set()
                return True
        return False

    def spawn(self, func, *args, **kwargs) -> Reply:
        """Asynchronously dispatch func(*args, **kwargs) and return a Reply."""
        reply = Reply((func, args, kwargs), self.execmodel)
        with self._running_lock:
            if self._shuttingdown:
                raise ValueError("pool is shutting down")
            self._running.add(reply)
            if not self._try_send_to_primary_thread(reply):
                self.execmodel.start(self._perform_spawn, (reply,))
        return reply

    def terminate(self, timeout: float | None = None) -> bool:
        """Trigger shutdown and wait for completion of all executions."""
        self.trigger_shutdown()
        return self.waitall(timeout=timeout)

    def waitall(self, timeout: float | None = None) -> bool:
        """Wait until all active spawns have finished executing."""
        with self._running_lock:
            if not self._running:
                return True
            # if a Reply still runs, we let run_and_release
            # signal us -- note that we are still holding the
            # _running_lock to avoid race conditions
            my_waitall_event = self.execmodel.Event()
            self._waitall_events.append(my_waitall_event)
        return my_waitall_event.wait(timeout=timeout)


sysex = (KeyboardInterrupt, SystemExit)


DEBUG = os.environ.get("EXECNET_DEBUG")
pid = os.getpid()
if DEBUG == "2":

    def trace(*msg: object) -> None:
        try:
            line = " ".join(map(str, msg))
            sys.stderr.write(f"[{pid}] {line}\n")
            sys.stderr.flush()
        except Exception:
            pass  # nothing we can do, likely interpreter-shutdown

elif DEBUG:
    import os
    import tempfile

    fn = os.path.join(tempfile.gettempdir(), "execnet-debug-%d" % pid)
    # sys.stderr.write("execnet-debug at %r" % (fn,))
    debugfile = open(fn, "w")

    def trace(*msg: object) -> None:
        try:
            line = " ".join(map(str, msg))
            debugfile.write(line + "\n")
            debugfile.flush()
        except Exception as exc:
            try:
                sys.stderr.write(f"[{pid}] exception during tracing: {exc!r}\n")
            except Exception:
                pass  # nothing we can do, likely interpreter-shutdown

else:
    notrace = trace = lambda *msg: None


class Popen2IO:
    error = (IOError, OSError, EOFError)

    def __init__(self, outfile, infile, execmodel: ExecModel) -> None:
        # we need raw byte streams
        self.outfile, self.infile = outfile, infile
        if sys.platform == "win32":
            import msvcrt

            try:
                msvcrt.setmode(infile.fileno(), os.O_BINARY)
                msvcrt.setmode(outfile.fileno(), os.O_BINARY)
            except (AttributeError, OSError):
                pass
        self._read = getattr(infile, "buffer", infile).read
        self._write = getattr(outfile, "buffer", outfile).write
        self.execmodel = execmodel

    def read(self, numbytes: int) -> bytes:
        """Read exactly 'numbytes' bytes from the pipe."""
        # a file in non-blocking mode may return less bytes, so we loop
        buf = b""
        while numbytes > len(buf):
            data = self._read(numbytes - len(buf))
            if not data:
                raise EOFError("expected %d bytes, got %d" % (numbytes, len(buf)))
            buf += data
        return buf

    def write(self, data: bytes) -> None:
        """Write out all data bytes."""
        assert isinstance(data, bytes)
        self._write(data)
        self.outfile.flush()

    def close_read(self) -> None:
        self.infile.close()

    def close_write(self) -> None:
        self.outfile.close()


class Message:
    """Encapsulates Messages and their wire protocol."""

    # message code -> name, handler
    _types: dict[int, tuple[str, Callable[[Message, BaseGateway], None]]] = {}

    def __init__(self, msgcode: int, channelid: int = 0, data: bytes = b"") -> None:
        self.msgcode = msgcode
        self.channelid = channelid
        self.data = data

    @staticmethod
    def from_io(io: ReadIO) -> Message:
        try:
            header = io.read(9)  # type 1, channel 4, payload 4
            if not header:
                raise EOFError("empty read")
        except EOFError as e:
            raise EOFError("couldn't load message header, " + e.args[0]) from None
        msgtype, channel, payload = struct.unpack("!bii", header)
        return Message(msgtype, channel, io.read(payload))

    def to_io(self, io: WriteIO) -> None:
        header = struct.pack("!bii", self.msgcode, self.channelid, len(self.data))
        io.write(header + self.data)

    def received(self, gateway: BaseGateway) -> None:
        handler = self._types[self.msgcode][1]
        handler(self, gateway)

    def __repr__(self) -> str:
        name = self._types[self.msgcode][0]
        return f"<Message {name} channel={self.channelid} lendata={len(self.data)}>"

    def _status(message: Message, gateway: BaseGateway) -> None:
        # we use the channelid to send back information
        # but don't instantiate a channel object
        d = {
            "numchannels": len(gateway._channelfactory._channels),
            # TODO(typing): Attribute `_execpool` is only on WorkerGateway.
            "numexecuting": gateway._execpool.active_count(),  # type: ignore[attr-defined]
            "execmodel": gateway.execmodel.backend,
        }
        gateway._send(Message.CHANNEL_DATA, message.channelid, dumps_internal(d))
        gateway._send(Message.CHANNEL_CLOSE, message.channelid)

    STATUS = 0
    _types[STATUS] = ("STATUS", _status)

    def _reconfigure(message: Message, gateway: BaseGateway) -> None:
        data = loads_internal(message.data, gateway)
        assert isinstance(data, tuple)
        strconfig: tuple[bool, bool] = data
        if message.channelid == 0:
            gateway._strconfig = strconfig
        else:
            gateway._channelfactory.new(message.channelid)._strconfig = strconfig

    RECONFIGURE = 1
    _types[RECONFIGURE] = ("RECONFIGURE", _reconfigure)

    def _gateway_terminate(message: Message, gateway: BaseGateway) -> None:
        raise GatewayReceivedTerminate(gateway)

    GATEWAY_TERMINATE = 2
    _types[GATEWAY_TERMINATE] = ("GATEWAY_TERMINATE", _gateway_terminate)

    def _channel_exec(message: Message, gateway: BaseGateway) -> None:
        channel = gateway._channelfactory.new(message.channelid)
        gateway._local_schedulexec(channel=channel, sourcetask=message.data)

    CHANNEL_EXEC = 3
    _types[CHANNEL_EXEC] = ("CHANNEL_EXEC", _channel_exec)

    def _channel_data(message: Message, gateway: BaseGateway) -> None:
        gateway._channelfactory._local_receive(message.channelid, message.data)

    CHANNEL_DATA = 4
    _types[CHANNEL_DATA] = ("CHANNEL_DATA", _channel_data)

    def _channel_close(message: Message, gateway: BaseGateway) -> None:
        gateway._channelfactory._local_close(message.channelid)

    CHANNEL_CLOSE = 5
    _types[CHANNEL_CLOSE] = ("CHANNEL_CLOSE", _channel_close)

    def _channel_close_error(message: Message, gateway: BaseGateway) -> None:
        error_message = loads_internal(message.data)
        assert isinstance(error_message, str)
        remote_error = RemoteError(error_message)
        gateway._channelfactory._local_close(message.channelid, remote_error)

    CHANNEL_CLOSE_ERROR = 6
    _types[CHANNEL_CLOSE_ERROR] = ("CHANNEL_CLOSE_ERROR", _channel_close_error)

    def _channel_last_message(message: Message, gateway: BaseGateway) -> None:
        gateway._channelfactory._local_close(message.channelid, sendonly=True)

    CHANNEL_LAST_MESSAGE = 7
    _types[CHANNEL_LAST_MESSAGE] = ("CHANNEL_LAST_MESSAGE", _channel_last_message)


class GatewayReceivedTerminate(Exception):
    """Receiverthread got termination message."""


def geterrortext(
    exc: BaseException,
    format_exception=traceback.format_exception,
    sysex: tuple[type[BaseException], ...] = sysex,
) -> str:
    try:
        # In py310, can change this to:
        # l = format_exception(exc)
        l = format_exception(type(exc), exc, exc.__traceback__)
        errortext = "".join(l)
    except sysex:
        raise
    except BaseException:
        errortext = f"{type(exc).__name__}: {exc}"
    return errortext


class RemoteError(Exception):
    """Exception containing a stringified error from the other side."""

    def __init__(self, formatted: str) -> None:
        super().__init__()
        self.formatted = formatted

    def __str__(self) -> str:
        return self.formatted

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}: {self.formatted}"

    def warn(self) -> None:
        if self.formatted != INTERRUPT_TEXT:
            # XXX do this better
            sys.stderr.write(f"[{os.getpid()}] Warning: unhandled {self!r}\n")


class TimeoutError(IOError):
    """Exception indicating that a timeout was reached."""


NO_ENDMARKER_WANTED = object()


class Channel:
    """Communication channel between two Python Interpreter execution points."""

    RemoteError = RemoteError
    TimeoutError = TimeoutError
    _INTERNALWAKEUP = 1000
    _executing = False

    def __init__(self, gateway: BaseGateway, id: int) -> None:
        """:private:"""
        assert isinstance(id, int)
        assert not isinstance(gateway, type)
        self.gateway = gateway
        # XXX: defaults copied from Unserializer
        self._strconfig = getattr(gateway, "_strconfig", (True, False))
        self.id = id
        self._items = self.gateway.execmodel.queue.Queue()
        self._closed = False
        self._receiveclosed = self.gateway.execmodel.Event()
        self._remoteerrors: list[RemoteError] = []

    def _trace(self, *msg: object) -> None:
        self.gateway._trace(self.id, *msg)

    def setcallback(
        self,
        callback: Callable[[Any], Any],
        endmarker: object = NO_ENDMARKER_WANTED,
    ) -> None:
        """Set a callback function for receiving items.

        All already-queued items will immediately trigger the callback.
        Afterwards the callback will execute in the receiver thread
        for each received data item and calls to ``receive()`` will
        raise an error.
        If an endmarker is specified the callback will eventually
        be called with the endmarker when the channel closes.
        """
        _callbacks = self.gateway._channelfactory._callbacks
        with self.gateway._receivelock:
            if self._items is None:
                raise OSError(f"{self!r} has callback already registered")
            items = self._items
            self._items = None
            while 1:
                try:
                    olditem = items.get(block=False)
                except self.gateway.execmodel.queue.Empty:
                    if not (self._closed or self._receiveclosed.is_set()):
                        _callbacks[self.id] = (callback, endmarker, self._strconfig)
                    break
                else:
                    if olditem is ENDMARKER:
                        items.put(olditem)  # for other receivers
                        if endmarker is not NO_ENDMARKER_WANTED:
                            callback(endmarker)
                        break
                    else:
                        callback(olditem)

    def __repr__(self) -> str:
        flag = self.isclosed() and "closed" or "open"
        return "<Channel id=%d %s>" % (self.id, flag)

    def __del__(self) -> None:
        if self.gateway is None:  # can be None in tests
            return  # type: ignore[unreachable]

        self._trace("channel.__del__")
        # no multithreading issues here, because we have the last ref to 'self'
        if self._closed:
            # state transition "closed" --> "deleted"
            for error in self._remoteerrors:
                error.warn()
        elif self._receiveclosed.is_set():
            # state transition "sendonly" --> "deleted"
            # the remote channel is already in "deleted" state, nothing to do
            pass
        else:
            # state transition "opened" --> "deleted"
            # check if we are in the middle of interpreter shutdown
            # in which case the process will go away and we probably
            # don't need to try to send a closing or last message
            # (and often it won't work anymore to send things out)
            if Message is not None:
                if self._items is None:  # has_callback
                    msgcode = Message.CHANNEL_LAST_MESSAGE
                else:
                    msgcode = Message.CHANNEL_CLOSE
                try:
                    self.gateway._send(msgcode, self.id)
                except (OSError, ValueError):  # ignore problems with sending
                    pass

    def _getremoteerror(self):
        try:
            return self._remoteerrors.pop(0)
        except IndexError:
            try:
                return self.gateway._error
            except AttributeError:
                pass
            return None

    #
    # public API for channel objects
    #
    def isclosed(self) -> bool:
        """Return True if the channel is closed.

        A closed channel may still hold items.
        """
        return self._closed

    @overload
    def makefile(self, mode: Literal["r"], proxyclose: bool = ...) -> ChannelFileRead:
        pass

    @overload
    def makefile(
        self,
        mode: Literal["w"] = ...,
        proxyclose: bool = ...,
    ) -> ChannelFileWrite:
        pass

    def makefile(
        self,
        mode: Literal["r", "w"] = "w",
        proxyclose: bool = False,
    ) -> ChannelFileWrite | ChannelFileRead:
        """Return a file-like object.

        mode can be 'w' or 'r' for writeable/readable files.
        If proxyclose is true, file.close() will also close the channel.
        """
        if mode == "w":
            return ChannelFileWrite(channel=self, proxyclose=proxyclose)
        elif mode == "r":
            return ChannelFileRead(channel=self, proxyclose=proxyclose)
        raise ValueError(f"mode {mode!r} not available")

    def close(self, error=None) -> None:
        """Close down this channel with an optional error message.

        Note that closing of a channel tied to remote_exec happens
        automatically at the end of execution and cannot
        be done explicitly.
        """
        if self._executing:
            raise OSError("cannot explicitly close channel within remote_exec")
        if self._closed:
            self.gateway._trace(self, "ignoring redundant call to close()")
        if not self._closed:
            # state transition "opened/sendonly" --> "closed"
            # threads warning: the channel might be closed under our feet,
            # but it's never damaging to send too many CHANNEL_CLOSE messages
            # however, if the other side triggered a close already, we
            # do not send back a closed message.
            if not self._receiveclosed.is_set():
                put = self.gateway._send
                if error is not None:
                    put(Message.CHANNEL_CLOSE_ERROR, self.id, dumps_internal(error))
                else:
                    put(Message.CHANNEL_CLOSE, self.id)
                self._trace("sent channel close message")
            if isinstance(error, RemoteError):
                self._remoteerrors.append(error)
            self._closed = True  # --> "closed"
            self._receiveclosed.set()
            queue = self._items
            if queue is not None:
                queue.put(ENDMARKER)
            self.gateway._channelfactory._no_longer_opened(self.id)

    def waitclose(self, timeout: float | None = None) -> None:
        """Wait until this channel is closed (or the remote side
        otherwise signalled that no more data was being sent).

        The channel may still hold receiveable items, but not receive
        any more after waitclose() has returned.

        Exceptions from executing code on the other side are reraised as local
        channel.RemoteErrors.

        EOFError is raised if the reading-connection was prematurely closed,
        which often indicates a dying process.

        self.TimeoutError is raised after the specified number of seconds
        (default is None, i.e. wait indefinitely).
        """
        # wait for non-"opened" state
        self._receiveclosed.wait(timeout=timeout)
        if not self._receiveclosed.is_set():
            raise self.TimeoutError("Timeout after %r seconds" % timeout)
        error = self._getremoteerror()
        if error:
            raise error

    def send(self, item: object) -> None:
        """Sends the given item to the other side of the channel,
        possibly blocking if the sender queue is full.

        The item must be a simple Python type and will be
        copied to the other side by value.

        OSError is raised if the write pipe was prematurely closed.
        """
        if self.isclosed():
            raise OSError(f"cannot send to {self!r}")
        self.gateway._send(Message.CHANNEL_DATA, self.id, dumps_internal(item))

    def receive(self, timeout: float | None = None) -> Any:
        """Receive a data item that was sent from the other side.

        timeout: None [default] blocked waiting. A positive number
        indicates the number of seconds after which a channel.TimeoutError
        exception will be raised if no item was received.

        Note that exceptions from the remotely executing code will be
        reraised as channel.RemoteError exceptions containing
        a textual representation of the remote traceback.
        """
        itemqueue = self._items
        if itemqueue is None:
            raise OSError("cannot receive(), channel has receiver callback")
        try:
            x = itemqueue.get(timeout=timeout)
        except self.gateway.execmodel.queue.Empty:
            raise self.TimeoutError("no item after %r seconds" % timeout) from None
        if x is ENDMARKER:
            itemqueue.put(x)  # for other receivers
            raise self._getremoteerror() or EOFError()
        else:
            return x

    def __iter__(self) -> Iterator[Any]:
        return self

    def next(self) -> Any:
        try:
            return self.receive()
        except EOFError:
            raise StopIteration from None

    __next__ = next

    def reconfigure(
        self, py2str_as_py3str: bool = True, py3str_as_py2str: bool = False
    ) -> None:
        """Set the string coercion for this channel.

        The default is to try to convert py2 str as py3 str,
        but not to try and convert py3 str to py2 str
        """
        self._strconfig = (py2str_as_py3str, py3str_as_py2str)
        data = dumps_internal(self._strconfig)
        self.gateway._send(Message.RECONFIGURE, self.id, data=data)


ENDMARKER = object()
INTERRUPT_TEXT = "keyboard-interrupted"
MAIN_THREAD_ONLY_DEADLOCK_TEXT = (
    "concurrent remote_exec would cause deadlock for main_thread_only execmodel"
)


class ChannelFactory:
    def __init__(self, gateway: BaseGateway, startcount: int = 1) -> None:
        self._channels: weakref.WeakValueDictionary[int, Channel] = (
            weakref.WeakValueDictionary()
        )
        # Channel ID => (callback, end marker, strconfig)
        self._callbacks: dict[
            int, tuple[Callable[[Any], Any], object, tuple[bool, bool]]
        ] = {}
        self._writelock = gateway.execmodel.Lock()
        self.gateway = gateway
        self.count = startcount
        self.finished = False
        self._list = list  # needed during interp-shutdown

    def new(self, id: int | None = None) -> Channel:
        """Create a new Channel with 'id' (or create new id if None)."""
        with self._writelock:
            if self.finished:
                raise OSError(f"connection already closed: {self.gateway}")
            if id is None:
                id = self.count
                self.count += 2
            try:
                channel = self._channels[id]
            except KeyError:
                channel = self._channels[id] = Channel(self.gateway, id)
            return channel

    def channels(self) -> list[Channel]:
        return self._list(self._channels.values())

    #
    # internal methods, called from the receiver thread
    #
    def _no_longer_opened(self, id: int) -> None:
        try:
            del self._channels[id]
        except KeyError:
            pass
        try:
            callback, endmarker, strconfig = self._callbacks.pop(id)
        except KeyError:
            pass
        else:
            if endmarker is not NO_ENDMARKER_WANTED:
                callback(endmarker)

    def _local_close(self, id: int, remoteerror=None, sendonly: bool = False) -> None:
        channel = self._channels.get(id)
        if channel is None:
            # channel already in "deleted" state
            if remoteerror:
                remoteerror.warn()
            self._no_longer_opened(id)
        else:
            # state transition to "closed" state
            if remoteerror:
                channel._remoteerrors.append(remoteerror)
            queue = channel._items
            if queue is not None:
                queue.put(ENDMARKER)
            self._no_longer_opened(id)
            if not sendonly:  # otherwise #--> "sendonly"
                channel._closed = True  # --> "closed"
            channel._receiveclosed.set()

    def _local_receive(self, id: int, data) -> None:
        # executes in receiver thread
        channel = self._channels.get(id)
        try:
            callback, endmarker, strconfig = self._callbacks[id]
        except KeyError:
            queue = channel._items if channel is not None else None
            if queue is None:
                pass  # drop data
            else:
                item = loads_internal(data, channel)
                queue.put(item)
        else:
            try:
                data = loads_internal(data, channel, strconfig)
                callback(data)  # even if channel may be already closed
            except Exception as exc:
                self.gateway._trace("exception during callback: %s" % exc)
                errortext = self.gateway._geterrortext(exc)
                self.gateway._send(
                    Message.CHANNEL_CLOSE_ERROR, id, dumps_internal(errortext)
                )
                self._local_close(id, errortext)

    def _finished_receiving(self) -> None:
        with self._writelock:
            self.finished = True
        for id in self._list(self._channels):
            self._local_close(id, sendonly=True)
        for id in self._list(self._callbacks):
            self._no_longer_opened(id)


class ChannelFile:
    def __init__(self, channel: Channel, proxyclose: bool = True) -> None:
        self.channel = channel
        self._proxyclose = proxyclose

    def isatty(self) -> bool:
        return False

    def close(self) -> None:
        if self._proxyclose:
            self.channel.close()

    def __repr__(self) -> str:
        state = self.channel.isclosed() and "closed" or "open"
        return "<ChannelFile %d %s>" % (self.channel.id, state)


class ChannelFileWrite(ChannelFile):
    def write(self, out: bytes) -> None:
        self.channel.send(out)

    def flush(self) -> None:
        pass


class ChannelFileRead(ChannelFile):
    def __init__(self, channel: Channel, proxyclose: bool = True) -> None:
        super().__init__(channel, proxyclose)
        self._buffer: str | None = None

    def read(self, n: int) -> str:
        try:
            if self._buffer is None:
                self._buffer = cast(str, self.channel.receive())
            while len(self._buffer) < n:
                self._buffer += cast(str, self.channel.receive())
        except EOFError:
            self.close()
        if self._buffer is None:
            ret = ""
        else:
            ret = self._buffer[:n]
            self._buffer = self._buffer[n:]
        return ret

    def readline(self) -> str:
        if self._buffer is not None:
            i = self._buffer.find("\n")
            if i != -1:
                return self.read(i + 1)
            line = self.read(len(self._buffer) + 1)
        else:
            line = self.read(1)
        while line and line[-1] != "\n":
            c = self.read(1)
            if not c:
                break
            line += c
        return line


class BaseGateway:
    _sysex = sysex
    id = "<worker>"

    def __init__(self, io: IO, id, _startcount: int = 2) -> None:
        self.execmodel = io.execmodel
        self._io = io
        self.id = id
        self._strconfig = (Unserializer.py2str_as_py3str, Unserializer.py3str_as_py2str)
        self._channelfactory = ChannelFactory(self, _startcount)
        self._receivelock = self.execmodel.RLock()
        # globals may be NONE at process-termination
        self.__trace = trace
        self._geterrortext = geterrortext
        self._receivepool = WorkerPool(self.execmodel)

    def _trace(self, *msg: object) -> None:
        self.__trace(self.id, *msg)

    def _initreceive(self) -> None:
        self._receivepool.spawn(self._thread_receiver)

    def _thread_receiver(self) -> None:
        def log(*msg: object) -> None:
            self._trace("[receiver-thread]", *msg)

        log("RECEIVERTHREAD: starting to run")
        io = self._io
        try:
            while 1:
                msg = Message.from_io(io)
                log("received", msg)
                with self._receivelock:
                    msg.received(self)
                    del msg
        except (KeyboardInterrupt, GatewayReceivedTerminate):
            pass
        except EOFError as exc:
            log("EOF without prior gateway termination message")
            self._error = exc
        except Exception as exc:
            log(self._geterrortext(exc))
        log("finishing receiving thread")
        # wake up and terminate any execution waiting to receive
        self._channelfactory._finished_receiving()
        log("terminating execution")
        self._terminate_execution()
        log("closing read")
        self._io.close_read()
        log("closing write")
        self._io.close_write()
        log("terminating our receive pseudo pool")
        self._receivepool.trigger_shutdown()

    def _terminate_execution(self) -> None:
        pass

    def _send(self, msgcode: int, channelid: int = 0, data: bytes = b"") -> None:
        message = Message(msgcode, channelid, data)
        try:
            message.to_io(self._io)
            self._trace("sent", message)
        except (OSError, ValueError) as e:
            self._trace("failed to send", message, e)
            # ValueError might be because the IO is already closed
            raise OSError("cannot send (already closed?)") from e

    def _local_schedulexec(self, channel: Channel, sourcetask: bytes) -> None:
        channel.close("execution disallowed")

    # _____________________________________________________________________
    #
    # High Level Interface
    # _____________________________________________________________________
    #
    def newchannel(self) -> Channel:
        """Return a new independent channel."""
        return self._channelfactory.new()

    def join(self, timeout: float | None = None) -> None:
        """Wait for receiverthread to terminate."""
        self._trace("waiting for receiver thread to finish")
        self._receivepool.waitall(timeout)


class WorkerGateway(BaseGateway):
    def _local_schedulexec(self, channel: Channel, sourcetask: bytes) -> None:
        if self._execpool.execmodel.backend == "main_thread_only":
            assert self._executetask_complete is not None
            # It's necessary to wait for a short time in order to ensure
            # that we do not report a false-positive deadlock error, since
            # channel close does not elicit a response that would provide
            # a guarantee to remote_exec callers that the previous task
            # has released the main thread. If the timeout expires then it
            # should be practically impossible to report a false-positive.
            if not self._executetask_complete.wait(timeout=1):
                channel.close(MAIN_THREAD_ONLY_DEADLOCK_TEXT)
                return
            # It's only safe to clear here because the above wait proves
            # that there is not a previous task about to set it again.
            self._executetask_complete.clear()

        sourcetask_ = loads_internal(sourcetask)
        self._execpool.spawn(self.executetask, (channel, sourcetask_))

    def _terminate_execution(self) -> None:
        # called from receiverthread
        self._trace("shutting down execution pool")
        self._execpool.trigger_shutdown()
        if not self._execpool.waitall(5.0):
            self._trace("execution ongoing after 5 secs," " trying interrupt_main")
            # We try hard to terminate execution based on the assumption
            # that there is only one gateway object running per-process.
            if sys.platform != "win32":
                self._trace("sending ourselves a SIGINT")
                os.kill(os.getpid(), 2)  # send ourselves a SIGINT
            elif interrupt_main is not None:
                self._trace("calling interrupt_main()")
                interrupt_main()
            if not self._execpool.waitall(10.0):
                self._trace(
                    "execution did not finish in another 10 secs, " "calling os._exit()"
                )
                os._exit(1)

    def serve(self) -> None:
        def trace(msg: str) -> None:
            self._trace("[serve] " + msg)

        hasprimary = self.execmodel.backend in ("thread", "main_thread_only")
        self._execpool = WorkerPool(self.execmodel, hasprimary=hasprimary)
        self._executetask_complete = None
        if self.execmodel.backend == "main_thread_only":
            self._executetask_complete = self.execmodel.Event()
            # Initialize state to indicate that there is no previous task
            # executing so that we don't need a separate flag to track this.
            self._executetask_complete.set()
        trace("spawning receiver thread")
        self._initreceive()
        try:
            if hasprimary:
                # this will return when we are in shutdown
                trace("integrating as primary thread")
                self._execpool.integrate_as_primary_thread()
            trace("joining receiver thread")
            self.join()
        except KeyboardInterrupt:
            # in the worker we can't really do anything sensible
            trace("swallowing keyboardinterrupt, serve finished")

    def executetask(
        self,
        item: tuple[Channel, tuple[str, str | None, str | None, dict[str, object]]],
    ) -> None:
        try:
            channel, (source, file_name, call_name, kwargs) = item
            loc: dict[str, Any] = {"channel": channel, "__name__": "__channelexec__"}
            self._trace(f"execution starts[{channel.id}]: {repr(source)[:50]}")
            channel._executing = True
            try:
                co = compile(source + "\n", file_name or "<remote exec>", "exec")
                exec(co, loc)
                if call_name:
                    self._trace("calling %s(**%60r)" % (call_name, kwargs))
                    function = loc[call_name]
                    function(channel, **kwargs)
            finally:
                channel._executing = False
                self._trace("execution finished")
        except KeyboardInterrupt:
            channel.close(INTERRUPT_TEXT)
            raise
        except BaseException as exc:
            if not isinstance(exc, EOFError):
                if not channel.gateway._channelfactory.finished:
                    self._trace(f"got exception: {exc!r}")
                    errortext = self._geterrortext(exc)
                    channel.close(errortext)
                    return
            self._trace("ignoring EOFError because receiving finished")
        channel.close()
        if self._executetask_complete is not None:
            # Indicate that this task has finished executing, meaning
            # that there is no possibility of it triggering a deadlock
            # for the next spawn call.
            self._executetask_complete.set()


#
# Cross-Python pickling code, tested from test_serializer.py
#


class DataFormatError(Exception):
    pass


class DumpError(DataFormatError):
    """Error while serializing an object."""


class LoadError(DataFormatError):
    """Error while unserializing an object."""


def bchr(n: int) -> bytes:
    return bytes([n])


DUMPFORMAT_VERSION = bchr(2)

FOUR_BYTE_INT_MAX = 2147483647

FLOAT_FORMAT = "!d"
FLOAT_FORMAT_SIZE = struct.calcsize(FLOAT_FORMAT)
COMPLEX_FORMAT = "!dd"
COMPLEX_FORMAT_SIZE = struct.calcsize(COMPLEX_FORMAT)


class _Stop(Exception):
    pass


class opcode:
    """Container for name -> num mappings."""

    BUILDTUPLE = b"@"
    BYTES = b"A"
    CHANNEL = b"B"
    FALSE = b"C"
    FLOAT = b"D"
    FROZENSET = b"E"
    INT = b"F"
    LONG = b"G"
    LONGINT = b"H"
    LONGLONG = b"I"
    NEWDICT = b"J"
    NEWLIST = b"K"
    NONE = b"L"
    PY2STRING = b"M"
    PY3STRING = b"N"
    SET = b"O"
    SETITEM = b"P"
    STOP = b"Q"
    TRUE = b"R"
    UNICODE = b"S"
    COMPLEX = b"T"


class Unserializer:
    num2func: dict[bytes, Callable[[Unserializer], None]] = {}
    py2str_as_py3str = True  # True
    py3str_as_py2str = False  # false means py2 will get unicode

    def __init__(
        self,
        stream: ReadIO,
        channel_or_gateway: Channel | BaseGateway | None = None,
        strconfig: tuple[bool, bool] | None = None,
    ) -> None:
        if isinstance(channel_or_gateway, Channel):
            gw: BaseGateway | None = channel_or_gateway.gateway
        else:
            gw = channel_or_gateway
        if channel_or_gateway is not None:
            strconfig = channel_or_gateway._strconfig
        if strconfig:
            self.py2str_as_py3str, self.py3str_as_py2str = strconfig
        self.stream = stream
        if gw is None:
            self.channelfactory = None
        else:
            self.channelfactory = gw._channelfactory

    def load(self, versioned: bool = False) -> Any:
        if versioned:
            ver = self.stream.read(1)
            if ver != DUMPFORMAT_VERSION:
                raise LoadError("wrong dumpformat version %r" % ver)
        self.stack: list[object] = []
        try:
            while True:
                opcode = self.stream.read(1)
                if not opcode:
                    raise EOFError
                try:
                    loader = self.num2func[opcode]
                except KeyError:
                    raise LoadError(
                        f"unknown opcode {opcode!r} - wire protocol corruption?"
                    ) from None
                loader(self)
        except _Stop:
            if len(self.stack) != 1:
                raise LoadError("internal unserialization error") from None
            return self.stack.pop(0)
        else:
            raise LoadError("didn't get STOP")

    def load_none(self) -> None:
        self.stack.append(None)

    num2func[opcode.NONE] = load_none

    def load_true(self) -> None:
        self.stack.append(True)

    num2func[opcode.TRUE] = load_true

    def load_false(self) -> None:
        self.stack.append(False)

    num2func[opcode.FALSE] = load_false

    def load_int(self) -> None:
        i = self._read_int4()
        self.stack.append(i)

    num2func[opcode.INT] = load_int

    def load_longint(self) -> None:
        s = self._read_byte_string()
        self.stack.append(int(s))

    num2func[opcode.LONGINT] = load_longint

    load_long = load_int
    num2func[opcode.LONG] = load_long
    load_longlong = load_longint
    num2func[opcode.LONGLONG] = load_longlong

    def load_float(self) -> None:
        binary = self.stream.read(FLOAT_FORMAT_SIZE)
        self.stack.append(struct.unpack(FLOAT_FORMAT, binary)[0])

    num2func[opcode.FLOAT] = load_float

    def load_complex(self) -> None:
        binary = self.stream.read(COMPLEX_FORMAT_SIZE)
        self.stack.append(complex(*struct.unpack(COMPLEX_FORMAT, binary)))

    num2func[opcode.COMPLEX] = load_complex

    def _read_int4(self) -> int:
        value: int = struct.unpack("!i", self.stream.read(4))[0]
        return value

    def _read_byte_string(self) -> bytes:
        length = self._read_int4()
        as_bytes = self.stream.read(length)
        return as_bytes

    def load_py3string(self) -> None:
        as_bytes = self._read_byte_string()
        if self.py3str_as_py2str:
            # XXX Should we try to decode into latin-1?
            self.stack.append(as_bytes)
        else:
            self.stack.append(as_bytes.decode("utf-8"))

    num2func[opcode.PY3STRING] = load_py3string

    def load_py2string(self) -> None:
        as_bytes = self._read_byte_string()
        if self.py2str_as_py3str:
            s: bytes | str = as_bytes.decode("latin-1")
        else:
            s = as_bytes
        self.stack.append(s)

    num2func[opcode.PY2STRING] = load_py2string

    def load_bytes(self) -> None:
        s = self._read_byte_string()
        self.stack.append(s)

    num2func[opcode.BYTES] = load_bytes

    def load_unicode(self) -> None:
        self.stack.append(self._read_byte_string().decode("utf-8"))

    num2func[opcode.UNICODE] = load_unicode

    def load_newlist(self) -> None:
        length = self._read_int4()
        self.stack.append([None] * length)

    num2func[opcode.NEWLIST] = load_newlist

    def load_setitem(self) -> None:
        if len(self.stack) < 3:
            raise LoadError("not enough items for setitem")
        value = self.stack.pop()
        key = self.stack.pop()
        self.stack[-1][key] = value  # type: ignore[index]

    num2func[opcode.SETITEM] = load_setitem

    def load_newdict(self) -> None:
        self.stack.append({})

    num2func[opcode.NEWDICT] = load_newdict

    def _load_collection(self, type_: type) -> None:
        length = self._read_int4()
        if length:
            res = type_(self.stack[-length:])
            del self.stack[-length:]
            self.stack.append(res)
        else:
            self.stack.append(type_())

    def load_buildtuple(self) -> None:
        self._load_collection(tuple)

    num2func[opcode.BUILDTUPLE] = load_buildtuple

    def load_set(self) -> None:
        self._load_collection(set)

    num2func[opcode.SET] = load_set

    def load_frozenset(self) -> None:
        self._load_collection(frozenset)

    num2func[opcode.FROZENSET] = load_frozenset

    def load_stop(self) -> None:
        raise _Stop

    num2func[opcode.STOP] = load_stop

    def load_channel(self) -> None:
        id = self._read_int4()
        assert self.channelfactory is not None
        newchannel = self.channelfactory.new(id)
        self.stack.append(newchannel)

    num2func[opcode.CHANNEL] = load_channel


def dumps(obj: object) -> bytes:
    """Serialize the given obj to a bytestring.

    The obj and all contained objects must be of a builtin
    Python type (so nested dicts, sets, etc. are all OK but
    not user-level instances).
    """
    return _Serializer().save(obj, versioned=True)  # type: ignore[return-value]


def dump(byteio, obj: object) -> None:
    """write a serialized bytestring of the given obj to the given stream."""
    _Serializer(write=byteio.write).save(obj, versioned=True)


def loads(
    bytestring: bytes, py2str_as_py3str: bool = False, py3str_as_py2str: bool = False
) -> Any:
    """Deserialize the given bytestring to an object.

    py2str_as_py3str: If true then string (str) objects previously
                      dumped on Python2 will be loaded as Python3
                      strings which really are text objects.
    py3str_as_py2str: If true then string (str) objects previously
                      dumped on Python3 will be loaded as Python2
                      strings instead of unicode objects.

    If the bytestring was dumped with an incompatible protocol
    version or if the bytestring is corrupted, the
    ``execnet.DataFormatError`` will be raised.
    """
    io = BytesIO(bytestring)
    return load(
        io, py2str_as_py3str=py2str_as_py3str, py3str_as_py2str=py3str_as_py2str
    )


def load(
    io: ReadIO, py2str_as_py3str: bool = False, py3str_as_py2str: bool = False
) -> Any:
    """Derserialize an object form the specified stream.

    Behaviour and parameters are otherwise the same as with ``loads``
    """
    strconfig = (py2str_as_py3str, py3str_as_py2str)
    return Unserializer(io, strconfig=strconfig).load(versioned=True)


def loads_internal(
    bytestring: bytes,
    channelfactory=None,
    strconfig: tuple[bool, bool] | None = None,
) -> Any:
    io = BytesIO(bytestring)
    return Unserializer(io, channelfactory, strconfig).load()


def dumps_internal(obj: object) -> bytes:
    return _Serializer().save(obj)  # type: ignore[return-value]


class _Serializer:
    _dispatch: dict[type, Callable[[_Serializer, object], None]] = {}

    def __init__(self, write: Callable[[bytes], None] | None = None) -> None:
        if write is None:
            self._streamlist: list[bytes] = []
            write = self._streamlist.append
        self._write = write

    def save(self, obj: object, versioned: bool = False) -> bytes | None:
        # calling here is not re-entrant but multiple instances
        # may write to the same stream because of the common platform
        # atomic-write guarantee (concurrent writes each happen atomically)
        if versioned:
            self._write(DUMPFORMAT_VERSION)
        self._save(obj)
        self._write(opcode.STOP)
        try:
            streamlist = self._streamlist
        except AttributeError:
            return None
        return b"".join(streamlist)

    def _save(self, obj: object) -> None:
        tp = type(obj)
        try:
            dispatch = self._dispatch[tp]
        except KeyError:
            methodname = "save_" + tp.__name__
            meth: Callable[[_Serializer, object], None] | None = getattr(
                self.__class__, methodname, None
            )
            if meth is None:
                raise DumpError(f"can't serialize {tp}") from None
            dispatch = self._dispatch[tp] = meth
        dispatch(self, obj)

    def save_NoneType(self, non: None) -> None:
        self._write(opcode.NONE)

    def save_bool(self, boolean: bool) -> None:
        if boolean:
            self._write(opcode.TRUE)
        else:
            self._write(opcode.FALSE)

    def save_bytes(self, bytes_: bytes) -> None:
        self._write(opcode.BYTES)
        self._write_byte_sequence(bytes_)

    def save_str(self, s: str) -> None:
        self._write(opcode.PY3STRING)
        self._write_unicode_string(s)

    def _write_unicode_string(self, s: str) -> None:
        try:
            as_bytes = s.encode("utf-8")
        except UnicodeEncodeError as e:
            raise DumpError("strings must be utf-8 encodable") from e
        self._write_byte_sequence(as_bytes)

    def _write_byte_sequence(self, bytes_: bytes) -> None:
        self._write_int4(len(bytes_), "string is too long")
        self._write(bytes_)

    def _save_integral(self, i: int, short_op: bytes, long_op: bytes) -> None:
        if i <= FOUR_BYTE_INT_MAX:
            self._write(short_op)
            self._write_int4(i)
        else:
            self._write(long_op)
            self._write_byte_sequence(str(i).rstrip("L").encode("ascii"))

    def save_int(self, i: int) -> None:
        self._save_integral(i, opcode.INT, opcode.LONGINT)

    def save_long(self, l: int) -> None:
        self._save_integral(l, opcode.LONG, opcode.LONGLONG)

    def save_float(self, flt: float) -> None:
        self._write(opcode.FLOAT)
        self._write(struct.pack(FLOAT_FORMAT, flt))

    def save_complex(self, cpx: complex) -> None:
        self._write(opcode.COMPLEX)
        self._write(struct.pack(COMPLEX_FORMAT, cpx.real, cpx.imag))

    def _write_int4(
        self, i: int, error: str = "int must be less than %i" % (FOUR_BYTE_INT_MAX,)
    ) -> None:
        if i > FOUR_BYTE_INT_MAX:
            raise DumpError(error)
        self._write(struct.pack("!i", i))

    def save_list(self, L: list[object]) -> None:
        self._write(opcode.NEWLIST)
        self._write_int4(len(L), "list is too long")
        for i, item in enumerate(L):
            self._write_setitem(i, item)

    def _write_setitem(self, key: object, value: object) -> None:
        self._save(key)
        self._save(value)
        self._write(opcode.SETITEM)

    def save_dict(self, d: dict[object, object]) -> None:
        self._write(opcode.NEWDICT)
        for key, value in d.items():
            self._write_setitem(key, value)

    def save_tuple(self, tup: tuple[object, ...]) -> None:
        for item in tup:
            self._save(item)
        self._write(opcode.BUILDTUPLE)
        self._write_int4(len(tup), "tuple is too long")

    def _write_set(self, s: set[object] | frozenset[object], op: bytes) -> None:
        for item in s:
            self._save(item)
        self._write(op)
        self._write_int4(len(s), "set is too long")

    def save_set(self, s: set[object]) -> None:
        self._write_set(s, opcode.SET)

    def save_frozenset(self, s: frozenset[object]) -> None:
        self._write_set(s, opcode.FROZENSET)

    def save_Channel(self, channel: Channel) -> None:
        self._write(opcode.CHANNEL)
        self._write_int4(channel.id)


def init_popen_io(execmodel: ExecModel) -> Popen2IO:
    if not hasattr(os, "dup"):  # jython
        io = Popen2IO(sys.stdout, sys.stdin, execmodel)
        import tempfile

        sys.stdin = tempfile.TemporaryFile("r")
        sys.stdout = tempfile.TemporaryFile("w")
    else:
        try:
            devnull = os.devnull
        except AttributeError:
            if os.name == "nt":
                devnull = "NUL"
            else:
                devnull = "/dev/null"
        # stdin
        stdin = execmodel.fdopen(os.dup(0), "r", 1)
        fd = os.open(devnull, os.O_RDONLY)
        os.dup2(fd, 0)
        os.close(fd)

        # stdout
        stdout = execmodel.fdopen(os.dup(1), "w", 1)
        fd = os.open(devnull, os.O_WRONLY)
        os.dup2(fd, 1)

        # stderr for win32
        if os.name == "nt":
            sys.stderr = execmodel.fdopen(os.dup(2), "w", 1)
            os.dup2(fd, 2)
        os.close(fd)
        io = Popen2IO(stdout, stdin, execmodel)
        # Use closefd=False since 0 and 1 are shared with
        # sys.__stdin__ and sys.__stdout__.
        sys.stdin = execmodel.fdopen(0, "r", 1, closefd=False)
        sys.stdout = execmodel.fdopen(1, "w", 1, closefd=False)
    return io


def serve(io: IO, id) -> None:
    trace(f"creating workergateway on {io!r}")
    WorkerGateway(io=io, id=id, _startcount=2).serve()
