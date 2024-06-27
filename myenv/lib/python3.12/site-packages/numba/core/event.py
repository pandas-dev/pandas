"""
The ``numba.core.event`` module provides a simple event system for applications
to register callbacks to listen to specific compiler events.

The following events are built in:

- ``"numba:compile"`` is broadcast when a dispatcher is compiling. Events of
  this kind have ``data`` defined to be a ``dict`` with the following
  key-values:

  - ``"dispatcher"``: the dispatcher object that is compiling.
  - ``"args"``: the argument types.
  - ``"return_type"``: the return type.

- ``"numba:compiler_lock"`` is broadcast when the internal compiler-lock is
  acquired. This is mostly used internally to measure time spent with the lock
  acquired.

- ``"numba:llvm_lock"`` is broadcast when the internal LLVM-lock is acquired.
  This is used internally to measure time spent with the lock acquired.

- ``"numba:run_pass"`` is broadcast when a compiler pass is running.

    - ``"name"``: pass name.
    - ``"qualname"``: qualified name of the function being compiled.
    - ``"module"``: module name of the function being compiled.
    - ``"flags"``: compilation flags.
    - ``"args"``: argument types.
    - ``"return_type"`` return type.

Applications can register callbacks that are listening for specific events using
``register(kind: str, listener: Listener)``, where ``listener`` is an instance
of ``Listener`` that defines custom actions on occurrence of the specific event.
"""

import os
import json
import atexit
import abc
import enum
import time
import threading
from timeit import default_timer as timer
from contextlib import contextmanager, ExitStack
from collections import defaultdict

from numba.core import config


class EventStatus(enum.Enum):
    """Status of an event.
    """
    START = enum.auto()
    END = enum.auto()


# Builtin event kinds.
_builtin_kinds = frozenset([
    "numba:compiler_lock",
    "numba:compile",
    "numba:llvm_lock",
    "numba:run_pass",
])


def _guard_kind(kind):
    """Guard to ensure that an event kind is valid.

    All event kinds with a "numba:" prefix must be defined in the pre-defined
    ``numba.core.event._builtin_kinds``.
    Custom event kinds are allowed by not using the above prefix.

    Parameters
    ----------
    kind : str

    Return
    ------
    res : str
    """
    if kind.startswith("numba:") and kind not in _builtin_kinds:
        msg = (f"{kind} is not a valid event kind, "
               "it starts with the reserved prefix 'numba:'")
        raise ValueError(msg)
    return kind


class Event:
    """An event.

    Parameters
    ----------
    kind : str
    status : EventStatus
    data : any; optional
        Additional data for the event.
    exc_details : 3-tuple; optional
        Same 3-tuple for ``__exit__``.
    """
    def __init__(self, kind, status, data=None, exc_details=None):
        self._kind = _guard_kind(kind)
        self._status = status
        self._data = data
        self._exc_details = (None
                             if exc_details is None or exc_details[0] is None
                             else exc_details)

    @property
    def kind(self):
        """Event kind

        Returns
        -------
        res : str
        """
        return self._kind

    @property
    def status(self):
        """Event status

        Returns
        -------
        res : EventStatus
        """
        return self._status

    @property
    def data(self):
        """Event data

        Returns
        -------
        res : object
        """
        return self._data

    @property
    def is_start(self):
        """Is it a *START* event?

        Returns
        -------
        res : bool
        """
        return self._status == EventStatus.START

    @property
    def is_end(self):
        """Is it an *END* event?

        Returns
        -------
        res : bool
        """
        return self._status == EventStatus.END

    @property
    def is_failed(self):
        """Is the event carrying an exception?

        This is used for *END* event. This method will never return ``True``
        in a *START* event.

        Returns
        -------
        res : bool
        """
        return self._exc_details is None

    def __str__(self):
        data = (f"{type(self.data).__qualname__}"
                if self.data is not None else "None")
        return f"Event({self._kind}, {self._status}, data: {data})"

    __repr__ = __str__


_registered = defaultdict(list)


def register(kind, listener):
    """Register a listener for a given event kind.

    Parameters
    ----------
    kind : str
    listener : Listener
    """
    assert isinstance(listener, Listener)
    kind = _guard_kind(kind)
    _registered[kind].append(listener)


def unregister(kind, listener):
    """Unregister a listener for a given event kind.

    Parameters
    ----------
    kind : str
    listener : Listener
    """
    assert isinstance(listener, Listener)
    kind = _guard_kind(kind)
    lst = _registered[kind]
    lst.remove(listener)


def broadcast(event):
    """Broadcast an event to all registered listeners.

    Parameters
    ----------
    event : Event
    """
    for listener in _registered[event.kind]:
        listener.notify(event)


class Listener(abc.ABC):
    """Base class for all event listeners.
    """
    @abc.abstractmethod
    def on_start(self, event):
        """Called when there is a *START* event.

        Parameters
        ----------
        event : Event
        """
        pass

    @abc.abstractmethod
    def on_end(self, event):
        """Called when there is a *END* event.

        Parameters
        ----------
        event : Event
        """
        pass

    def notify(self, event):
        """Notify this Listener with the given Event.

        Parameters
        ----------
        event : Event
        """
        if event.is_start:
            self.on_start(event)
        elif event.is_end:
            self.on_end(event)
        else:
            raise AssertionError("unreachable")


class TimingListener(Listener):
    """A listener that measures the total time spent between *START* and
    *END* events during the time this listener is active.
    """
    def __init__(self):
        self._depth = 0

    def on_start(self, event):
        if self._depth == 0:
            self._ts = timer()
        self._depth += 1

    def on_end(self, event):
        self._depth -= 1
        if self._depth == 0:
            last = getattr(self, "_duration", 0)
            self._duration = (timer() - self._ts) + last

    @property
    def done(self):
        """Returns a ``bool`` indicating whether a measurement has been made.

        When this returns ``False``, the matching event has never fired.
        If and only if this returns ``True``, ``.duration`` can be read without
        error.
        """
        return hasattr(self, "_duration")

    @property
    def duration(self):
        """Returns the measured duration.

        This may raise ``AttributeError``. Users can use ``.done`` to check
        that a measurement has been made.
        """
        return self._duration


class RecordingListener(Listener):
    """A listener that records all events and stores them in the ``.buffer``
    attribute as a list of 2-tuple ``(float, Event)``, where the first element
    is the time the event occurred as returned by ``time.time()`` and the second
    element is the event.
    """
    def __init__(self):
        self.buffer = []

    def on_start(self, event):
        self.buffer.append((time.time(), event))

    def on_end(self, event):
        self.buffer.append((time.time(), event))


@contextmanager
def install_listener(kind, listener):
    """Install a listener for event "kind" temporarily within the duration of
    the context.

    Returns
    -------
    res : Listener
        The *listener* provided.

    Examples
    --------

    >>> with install_listener("numba:compile", listener):
    >>>     some_code()  # listener will be active here.
    >>> other_code()     # listener will be unregistered by this point.

    """
    register(kind, listener)
    try:
        yield listener
    finally:
        unregister(kind, listener)


@contextmanager
def install_timer(kind, callback):
    """Install a TimingListener temporarily to measure the duration of
    an event.

    If the context completes successfully, the *callback* function is executed.
    The *callback* function is expected to take a float argument for the
    duration in seconds.

    Returns
    -------
    res : TimingListener

    Examples
    --------

    This is equivalent to:

    >>> with install_listener(kind, TimingListener()) as res:
    >>>    ...
    """
    tl = TimingListener()
    with install_listener(kind, tl):
        yield tl

    if tl.done:
        callback(tl.duration)


@contextmanager
def install_recorder(kind):
    """Install a RecordingListener temporarily to record all events.

    Once the context is closed, users can use ``RecordingListener.buffer``
    to access the recorded events.

    Returns
    -------
    res : RecordingListener

    Examples
    --------

    This is equivalent to:

    >>> with install_listener(kind, RecordingListener()) as res:
    >>>    ...
    """
    rl = RecordingListener()
    with install_listener(kind, rl):
        yield rl


def start_event(kind, data=None):
    """Trigger the start of an event of *kind* with *data*.

    Parameters
    ----------
    kind : str
        Event kind.
    data : any; optional
        Extra event data.
    """
    evt = Event(kind=kind, status=EventStatus.START, data=data)
    broadcast(evt)


def end_event(kind, data=None, exc_details=None):
    """Trigger the end of an event of *kind*, *exc_details*.

    Parameters
    ----------
    kind : str
        Event kind.
    data : any; optional
        Extra event data.
    exc_details : 3-tuple; optional
        Same 3-tuple for ``__exit__``. Or, ``None`` if no error.
    """
    evt = Event(
        kind=kind, status=EventStatus.END, data=data, exc_details=exc_details,
    )
    broadcast(evt)


@contextmanager
def trigger_event(kind, data=None):
    """A context manager to trigger the start and end events of *kind* with
    *data*. The start event is triggered when entering the context.
    The end event is triggered when exiting the context.

    Parameters
    ----------
    kind : str
        Event kind.
    data : any; optional
        Extra event data.
    """
    with ExitStack() as scope:
        @scope.push
        def on_exit(*exc_details):
            end_event(kind, data=data, exc_details=exc_details)

        start_event(kind, data=data)
        yield


def _prepare_chrome_trace_data(listener: RecordingListener):
    """Prepare events in `listener` for serializing as chrome trace data.
    """
    # The spec for the trace event format can be found at:
    # https://docs.google.com/document/d/1CvAClvFfyA5R-PhYUmn5OOQtYMH4h6I0nSsKchNAySU/edit   # noqa
    # This code only uses the JSON Array Format for simplicity.
    pid = os.getpid()
    tid = threading.get_native_id()
    evs = []
    for ts, rec in listener.buffer:
        data = rec.data
        cat = str(rec.kind)
        ts_scaled = ts * 1_000_000   # scale to microseconds
        ph = 'B' if rec.is_start else 'E'
        name = data['name']
        args = data
        ev = dict(
            cat=cat, pid=pid, tid=tid, ts=ts_scaled, ph=ph, name=name,
            args=args,
        )
        evs.append(ev)
    return evs


def _setup_chrome_trace_exit_handler():
    """Setup a RecordingListener and an exit handler to write the captured
    events to file.
    """
    listener = RecordingListener()
    register("numba:run_pass", listener)
    filename = config.CHROME_TRACE

    @atexit.register
    def _write_chrome_trace():
        # The following output file is not multi-process safe.
        evs = _prepare_chrome_trace_data(listener)
        with open(filename, "w") as out:
            json.dump(evs, out)


if config.CHROME_TRACE:
    _setup_chrome_trace_exit_handler()
