"""
This module houses a utility class that acts as an idempotent context manager,
ensuring that, within an active context, the given callback is called exactly once.
"""
from collections import defaultdict
import contextvars
from typing import Callable, Hashable, Iterator, MutableMapping, Optional, Union


class _ContextMapping(MutableMapping):
    """
    :class:`~contextvars.Context`-aware mapping.

    Ensures that modifications to the mapping aren't leaked across contexts.
    """

    # We need to copy the internal dictionary when modifying it because contextvars
    # only makes shallow copies of Context objects, so modifications to the same dict
    # would leak across contexts.

    _NO_DEFAULT_VALUE = object()

    def __init__(self, default_value=_NO_DEFAULT_VALUE):
        d: dict = (
            defaultdict(lambda: default_value)
            if default_value is not self._NO_DEFAULT_VALUE
            else dict()
        )
        # yes, we're creating a contextvar inside a closure, but it doesn't matter
        # because objects of this class will only be created at module level
        self._dict_var: contextvars.ContextVar[dict] = contextvars.ContextVar(
            "_ContextMapping<{}>._dict_var".format(id(self))
        )
        self._dict_var.set(d)

    def __setitem__(self, k, v) -> None:
        d = self._dict_var.get().copy()
        d[k] = v
        self._dict_var.set(d)

    def __delitem__(self, k) -> None:
        d = self._dict_var.get().copy()
        del d[k]
        self._dict_var.set(d)

    def __getitem__(self, k):
        return self._dict_var.get()[k]

    def __len__(self) -> int:
        return len(self._dict_var.get())

    def __iter__(self) -> Iterator:
        return iter(self._dict_var.get())


class CallOnceContextManager:
    """
    Ensures the given callback is called exactly once within an active context.

    The first (outermost) context manager, when entered, calls the given callback,
    and any nested context managers (including those that are further down the call
    stack) do nothing while the first one is still active. This ensures that only the
    outermost active context in the current call stack is the one that called the
    callback upon being entered, and thus the callback is guaranteed to (have) be(en)
    called exactly once within any active context.


    It is concurrency-safe (through usage of contextvars).

    Notes
    -----
    To distinguish between context managers with different callbacks ( which
    shouldn't interfere with each other, as they are unrelated), unless an explicit
    ``key`` is given, the :func:`hash` of the callback is used. This is important to
    keep in mind, as, for example, using ``lambda``s won't work as expected:

    >>> with call_once(lambda: print("callback called")):
    ...     with call_once(lambda: print("callback called")):
    ...         pass  # "callback called" printed twice

    So, to ensure it works as intended, either pass an explicit key (preferred) or
    always use identical callback objects (i.e. with the same object ID):

    >>> with call_once(lambda: print("called"), key='foo'):
    ...     with call_once(lambda: print("called"), key='foo'):
    ...         pass  # "callback called" printed once

    The function :func:`call_once` is an alias to this class.

    Examples
    --------
    First, let's define some callback functions:

    >>> cb1, cb2 = map(lambda n: lambda: print(n, "called"), ["callback1", "callback2"])

    A basic illustrative example:

    >>> with call_once(cb1):
    ...    with call_once(cb1):
    ...       pass
    callback1 called

    Context managers for different callbacks don't interfere with each other:
    >>> with call_once(cb1):
    ...     with call_once(cb2):
    ...         with call_once(cb1):
    ...             with call_once(cb2):
    ...                 pass
    callback1 called
    callback2 called

    Example across stack frames:

    >>> def f():
    ...     print("f")
    ...     with call_once(cb1):
    ...         print("f inside with")
    ...
    >>> def g():
    ...     print("g")
    ...     with call_once(cb1):
    ...         print("g inside with")
    ...         f()
    ...

    If we call ``f`` directly, the callback will be called in ``f``:

    >>> f()
    f
    callback called
    f inside with

    But if we call ``g``, the callback will only be called in ``g``, which is where
    the outermost context manager is:

    >>> g()
    g
    callback called
    g inside with
    f
    f inside with

    Concurrency-safe (the state of the class doesn't leak across asyncio tasks or
    different threads):

    >>> import asyncio
    >>> callback = lambda: print("callback called")
    >>> async def task(n: int):
    ...     print("task {} started".format(n))
    ...     with call_once(callback):
    ...         print("task {} sleeping".format(n))
    ...         await asyncio.sleep(1)
    ...         print("task {} exiting with".format(n))
    ...
    >>> async def main():
    ...     await asyncio.gather(task(1), task(2))
    ...
    >>> asyncio.run(main())
    task 1 started
    callback called
    task 1 sleeping
    task 2 started
    callback called
    task 2 sleeping
    task 1 exiting with
    task 2 exiting with
    """

    # mapping of callbacks to the corresponding outermost context manager instances
    _outermost_instances: MutableMapping[
        Union[Callable, Hashable], Optional["CallOnceContextManager"]
    ] = _ContextMapping(default_value=None)

    def __init__(self, callback: Callable[[], None], key: Optional[Hashable] = None):
        """
        Parameters
        ----------
        callback : callable taking no arguments and returning ``None``
            The callback to run.
        key : hashable object, default ``callback``
            Key object to uniquely identify the callback and to distinguish it from
            other callbacks. If not specified, the ``callback`` object itself is used.
        """
        self.key = key if key is not None else callback
        self.callback = callback

    def __enter__(self):
        outermost = self._outermost_instances[self.key]
        if outermost is None:
            self._outermost_instances[self.key] = self
            self.callback()
            return self
        else:
            return outermost

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._outermost_instances[self.key] is self:
            self._outermost_instances[self.key] = None


def call_once(callback: Callable[[], None], key: Optional[Hashable] = None):
    """
    Alias for :class:`CallOnceContextManager`.
    """
    return CallOnceContextManager(callback, key=key)
