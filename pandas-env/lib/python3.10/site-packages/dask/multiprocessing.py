from __future__ import annotations

import copyreg
import multiprocessing
import multiprocessing.pool
import os
import pickle
import sys
import traceback
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from warnings import warn

import cloudpickle

from dask import config
from dask.local import MultiprocessingPoolExecutor, get_async, reraise
from dask.optimization import cull, fuse
from dask.system import CPU_COUNT
from dask.typing import Key
from dask.utils import ensure_dict


def _reduce_method_descriptor(m):
    return getattr, (m.__objclass__, m.__name__)


# type(set.union) is used as a proxy to <class 'method_descriptor'>
copyreg.pickle(type(set.union), _reduce_method_descriptor)

_dumps = partial(cloudpickle.dumps, protocol=pickle.HIGHEST_PROTOCOL)
_loads = cloudpickle.loads


def _process_get_id():
    return multiprocessing.current_process().ident


# -- Remote Exception Handling --
# By default, tracebacks can't be serialized using pickle. However, the
# `tblib` library can enable support for this. Since we don't mandate
# that tblib is installed, we do the following:
#
# - If tblib is installed, use it to serialize the traceback and reraise
#   in the scheduler process
# - Otherwise, use a ``RemoteException`` class to contain a serialized
#   version of the formatted traceback, which will then print in the
#   scheduler process.
#
# To enable testing of the ``RemoteException`` class even when tblib is
# installed, we don't wrap the class in the try block below
class RemoteException(Exception):
    """Remote Exception

    Contains the exception and traceback from a remotely run task
    """

    def __init__(self, exception, traceback):
        self.exception = exception
        self.traceback = traceback

    def __str__(self):
        return str(self.exception) + "\n\nTraceback\n---------\n" + self.traceback

    def __dir__(self):
        return sorted(set(dir(type(self)) + list(self.__dict__) + dir(self.exception)))

    def __getattr__(self, key):
        try:
            return object.__getattribute__(self, key)
        except AttributeError:
            return getattr(self.exception, key)


exceptions: dict[type[Exception], type[Exception]] = {}


def remote_exception(exc: Exception, tb) -> Exception:
    """Metaclass that wraps exception type in RemoteException"""
    if type(exc) in exceptions:
        typ = exceptions[type(exc)]
        return typ(exc, tb)
    else:
        try:
            typ = type(
                exc.__class__.__name__,
                (RemoteException, type(exc)),
                {"exception_type": type(exc)},
            )
            exceptions[type(exc)] = typ
            return typ(exc, tb)
        except TypeError:
            return exc


try:
    import tblib.pickling_support

    tblib.pickling_support.install()

    def _pack_traceback(tb):
        return tb

except ImportError:

    def _pack_traceback(tb):
        return "".join(traceback.format_tb(tb))

    def reraise(exc, tb=None):
        exc = remote_exception(exc, tb)
        raise exc


def pack_exception(e, dumps):
    exc_type, exc_value, exc_traceback = sys.exc_info()
    tb = _pack_traceback(exc_traceback)
    try:
        result = dumps((e, tb))
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = _pack_traceback(exc_traceback)
        result = dumps((e, tb))
    return result


_CONTEXT_UNSUPPORTED = """\
The 'multiprocessing.context' configuration option will be ignored on Python 2
and on Windows, because they each only support a single context.
"""


def get_context():
    """Return the current multiprocessing context."""
    # fork context does fork()-without-exec(), which can lead to deadlocks,
    # so default to "spawn".
    context_name = config.get("multiprocessing.context", "spawn")
    if sys.platform == "win32":
        if context_name != "spawn":
            # Only spawn is supported on Win32, can't change it:
            warn(_CONTEXT_UNSUPPORTED, UserWarning)
        return multiprocessing
    else:
        return multiprocessing.get_context(context_name)


def get(
    dsk: Mapping,
    keys: Sequence[Key] | Key,
    num_workers=None,
    func_loads=None,
    func_dumps=None,
    optimize_graph=True,
    pool=None,
    initializer=None,
    chunksize=None,
    **kwargs,
):
    """Multiprocessed get function appropriate for Bags

    Parameters
    ----------
    dsk : dict
        dask graph
    keys : object or list
        Desired results from graph
    num_workers : int
        Number of worker processes (defaults to number of cores)
    func_dumps : function
        Function to use for function serialization (defaults to cloudpickle.dumps)
    func_loads : function
        Function to use for function deserialization (defaults to cloudpickle.loads)
    optimize_graph : bool
        If True [default], `fuse` is applied to the graph before computation.
    pool : Executor or Pool
        Some sort of `Executor` or `Pool` to use
    initializer: function
        Ignored if ``pool`` has been set.
        Function to initialize a worker process before running any tasks in it.
    chunksize: int, optional
        Size of chunks to use when dispatching work.
        Defaults to 6 as some batching is helpful.
        If -1, will be computed to evenly divide ready work across workers.
    """
    chunksize = chunksize or config.get("chunksize", 6)
    pool = pool or config.get("pool", None)
    initializer = initializer or config.get("multiprocessing.initializer", None)
    num_workers = num_workers or config.get("num_workers", None) or CPU_COUNT
    if pool is None:
        # In order to get consistent hashing in subprocesses, we need to set a
        # consistent seed for the Python hash algorithm. Unfortunately, there
        # is no way to specify environment variables only for the Pool
        # processes, so we have to rely on environment variables being
        # inherited.
        if os.environ.get("PYTHONHASHSEED") in (None, "0"):
            # This number is arbitrary; it was chosen to commemorate
            # https://github.com/dask/dask/issues/6640.
            os.environ["PYTHONHASHSEED"] = "6640"
        context = get_context()
        initializer = partial(initialize_worker_process, user_initializer=initializer)
        pool = ProcessPoolExecutor(
            num_workers, mp_context=context, initializer=initializer
        )
        cleanup = True
    else:
        if initializer is not None:
            warn(
                "The ``initializer`` argument is ignored when ``pool`` is provided. "
                "The user should configure ``pool`` with the needed ``initializer`` "
                "on creation."
            )
        if isinstance(pool, multiprocessing.pool.Pool):
            pool = MultiprocessingPoolExecutor(pool)
        cleanup = False

    # Optimize Dask
    dsk = ensure_dict(dsk)
    dsk2, dependencies = cull(dsk, keys)
    if optimize_graph:
        dsk3, dependencies = fuse(dsk2, keys, dependencies)
    else:
        dsk3 = dsk2

    # We specify marshalling functions in order to catch serialization
    # errors and report them to the user.
    loads = func_loads or config.get("func_loads", None) or _loads
    dumps = func_dumps or config.get("func_dumps", None) or _dumps

    # Note former versions used a multiprocessing Manager to share
    # a Queue between parent and workers, but this is fragile on Windows
    # (issue #1652).
    try:
        # Run
        result = get_async(
            pool.submit,
            pool._max_workers,
            dsk3,
            keys,
            get_id=_process_get_id,
            dumps=dumps,
            loads=loads,
            pack_exception=pack_exception,
            raise_exception=reraise,
            chunksize=chunksize,
            **kwargs,
        )
    finally:
        if cleanup:
            pool.shutdown()
    return result


def default_initializer():
    # If Numpy is already imported, presumably its random state was
    # inherited from the parent => re-seed it.
    np = sys.modules.get("numpy")
    if np is not None:
        np.random.seed()


def initialize_worker_process(user_initializer=None):
    """
    Initialize a worker process before running any tasks in it.
    """
    default_initializer()

    if user_initializer is not None:
        user_initializer()
