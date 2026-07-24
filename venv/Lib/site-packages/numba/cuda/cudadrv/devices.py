"""
Expose each GPU devices directly.

This module implements a API that is like the "CUDA runtime" context manager
for managing CUDA context stack and clean up.  It relies on thread-local globals
to separate the context stack management of each thread. Contexts are also
shareable among threads.  Only the main thread can destroy Contexts.

Note:
- This module must be imported by the main-thread.

"""
import functools
import threading
from contextlib import contextmanager

from .driver import driver, USE_NV_BINDING


class _DeviceList(object):
    def __getattr__(self, attr):
        # First time looking at "lst" attribute.
        if attr == "lst":
            # Device list is not initialized.
            # Query all CUDA devices.
            numdev = driver.get_device_count()
            gpus = [_DeviceContextManager(driver.get_device(devid))
                    for devid in range(numdev)]
            # Define "lst" to avoid re-initialization
            self.lst = gpus
            return gpus

        # Other attributes
        return super(_DeviceList, self).__getattr__(attr)

    def __getitem__(self, devnum):
        '''
        Returns the context manager for device *devnum*.
        '''
        return self.lst[devnum]

    def __str__(self):
        return ', '.join([str(d) for d in self.lst])

    def __iter__(self):
        return iter(self.lst)

    def __len__(self):
        return len(self.lst)

    @property
    def current(self):
        """Returns the active device or None if there's no active device
        """
        with driver.get_active_context() as ac:
            devnum = ac.devnum
            if devnum is not None:
                return self[devnum]


class _DeviceContextManager(object):
    """
    Provides a context manager for executing in the context of the chosen
    device. The normal use of instances of this type is from
    ``numba.cuda.gpus``. For example, to execute on device 2::

       with numba.cuda.gpus[2]:
           d_a = numba.cuda.to_device(a)

    to copy the array *a* onto device 2, referred to by *d_a*.
    """

    def __init__(self, device):
        self._device = device

    def __getattr__(self, item):
        return getattr(self._device, item)

    def __enter__(self):
        _runtime.get_or_create_context(self._device.id)

    def __exit__(self, exc_type, exc_val, exc_tb):
        # this will verify that we are popping the right device context.
        self._device.get_primary_context().pop()

    def __str__(self):
        return "<Managed Device {self.id}>".format(self=self)


class _Runtime(object):
    """Emulate the CUDA runtime context management.

    It owns all Devices and Contexts.
    Keeps at most one Context per Device
    """

    def __init__(self):
        self.gpus = _DeviceList()

        # For caching the attached CUDA Context
        self._tls = threading.local()

        # Remember the main thread
        # Only the main thread can *actually* destroy
        self._mainthread = threading.current_thread()

        # Avoid mutation of runtime state in multithreaded programs
        self._lock = threading.RLock()

    @contextmanager
    def ensure_context(self):
        """Ensure a CUDA context is available inside the context.

        On entrance, queries the CUDA driver for an active CUDA context and
        attaches it in TLS for subsequent calls so they do not need to query
        the CUDA driver again.  On exit, detach the CUDA context from the TLS.

        This will allow us to pickup thirdparty activated CUDA context in
        any top-level Numba CUDA API.
        """
        with driver.get_active_context():
            oldctx = self._get_attached_context()
            newctx = self.get_or_create_context(None)
            self._set_attached_context(newctx)
            try:
                yield
            finally:
                self._set_attached_context(oldctx)

    def get_or_create_context(self, devnum):
        """Returns the primary context and push+create it if needed
        for *devnum*.  If *devnum* is None, use the active CUDA context (must
        be primary) or create a new one with ``devnum=0``.
        """
        if devnum is None:
            attached_ctx = self._get_attached_context()
            if attached_ctx is None:
                return self._get_or_create_context_uncached(devnum)
            else:
                return attached_ctx
        else:
            if USE_NV_BINDING:
                devnum = int(devnum)
            return self._activate_context_for(devnum)

    def _get_or_create_context_uncached(self, devnum):
        """See also ``get_or_create_context(devnum)``.
        This version does not read the cache.
        """
        with self._lock:
            # Try to get the active context in the CUDA stack or
            # activate GPU-0 with the primary context
            with driver.get_active_context() as ac:
                if not ac:
                    return self._activate_context_for(0)
                else:
                    # Get primary context for the active device
                    ctx = self.gpus[ac.devnum].get_primary_context()
                    # Is active context the primary context?
                    if USE_NV_BINDING:
                        ctx_handle = int(ctx.handle)
                        ac_ctx_handle = int(ac.context_handle)
                    else:
                        ctx_handle = ctx.handle.value
                        ac_ctx_handle = ac.context_handle.value
                    if ctx_handle != ac_ctx_handle:
                        msg = ('Numba cannot operate on non-primary'
                               ' CUDA context {:x}')
                        raise RuntimeError(msg.format(ac_ctx_handle))
                    # Ensure the context is ready
                    ctx.prepare_for_use()
                return ctx

    def _activate_context_for(self, devnum):
        with self._lock:
            gpu = self.gpus[devnum]
            newctx = gpu.get_primary_context()
            # Detect unexpected context switch
            cached_ctx = self._get_attached_context()
            if cached_ctx is not None and cached_ctx is not newctx:
                raise RuntimeError('Cannot switch CUDA-context.')
            newctx.push()
            return newctx

    def _get_attached_context(self):
        return getattr(self._tls, 'attached_context', None)

    def _set_attached_context(self, ctx):
        self._tls.attached_context = ctx

    def reset(self):
        """Clear all contexts in the thread.  Destroy the context if and only
        if we are in the main thread.
        """
        # Pop all active context.
        while driver.pop_active_context() is not None:
            pass

        # If it is the main thread
        if threading.current_thread() == self._mainthread:
            self._destroy_all_contexts()

    def _destroy_all_contexts(self):
        # Reset all devices
        for gpu in self.gpus:
            gpu.reset()


_runtime = _Runtime()

# ================================ PUBLIC API ================================

gpus = _runtime.gpus


def get_context(devnum=None):
    """Get the current device or use a device by device number, and
    return the CUDA context.
    """
    return _runtime.get_or_create_context(devnum)


def require_context(fn):
    """
    A decorator that ensures a CUDA context is available when *fn* is executed.

    Note: The function *fn* cannot switch CUDA-context.
    """
    @functools.wraps(fn)
    def _require_cuda_context(*args, **kws):
        with _runtime.ensure_context():
            return fn(*args, **kws)

    return _require_cuda_context


def reset():
    """Reset the CUDA subsystem for the current thread.

    In the main thread:
    This removes all CUDA contexts.  Only use this at shutdown or for
    cleaning up between tests.

    In non-main threads:
    This clear the CUDA context stack only.

    """
    _runtime.reset()
