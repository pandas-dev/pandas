"""
API that are reported to numba.cuda
"""


import contextlib
import os

import numpy as np

from .cudadrv import devicearray, devices, driver
from numba.core import config
from numba.cuda.api_util import prepare_shape_strides_dtype

# NDarray device helper

require_context = devices.require_context
current_context = devices.get_context
gpus = devices.gpus


@require_context
def from_cuda_array_interface(desc, owner=None, sync=True):
    """Create a DeviceNDArray from a cuda-array-interface description.
    The ``owner`` is the owner of the underlying memory.
    The resulting DeviceNDArray will acquire a reference from it.

    If ``sync`` is ``True``, then the imported stream (if present) will be
    synchronized.
    """
    version = desc.get('version')
    # Mask introduced in version 1
    if 1 <= version:
        mask = desc.get('mask')
        # Would ideally be better to detect if the mask is all valid
        if mask is not None:
            raise NotImplementedError('Masked arrays are not supported')

    shape = desc['shape']
    strides = desc.get('strides')
    dtype = np.dtype(desc['typestr'])

    shape, strides, dtype = prepare_shape_strides_dtype(
        shape, strides, dtype, order='C')
    size = driver.memory_size_from_info(shape, strides, dtype.itemsize)

    devptr = driver.get_devptr_for_active_ctx(desc['data'][0])
    data = driver.MemoryPointer(
        current_context(), devptr, size=size, owner=owner)
    stream_ptr = desc.get('stream', None)
    if stream_ptr is not None:
        stream = external_stream(stream_ptr)
        if sync and config.CUDA_ARRAY_INTERFACE_SYNC:
            stream.synchronize()
    else:
        stream = 0 # No "Numba default stream", not the CUDA default stream
    da = devicearray.DeviceNDArray(shape=shape, strides=strides,
                                   dtype=dtype, gpu_data=data,
                                   stream=stream)
    return da


def as_cuda_array(obj, sync=True):
    """Create a DeviceNDArray from any object that implements
    the :ref:`cuda array interface <cuda-array-interface>`.

    A view of the underlying GPU buffer is created.  No copying of the data
    is done.  The resulting DeviceNDArray will acquire a reference from `obj`.

    If ``sync`` is ``True``, then the imported stream (if present) will be
    synchronized.
    """
    if not is_cuda_array(obj):
        raise TypeError("*obj* doesn't implement the cuda array interface.")
    else:
        return from_cuda_array_interface(obj.__cuda_array_interface__,
                                         owner=obj, sync=sync)


def is_cuda_array(obj):
    """Test if the object has defined the `__cuda_array_interface__` attribute.

    Does not verify the validity of the interface.
    """
    return hasattr(obj, '__cuda_array_interface__')


def is_float16_supported():
    """Whether 16-bit floats are supported.

    float16 is always supported in current versions of Numba - returns True.
    """
    return True


@require_context
def to_device(obj, stream=0, copy=True, to=None):
    """to_device(obj, stream=0, copy=True, to=None)

    Allocate and transfer a numpy ndarray or structured scalar to the device.

    To copy host->device a numpy array::

        ary = np.arange(10)
        d_ary = cuda.to_device(ary)

    To enqueue the transfer to a stream::

        stream = cuda.stream()
        d_ary = cuda.to_device(ary, stream=stream)

    The resulting ``d_ary`` is a ``DeviceNDArray``.

    To copy device->host::

        hary = d_ary.copy_to_host()

    To copy device->host to an existing array::

        ary = np.empty(shape=d_ary.shape, dtype=d_ary.dtype)
        d_ary.copy_to_host(ary)

    To enqueue the transfer to a stream::

        hary = d_ary.copy_to_host(stream=stream)
    """
    if to is None:
        to, new = devicearray.auto_device(obj, stream=stream, copy=copy,
                                          user_explicit=True)
        return to
    if copy:
        to.copy_to_device(obj, stream=stream)
    return to


@require_context
def device_array(shape, dtype=np.float64, strides=None, order='C', stream=0):
    """device_array(shape, dtype=np.float64, strides=None, order='C', stream=0)

    Allocate an empty device ndarray. Similar to :meth:`numpy.empty`.
    """
    shape, strides, dtype = prepare_shape_strides_dtype(shape, strides, dtype,
                                                        order)
    return devicearray.DeviceNDArray(shape=shape, strides=strides, dtype=dtype,
                                     stream=stream)


@require_context
def managed_array(shape, dtype=np.float64, strides=None, order='C', stream=0,
                  attach_global=True):
    """managed_array(shape, dtype=np.float64, strides=None, order='C', stream=0,
                     attach_global=True)

    Allocate a np.ndarray with a buffer that is managed.
    Similar to np.empty().

    Managed memory is supported on Linux / x86 and PowerPC, and is considered
    experimental on Windows and Linux / AArch64.

    :param attach_global: A flag indicating whether to attach globally. Global
                          attachment implies that the memory is accessible from
                          any stream on any device. If ``False``, attachment is
                          *host*, and memory is only accessible by devices
                          with Compute Capability 6.0 and later.
    """
    shape, strides, dtype = prepare_shape_strides_dtype(shape, strides, dtype,
                                                        order)
    bytesize = driver.memory_size_from_info(shape, strides, dtype.itemsize)
    buffer = current_context().memallocmanaged(bytesize,
                                               attach_global=attach_global)
    npary = np.ndarray(shape=shape, strides=strides, dtype=dtype, order=order,
                       buffer=buffer)
    managedview = np.ndarray.view(npary, type=devicearray.ManagedNDArray)
    managedview.device_setup(buffer, stream=stream)
    return managedview


@require_context
def pinned_array(shape, dtype=np.float64, strides=None, order='C'):
    """pinned_array(shape, dtype=np.float64, strides=None, order='C')

    Allocate an :class:`ndarray <numpy.ndarray>` with a buffer that is pinned
    (pagelocked).  Similar to :func:`np.empty() <numpy.empty>`.
    """
    shape, strides, dtype = prepare_shape_strides_dtype(shape, strides, dtype,
                                                        order)
    bytesize = driver.memory_size_from_info(shape, strides,
                                            dtype.itemsize)
    buffer = current_context().memhostalloc(bytesize)
    return np.ndarray(shape=shape, strides=strides, dtype=dtype, order=order,
                      buffer=buffer)


@require_context
def mapped_array(shape, dtype=np.float64, strides=None, order='C', stream=0,
                 portable=False, wc=False):
    """mapped_array(shape, dtype=np.float64, strides=None, order='C', stream=0,
                    portable=False, wc=False)

    Allocate a mapped ndarray with a buffer that is pinned and mapped on
    to the device. Similar to np.empty()

    :param portable: a boolean flag to allow the allocated device memory to be
              usable in multiple devices.
    :param wc: a boolean flag to enable writecombined allocation which is faster
        to write by the host and to read by the device, but slower to
        write by the host and slower to write by the device.
    """
    shape, strides, dtype = prepare_shape_strides_dtype(shape, strides, dtype,
                                                        order)
    bytesize = driver.memory_size_from_info(shape, strides, dtype.itemsize)
    buffer = current_context().memhostalloc(bytesize, mapped=True)
    npary = np.ndarray(shape=shape, strides=strides, dtype=dtype, order=order,
                       buffer=buffer)
    mappedview = np.ndarray.view(npary, type=devicearray.MappedNDArray)
    mappedview.device_setup(buffer, stream=stream)
    return mappedview


@contextlib.contextmanager
@require_context
def open_ipc_array(handle, shape, dtype, strides=None, offset=0):
    """
    A context manager that opens a IPC *handle* (*CUipcMemHandle*) that is
    represented as a sequence of bytes (e.g. *bytes*, tuple of int)
    and represent it as an array of the given *shape*, *strides* and *dtype*.
    The *strides* can be omitted.  In that case, it is assumed to be a 1D
    C contiguous array.

    Yields a device array.

    The IPC handle is closed automatically when context manager exits.
    """
    dtype = np.dtype(dtype)
    # compute size
    size = np.prod(shape) * dtype.itemsize
    # manually recreate the IPC mem handle
    if driver.USE_NV_BINDING:
        driver_handle = driver.binding.CUipcMemHandle()
        driver_handle.reserved = handle
    else:
        driver_handle = driver.drvapi.cu_ipc_mem_handle(*handle)
    # use *IpcHandle* to open the IPC memory
    ipchandle = driver.IpcHandle(None, driver_handle, size, offset=offset)
    yield ipchandle.open_array(current_context(), shape=shape,
                               strides=strides, dtype=dtype)
    ipchandle.close()


def synchronize():
    "Synchronize the current context."
    return current_context().synchronize()


def _contiguous_strides_like_array(ary):
    """
    Given an array, compute strides for a new contiguous array of the same
    shape.
    """
    # Don't recompute strides if the default strides will be sufficient to
    # create a contiguous array.
    if ary.flags['C_CONTIGUOUS'] or ary.flags['F_CONTIGUOUS'] or ary.ndim <= 1:
        return None

    # Otherwise, we need to compute new strides using an algorithm adapted from
    # NumPy v1.17.4's PyArray_NewLikeArrayWithShape in
    # core/src/multiarray/ctors.c. We permute the strides in ascending order
    # then compute the stride for the dimensions with the same permutation.

    # Stride permutation. E.g. a stride array (4, -2, 12) becomes
    # [(1, -2), (0, 4), (2, 12)]
    strideperm = [ x for x in enumerate(ary.strides) ]
    strideperm.sort(key=lambda x: x[1])

    # Compute new strides using permutation
    strides = [0] * len(ary.strides)
    stride = ary.dtype.itemsize
    for i_perm, _ in strideperm:
        strides[i_perm] = stride
        stride *= ary.shape[i_perm]
    return tuple(strides)


def _order_like_array(ary):
    if ary.flags['F_CONTIGUOUS'] and not ary.flags['C_CONTIGUOUS']:
        return 'F'
    else:
        return 'C'


def device_array_like(ary, stream=0):
    """
    Call :func:`device_array() <numba.cuda.device_array>` with information from
    the array.
    """
    strides = _contiguous_strides_like_array(ary)
    order = _order_like_array(ary)
    return device_array(shape=ary.shape, dtype=ary.dtype, strides=strides,
                        order=order, stream=stream)


def mapped_array_like(ary, stream=0, portable=False, wc=False):
    """
    Call :func:`mapped_array() <numba.cuda.mapped_array>` with the information
    from the array.
    """
    strides = _contiguous_strides_like_array(ary)
    order = _order_like_array(ary)
    return mapped_array(shape=ary.shape, dtype=ary.dtype, strides=strides,
                        order=order, stream=stream, portable=portable, wc=wc)


def pinned_array_like(ary):
    """
    Call :func:`pinned_array() <numba.cuda.pinned_array>` with the information
    from the array.
    """
    strides = _contiguous_strides_like_array(ary)
    order = _order_like_array(ary)
    return pinned_array(shape=ary.shape, dtype=ary.dtype, strides=strides,
                        order=order)


# Stream helper
@require_context
def stream():
    """
    Create a CUDA stream that represents a command queue for the device.
    """
    return current_context().create_stream()


@require_context
def default_stream():
    """
    Get the default CUDA stream. CUDA semantics in general are that the default
    stream is either the legacy default stream or the per-thread default stream
    depending on which CUDA APIs are in use. In Numba, the APIs for the legacy
    default stream are always the ones in use, but an option to use APIs for
    the per-thread default stream may be provided in future.
    """
    return current_context().get_default_stream()


@require_context
def legacy_default_stream():
    """
    Get the legacy default CUDA stream.
    """
    return current_context().get_legacy_default_stream()


@require_context
def per_thread_default_stream():
    """
    Get the per-thread default CUDA stream.
    """
    return current_context().get_per_thread_default_stream()


@require_context
def external_stream(ptr):
    """Create a Numba stream object for a stream allocated outside Numba.

    :param ptr: Pointer to the external stream to wrap in a Numba Stream
    :type ptr: int
    """
    return current_context().create_external_stream(ptr)


# Page lock
@require_context
@contextlib.contextmanager
def pinned(*arylist):
    """A context manager for temporary pinning a sequence of host ndarrays.
    """
    pmlist = []
    for ary in arylist:
        pm = current_context().mempin(ary, driver.host_pointer(ary),
                                      driver.host_memory_size(ary),
                                      mapped=False)
        pmlist.append(pm)
    yield


@require_context
@contextlib.contextmanager
def mapped(*arylist, **kws):
    """A context manager for temporarily mapping a sequence of host ndarrays.
    """
    assert not kws or 'stream' in kws, "Only accept 'stream' as keyword."
    stream = kws.get('stream', 0)
    pmlist = []
    devarylist = []
    for ary in arylist:
        pm = current_context().mempin(ary, driver.host_pointer(ary),
                                      driver.host_memory_size(ary),
                                      mapped=True)
        pmlist.append(pm)
        devary = devicearray.from_array_like(ary, gpu_data=pm, stream=stream)
        devarylist.append(devary)
    try:
        if len(devarylist) == 1:
            yield devarylist[0]
        else:
            yield devarylist
    finally:
        # When exiting from `with cuda.mapped(*arrs) as mapped_arrs:`, the name
        # `mapped_arrs` stays in scope, blocking automatic unmapping based on
        # reference count. We therefore invoke the finalizer manually.
        for pm in pmlist:
            pm.free()


def event(timing=True):
    """
    Create a CUDA event. Timing data is only recorded by the event if it is
    created with ``timing=True``.
    """
    evt = current_context().create_event(timing=timing)
    return evt


event_elapsed_time = driver.event_elapsed_time


# Device selection

def select_device(device_id):
    """
    Make the context associated with device *device_id* the current context.

    Returns a Device instance.

    Raises exception on error.
    """
    context = devices.get_context(device_id)
    return context.device


def get_current_device():
    "Get current device associated with the current thread"
    return current_context().device


def list_devices():
    "Return a list of all detected devices"
    return devices.gpus


def close():
    """
    Explicitly clears all contexts in the current thread, and destroys all
    contexts if the current thread is the main thread.
    """
    devices.reset()


def _auto_device(ary, stream=0, copy=True):
    return devicearray.auto_device(ary, stream=stream, copy=copy)


def detect():
    """
    Detect supported CUDA hardware and print a summary of the detected hardware.

    Returns a boolean indicating whether any supported devices were detected.
    """
    devlist = list_devices()
    print('Found %d CUDA devices' % len(devlist))
    supported_count = 0
    for dev in devlist:
        attrs = []
        cc = dev.compute_capability
        kernel_timeout = dev.KERNEL_EXEC_TIMEOUT
        tcc = dev.TCC_DRIVER
        fp32_to_fp64_ratio = dev.SINGLE_TO_DOUBLE_PRECISION_PERF_RATIO
        attrs += [('Compute Capability', '%d.%d' % cc)]
        attrs += [('PCI Device ID', dev.PCI_DEVICE_ID)]
        attrs += [('PCI Bus ID', dev.PCI_BUS_ID)]
        attrs += [('UUID', dev.uuid)]
        attrs += [('Watchdog', 'Enabled' if kernel_timeout else 'Disabled')]
        if os.name == "nt":
            attrs += [('Compute Mode', 'TCC' if tcc else 'WDDM')]
        attrs += [('FP32/FP64 Performance Ratio', fp32_to_fp64_ratio)]
        if cc < (3, 5):
            support = '[NOT SUPPORTED: CC < 3.5]'
        elif cc < (5, 0):
            support = '[SUPPORTED (DEPRECATED)]'
            supported_count += 1
        else:
            support = '[SUPPORTED]'
            supported_count += 1

        print('id %d    %20s %40s' % (dev.id, dev.name, support))
        for key, val in attrs:
            print('%40s: %s' % (key, val))

    print('Summary:')
    print('\t%d/%d devices are supported' % (supported_count, len(devlist)))
    return supported_count > 0


@contextlib.contextmanager
def defer_cleanup():
    """
    Temporarily disable memory deallocation.
    Use this to prevent resource deallocation breaking asynchronous execution.

    For example::

        with defer_cleanup():
            # all cleanup is deferred in here
            do_speed_critical_code()
        # cleanup can occur here

    Note: this context manager can be nested.
    """
    with current_context().defer_cleanup():
        yield


profiling = require_context(driver.profiling)
profile_start = require_context(driver.profile_start)
profile_stop = require_context(driver.profile_stop)
