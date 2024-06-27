"""
A CUDA ND Array is recognized by checking the __cuda_memory__ attribute
on the object.  If it exists and evaluate to True, it must define shape,
strides, dtype and size attributes similar to a NumPy ndarray.
"""

import math
import functools
import operator
import copy
from ctypes import c_void_p

import numpy as np

import numba
from numba import _devicearray
from numba.cuda.cudadrv import devices
from numba.cuda.cudadrv import driver as _driver
from numba.core import types, config
from numba.np.unsafe.ndarray import to_fixed_tuple
from numba.np.numpy_support import numpy_version
from numba.misc import dummyarray
from numba.np import numpy_support
from numba.cuda.api_util import prepare_shape_strides_dtype
from numba.core.errors import NumbaPerformanceWarning
from warnings import warn

try:
    lru_cache = getattr(functools, 'lru_cache')(None)
except AttributeError:
    # Python 3.1 or lower
    def lru_cache(func):
        return func


def is_cuda_ndarray(obj):
    "Check if an object is a CUDA ndarray"
    return getattr(obj, '__cuda_ndarray__', False)


def verify_cuda_ndarray_interface(obj):
    "Verify the CUDA ndarray interface for an obj"
    require_cuda_ndarray(obj)

    def requires_attr(attr, typ):
        if not hasattr(obj, attr):
            raise AttributeError(attr)
        if not isinstance(getattr(obj, attr), typ):
            raise AttributeError('%s must be of type %s' % (attr, typ))

    requires_attr('shape', tuple)
    requires_attr('strides', tuple)
    requires_attr('dtype', np.dtype)
    requires_attr('size', int)


def require_cuda_ndarray(obj):
    "Raises ValueError is is_cuda_ndarray(obj) evaluates False"
    if not is_cuda_ndarray(obj):
        raise ValueError('require an cuda ndarray object')


class DeviceNDArrayBase(_devicearray.DeviceArray):
    """A on GPU NDArray representation
    """
    __cuda_memory__ = True
    __cuda_ndarray__ = True     # There must be gpu_data attribute

    def __init__(self, shape, strides, dtype, stream=0, gpu_data=None):
        """
        Args
        ----

        shape
            array shape.
        strides
            array strides.
        dtype
            data type as np.dtype coercible object.
        stream
            cuda stream.
        gpu_data
            user provided device memory for the ndarray data buffer
        """
        if isinstance(shape, int):
            shape = (shape,)
        if isinstance(strides, int):
            strides = (strides,)
        dtype = np.dtype(dtype)
        self.ndim = len(shape)
        if len(strides) != self.ndim:
            raise ValueError('strides not match ndim')
        self._dummy = dummyarray.Array.from_desc(0, shape, strides,
                                                 dtype.itemsize)
        self.shape = tuple(shape)
        self.strides = tuple(strides)
        self.dtype = dtype
        self.size = int(functools.reduce(operator.mul, self.shape, 1))
        # prepare gpu memory
        if self.size > 0:
            if gpu_data is None:
                self.alloc_size = _driver.memory_size_from_info(
                    self.shape, self.strides, self.dtype.itemsize)
                gpu_data = devices.get_context().memalloc(self.alloc_size)
            else:
                self.alloc_size = _driver.device_memory_size(gpu_data)
        else:
            # Make NULL pointer for empty allocation
            if _driver.USE_NV_BINDING:
                null = _driver.binding.CUdeviceptr(0)
            else:
                null = c_void_p(0)
            gpu_data = _driver.MemoryPointer(context=devices.get_context(),
                                             pointer=null, size=0)
            self.alloc_size = 0

        self.gpu_data = gpu_data
        self.stream = stream

    @property
    def __cuda_array_interface__(self):
        if _driver.USE_NV_BINDING:
            if self.device_ctypes_pointer is not None:
                ptr = int(self.device_ctypes_pointer)
            else:
                ptr = 0
        else:
            if self.device_ctypes_pointer.value is not None:
                ptr = self.device_ctypes_pointer.value
            else:
                ptr = 0

        return {
            'shape': tuple(self.shape),
            'strides': None if is_contiguous(self) else tuple(self.strides),
            'data': (ptr, False),
            'typestr': self.dtype.str,
            'stream': int(self.stream) if self.stream != 0 else None,
            'version': 3,
        }

    def bind(self, stream=0):
        """Bind a CUDA stream to this object so that all subsequent operation
        on this array defaults to the given stream.
        """
        clone = copy.copy(self)
        clone.stream = stream
        return clone

    @property
    def T(self):
        return self.transpose()

    def transpose(self, axes=None):
        if axes and tuple(axes) == tuple(range(self.ndim)):
            return self
        elif self.ndim != 2:
            msg = "transposing a non-2D DeviceNDArray isn't supported"
            raise NotImplementedError(msg)
        elif axes is not None and set(axes) != set(range(self.ndim)):
            raise ValueError("invalid axes list %r" % (axes,))
        else:
            from numba.cuda.kernels.transpose import transpose
            return transpose(self)

    def _default_stream(self, stream):
        return self.stream if not stream else stream

    @property
    def _numba_type_(self):
        """
        Magic attribute expected by Numba to get the numba type that
        represents this object.
        """
        # Typing considerations:
        #
        # 1. The preference is to use 'C' or 'F' layout since this enables
        # hardcoding stride values into compiled kernels, which is more
        # efficient than storing a passed-in value in a register.
        #
        # 2. If an array is both C- and F-contiguous, prefer 'C' layout as it's
        # the more likely / common case.
        #
        # 3. If an array is broadcast then it must be typed as 'A' - using 'C'
        # or 'F' does not apply for broadcast arrays, because the strides, some
        # of which will be 0, will not match those hardcoded in for 'C' or 'F'
        # layouts.

        broadcast = 0 in self.strides
        if self.flags['C_CONTIGUOUS'] and not broadcast:
            layout = 'C'
        elif self.flags['F_CONTIGUOUS'] and not broadcast:
            layout = 'F'
        else:
            layout = 'A'

        dtype = numpy_support.from_dtype(self.dtype)
        return types.Array(dtype, self.ndim, layout)

    @property
    def device_ctypes_pointer(self):
        """Returns the ctypes pointer to the GPU data buffer
        """
        if self.gpu_data is None:
            if _driver.USE_NV_BINDING:
                return _driver.binding.CUdeviceptr(0)
            else:
                return c_void_p(0)
        else:
            return self.gpu_data.device_ctypes_pointer

    @devices.require_context
    def copy_to_device(self, ary, stream=0):
        """Copy `ary` to `self`.

        If `ary` is a CUDA memory, perform a device-to-device transfer.
        Otherwise, perform a a host-to-device transfer.
        """
        if ary.size == 0:
            # Nothing to do
            return

        sentry_contiguous(self)
        stream = self._default_stream(stream)

        self_core, ary_core = array_core(self), array_core(ary)
        if _driver.is_device_memory(ary):
            sentry_contiguous(ary)
            check_array_compatibility(self_core, ary_core)
            _driver.device_to_device(self, ary, self.alloc_size, stream=stream)
        else:
            # Ensure same contiguity. Only makes a host-side copy if necessary
            # (i.e., in order to materialize a writable strided view)
            ary_core = np.array(
                ary_core,
                order='C' if self_core.flags['C_CONTIGUOUS'] else 'F',
                subok=True,
                copy=(not ary_core.flags['WRITEABLE'])
                if numpy_version < (2, 0) else None)
            check_array_compatibility(self_core, ary_core)
            _driver.host_to_device(self, ary_core, self.alloc_size,
                                   stream=stream)

    @devices.require_context
    def copy_to_host(self, ary=None, stream=0):
        """Copy ``self`` to ``ary`` or create a new Numpy ndarray
        if ``ary`` is ``None``.

        If a CUDA ``stream`` is given, then the transfer will be made
        asynchronously as part as the given stream.  Otherwise, the transfer is
        synchronous: the function returns after the copy is finished.

        Always returns the host array.

        Example::

            import numpy as np
            from numba import cuda

            arr = np.arange(1000)
            d_arr = cuda.to_device(arr)

            my_kernel[100, 100](d_arr)

            result_array = d_arr.copy_to_host()
        """
        if any(s < 0 for s in self.strides):
            msg = 'D->H copy not implemented for negative strides: {}'
            raise NotImplementedError(msg.format(self.strides))
        assert self.alloc_size >= 0, "Negative memory size"
        stream = self._default_stream(stream)
        if ary is None:
            hostary = np.empty(shape=self.alloc_size, dtype=np.byte)
        else:
            check_array_compatibility(self, ary)
            hostary = ary

        if self.alloc_size != 0:
            _driver.device_to_host(hostary, self, self.alloc_size,
                                   stream=stream)

        if ary is None:
            if self.size == 0:
                hostary = np.ndarray(shape=self.shape, dtype=self.dtype,
                                     buffer=hostary)
            else:
                hostary = np.ndarray(shape=self.shape, dtype=self.dtype,
                                     strides=self.strides, buffer=hostary)
        return hostary

    def split(self, section, stream=0):
        """Split the array into equal partition of the `section` size.
        If the array cannot be equally divided, the last section will be
        smaller.
        """
        stream = self._default_stream(stream)
        if self.ndim != 1:
            raise ValueError("only support 1d array")
        if self.strides[0] != self.dtype.itemsize:
            raise ValueError("only support unit stride")
        nsect = int(math.ceil(float(self.size) / section))
        strides = self.strides
        itemsize = self.dtype.itemsize
        for i in range(nsect):
            begin = i * section
            end = min(begin + section, self.size)
            shape = (end - begin,)
            gpu_data = self.gpu_data.view(begin * itemsize, end * itemsize)
            yield DeviceNDArray(shape, strides, dtype=self.dtype, stream=stream,
                                gpu_data=gpu_data)

    def as_cuda_arg(self):
        """Returns a device memory object that is used as the argument.
        """
        return self.gpu_data

    def get_ipc_handle(self):
        """
        Returns a *IpcArrayHandle* object that is safe to serialize and transfer
        to another process to share the local allocation.

        Note: this feature is only available on Linux.
        """
        ipch = devices.get_context().get_ipc_handle(self.gpu_data)
        desc = dict(shape=self.shape, strides=self.strides, dtype=self.dtype)
        return IpcArrayHandle(ipc_handle=ipch, array_desc=desc)

    def squeeze(self, axis=None, stream=0):
        """
        Remove axes of size one from the array shape.

        Parameters
        ----------
        axis : None or int or tuple of ints, optional
            Subset of dimensions to remove. A `ValueError` is raised if an axis
            with size greater than one is selected. If `None`, all axes with
            size one are removed.
        stream : cuda stream or 0, optional
            Default stream for the returned view of the array.

        Returns
        -------
        DeviceNDArray
            Squeezed view into the array.

        """
        new_dummy, _ = self._dummy.squeeze(axis=axis)
        return DeviceNDArray(
            shape=new_dummy.shape,
            strides=new_dummy.strides,
            dtype=self.dtype,
            stream=self._default_stream(stream),
            gpu_data=self.gpu_data,
        )

    def view(self, dtype):
        """Returns a new object by reinterpretting the dtype without making a
        copy of the data.
        """
        dtype = np.dtype(dtype)
        shape = list(self.shape)
        strides = list(self.strides)

        if self.dtype.itemsize != dtype.itemsize:
            if not self.is_c_contiguous():
                raise ValueError(
                    "To change to a dtype of a different size,"
                    " the array must be C-contiguous"
                )

            shape[-1], rem = divmod(
                shape[-1] * self.dtype.itemsize,
                dtype.itemsize
            )

            if rem != 0:
                raise ValueError(
                    "When changing to a larger dtype,"
                    " its size must be a divisor of the total size in bytes"
                    " of the last axis of the array."
                )

            strides[-1] = dtype.itemsize

        return DeviceNDArray(
            shape=shape,
            strides=strides,
            dtype=dtype,
            stream=self.stream,
            gpu_data=self.gpu_data,
        )

    @property
    def nbytes(self):
        # Note: not using `alloc_size`.  `alloc_size` reports memory
        # consumption of the allocation, not the size of the array
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.nbytes.html
        return self.dtype.itemsize * self.size


class DeviceRecord(DeviceNDArrayBase):
    '''
    An on-GPU record type
    '''
    def __init__(self, dtype, stream=0, gpu_data=None):
        shape = ()
        strides = ()
        super(DeviceRecord, self).__init__(shape, strides, dtype, stream,
                                           gpu_data)

    @property
    def flags(self):
        """
        For `numpy.ndarray` compatibility. Ideally this would return a
        `np.core.multiarray.flagsobj`, but that needs to be constructed
        with an existing `numpy.ndarray` (as the C- and F- contiguous flags
        aren't writeable).
        """
        return dict(self._dummy.flags) # defensive copy

    @property
    def _numba_type_(self):
        """
        Magic attribute expected by Numba to get the numba type that
        represents this object.
        """
        return numpy_support.from_dtype(self.dtype)

    @devices.require_context
    def __getitem__(self, item):
        return self._do_getitem(item)

    @devices.require_context
    def getitem(self, item, stream=0):
        """Do `__getitem__(item)` with CUDA stream
        """
        return self._do_getitem(item, stream)

    def _do_getitem(self, item, stream=0):
        stream = self._default_stream(stream)
        typ, offset = self.dtype.fields[item]
        newdata = self.gpu_data.view(offset)

        if typ.shape == ():
            if typ.names is not None:
                return DeviceRecord(dtype=typ, stream=stream,
                                    gpu_data=newdata)
            else:
                hostary = np.empty(1, dtype=typ)
                _driver.device_to_host(dst=hostary, src=newdata,
                                       size=typ.itemsize,
                                       stream=stream)
            return hostary[0]
        else:
            shape, strides, dtype = \
                prepare_shape_strides_dtype(typ.shape,
                                            None,
                                            typ.subdtype[0], 'C')
            return DeviceNDArray(shape=shape, strides=strides,
                                 dtype=dtype, gpu_data=newdata,
                                 stream=stream)

    @devices.require_context
    def __setitem__(self, key, value):
        return self._do_setitem(key, value)

    @devices.require_context
    def setitem(self, key, value, stream=0):
        """Do `__setitem__(key, value)` with CUDA stream
        """
        return self._do_setitem(key, value, stream=stream)

    def _do_setitem(self, key, value, stream=0):

        stream = self._default_stream(stream)

        # If the record didn't have a default stream, and the user didn't
        # provide a stream, then we will use the default stream for the
        # assignment kernel and synchronize on it.
        synchronous = not stream
        if synchronous:
            ctx = devices.get_context()
            stream = ctx.get_default_stream()

        # (1) prepare LHS

        typ, offset = self.dtype.fields[key]
        newdata = self.gpu_data.view(offset)

        lhs = type(self)(dtype=typ, stream=stream, gpu_data=newdata)

        # (2) prepare RHS

        rhs, _ = auto_device(lhs.dtype.type(value), stream=stream)

        # (3) do the copy

        _driver.device_to_device(lhs, rhs, rhs.dtype.itemsize, stream)

        if synchronous:
            stream.synchronize()


@lru_cache
def _assign_kernel(ndim):
    """
    A separate method so we don't need to compile code every assignment (!).

    :param ndim: We need to have static array sizes for cuda.local.array, so
        bake in the number of dimensions into the kernel
    """
    from numba import cuda  # circular!

    if ndim == 0:
        # the (2, ndim) allocation below is not yet supported, so avoid it
        @cuda.jit
        def kernel(lhs, rhs):
            lhs[()] = rhs[()]
        return kernel

    @cuda.jit
    def kernel(lhs, rhs):
        location = cuda.grid(1)

        n_elements = 1
        for i in range(lhs.ndim):
            n_elements *= lhs.shape[i]
        if location >= n_elements:
            # bake n_elements into the kernel, better than passing it in
            # as another argument.
            return

        # [0, :] is the to-index (into `lhs`)
        # [1, :] is the from-index (into `rhs`)
        idx = cuda.local.array(
            shape=(2, ndim),
            dtype=types.int64)

        for i in range(ndim - 1, -1, -1):
            idx[0, i] = location % lhs.shape[i]
            idx[1, i] = (location % lhs.shape[i]) * (rhs.shape[i] > 1)
            location //= lhs.shape[i]

        lhs[to_fixed_tuple(idx[0], ndim)] = rhs[to_fixed_tuple(idx[1], ndim)]
    return kernel


class DeviceNDArray(DeviceNDArrayBase):
    '''
    An on-GPU array type
    '''
    def is_f_contiguous(self):
        '''
        Return true if the array is Fortran-contiguous.
        '''
        return self._dummy.is_f_contig

    @property
    def flags(self):
        """
        For `numpy.ndarray` compatibility. Ideally this would return a
        `np.core.multiarray.flagsobj`, but that needs to be constructed
        with an existing `numpy.ndarray` (as the C- and F- contiguous flags
        aren't writeable).
        """
        return dict(self._dummy.flags) # defensive copy

    def is_c_contiguous(self):
        '''
        Return true if the array is C-contiguous.
        '''
        return self._dummy.is_c_contig

    def __array__(self, dtype=None):
        """
        :return: an `numpy.ndarray`, so copies to the host.
        """
        if dtype:
            return self.copy_to_host().__array__(dtype)
        else:
            return self.copy_to_host().__array__()

    def __len__(self):
        return self.shape[0]

    def reshape(self, *newshape, **kws):
        """
        Reshape the array without changing its contents, similarly to
        :meth:`numpy.ndarray.reshape`. Example::

            d_arr = d_arr.reshape(20, 50, order='F')
        """
        if len(newshape) == 1 and isinstance(newshape[0], (tuple, list)):
            newshape = newshape[0]

        cls = type(self)
        if newshape == self.shape:
            # nothing to do
            return cls(shape=self.shape, strides=self.strides,
                       dtype=self.dtype, gpu_data=self.gpu_data)

        newarr, extents = self._dummy.reshape(*newshape, **kws)

        if extents == [self._dummy.extent]:
            return cls(shape=newarr.shape, strides=newarr.strides,
                       dtype=self.dtype, gpu_data=self.gpu_data)
        else:
            raise NotImplementedError("operation requires copying")

    def ravel(self, order='C', stream=0):
        '''
        Flattens a contiguous array without changing its contents, similar to
        :meth:`numpy.ndarray.ravel`. If the array is not contiguous, raises an
        exception.
        '''
        stream = self._default_stream(stream)
        cls = type(self)
        newarr, extents = self._dummy.ravel(order=order)

        if extents == [self._dummy.extent]:
            return cls(shape=newarr.shape, strides=newarr.strides,
                       dtype=self.dtype, gpu_data=self.gpu_data,
                       stream=stream)

        else:
            raise NotImplementedError("operation requires copying")

    @devices.require_context
    def __getitem__(self, item):
        return self._do_getitem(item)

    @devices.require_context
    def getitem(self, item, stream=0):
        """Do `__getitem__(item)` with CUDA stream
        """
        return self._do_getitem(item, stream)

    def _do_getitem(self, item, stream=0):
        stream = self._default_stream(stream)

        arr = self._dummy.__getitem__(item)
        extents = list(arr.iter_contiguous_extent())
        cls = type(self)
        if len(extents) == 1:
            newdata = self.gpu_data.view(*extents[0])

            if not arr.is_array:
                # Check for structured array type (record)
                if self.dtype.names is not None:
                    return DeviceRecord(dtype=self.dtype, stream=stream,
                                        gpu_data=newdata)
                else:
                    # Element indexing
                    hostary = np.empty(1, dtype=self.dtype)
                    _driver.device_to_host(dst=hostary, src=newdata,
                                           size=self._dummy.itemsize,
                                           stream=stream)
                return hostary[0]
            else:
                return cls(shape=arr.shape, strides=arr.strides,
                           dtype=self.dtype, gpu_data=newdata, stream=stream)
        else:
            newdata = self.gpu_data.view(*arr.extent)
            return cls(shape=arr.shape, strides=arr.strides,
                       dtype=self.dtype, gpu_data=newdata, stream=stream)

    @devices.require_context
    def __setitem__(self, key, value):
        return self._do_setitem(key, value)

    @devices.require_context
    def setitem(self, key, value, stream=0):
        """Do `__setitem__(key, value)` with CUDA stream
        """
        return self._do_setitem(key, value, stream=stream)

    def _do_setitem(self, key, value, stream=0):

        stream = self._default_stream(stream)

        # If the array didn't have a default stream, and the user didn't provide
        # a stream, then we will use the default stream for the assignment
        # kernel and synchronize on it.
        synchronous = not stream
        if synchronous:
            ctx = devices.get_context()
            stream = ctx.get_default_stream()

        # (1) prepare LHS

        arr = self._dummy.__getitem__(key)
        newdata = self.gpu_data.view(*arr.extent)

        if isinstance(arr, dummyarray.Element):
            # convert to a 0d array
            shape = ()
            strides = ()
        else:
            shape = arr.shape
            strides = arr.strides

        lhs = type(self)(
            shape=shape,
            strides=strides,
            dtype=self.dtype,
            gpu_data=newdata,
            stream=stream)

        # (2) prepare RHS

        rhs, _ = auto_device(value, stream=stream, user_explicit=True)
        if rhs.ndim > lhs.ndim:
            raise ValueError("Can't assign %s-D array to %s-D self" % (
                rhs.ndim,
                lhs.ndim))
        rhs_shape = np.ones(lhs.ndim, dtype=np.int64)
        # negative indices would not work if rhs.ndim == 0
        rhs_shape[lhs.ndim - rhs.ndim:] = rhs.shape
        rhs = rhs.reshape(*rhs_shape)
        for i, (l, r) in enumerate(zip(lhs.shape, rhs.shape)):
            if r != 1 and l != r:
                raise ValueError("Can't copy sequence with size %d to array "
                                 "axis %d with dimension %d" % ( r, i, l))

        # (3) do the copy

        n_elements = functools.reduce(operator.mul, lhs.shape, 1)
        _assign_kernel(lhs.ndim).forall(n_elements, stream=stream)(lhs, rhs)
        if synchronous:
            stream.synchronize()


class IpcArrayHandle(object):
    """
    An IPC array handle that can be serialized and transfer to another process
    in the same machine for share a GPU allocation.

    On the destination process, use the *.open()* method to creates a new
    *DeviceNDArray* object that shares the allocation from the original process.
    To release the resources, call the *.close()* method.  After that, the
    destination can no longer use the shared array object.  (Note: the
    underlying weakref to the resource is now dead.)

    This object implements the context-manager interface that calls the
    *.open()* and *.close()* method automatically::

        with the_ipc_array_handle as ipc_array:
            # use ipc_array here as a normal gpu array object
            some_code(ipc_array)
        # ipc_array is dead at this point
    """
    def __init__(self, ipc_handle, array_desc):
        self._array_desc = array_desc
        self._ipc_handle = ipc_handle

    def open(self):
        """
        Returns a new *DeviceNDArray* that shares the allocation from the
        original process.  Must not be used on the original process.
        """
        dptr = self._ipc_handle.open(devices.get_context())
        return DeviceNDArray(gpu_data=dptr, **self._array_desc)

    def close(self):
        """
        Closes the IPC handle to the array.
        """
        self._ipc_handle.close()

    def __enter__(self):
        return self.open()

    def __exit__(self, type, value, traceback):
        self.close()


class MappedNDArray(DeviceNDArrayBase, np.ndarray):
    """
    A host array that uses CUDA mapped memory.
    """

    def device_setup(self, gpu_data, stream=0):
        self.gpu_data = gpu_data
        self.stream = stream


class ManagedNDArray(DeviceNDArrayBase, np.ndarray):
    """
    A host array that uses CUDA managed memory.
    """

    def device_setup(self, gpu_data, stream=0):
        self.gpu_data = gpu_data
        self.stream = stream


def from_array_like(ary, stream=0, gpu_data=None):
    "Create a DeviceNDArray object that is like ary."
    return DeviceNDArray(ary.shape, ary.strides, ary.dtype, stream=stream,
                         gpu_data=gpu_data)


def from_record_like(rec, stream=0, gpu_data=None):
    "Create a DeviceRecord object that is like rec."
    return DeviceRecord(rec.dtype, stream=stream, gpu_data=gpu_data)


def array_core(ary):
    """
    Extract the repeated core of a broadcast array.

    Broadcast arrays are by definition non-contiguous due to repeated
    dimensions, i.e., dimensions with stride 0. In order to ascertain memory
    contiguity and copy the underlying data from such arrays, we must create
    a view without the repeated dimensions.

    """
    if not ary.strides or not ary.size:
        return ary
    core_index = []
    for stride in ary.strides:
        core_index.append(0 if stride == 0 else slice(None))
    return ary[tuple(core_index)]


def is_contiguous(ary):
    """
    Returns True iff `ary` is C-style contiguous while ignoring
    broadcasted and 1-sized dimensions.
    As opposed to array_core(), it does not call require_context(),
    which can be quite expensive.
    """
    size = ary.dtype.itemsize
    for shape, stride in zip(reversed(ary.shape), reversed(ary.strides)):
        if shape > 1 and stride != 0:
            if size != stride:
                return False
            size *= shape
    return True


errmsg_contiguous_buffer = ("Array contains non-contiguous buffer and cannot "
                            "be transferred as a single memory region. Please "
                            "ensure contiguous buffer with numpy "
                            ".ascontiguousarray()")


def sentry_contiguous(ary):
    core = array_core(ary)
    if not core.flags['C_CONTIGUOUS'] and not core.flags['F_CONTIGUOUS']:
        raise ValueError(errmsg_contiguous_buffer)


def auto_device(obj, stream=0, copy=True, user_explicit=False):
    """
    Create a DeviceRecord or DeviceArray like obj and optionally copy data from
    host to device. If obj already represents device memory, it is returned and
    no copy is made.
    """
    if _driver.is_device_memory(obj):
        return obj, False
    elif hasattr(obj, '__cuda_array_interface__'):
        return numba.cuda.as_cuda_array(obj), False
    else:
        if isinstance(obj, np.void):
            devobj = from_record_like(obj, stream=stream)
        else:
            # This allows you to pass non-array objects like constants and
            # objects implementing the array interface
            # https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.interface.html
            # into this function (with no overhead -- copies -- for `obj`s
            # that are already `ndarray`s.
            obj = np.array(
                obj,
                copy=False if numpy_version < (2, 0) else None,
                subok=True)
            sentry_contiguous(obj)
            devobj = from_array_like(obj, stream=stream)
        if copy:
            if config.CUDA_WARN_ON_IMPLICIT_COPY:
                if (
                    not user_explicit and
                    (not isinstance(obj, DeviceNDArray)
                     and isinstance(obj, np.ndarray))
                ):
                    msg = ("Host array used in CUDA kernel will incur "
                           "copy overhead to/from device.")
                    warn(NumbaPerformanceWarning(msg))
            devobj.copy_to_device(obj, stream=stream)
        return devobj, True


def check_array_compatibility(ary1, ary2):
    ary1sq, ary2sq = ary1.squeeze(), ary2.squeeze()
    if ary1.dtype != ary2.dtype:
        raise TypeError('incompatible dtype: %s vs. %s' %
                        (ary1.dtype, ary2.dtype))
    if ary1sq.shape != ary2sq.shape:
        raise ValueError('incompatible shape: %s vs. %s' %
                         (ary1.shape, ary2.shape))
    # We check strides only if the size is nonzero, because strides are
    # irrelevant (and can differ) for zero-length copies.
    if ary1.size and ary1sq.strides != ary2sq.strides:
        raise ValueError('incompatible strides: %s vs. %s' %
                         (ary1.strides, ary2.strides))
