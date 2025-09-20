'''
The Device Array API is not implemented in the simulator. This module provides
stubs to allow tests to import correctly.
'''
from contextlib import contextmanager
from numba.np.numpy_support import numpy_version

import numpy as np


DeviceRecord = None
from_record_like = None


errmsg_contiguous_buffer = ("Array contains non-contiguous buffer and cannot "
                            "be transferred as a single memory region. Please "
                            "ensure contiguous buffer with numpy "
                            ".ascontiguousarray()")


class FakeShape(tuple):
    '''
    The FakeShape class is used to provide a shape which does not allow negative
    indexing, similar to the shape in CUDA Python. (Numpy shape arrays allow
    negative indexing)
    '''

    def __getitem__(self, k):
        if isinstance(k, int) and k < 0:
            raise IndexError('tuple index out of range')
        return super(FakeShape, self).__getitem__(k)


class FakeWithinKernelCUDAArray(object):
    '''
    Created to emulate the behavior of arrays within kernels, where either
    array.item or array['item'] is valid (that is, give all structured
    arrays `numpy.recarray`-like semantics). This behaviour does not follow
    the semantics of Python and NumPy with non-jitted code, and will be
    deprecated and removed.
    '''

    def __init__(self, item):
        assert isinstance(item, FakeCUDAArray)
        self.__dict__['_item'] = item

    def __wrap_if_fake(self, item):
        if isinstance(item, FakeCUDAArray):
            return FakeWithinKernelCUDAArray(item)
        else:
            return item

    def __getattr__(self, attrname):
        try:
            if attrname in dir(self._item._ary):  # For e.g. array size.
                return self.__wrap_if_fake(getattr(self._item._ary, attrname))
            else:
                return self.__wrap_if_fake(self._item.__getitem__(attrname))
        except Exception as e:
            if not isinstance(e, AttributeError):
                raise AttributeError(attrname) from e

    def __setattr__(self, nm, val):
        self._item.__setitem__(nm, val)

    def __getitem__(self, idx):
        return self.__wrap_if_fake(self._item.__getitem__(idx))

    def __setitem__(self, idx, val):
        self._item.__setitem__(idx, val)

    def __len__(self):
        return len(self._item)

    def __array_ufunc__(self, ufunc, method, *args, **kwargs):
        # ufuncs can only be called directly on instances of numpy.ndarray (not
        # things that implement its interfaces, like the FakeCUDAArray or
        # FakeWithinKernelCUDAArray). For other objects, __array_ufunc__ is
        # called when they are arguments to ufuncs, to provide an opportunity
        # to somehow implement the ufunc. Since the FakeWithinKernelCUDAArray
        # is just a thin wrapper over an ndarray, we can implement all ufuncs
        # by passing the underlying ndarrays to a call to the intended ufunc.
        call = getattr(ufunc, method)

        def convert_fakes(obj):
            if isinstance(obj, FakeWithinKernelCUDAArray):
                obj = obj._item._ary

            return obj

        out = kwargs.get('out')
        if out:
            kwargs['out'] = tuple(convert_fakes(o) for o in out)
        args = tuple(convert_fakes(a) for a in args)
        return call(*args, **kwargs)


class FakeCUDAArray(object):
    '''
    Implements the interface of a DeviceArray/DeviceRecord, but mostly just
    wraps a NumPy array.
    '''

    __cuda_ndarray__ = True  # There must be gpu_data attribute

    def __init__(self, ary, stream=0):
        self._ary = ary
        self.stream = stream

    @property
    def alloc_size(self):
        return self._ary.nbytes

    @property
    def nbytes(self):
        # return nbytes -- FakeCUDAArray is a wrapper around NumPy
        return self._ary.nbytes

    def __getattr__(self, attrname):
        try:
            attr = getattr(self._ary, attrname)
            return attr
        except AttributeError as e:
            msg = "Wrapped array has no attribute '%s'" % attrname
            raise AttributeError(msg) from e

    def bind(self, stream=0):
        return FakeCUDAArray(self._ary, stream)

    @property
    def T(self):
        return self.transpose()

    def transpose(self, axes=None):
        return FakeCUDAArray(np.transpose(self._ary, axes=axes))

    def __getitem__(self, idx):
        ret = self._ary.__getitem__(idx)
        if type(ret) not in [np.ndarray, np.void]:
            return ret
        else:
            return FakeCUDAArray(ret, stream=self.stream)

    def __setitem__(self, idx, val):
        return self._ary.__setitem__(idx, val)

    def copy_to_host(self, ary=None, stream=0):
        if ary is None:
            ary = np.empty_like(self._ary)
        else:
            check_array_compatibility(self, ary)
        np.copyto(ary, self._ary)
        return ary

    def copy_to_device(self, ary, stream=0):
        '''
        Copy from the provided array into this array.

        This may be less forgiving than the CUDA Python implementation, which
        will copy data up to the length of the smallest of the two arrays,
        whereas this expects the size of the arrays to be equal.
        '''
        sentry_contiguous(self)
        self_core, ary_core = array_core(self), array_core(ary)
        if isinstance(ary, FakeCUDAArray):
            sentry_contiguous(ary)
            check_array_compatibility(self_core, ary_core)
        else:
            ary_core = np.array(
                ary_core,
                order='C' if self_core.flags['C_CONTIGUOUS'] else 'F',
                subok=True,
                copy=False if numpy_version < (2, 0) else None)
            check_array_compatibility(self_core, ary_core)
        np.copyto(self_core._ary, ary_core)

    @property
    def shape(self):
        return FakeShape(self._ary.shape)

    def ravel(self, *args, **kwargs):
        return FakeCUDAArray(self._ary.ravel(*args, **kwargs))

    def reshape(self, *args, **kwargs):
        return FakeCUDAArray(self._ary.reshape(*args, **kwargs))

    def view(self, *args, **kwargs):
        return FakeCUDAArray(self._ary.view(*args, **kwargs))

    def is_c_contiguous(self):
        return self._ary.flags.c_contiguous

    def is_f_contiguous(self):
        return self._ary.flags.f_contiguous

    def __str__(self):
        return str(self._ary)

    def __repr__(self):
        return repr(self._ary)

    def __len__(self):
        return len(self._ary)

    # TODO: Add inplace, bitwise, unary magic methods
    #  (or maybe inherit this class from numpy)?
    def __eq__(self, other):
        return FakeCUDAArray(self._ary == other)

    def __ne__(self, other):
        return FakeCUDAArray(self._ary != other)

    def __lt__(self, other):
        return FakeCUDAArray(self._ary < other)

    def __le__(self, other):
        return FakeCUDAArray(self._ary <= other)

    def __gt__(self, other):
        return FakeCUDAArray(self._ary > other)

    def __ge__(self, other):
        return FakeCUDAArray(self._ary >= other)

    def __add__(self, other):
        return FakeCUDAArray(self._ary + other)

    def __sub__(self, other):
        return FakeCUDAArray(self._ary - other)

    def __mul__(self, other):
        return FakeCUDAArray(self._ary * other)

    def __floordiv__(self, other):
        return FakeCUDAArray(self._ary // other)

    def __truediv__(self, other):
        return FakeCUDAArray(self._ary / other)

    def __mod__(self, other):
        return FakeCUDAArray(self._ary % other)

    def __pow__(self, other):
        return FakeCUDAArray(self._ary ** other)

    def split(self, section, stream=0):
        return [
            FakeCUDAArray(a)
            for a in np.split(self._ary, range(section, len(self), section))
        ]


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


def sentry_contiguous(ary):
    core = array_core(ary)
    if not core.flags['C_CONTIGUOUS'] and not core.flags['F_CONTIGUOUS']:
        raise ValueError(errmsg_contiguous_buffer)


def check_array_compatibility(ary1, ary2):
    ary1sq, ary2sq = ary1.squeeze(), ary2.squeeze()
    if ary1.dtype != ary2.dtype:
        raise TypeError('incompatible dtype: %s vs. %s' %
                        (ary1.dtype, ary2.dtype))
    if ary1sq.shape != ary2sq.shape:
        raise ValueError('incompatible shape: %s vs. %s' %
                         (ary1.shape, ary2.shape))
    if ary1sq.strides != ary2sq.strides:
        raise ValueError('incompatible strides: %s vs. %s' %
                         (ary1.strides, ary2.strides))


def to_device(ary, stream=0, copy=True, to=None):
    ary = np.array(ary,
                   copy=False if numpy_version < (2, 0) else None,
                   subok=True)
    sentry_contiguous(ary)
    if to is None:
        buffer_dtype = np.int64 if ary.dtype.char in 'Mm' else ary.dtype
        return FakeCUDAArray(
            np.ndarray(
                buffer=np.copy(array_core(ary)).view(buffer_dtype),
                dtype=ary.dtype,
                shape=ary.shape,
                strides=ary.strides,
            ).view(type=type(ary)),
        )
    else:
        to.copy_to_device(ary, stream=stream)


@contextmanager
def pinned(arg):
    yield


def mapped_array(*args, **kwargs):
    for unused_arg in ('portable', 'wc'):
        if unused_arg in kwargs:
            kwargs.pop(unused_arg)
    return device_array(*args, **kwargs)


def pinned_array(shape, dtype=np.float64, strides=None, order='C'):
    return np.ndarray(shape=shape, strides=strides, dtype=dtype, order=order)


def managed_array(shape, dtype=np.float64, strides=None, order='C'):
    return np.ndarray(shape=shape, strides=strides, dtype=dtype, order=order)


def device_array(*args, **kwargs):
    stream = kwargs.pop('stream') if 'stream' in kwargs else 0
    return FakeCUDAArray(np.ndarray(*args, **kwargs), stream=stream)


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
    strides = _contiguous_strides_like_array(ary)
    order = _order_like_array(ary)
    return device_array(shape=ary.shape, dtype=ary.dtype, strides=strides,
                        order=order)


def pinned_array_like(ary):
    strides = _contiguous_strides_like_array(ary)
    order = _order_like_array(ary)
    return pinned_array(shape=ary.shape, dtype=ary.dtype, strides=strides,
                        order=order)


def auto_device(ary, stream=0, copy=True):
    if isinstance(ary, FakeCUDAArray):
        return ary, False

    if not isinstance(ary, np.void):
        ary = np.array(
            ary,
            copy=False if numpy_version < (2, 0) else None,
            subok=True)
    return to_device(ary, stream, copy), True


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
