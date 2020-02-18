from copy import copy
from distutils.version import LooseVersion

from cython import Py_ssize_t
from cpython.ref cimport Py_INCREF

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray,
                    int64_t,
                    PyArray_SETITEM,
                    PyArray_ITER_NEXT, PyArray_ITER_DATA, PyArray_IterNew,
                    flatiter)
cnp.import_array()

cimport pandas._libs.util as util
from pandas._libs.lib import maybe_convert_objects, is_scalar


cdef _check_result_array(object obj, Py_ssize_t cnt):

    if (util.is_array(obj) or
            (isinstance(obj, list) and len(obj) == cnt) or
            getattr(obj, 'shape', None) == (cnt,)):
        raise ValueError('Function does not reduce')


cdef class Reducer:
    """
    Performs generic reduction operation on a C or Fortran-contiguous ndarray
    while avoiding ndarray construction overhead
    """
    cdef:
        Py_ssize_t increment, chunksize, nresults
        object dummy, f, labels, typ, ityp, index
        ndarray arr

    def __init__(self, ndarray arr, object f, axis=1, dummy=None, labels=None):
        n, k = (<object>arr).shape

        if axis == 0:
            if not arr.flags.f_contiguous:
                arr = arr.copy('F')

            self.nresults = k
            self.chunksize = n
            self.increment = n * arr.dtype.itemsize
        else:
            if not arr.flags.c_contiguous:
                arr = arr.copy('C')

            self.nresults = n
            self.chunksize = k
            self.increment = k * arr.dtype.itemsize

        self.f = f
        self.arr = arr
        self.labels = labels
        self.dummy, self.typ, self.index, self.ityp = self._check_dummy(
            dummy=dummy)

    cdef _check_dummy(self, dummy=None):
        cdef:
            object index = None, typ = None, ityp = None

        if dummy is None:
            dummy = np.empty(self.chunksize, dtype=self.arr.dtype)

            # our ref is stolen later since we are creating this array
            # in cython, so increment first
            Py_INCREF(dummy)

        else:

            # we passed a Series
            typ = type(dummy)
            index = dummy.index
            dummy = dummy.values

            if dummy.dtype != self.arr.dtype:
                raise ValueError('Dummy array must be same dtype')
            if len(dummy) != self.chunksize:
                raise ValueError(f'Dummy array must be length {self.chunksize}')

        return dummy, typ, index, ityp

    def get_result(self):
        cdef:
            char* dummy_buf
            ndarray arr, result, chunk
            Py_ssize_t i
            flatiter it
            object res, name, labels
            object cached_typ = None

        arr = self.arr
        chunk = self.dummy
        dummy_buf = chunk.data
        chunk.data = arr.data
        labels = self.labels

        result = np.empty(self.nresults, dtype='O')
        it = <flatiter>PyArray_IterNew(result)

        try:
            for i in range(self.nresults):

                # create the cached type
                # each time just reassign the data
                if i == 0:

                    if self.typ is not None:
                        # In this case, we also have self.index
                        name = labels[i]
                        cached_typ = self.typ(
                            chunk, index=self.index, name=name, dtype=arr.dtype)

                # use the cached_typ if possible
                if cached_typ is not None:
                    # In this case, we also have non-None labels
                    name = labels[i]

                    object.__setattr__(
                        cached_typ._data._block, 'values', chunk)
                    object.__setattr__(cached_typ, 'name', name)
                    res = self.f(cached_typ)
                else:
                    res = self.f(chunk)

                # TODO: reason for not squeezing here?
                res = _extract_result(res, squeeze=False)
                if i == 0:
                    # On the first pass, we check the output shape to see
                    #  if this looks like a reduction.
                    _check_result_array(res, len(self.dummy))

                PyArray_SETITEM(result, PyArray_ITER_DATA(it), res)
                chunk.data = chunk.data + self.increment
                PyArray_ITER_NEXT(it)
        finally:
            # so we don't free the wrong memory
            chunk.data = dummy_buf

        result = maybe_convert_objects(result)
        return result


cdef inline _extract_result(object res, bint squeeze=True):
    """ extract the result object, it might be a 0-dim ndarray
        or a len-1 0-dim, or a scalar """
    if hasattr(res, 'values') and util.is_array(res.values):
        res = res.values
    if util.is_array(res):
        if res.ndim == 0:
            res = res.item()
        elif squeeze and res.ndim == 1 and len(res) == 1:
            res = res[0]
    return res


cdef class Slider:
    """
    Only handles contiguous data for now
    """
    cdef:
        ndarray values, buf
        Py_ssize_t stride, orig_len, orig_stride
        char *orig_data

    def __init__(self, ndarray values, ndarray buf):
        assert values.ndim == 1
        assert values.dtype == buf.dtype

        if not values.flags.contiguous:
            values = values.copy()

        self.values = values
        self.buf = buf
        self.stride = values.strides[0]

        self.orig_data = self.buf.data
        self.orig_len = self.buf.shape[0]
        self.orig_stride = self.buf.strides[0]

        self.buf.data = self.values.data
        self.buf.strides[0] = self.stride

    cdef advance(self, Py_ssize_t k):
        self.buf.data = <char*>self.buf.data + self.stride * k

    cdef move(self, int start, int end):
        """
        For slicing
        """
        self.buf.data = self.values.data + self.stride * start
        self.buf.shape[0] = end - start

    cdef set_length(self, Py_ssize_t length):
        self.buf.shape[0] = length

    cdef reset(self):

        self.buf.shape[0] = self.orig_len
        self.buf.data = self.orig_data
        self.buf.strides[0] = self.orig_stride


class InvalidApply(Exception):
    pass


def apply_frame_axis0(object frame, object f, object names,
                      const int64_t[:] starts, const int64_t[:] ends):
    cdef:
        BlockSlider slider
        Py_ssize_t i, n = len(starts)
        list results
        object piece
        dict item_cache

    # We have already checked that we don't have a MultiIndex before calling
    assert frame.index.nlevels == 1

    results = []

    slider = BlockSlider(frame)

    mutated = False
    item_cache = slider.dummy._item_cache
    try:
        for i in range(n):
            slider.move(starts[i], ends[i])

            item_cache.clear()  # ugh
            chunk = slider.dummy
            object.__setattr__(chunk, 'name', names[i])

            try:
                piece = f(chunk)
            except Exception:
                # We can't be more specific without knowing something about `f`
                raise InvalidApply('Let this error raise above us')

            # Need to infer if low level index slider will cause segfaults
            require_slow_apply = i == 0 and piece is chunk
            try:
                if piece.index is not chunk.index:
                    mutated = True
            except AttributeError:
                # `piece` might not have an index, could be e.g. an int
                pass

            if not is_scalar(piece):
                # Need to copy data to avoid appending references
                try:
                    piece = piece.copy(deep="all")
                except (TypeError, AttributeError):
                    piece = copy(piece)

            results.append(piece)

            # If the data was modified inplace we need to
            # take the slow path to not risk segfaults
            # we have already computed the first piece
            if require_slow_apply:
                break
    finally:
        slider.reset()

    return results, mutated


cdef class BlockSlider:
    """
    Only capable of sliding on axis=0
    """

    cdef public:
        object frame, dummy, index
        int nblocks
        Slider idx_slider
        list blocks

    cdef:
        char **base_ptrs

    def __init__(self, frame):
        self.frame = frame
        self.dummy = frame[:0]
        self.index = self.dummy.index

        self.blocks = [b.values for b in self.dummy._data.blocks]

        for x in self.blocks:
            util.set_array_not_contiguous(x)

        self.nblocks = len(self.blocks)
        # See the comment in indexes/base.py about _index_data.
        # We need this for EA-backed indexes that have a reference to a 1-d
        # ndarray like datetime / timedelta / period.
        self.idx_slider = Slider(
            self.frame.index._index_data, self.dummy.index._index_data)

        self.base_ptrs = <char**>malloc(sizeof(char*) * len(self.blocks))
        for i, block in enumerate(self.blocks):
            self.base_ptrs[i] = (<ndarray>block).data

    def __dealloc__(self):
        free(self.base_ptrs)

    cdef move(self, int start, int end):
        cdef:
            ndarray arr
            Py_ssize_t i

        # move blocks
        for i in range(self.nblocks):
            arr = self.blocks[i]

            # axis=1 is the frame's axis=0
            arr.data = self.base_ptrs[i] + arr.strides[1] * start
            arr.shape[1] = end - start

        # move and set the index
        self.idx_slider.move(start, end)

        object.__setattr__(self.index, '_index_data', self.idx_slider.buf)
        self.index._engine.clear_mapping()

    cdef reset(self):
        cdef:
            ndarray arr
            Py_ssize_t i

        # reset blocks
        for i in range(self.nblocks):
            arr = self.blocks[i]

            # axis=1 is the frame's axis=0
            arr.data = self.base_ptrs[i]
            arr.shape[1] = 0


def compute_reduction(arr: np.ndarray, f, axis: int = 0, dummy=None, labels=None):
    """

    Parameters
    -----------
    arr : np.ndarray
    f : function
    axis : integer axis
    dummy : type of reduced output (series)
    labels : Index or None
    """

    # We either have both dummy and labels, or neither of them
    if (labels is None) ^ (dummy is None):
        raise ValueError("Must pass either dummy and labels, or neither")

    if labels is not None:
        # Caller is responsible for ensuring we don't have MultiIndex
        assert labels.nlevels == 1

        # pass as an ndarray/ExtensionArray
        labels = labels._values

    reducer = Reducer(arr, f, axis=axis, dummy=dummy, labels=labels)
    return reducer.get_result()
