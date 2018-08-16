# -*- coding: utf-8 -*-
from distutils.version import LooseVersion

from cython cimport Py_ssize_t
from cpython cimport Py_INCREF

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray,
                    int64_t,
                    PyArray_SETITEM,
                    PyArray_ITER_NEXT, PyArray_ITER_DATA, PyArray_IterNew,
                    flatiter)
cnp.import_array()

cimport util
from lib import maybe_convert_objects

is_numpy_prior_1_6_2 = LooseVersion(np.__version__) < '1.6.2'


cdef _get_result_array(object obj, Py_ssize_t size, Py_ssize_t cnt):

    if (util.is_array(obj) or
            isinstance(obj, list) and len(obj) == cnt or
            getattr(obj, 'shape', None) == (cnt,)):
        raise ValueError('function does not reduce')

    return np.empty(size, dtype='O')


cdef class Reducer:
    """
    Performs generic reduction operation on a C or Fortran-contiguous ndarray
    while avoiding ndarray construction overhead
    """
    cdef:
        Py_ssize_t increment, chunksize, nresults
        object arr, dummy, f, labels, typ, ityp, index

    def __init__(self, object arr, object f, axis=1, dummy=None,
                 labels=None):
        n, k = arr.shape

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

    def _check_dummy(self, dummy=None):
        cdef object index=None, typ=None, ityp=None

        if dummy is None:
            dummy = np.empty(self.chunksize, dtype=self.arr.dtype)

            # our ref is stolen later since we are creating this array
            # in cython, so increment first
            Py_INCREF(dummy)

        else:

            # we passed a series-like
            if hasattr(dummy, 'values'):

                typ = type(dummy)
                index = getattr(dummy, 'index', None)
                dummy = dummy.values

            if dummy.dtype != self.arr.dtype:
                raise ValueError('Dummy array must be same dtype')
            if len(dummy) != self.chunksize:
                raise ValueError('Dummy array must be length %d' %
                                 self.chunksize)

        return dummy, typ, index, ityp

    def get_result(self):
        cdef:
            char* dummy_buf
            ndarray arr, result, chunk
            Py_ssize_t i, incr
            flatiter it
            bint has_labels
            object res, name, labels, index
            object cached_typ=None

        arr = self.arr
        chunk = self.dummy
        dummy_buf = chunk.data
        chunk.data = arr.data
        labels = self.labels
        has_labels = labels is not None
        has_index = self.index is not None
        incr = self.increment

        try:
            for i in range(self.nresults):

                if has_labels:
                    name = util.get_value_at(labels, i)
                else:
                    name = None

                # create the cached type
                # each time just reassign the data
                if i == 0:

                    if self.typ is not None:

                        # recreate with the index if supplied
                        if has_index:

                            cached_typ = self.typ(
                                chunk, index=self.index, name=name)

                        else:

                            # use the passsed typ, sans index
                            cached_typ = self.typ(chunk, name=name)

                # use the cached_typ if possible
                if cached_typ is not None:

                    if has_index:
                        object.__setattr__(cached_typ, 'index', self.index)

                    object.__setattr__(
                        cached_typ._data._block, 'values', chunk)
                    object.__setattr__(cached_typ, 'name', name)
                    res = self.f(cached_typ)
                else:
                    res = self.f(chunk)

                if hasattr(res, 'values') and util.is_array(res.values):
                    res = res.values
                if i == 0:
                    result = _get_result_array(res,
                                               self.nresults,
                                               len(self.dummy))
                    it = <flatiter> PyArray_IterNew(result)

                PyArray_SETITEM(result, PyArray_ITER_DATA(it), res)
                chunk.data = chunk.data + self.increment
                PyArray_ITER_NEXT(it)
        except Exception, e:
            if hasattr(e, 'args'):
                e.args = e.args + (i,)
            raise
        finally:
            # so we don't free the wrong memory
            chunk.data = dummy_buf

        if result.dtype == np.object_:
            result = maybe_convert_objects(result)

        return result


cdef class SeriesBinGrouper:
    """
    Performs grouping operation according to bin edges, rather than labels
    """
    cdef:
        Py_ssize_t nresults, ngroups
        bint passed_dummy

    cdef public:
        object arr, index, dummy_arr, dummy_index
        object values, f, bins, typ, ityp, name

    def __init__(self, object series, object f, object bins, object dummy):
        n = len(series)

        self.bins = bins
        self.f = f

        values = series.values
        if not values.flags.c_contiguous:
            values = values.copy('C')
        self.arr = values
        self.typ = series._constructor
        self.ityp = series.index._constructor
        self.index = series.index.values
        self.name = getattr(series, 'name', None)

        self.dummy_arr, self.dummy_index = self._check_dummy(dummy)
        self.passed_dummy = dummy is not None

        # kludge for #1688
        if len(bins) > 0 and bins[-1] == len(series):
            self.ngroups = len(bins)
        else:
            self.ngroups = len(bins) + 1

    def _check_dummy(self, dummy=None):
        # both values and index must be an ndarray!

        if dummy is None:
            values = np.empty(0, dtype=self.arr.dtype)
            index = None
        else:
            values = dummy.values
            if values.dtype != self.arr.dtype:
                raise ValueError('Dummy array must be same dtype')
            if not values.flags.contiguous:
                values = values.copy()
            index = dummy.index.values
            if not index.flags.contiguous:
                index = index.copy()

        return values, index

    def get_result(self):
        cdef:
            ndarray arr, result
            ndarray[int64_t] counts
            Py_ssize_t i, n, group_size
            object res
            bint initialized = 0
            Slider vslider, islider
            object name, cached_typ=None, cached_ityp=None

        counts = np.zeros(self.ngroups, dtype=np.int64)

        if self.ngroups > 0:
            counts[0] = self.bins[0]
            for i in range(1, self.ngroups):
                if i == self.ngroups - 1:
                    counts[i] = len(self.arr) - self.bins[i - 1]
                else:
                    counts[i] = self.bins[i] - self.bins[i - 1]

        group_size = 0
        n = len(self.arr)
        name = self.name

        vslider = Slider(self.arr, self.dummy_arr)
        islider = Slider(self.index, self.dummy_index)

        try:
            for i in range(self.ngroups):
                group_size = counts[i]

                islider.set_length(group_size)
                vslider.set_length(group_size)

                if cached_typ is None:
                    cached_ityp = self.ityp(islider.buf)
                    cached_typ = self.typ(vslider.buf, index=cached_ityp,
                                          name=name)
                else:
                    object.__setattr__(cached_ityp, '_data', islider.buf)
                    cached_ityp._engine.clear_mapping()
                    object.__setattr__(
                        cached_typ._data._block, 'values', vslider.buf)
                    object.__setattr__(cached_typ, '_index', cached_ityp)
                    object.__setattr__(cached_typ, 'name', name)

                cached_ityp._engine.clear_mapping()
                res = self.f(cached_typ)
                res = _extract_result(res)
                if not initialized:
                    initialized = 1
                    result = _get_result_array(res,
                                               self.ngroups,
                                               len(self.dummy_arr))

                util.assign_value_1d(result, i, res)

                islider.advance(group_size)
                vslider.advance(group_size)

        except:
            raise
        finally:
            # so we don't free the wrong memory
            islider.reset()
            vslider.reset()

        if result.dtype == np.object_:
            result = maybe_convert_objects(result)

        return result, counts


cdef class SeriesGrouper:
    """
    Performs generic grouping operation while avoiding ndarray construction
    overhead
    """
    cdef:
        Py_ssize_t nresults, ngroups
        bint passed_dummy

    cdef public:
        object arr, index, dummy_arr, dummy_index
        object f, labels, values, typ, ityp, name

    def __init__(self, object series, object f, object labels,
                 Py_ssize_t ngroups, object dummy):
        n = len(series)

        self.labels = labels
        self.f = f

        values = series.values
        if not values.flags.c_contiguous:
            values = values.copy('C')
        self.arr = values
        self.typ = series._constructor
        self.ityp = series.index._constructor
        self.index = series.index.values
        self.name = getattr(series, 'name', None)

        self.dummy_arr, self.dummy_index = self._check_dummy(dummy)
        self.passed_dummy = dummy is not None
        self.ngroups = ngroups

    def _check_dummy(self, dummy=None):
        # both values and index must be an ndarray!

        if dummy is None:
            values = np.empty(0, dtype=self.arr.dtype)
            index = None
        else:
            values = dummy.values
            if dummy.dtype != self.arr.dtype:
                raise ValueError('Dummy array must be same dtype')
            if not values.flags.contiguous:
                values = values.copy()
            index = dummy.index.values
            if not index.flags.contiguous:
                index = index.copy()

        return values, index

    def get_result(self):
        cdef:
            ndarray arr, result
            ndarray[int64_t] labels, counts
            Py_ssize_t i, n, group_size, lab
            object res
            bint initialized = 0
            Slider vslider, islider
            object name, cached_typ=None, cached_ityp=None

        labels = self.labels
        counts = np.zeros(self.ngroups, dtype=np.int64)
        group_size = 0
        n = len(self.arr)
        name = self.name

        vslider = Slider(self.arr, self.dummy_arr)
        islider = Slider(self.index, self.dummy_index)

        try:
            for i in range(n):
                group_size += 1

                lab = labels[i]

                if i == n - 1 or lab != labels[i + 1]:
                    if lab == -1:
                        islider.advance(group_size)
                        vslider.advance(group_size)
                        group_size = 0
                        continue

                    islider.set_length(group_size)
                    vslider.set_length(group_size)

                    if cached_typ is None:
                        cached_ityp = self.ityp(islider.buf)
                        cached_typ = self.typ(vslider.buf, index=cached_ityp,
                                              name=name)
                    else:
                        object.__setattr__(cached_ityp, '_data', islider.buf)
                        cached_ityp._engine.clear_mapping()
                        object.__setattr__(
                            cached_typ._data._block, 'values', vslider.buf)
                        object.__setattr__(cached_typ, '_index', cached_ityp)
                        object.__setattr__(cached_typ, 'name', name)

                    cached_ityp._engine.clear_mapping()
                    res = self.f(cached_typ)
                    res = _extract_result(res)
                    if not initialized:
                        initialized = 1
                        result = _get_result_array(res,
                                                   self.ngroups,
                                                   len(self.dummy_arr))

                    util.assign_value_1d(result, lab, res)
                    counts[lab] = group_size
                    islider.advance(group_size)
                    vslider.advance(group_size)

                    group_size = 0

        except:
            raise
        finally:
            # so we don't free the wrong memory
            islider.reset()
            vslider.reset()

        if result.dtype == np.object_:
            result = maybe_convert_objects(result)

        return result, counts


cdef inline _extract_result(object res):
    """ extract the result object, it might be a 0-dim ndarray
        or a len-1 0-dim, or a scalar """
    if hasattr(res, 'values') and util.is_array(res.values):
        res = res.values
    if not np.isscalar(res):
        if util.is_array(res):
            if res.ndim == 0:
                res = res.item()
            elif res.ndim == 1 and len(res) == 1:
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

    def __init__(self, object values, object buf):
        assert(values.ndim == 1)

        if not values.flags.contiguous:
            values = values.copy()

        assert(values.dtype == buf.dtype)
        self.values = values
        self.buf = buf
        self.stride = values.strides[0]

        self.orig_data = self.buf.data
        self.orig_len = self.buf.shape[0]
        self.orig_stride = self.buf.strides[0]

        self.buf.data = self.values.data
        self.buf.strides[0] = self.stride

    cpdef advance(self, Py_ssize_t k):
        self.buf.data = <char*> self.buf.data + self.stride * k

    cdef move(self, int start, int end):
        """
        For slicing
        """
        self.buf.data = self.values.data + self.stride * start
        self.buf.shape[0] = end - start

    cpdef set_length(self, Py_ssize_t length):
        self.buf.shape[0] = length

    cpdef reset(self):

        self.buf.shape[0] = self.orig_len
        self.buf.data = self.orig_data
        self.buf.strides[0] = self.orig_stride


class InvalidApply(Exception):
    pass


def apply_frame_axis0(object frame, object f, object names,
                      ndarray[int64_t] starts, ndarray[int64_t] ends):
    cdef:
        BlockSlider slider
        Py_ssize_t i, n = len(starts)
        list results
        object piece
        dict item_cache

    if frame.index._has_complex_internals:
        raise InvalidApply('Cannot modify frame index internals')

    results = []

    # Need to infer if our low-level mucking is going to cause a segfault
    if n > 0:
        chunk = frame.iloc[starts[0]:ends[0]]
        object.__setattr__(chunk, 'name', names[0])
        try:
            result = f(chunk)
            if result is chunk:
                raise InvalidApply('Function unsafe for fast apply')
        except:
            raise InvalidApply('Let this error raise above us')

    slider = BlockSlider(frame)

    mutated = False
    item_cache = slider.dummy._item_cache
    try:
        for i in range(n):
            slider.move(starts[i], ends[i])

            item_cache.clear()  # ugh

            object.__setattr__(slider.dummy, 'name', names[i])
            piece = f(slider.dummy)

            # I'm paying the price for index-sharing, ugh
            try:
                if piece.index is slider.dummy.index:
                    piece = piece.copy(deep='all')
                else:
                    mutated = True
            except AttributeError:
                pass

            results.append(piece)
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
        self.idx_slider = Slider(
            self.frame.index.values, self.dummy.index.values)

        self.base_ptrs = <char**> malloc(sizeof(char*) * len(self.blocks))
        for i, block in enumerate(self.blocks):
            self.base_ptrs[i] = (<ndarray> block).data

    def __dealloc__(self):
        free(self.base_ptrs)

    cpdef move(self, int start, int end):
        cdef:
            ndarray arr
            object index

        # move blocks
        for i in range(self.nblocks):
            arr = self.blocks[i]

            # axis=1 is the frame's axis=0
            arr.data = self.base_ptrs[i] + arr.strides[1] * start
            arr.shape[1] = end - start

        # move and set the index
        self.idx_slider.move(start, end)
        object.__setattr__(self.index, '_data', self.idx_slider.buf)
        self.index._engine.clear_mapping()

    cdef reset(self):
        cdef:
            ndarray arr

        # reset blocks
        for i in range(self.nblocks):
            arr = self.blocks[i]

            # axis=1 is the frame's axis=0
            arr.data = self.base_ptrs[i]
            arr.shape[1] = 0


def reduce(arr, f, axis=0, dummy=None, labels=None):
    """

    Parameters
    -----------
    arr : NDFrame object
    f : function
    axis : integer axis
    dummy : type of reduced output (series)
    labels : Index or None
    """

    if labels is not None:
        if labels._has_complex_internals:
            raise Exception('Cannot use shortcut')

        # pass as an ndarray
        if hasattr(labels, 'values'):
            labels = labels.values

    reducer = Reducer(arr, f, axis=axis, dummy=dummy, labels=labels)
    return reducer.get_result()
