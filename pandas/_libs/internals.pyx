import cython
from collections import defaultdict
from cython import Py_ssize_t

from cpython.slice cimport PySlice_GetIndicesEx

cdef extern from "Python.h":
    Py_ssize_t PY_SSIZE_T_MAX

import numpy as np
cimport numpy as cnp
from numpy cimport NPY_INT64, int64_t
cnp.import_array()

from pandas._libs.algos import ensure_int64


cdef class BlockPlacement:
    # __slots__ = '_as_slice', '_as_array', '_len'
    cdef:
        slice _as_slice
        object _as_array

        bint _has_slice, _has_array, _is_known_slice_like

    def __init__(self, val):
        cdef:
            slice slc

        self._as_slice = None
        self._as_array = None
        self._has_slice = False
        self._has_array = False

        if isinstance(val, slice):
            slc = slice_canonize(val)

            if slc.start != slc.stop:
                self._as_slice = slc
                self._has_slice = True
            else:
                arr = np.empty(0, dtype=np.int64)
                self._as_array = arr
                self._has_array = True
        else:
            # Cython memoryview interface requires ndarray to be writeable.
            arr = np.require(val, dtype=np.int64, requirements='W')
            assert arr.ndim == 1
            self._as_array = arr
            self._has_array = True

    def __str__(self) -> str:
        cdef:
            slice s = self._ensure_has_slice()
        if s is not None:
            v = self._as_slice
        else:
            v = self._as_array

        return f'{type(self).__name__}({v})'

    def __repr__(self) -> str:
        return str(self)

    def __len__(self) -> int:
        cdef:
            slice s = self._ensure_has_slice()
        if s is not None:
            return slice_len(s)
        else:
            return len(self._as_array)

    def __iter__(self):
        cdef:
            slice s = self._ensure_has_slice()
            Py_ssize_t start, stop, step, _
        if s is not None:
            start, stop, step, _ = slice_get_indices_ex(s)
            return iter(range(start, stop, step))
        else:
            return iter(self._as_array)

    @property
    def as_slice(self) -> slice:
        cdef:
            slice s = self._ensure_has_slice()
        if s is None:
            raise TypeError('Not slice-like')
        else:
            return s

    @property
    def indexer(self):
        cdef:
            slice s = self._ensure_has_slice()
        if s is not None:
            return s
        else:
            return self._as_array

    def isin(self, arr):
        from pandas.core.indexes.api import Int64Index
        return Int64Index(self.as_array, copy=False).isin(arr)

    @property
    def as_array(self):
        cdef:
            Py_ssize_t start, stop, end, _
        if not self._has_array:
            start, stop, step, _ = slice_get_indices_ex(self._as_slice)
            # NOTE: this is the C-optimized equivalent of
            #  np.arange(start, stop, step, dtype=np.int64)
            self._as_array = cnp.PyArray_Arange(start, stop, step, NPY_INT64)
            self._has_array = True
        return self._as_array

    @property
    def is_slice_like(self) -> bool:
        cdef:
            slice s = self._ensure_has_slice()
        return s is not None

    def __getitem__(self, loc):
        cdef:
            slice s = self._ensure_has_slice()
        if s is not None:
            val = slice_getitem(s, loc)
        else:
            val = self._as_array[loc]

        if not isinstance(val, slice) and val.ndim == 0:
            return val

        return BlockPlacement(val)

    def delete(self, loc):
        return BlockPlacement(np.delete(self.as_array, loc, axis=0))

    def append(self, others):
        if len(others) == 0:
            return self

        return BlockPlacement(np.concatenate([self.as_array] +
                                             [o.as_array for o in others]))

    cdef iadd(self, other):
        cdef:
            slice s = self._ensure_has_slice()
            Py_ssize_t other_int, start, stop, step, l

        if isinstance(other, int) and s is not None:
            other_int = <Py_ssize_t>other

            if other_int == 0:
                # BlockPlacement is treated as immutable
                return self

            start, stop, step, l = slice_get_indices_ex(s)
            start += other_int
            stop += other_int

            if ((step > 0 and start < 0) or
                    (step < 0 and stop < step)):
                raise ValueError("iadd causes length change")

            if stop < 0:
                val = slice(start, None, step)
            else:
                val = slice(start, stop, step)

            return BlockPlacement(val)
        else:
            newarr = self.as_array + other
            if (newarr < 0).any():
                raise ValueError("iadd causes length change")

            val = newarr
            return BlockPlacement(val)

    def add(self, other):
        return self.iadd(other)

    def sub(self, other):
        return self.add(-other)

    cdef slice _ensure_has_slice(self):
        if not self._has_slice:
            self._as_slice = indexer_as_slice(self._as_array)
            self._has_slice = True
        return self._as_slice


cdef slice slice_canonize(slice s):
    """
    Convert slice to canonical bounded form.
    """
    cdef:
        Py_ssize_t start = 0, stop = 0, step = 1, length

    if s.step is None:
        step = 1
    else:
        step = <Py_ssize_t>s.step
        if step == 0:
            raise ValueError("slice step cannot be zero")

    if step > 0:
        if s.stop is None:
            raise ValueError("unbounded slice")

        stop = <Py_ssize_t>s.stop
        if s.start is None:
            start = 0
        else:
            start = <Py_ssize_t>s.start
            if start > stop:
                start = stop
    elif step < 0:
        if s.start is None:
            raise ValueError("unbounded slice")

        start = <Py_ssize_t>s.start
        if s.stop is None:
            stop = -1
        else:
            stop = <Py_ssize_t>s.stop
            if stop > start:
                stop = start

    if start < 0 or (stop < 0 and s.stop is not None):
        raise ValueError("unbounded slice")

    if stop < 0:
        return slice(start, None, step)
    else:
        return slice(start, stop, step)


cpdef Py_ssize_t slice_len(
        slice slc, Py_ssize_t objlen=PY_SSIZE_T_MAX) except -1:
    """
    Get length of a bounded slice.

    The slice must not have any "open" bounds that would create dependency on
    container size, i.e.:
    - if ``s.step is None or s.step > 0``, ``s.stop`` is not ``None``
    - if ``s.step < 0``, ``s.start`` is not ``None``

    Otherwise, the result is unreliable.
    """
    cdef:
        Py_ssize_t start, stop, step, length

    if slc is None:
        raise TypeError("slc must be slice")

    PySlice_GetIndicesEx(slc, objlen,
                         &start, &stop, &step, &length)

    return length


cdef slice_get_indices_ex(slice slc, Py_ssize_t objlen=PY_SSIZE_T_MAX):
    """
    Get (start, stop, step, length) tuple for a slice.

    If `objlen` is not specified, slice must be bounded, otherwise the result
    will be wrong.
    """
    cdef:
        Py_ssize_t start, stop, step, length

    if slc is None:
        raise TypeError("slc should be a slice")

    PySlice_GetIndicesEx(slc, objlen,
                         &start, &stop, &step, &length)

    return start, stop, step, length


cdef slice_getitem(slice slc, ind):
    cdef:
        Py_ssize_t s_start, s_stop, s_step, s_len
        Py_ssize_t ind_start, ind_stop, ind_step, ind_len

    s_start, s_stop, s_step, s_len = slice_get_indices_ex(slc)

    if isinstance(ind, slice):
        ind_start, ind_stop, ind_step, ind_len = slice_get_indices_ex(ind, s_len)

        if ind_step > 0 and ind_len == s_len:
            # short-cut for no-op slice
            if ind_len == s_len:
                return slc

        if ind_step < 0:
            s_start = s_stop - s_step
            ind_step = -ind_step

        s_step *= ind_step
        s_stop = s_start + ind_stop * s_step
        s_start = s_start + ind_start * s_step

        if s_step < 0 and s_stop < 0:
            return slice(s_start, None, s_step)
        else:
            return slice(s_start, s_stop, s_step)

    else:
        # NOTE:
        # this is the C-optimized equivalent of
        # `np.arange(s_start, s_stop, s_step, dtype=np.int64)[ind]`
        return cnp.PyArray_Arange(s_start, s_stop, s_step, NPY_INT64)[ind]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef slice indexer_as_slice(int64_t[:] vals):
    cdef:
        Py_ssize_t i, n, start, stop
        int64_t d

    if vals is None:
        raise TypeError("vals must be ndarray")

    n = vals.shape[0]

    if n == 0 or vals[0] < 0:
        return None

    if n == 1:
        return slice(vals[0], vals[0] + 1, 1)

    if vals[1] < 0:
        return None

    # n > 2
    d = vals[1] - vals[0]

    if d == 0:
        return None

    for i in range(2, n):
        if vals[i] < 0 or vals[i] - vals[i - 1] != d:
            return None

    start = vals[0]
    stop = start + n * d
    if stop < 0 and d < 0:
        return slice(start, None, d)
    else:
        return slice(start, stop, d)


@cython.boundscheck(False)
@cython.wraparound(False)
def get_blkno_indexers(int64_t[:] blknos, bint group=True):
    """
    Enumerate contiguous runs of integers in ndarray.

    Iterate over elements of `blknos` yielding ``(blkno, slice(start, stop))``
    pairs for each contiguous run found.

    If `group` is True and there is more than one run for a certain blkno,
    ``(blkno, array)`` with an array containing positions of all elements equal
    to blkno.

    Returns
    -------
    iter : iterator of (int, slice or array)
    """
    # There's blkno in this function's name because it's used in block &
    # blockno handling.
    cdef:
        int64_t cur_blkno
        Py_ssize_t i, start, stop, n, diff

        object blkno
        object group_dict = defaultdict(list)
        int64_t[:] res_view

    n = blknos.shape[0]

    if n == 0:
        return

    start = 0
    cur_blkno = blknos[start]

    if group is False:
        for i in range(1, n):
            if blknos[i] != cur_blkno:
                yield cur_blkno, slice(start, i)

                start = i
                cur_blkno = blknos[i]

        yield cur_blkno, slice(start, n)
    else:
        for i in range(1, n):
            if blknos[i] != cur_blkno:
                group_dict[cur_blkno].append((start, i))

                start = i
                cur_blkno = blknos[i]

        group_dict[cur_blkno].append((start, n))

        for blkno, slices in group_dict.items():
            if len(slices) == 1:
                yield blkno, slice(slices[0][0], slices[0][1])
            else:
                tot_len = sum(stop - start for start, stop in slices)
                result = np.empty(tot_len, dtype=np.int64)
                res_view = result

                i = 0
                for start, stop in slices:
                    for diff in range(start, stop):
                        res_view[i] = diff
                        i += 1

                yield blkno, result


def get_blkno_placements(blknos, group: bool = True):
    """
    Parameters
    ----------
    blknos : array of int64
    group : bool, default True

    Returns
    -------
    iterator
        yield (BlockPlacement, blkno)
    """
    blknos = ensure_int64(blknos)

    for blkno, indexer in get_blkno_indexers(blknos, group):
        yield blkno, BlockPlacement(indexer)
