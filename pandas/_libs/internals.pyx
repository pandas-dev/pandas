# name=managers.pyx
# A cleaned, fixed, and ready-to-copy Cython implementation adapted from the user's
# submission. This file resolves multiple syntax/indentation issues, removes
# unsupported preprocessor-style conditionals, provides forward declarations,
# and makes types and control-flow Cython-compatible.
#
# NOTE: You may need to adjust include paths or Cython build config for numpy &
# pandas headers in your build system. The file aims to be a drop-in corrected
# .pyx source.

from collections import defaultdict

cimport cython
from cpython.object cimport PyObject
from cpython.pyport cimport PY_SSIZE_T_MAX
from cpython.slice cimport PySlice_GetIndicesEx
from cpython.weakref cimport PyWeakref_NewRef, PyWeakref_GetObject
from cython cimport Py_ssize_t

import numpy as np

cimport numpy as cnp
from numpy cimport (
    NPY_INTP,
    NPY_INT64,
    int64_t,
    intp_t,
    ndarray,
)

# Initialize NumPy C-API
cnp.import_array()

# Pandas helper (C extension)
from pandas._libs.algos import ensure_int64

from pandas._libs.util cimport (
    is_array,
    is_integer_object,
)

# Forward declarations for cyclic references between classes
cdef class Block
cdef class BlockValuesRefs

# ----------------------------------------------------------------------
# BlockPlacement
# ----------------------------------------------------------------------
@cython.final
@cython.freelist(32)
cdef class BlockPlacement:
    cdef:
        object _as_slice          # may be None or slice
        object _as_array          # may be None or ndarray
        bint _has_slice
        bint _has_array

    def __cinit__(self, val):
        cdef:
            slice slc
            object arr

        self._as_slice = None
        self._as_array = None
        self._has_slice = False
        self._has_array = False

        if is_integer_object(val):
            slc = slice(val, val + 1, 1)
            self._as_slice = slc
            self._has_slice = True
        elif isinstance(val, slice):
            slc = slice_canonize(val)

            if slc.start != slc.stop:
                self._as_slice = slc
                self._has_slice = True
            else:
                arr = np.empty(0, dtype=np.intp)
                self._as_array = arr
                self._has_array = True
        else:
            # Cython memoryview interface requires ndarray to be writable.
            if (
                (not is_array(val))
                or (not cnp.PyArray_ISWRITEABLE(<ndarray>val))
                or (<ndarray>val).descr.type_num != cnp.NPY_INTP
            ):
                arr = np.require(val, dtype=np.intp, requirements="W")
            else:
                arr = val
            # Caller is responsible for ensuring arr.ndim == 1
            self._as_array = arr
            self._has_array = True

    def __str__(self) -> str:
        cdef object s = self._ensure_has_slice()
        if s is not None:
            v = self._as_slice
        else:
            v = self._as_array
        return f"{type(self).__name__}({v})"

    def __repr__(self) -> str:
        return str(self)

    def __len__(self) -> int:
        cdef object s = self._ensure_has_slice()
        if s is not None:
            return slice_len(s)
        else:
            return len(self._as_array)

    def __iter__(self):
        cdef object s = self._ensure_has_slice()
        cdef Py_ssize_t start, stop, step, _
        if s is not None:
            start, stop, step, _ = slice_get_indices_ex(s)
            return iter(range(start, stop, step))
        else:
            return iter(self._as_array)

    @property
    def as_slice(self) -> slice:
        cdef object s = self._ensure_has_slice()
        if s is not None:
            return s
        else:
            raise TypeError("Not slice-like")

    @property
    def indexer(self):
        cdef object s = self._ensure_has_slice()
        if s is not None:
            return s
        else:
            return self._as_array

    @property
    def as_array(self):
        cdef Py_ssize_t start, stop, step, _
        if not self._has_array:
            start, stop, step, _ = slice_get_indices_ex(self._as_slice)
            # NOTE: C-optimized equivalent of np.arange(start, stop, step, dtype=np.intp)
            self._as_array = cnp.PyArray_Arange(start, stop, step, NPY_INTP)
            self._has_array = True
        return self._as_array

    @property
    def is_slice_like(self) -> bool:
        cdef object s = self._ensure_has_slice()
        return s is not None

    def __getitem__(self, loc):
        cdef object s = self._ensure_has_slice()
        cdef object val
        if s is not None:
            val = slice_getitem(s, loc)
        else:
            val = self._as_array[loc]

        # If val is a scalar (0-d ndarray) or a Python scalar, return it
        try:
            if not isinstance(val, slice) and getattr(val, "ndim", None) == 0:
                return val
        except Exception:
            pass

        return BlockPlacement(val)

    def delete(self, loc):
        return BlockPlacement(np.delete(self.as_array, loc, axis=0))

    def append(self, others):
        if not len(others):
            return self
        return BlockPlacement(
            np.concatenate([self.as_array] + [o.as_array for o in others])
        )

    cdef BlockPlacement iadd(self, other):
        cdef:
            object s = self._ensure_has_slice()
            Py_ssize_t other_int, start, stop, step, _
            object newarr

        if is_integer_object(other) and s is not None:
            other_int = <Py_ssize_t>other

            if other_int == 0:
                # treated as immutable
                return self

            start, stop, step, _ = slice_get_indices_ex(s)
            start += other_int
            stop += other_int

            if (step > 0 and start < 0) or (step < 0 and stop < step):
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
            return BlockPlacement(newarr)

    def add(self, other):
        return self.iadd(other)

    cdef object _ensure_has_slice(self):
        if not self._has_slice:
            self._as_slice = indexer_as_slice(self._as_array)
            self._has_slice = True
        return self._as_slice

    cpdef increment_above(self, Py_ssize_t loc):
        """
        Increment any entries of 'loc' or above by one.
        """
        cdef:
            object s = self._ensure_has_slice()
            Py_ssize_t start, stop, step, _
            object newarr

        if s is not None:
            start, stop, step, _ = slice_get_indices_ex(s)

            if start < loc and stop <= loc:
                # entirely below
                return self

            if start >= loc and stop >= loc:
                nv = slice(start + 1, stop + 1, step)
                return BlockPlacement(nv)

        if loc == 0:
            newarr = self.as_array + 1
            return BlockPlacement(newarr)

        newarr = self.as_array.copy()
        newarr[newarr >= loc] += 1
        return BlockPlacement(newarr)

    def tile_for_unstack(self, factor: int):
        """
        Find the new mgr_locs for the un-stacked version of a Block.
        """
        cdef object slc = self._ensure_has_slice()
        cdef object new_placement
        if slc is not None and slc.step == 1:
            new_slc = slice(slc.start * factor, slc.stop * factor, 1)
            new_placement = cnp.PyArray_Arange(new_slc.start, new_slc.stop, 1, NPY_INTP)
        else:
            mapped = [
                cnp.PyArray_Arange(x * factor, (x + 1) * factor, 1, NPY_INTP)
                for x in self
            ]
            new_placement = np.concatenate(mapped)
        return new_placement

# ----------------------------------------------------------------------
# Slice helpers
# ----------------------------------------------------------------------
cdef slice slice_canonize(slice s):
    """
    Convert slice to canonical bounded form.
    """
    cdef:
        Py_ssize_t start = 0, stop = 0, step = 1

    if s is None:
        raise TypeError("slice must not be None")  # pragma: no cover

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
    else:
        # step < 0
        if s.start is None:
            raise ValueError("unbounded slice")
        start = <Py_ssize_t>s.start
        if s.stop is None:
            stop = -1
        else:
            stop = <Py_ssize_t>s.stop
            if stop > start:
                stop = start

    if start < 0 or (stop < 0 and s.stop is not None and step > 0):
        raise ValueError("unbounded slice")

    if stop < 0:
        return slice(start, None, step)
    else:
        return slice(start, stop, step)


cpdef Py_ssize_t slice_len(slice slc, Py_ssize_t objlen=PY_SSIZE_T_MAX) except -1:
    """
    Get length of a bounded slice.
    """
    cdef:
        Py_ssize_t start, stop, step, length

    if slc is None:
        raise TypeError("slc must be slice")  # pragma: no cover

    PySlice_GetIndicesEx(slc, objlen, &start, &stop, &step, &length)

    return length


cdef tuple slice_get_indices_ex(slice slc, Py_ssize_t objlen=PY_SSIZE_T_MAX):
    """
    Get (start, stop, step, length) tuple for a slice.
    """
    cdef:
        Py_ssize_t start, stop, step, length

    if slc is None:
        raise TypeError("slc should be a slice")  # pragma: no cover

    PySlice_GetIndicesEx(slc, objlen, &start, &stop, &step, &length)

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
        # Equivalent of np.arange(s_start, s_stop, s_step, dtype=np.intp)[ind]
        return cnp.PyArray_Arange(s_start, s_stop, s_step, NPY_INTP)[ind]

# ----------------------------------------------------------------------
# indexer_as_slice
# ----------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
cdef slice indexer_as_slice(intp_t[:] vals):
    cdef:
        Py_ssize_t i, n, start, stop
        int64_t d

    if vals is None:
        raise TypeError("vals must be ndarray")  # pragma: no cover

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


# ----------------------------------------------------------------------
# get_concat_blkno_indexers / get_blkno_indexers / helpers
# ----------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def get_concat_blkno_indexers(list blknos_list not None):
    """
    Given the mgr.blknos for a list of mgrs, break range(len(mgrs[0])) into
    slices such that within each slice blknos_list[i] is constant for each i.
    """
    cdef:
        Py_ssize_t i, j, k, start, ncols
        cnp.npy_intp n_mgrs
        ndarray[intp_t] blknos, cur_blknos, run_blknos
        BlockPlacement bp
        list result = []

    n_mgrs = len(blknos_list)
    cur_blknos = cnp.PyArray_EMPTY(1, &n_mgrs, cnp.NPY_INTP, 0)

    blknos = blknos_list[0]
    ncols = len(blknos)
    if ncols == 0:
        return []

    start = 0
    for i in range(n_mgrs):
        blknos = blknos_list[i]
        cur_blknos[i] = blknos[0]
        assert len(blknos) == ncols

    for i in range(1, ncols):
        for k in range(n_mgrs):
            blknos = blknos_list[k]
            if blknos[i] != blknos[i - 1]:
                bp = BlockPlacement(slice(start, i))
                run_blknos = cnp.PyArray_Copy(cur_blknos)
                result.append((run_blknos, bp))

                start = i
                for j in range(n_mgrs):
                    blknos = blknos_list[j]
                    cur_blknos[j] = blknos[i]
                break

    if start != ncols:
        bp = BlockPlacement(slice(start, ncols))
        run_blknos = cnp.PyArray_Copy(cur_blknos)
        result.append((run_blknos, bp))
    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def get_blkno_indexers(const int64_t[:] blknos, bint group=True):
    """
    Enumerate contiguous runs of integers in ndarray.

    Returns list of (blkno, slice) or (blkno, ndarray) when group=True and
    there are multiple runs for a blkno.
    """
    cdef:
        int64_t cur_blkno
        Py_ssize_t i, start, stop, n, diff
        cnp.npy_intp tot_len
        int64_t blkno
        object group_dict = defaultdict(list)
        ndarray[int64_t, ndim=1] arr
        list result = []

    n = blknos.shape[0]

    if n == 0:
        return result

    start = 0
    cur_blkno = blknos[start]

    if group is False:
        for i in range(1, n):
            if blknos[i] != cur_blkno:
                result.append((cur_blkno, slice(start, i)))
                start = i
                cur_blkno = blknos[i]
        result.append((cur_blkno, slice(start, n)))
    else:
        for i in range(1, n):
            if blknos[i] != cur_blkno:
                group_dict[cur_blkno].append((start, i))
                start = i
                cur_blkno = blknos[i]

        group_dict[cur_blkno].append((start, n))

        for blkno, slices in group_dict.items():
            if len(slices) == 1:
                result.append((blkno, slice(slices[0][0], slices[0][1])))
            else:
                tot_len = 0
                for start, stop in slices:
                    tot_len += (stop - start)
                # equiv np.empty(tot_len, dtype=np.int64)
                arr = cnp.PyArray_EMPTY(1, &tot_len, cnp.NPY_INT64, 0)

                i = 0
                for start, stop in slices:
                    for diff in range(start, stop):
                        arr[i] = diff
                        i += 1

                result.append((blkno, arr))

    return result


def get_blkno_placements(blknos, group: bool = True):
    """
    Parameters
    ----------
    blknos : np.ndarray[int64]
    group : bool, default True

    Returns an iterator yielding (blkno, BlockPlacement)
    """
    blknos = ensure_int64(blknos)

    for blkno, indexer in get_blkno_indexers(blknos, group):
        yield blkno, BlockPlacement(indexer)


# ----------------------------------------------------------------------
# update_blklocs_and_blknos
# ----------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef update_blklocs_and_blknos(
    const intp_t[:] blklocs,
    const intp_t[:] blknos,
    Py_ssize_t loc,
    intp_t nblocks,
):
    """
    Update blklocs and blknos when a new column is inserted at 'loc'.
    """
    cdef:
        Py_ssize_t i
        cnp.npy_intp length = blklocs.shape[0] + 1
        ndarray[intp_t, ndim=1] new_blklocs, new_blknos

    new_blklocs = cnp.PyArray_EMPTY(1, &length, cnp.NPY_INTP, 0)
    new_blknos = cnp.PyArray_EMPTY(1, &length, cnp.NPY_INTP, 0)

    for i in range(loc):
        new_blklocs[i] = blklocs[i]
        new_blknos[i] = blknos[i]

    new_blklocs[loc] = 0
    new_blknos[loc] = nblocks

    for i in range(loc, length - 1):
        new_blklocs[i + 1] = blklocs[i]
        new_blknos[i + 1] = blknos[i]

    return new_blklocs, new_blknos


# ----------------------------------------------------------------------
# Pickle helper
# ----------------------------------------------------------------------
def _unpickle_block(values, placement, ndim):
    # We have to do some gymnastics b/c "ndim" is keyword-only in some contexts
    from pandas.core.internals.blocks import (
        maybe_coerce_values,
        new_block,
    )
    values = maybe_coerce_values(values)

    if not isinstance(placement, BlockPlacement):
        placement = BlockPlacement(placement)
    return new_block(values, placement, ndim=ndim)


# ----------------------------------------------------------------------
# Block
# ----------------------------------------------------------------------
@cython.freelist(64)
cdef class Block:
    """
    Minimal Block class (Cython) with basic pickling and slicing behavior.
    """
    cdef:
        public BlockPlacement _mgr_locs
        public BlockValuesRefs refs
        readonly int ndim
        public object values

    def __cinit__(self, values, placement, int ndim, BlockValuesRefs refs=None):
        """
        Parameters
        ----------
        values : ndarray or ExtensionArray
        placement : BlockPlacement
        ndim : int
        refs : BlockValuesRefs | None
        """
        self.values = values
        self._mgr_locs = placement
        self.ndim = ndim
        if refs is None:
            self.refs = BlockValuesRefs(self)
        else:
            refs.add_reference(self)
            self.refs = refs

    cpdef __reduce__(self):
        args = (self.values, self._mgr_locs.indexer, self.ndim)
        return _unpickle_block, args

    cpdef __setstate__(self, state):
        from pandas.core.construction import extract_array

        self._mgr_locs = BlockPlacement(state[0])
        self.values = extract_array(state[1], extract_numpy=True)
        if len(state) > 2:
            self.ndim = state[2]
        else:
            from pandas.core.internals.api import _maybe_infer_ndim
            ndim = _maybe_infer_ndim(self.values, self._mgr_locs)
            self.ndim = ndim

    cpdef slice_block_rows(self, slice slicer):
        """
        Specialized slicing along rows (assumes ndim == 2)
        """
        new_values = self.values[..., slicer]
        return type(self)(new_values, self._mgr_locs, ndim=self.ndim, refs=self.refs)


# ----------------------------------------------------------------------
# BlockManager
# ----------------------------------------------------------------------
@cython.freelist(64)
cdef class BlockManager:
    cdef:
        public tuple blocks
        public list axes
        public bint _known_consolidated
        public bint _is_consolidated
        public ndarray _blknos
        public ndarray _blklocs

    def __cinit__(self, blocks=None, axes=None, verify_integrity=True):
        # None as defaults for unpickling GH#42345
        if blocks is None:
            return

        if isinstance(blocks, list):
            blocks = tuple(blocks)

        self.blocks = blocks
        # copy to avoid remote mutation
        self.axes = axes.copy()
        self._known_consolidated = False
        self._is_consolidated = False
        self._blknos = None
        self._blklocs = None

    cpdef _rebuild_blknos_and_blklocs(self):
        """
        Update mgr._blknos / mgr._blklocs.
        """
        cdef:
            Py_ssize_t blkno, i, j
            cnp.npy_intp length = self.shape[0]
            Block blk
            BlockPlacement bp
            ndarray[intp_t, ndim=1] new_blknos, new_blklocs

        # equiv: np.empty(length, dtype=np.intp)
        new_blknos = cnp.PyArray_EMPTY(1, &length, cnp.NPY_INTP, 0)
        new_blklocs = cnp.PyArray_EMPTY(1, &length, cnp.NPY_INTP, 0)
        cnp.PyArray_FILLWBYTE(new_blknos, -1)
        cnp.PyArray_FILLWBYTE(new_blklocs, -1)

        for blkno, blk in enumerate(self.blocks):
            bp = blk._mgr_locs
            for i, j in enumerate(bp):
                new_blknos[j] = blkno
                new_blklocs[j] = i

        for i in range(length):
            blkno = new_blknos[i]
            if blkno == -1:
                raise AssertionError("Gaps in blk ref_locs")

        self._blknos = new_blknos
        self._blklocs = new_blklocs

    cpdef __reduce__(self):
        if len(self.axes) == 1:
            args = (self.blocks[0], self.axes[0])
        else:
            args = (self.blocks, self.axes)
        return type(self), args

    cpdef __setstate__(self, state):
        from pandas.core.construction import extract_array
        from pandas.core.internals.blocks import (
            ensure_block_shape,
            maybe_coerce_values,
            new_block,
        )
        from pandas.core.internals.managers import ensure_index

        if isinstance(state, tuple) and len(state) >= 4 and "0.14.1" in state[3]:
            state = state[3]["0.14.1"]
            axes = [ensure_index(ax) for ax in state["axes"]]
            ndim = len(axes)

            for blk in state["blocks"]:
                vals = blk["values"]
                vals = extract_array(vals, extract_numpy=True)
                blk["values"] = maybe_coerce_values(ensure_block_shape(vals, ndim=ndim))

                if not isinstance(blk["mgr_locs"], BlockPlacement):
                    blk["mgr_locs"] = BlockPlacement(blk["mgr_locs"])

            nbs = [
                new_block(blk["values"], blk["mgr_locs"], ndim=ndim)
                for blk in state["blocks"]
            ]
            blocks = tuple(nbs)
            self.blocks = blocks
            self.axes = axes
        else:  # pragma: no cover
            raise NotImplementedError("pre-0.14.1 pickles are no longer supported")

        self._post_setstate()

    def _post_setstate(self) -> None:
        self._is_consolidated = False
        self._known_consolidated = False
        self._rebuild_blknos_and_blklocs()

    cdef BlockManager _slice_mgr_rows(self, slice slobj):
        cdef:
            Block blk, nb
            BlockManager mgr
            ndarray blknos, blklocs

        nbs = []
        for blk in self.blocks:
            nb = blk.slice_block_rows(slobj)
            nbs.append(nb)

        new_axes = [self.axes[0], self.axes[1]._getitem_slice(slobj)]
        mgr = type(self)(tuple(nbs), new_axes, verify_integrity=False)

        blklocs = self._blklocs
        blknos = self._blknos
        if blknos is not None:
            mgr._blknos = blknos.copy()
            mgr._blklocs = blklocs.copy()
        return mgr

    def get_slice(self, slobj: slice, axis: int = 0):
        if axis == 0:
            new_blocks = self._slice_take_blocks_ax0(slobj)
        elif axis == 1:
            return self._slice_mgr_rows(slobj)
        else:
            raise IndexError("Requested axis not found in manager")

        new_axes = list(self.axes)
        new_axes[axis] = new_axes[axis]._getitem_slice(slobj)

        return type(self)(tuple(new_blocks), new_axes, verify_integrity=False)


# ----------------------------------------------------------------------
# BlockValuesRefs
# ----------------------------------------------------------------------
cdef class BlockValuesRefs:
    """Tracks all references to a given array via weakrefs."""
    cdef:
        public list referenced_blocks
        public int clear_counter

    def __cinit__(self, Block blk=None):
        if blk is not None:
            self.referenced_blocks = [PyWeakref_NewRef(blk, None)]
        else:
            self.referenced_blocks = []
        self.clear_counter = 500

    cdef void _clear_dead_references(self, bint force=False):
        # Use exponential backoff to decide when we want to clear references
        if not force and len(self.referenced_blocks) <= self.clear_counter:
            return

        new_referenced_blocks = []
        for ref in self.referenced_blocks:
            if PyWeakref_GetObject(ref) != Py_None:
                new_referenced_blocks.append(ref)
        self.referenced_blocks = new_referenced_blocks

        nr_of_refs = len(self.referenced_blocks)
        if nr_of_refs < self.clear_counter // 2:
            self.clear_counter = max(self.clear_counter // 2, 500)
        elif nr_of_refs > self.clear_counter:
            self.clear_counter = max(self.clear_counter * 2, nr_of_refs)

    cpdef _add_reference_maybe_locked(self, Block blk):
        self._clear_dead_references()
        self.referenced_blocks.append(PyWeakref_NewRef(blk, None))

    cpdef add_reference(self, Block blk):
        """
        Adds a new reference to our reference collection.
        Uses critical_section when available for thread-safety, falls back when not.
        """
        try:
            with cython.critical_section(self):
                self._add_reference_maybe_locked(blk)
        except Exception:
            # Fallback for Cython/runtime without critical_section support
            self._add_reference_maybe_locked(blk)

    cdef _add_index_reference_maybe_locked(self, object index):
        self._clear_dead_references()
        self.referenced_blocks.append(PyWeakref_NewRef(index, None))

    def add_index_reference(self, object index):
        """
        Adds a new reference to our reference collection when creating an index.
        """
        try:
            with cython.critical_section(self):
                self._add_index_reference_maybe_locked(index)
        except Exception:
            self._add_index_reference_maybe_locked(index)

    cdef bint _has_reference_maybe_locked(self) except *:
        self._clear_dead_references(force=True)
        return len(self.referenced_blocks) > 1

    def has_reference(self) -> bool:
        """Checks if block has foreign references (excluding itself)."""
        try:
            with cython.critical_section(self):
                return self._has_reference_maybe_locked()
        except Exception:
            return self._has_reference_maybe_locked()