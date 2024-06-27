"""Here is defined the Index class."""

import math
import operator
import os
import sys
import tempfile
import warnings

from pathlib import Path
from time import perf_counter as clock
from time import process_time as cpuclock

import numpy as np

from .idxutils import (calc_chunksize, calcoptlevels,
                       get_reduction_level, nextafter, inftype)

from . import indexesextension
from .node import NotLoggedMixin
from .atom import UIntAtom, Atom
from .earray import EArray
from .carray import CArray
from .leaf import Filters
from .indexes import CacheArray, LastRowArray, IndexArray
from .group import Group
from .path import join_path
from .exceptions import PerformanceWarning
from .utils import is_idx, idx2long, lazyattr
from .utilsextension import (nan_aware_gt, nan_aware_ge,
                             nan_aware_lt, nan_aware_le,
                             bisect_left, bisect_right)
from .lrucacheextension import ObjectCache

# default version for INDEX objects
# obversion = "1.0"    # Version of indexes in PyTables 1.x series
# obversion = "2.0"    # Version of indexes in PyTables Pro 2.0 series
obversion = "2.1"     # Version of indexes in PyTables Pro 2.1 and up series,
#                     # including the join 2.3 Std + Pro version

debug = False
# debug = True  # Uncomment this for printing sizes purposes
profile = False
# profile = True  # Uncomment for profiling
if profile:
    from .utils import show_stats

# The default method for sorting
# defsort = "quicksort"
# Changing to mergesort to fix #441
defsort = "mergesort"

# Default policy for automatically updating indexes after a table
# append operation, or automatically reindexing after an
# index-invalidating operation like removing or modifying table rows.
default_auto_index = True
# Keep in sync with ``Table.autoindex`` docstring.

# Default filters used to compress indexes.  This is quite fast and
# compression is pretty good.
# Remember to keep these defaults in sync with the docstrings and UG.
default_index_filters = Filters(complevel=1, complib='zlib',
                                shuffle=True, fletcher32=False)

# Deprecated API
defaultAutoIndex = default_auto_index
defaultIndexFilters = default_index_filters

# The list of types for which an optimised search in cython and C has
# been implemented. Always add here the name of a new optimised type.
opt_search_types = ("int8", "int16", "int32", "int64",
                    "uint8", "uint16", "uint32", "uint64",
                    "float32", "float64")

# The upper limit for uint32 ints
max32 = 2**32


def _table_column_pathname_of_index(indexpathname):
    names = indexpathname.split("/")
    for i, name in enumerate(names):
        if name.startswith('_i_'):
            break
    tablepathname = "/".join(names[:i]) + "/" + name[3:]
    colpathname = "/".join(names[i + 1:])
    return (tablepathname, colpathname)


class Index(NotLoggedMixin, Group, indexesextension.Index):
    """Represents the index of a column in a table.

    This class is used to keep the indexing information for columns in a Table
    dataset (see :ref:`TableClassDescr`). It is actually a descendant of the
    Group class (see :ref:`GroupClassDescr`), with some added functionality. An
    Index is always associated with one and only one column in the table.

    .. note::

        This class is mainly intended for internal use, but some of its
        documented attributes and methods may be interesting for the
        programmer.

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.

        .. versionchanged:: 3.0
           Renamed from *parentNode* to *parentnode*.

    name : str
        The name of this node in its parent group.
    atom : Atom
        An Atom object representing the shape and type of the atomic objects to
        be saved. Only scalar atoms are supported.
    title
        Sets a TITLE attribute of the Index entity.
    kind
        The desired kind for this index.  The 'full' kind specifies a complete
        track of the row position (64-bit), while the 'medium', 'light' or
        'ultralight' kinds only specify in which chunk the row is (using
        32-bit, 16-bit and 8-bit respectively).
    optlevel
        The desired optimization level for this index.
    filters : Filters
        An instance of the Filters class that provides information about the
        desired I/O filters to be applied during the life of this object.
    tmp_dir
        The directory for the temporary files.
    expectedrows
        Represents an user estimate about the number of row slices that will be
        added to the growable dimension in the IndexArray object.
    byteorder
        The byteorder of the index datasets *on-disk*.
    blocksizes
        The four main sizes of the compound blocks in index datasets (a low
        level parameter).

    """

    _c_classid = 'INDEX'

    @property
    def kind(self):
        """The kind of this index."""
        return {1: 'ultralight', 2: 'light',
                4: 'medium', 8: 'full'}[self.indsize]

    @property
    def filters(self):
        """Filter properties for this index - see Filters in
        :ref:`FiltersClassDescr`."""
        return self._v_filters

    @property
    def dirty(self):
        """Whether the index is dirty or not.
        Dirty indexes are out of sync with column data, so they exist but they
        are not usable.
        """

        # If there is no ``DIRTY`` attribute, index should be clean.
        return getattr(self._v_attrs, 'DIRTY', False)

    @dirty.setter
    def dirty(self, dirty):
        wasdirty, isdirty = self.dirty, bool(dirty)
        self._v_attrs.DIRTY = dirty
        # If an *actual* change in dirtiness happens,
        # notify the condition cache by setting or removing a nail.
        conditioncache = self.table._condition_cache
        if not wasdirty and isdirty:
            conditioncache.nail()
        if wasdirty and not isdirty:
            conditioncache.unnail()

    @property
    def column(self):
        """The Column (see :ref:`ColumnClassDescr`) instance for the indexed
        column."""

        tablepath, columnpath = _table_column_pathname_of_index(
            self._v_pathname)
        table = self._v_file._get_node(tablepath)
        column = table.cols._g_col(columnpath)
        return column

    @property
    def table(self):
        """Accessor for the `Table` object of this index."""
        tablepath, columnpath = _table_column_pathname_of_index(
            self._v_pathname)
        table = self._v_file._get_node(tablepath)
        return table

    @property
    def nblockssuperblock(self):
        """The number of blocks in a superblock."""
        return self.superblocksize // self.blocksize

    @property
    def nslicesblock(self):
        """The number of slices in a block."""
        return self.blocksize // self.slicesize

    @property
    def nchunkslice(self):
        """The number of chunks in a slice."""
        return self.slicesize // self.chunksize

    @property
    def nsuperblocks(self):
        """The total number of superblocks in index."""
        # Last row should not be considered as a superblock
        nelements = self.nelements - self.nelementsILR
        nblocks = nelements // self.superblocksize
        if nelements % self.blocksize > 0:
            nblocks += 1
        return nblocks

    @property
    def nblocks(self):
        """The total number of blocks in index."""
        # Last row should not be considered as a block
        nelements = self.nelements - self.nelementsILR
        nblocks = nelements // self.blocksize
        if nelements % self.blocksize > 0:
            nblocks += 1
        return nblocks

    @property
    def nslices(self):
        """The number of complete slices in index."""
        return self.nelements // self.slicesize

    @property
    def nchunks(self):
        """The number of complete chunks in index."""
        return self.nelements // self.chunksize

    @property
    def shape(self):
        """The shape of this index (in slices and elements)."""
        return (self.nrows, self.slicesize)

    @property
    def temp_required(self):
        """Whether a temporary file for indexes is required or not."""
        return (self.indsize > 1 and
                self.optlevel > 0 and
                self.table.nrows > self.slicesize)

    @property
    def want_complete_sort(self):
        """Whether we should try to build a completely sorted index or not."""
        return self.indsize == 8 and self.optlevel == 9

    @property
    def is_csi(self):
        """Whether the index is completely sorted or not.

        .. versionchanged:: 3.0
           The *is_CSI* property has been renamed into *is_csi*.

        """

        if self.nelements == 0:
            # An index with 0 indexed elements is not a CSI one (by definition)
            return False
        if self.indsize < 8:
            # An index that is not full cannot be completely sorted
            return False
        # Try with the 'is_csi' attribute
        if 'is_csi' in self._v_attrs:
            return self._v_attrs.is_csi
        # If not, then compute the overlaps manually
        # (the attribute 'is_csi' will be set there)
        self.compute_overlaps(self, None, False)
        return self.noverlaps == 0

    @lazyattr
    def nrowsinchunk(self):
        """The number of rows that fits in a *table* chunk."""

        return self.table.chunkshape[0]

    @lazyattr
    def lbucket(self):
        """Return the length of a bucket based index type."""

        # Avoid to set a too large lbucket size (mainly useful for tests)
        lbucket = min(self.nrowsinchunk, self.chunksize)
        if self.indsize == 1:
            # For ultra-light, we will never have to keep track of a
            # bucket outside of a slice.
            maxnb = 2**8
            if self.slicesize > maxnb * lbucket:
                lbucket = math.ceil(self.slicesize / maxnb)
        elif self.indsize == 2:
            # For light, we will never have to keep track of a
            # bucket outside of a block.
            maxnb = 2**16
            if self.blocksize > maxnb * lbucket:
                lbucket = math.ceil(self.blocksize / maxnb)
        else:
            # For medium and full indexes there should not be a need to
            # increase lbucket
            pass
        return lbucket

    def __init__(self, parentnode, name,
                 atom=None, title="",
                 kind=None,
                 optlevel=None,
                 filters=None,
                 tmp_dir=None,
                 expectedrows=0,
                 byteorder=None,
                 blocksizes=None,
                 new=True):

        self._v_version = None
        """The object version of this index."""
        self.optlevel = optlevel
        """The optimization level for this index."""
        self.tmp_dir = tmp_dir
        """The directory for the temporary files."""
        self.expectedrows = expectedrows
        """The expected number of items of index arrays."""
        if byteorder in ["little", "big"]:
            self.byteorder = byteorder
        else:
            self.byteorder = sys.byteorder
        """The byteorder of the index datasets."""
        if atom is not None:
            self.dtype = atom.dtype.base
            self.type = atom.type
            """The datatypes to be stored by the sorted index array."""
            # ############## Important note ###########################
            # The datatypes saved as index values are NumPy native
            # types, so we get rid of type metainfo like Time* or Enum*
            # that belongs to HDF5 types (actually, this metainfo is
            # not needed for sorting and looking-up purposes).
            # #########################################################
            indsize = {
                'ultralight': 1, 'light': 2, 'medium': 4, 'full': 8}[kind]
            assert indsize in (1, 2, 4, 8), "indsize should be 1, 2, 4 or 8!"
            self.indsize = indsize
            """The itemsize for the indices part of the index."""

        self.nrows = None
        """The total number of slices in the index."""
        self.nelements = None
        """The number of currently indexed rows for this column."""
        self.blocksizes = blocksizes
        """The four main sizes of the compound blocks (if specified)."""
        self.dirtycache = True
        """Dirty cache (for ranges, bounds & sorted) flag."""
        self.superblocksize = None
        """Size of the superblock for this index."""
        self.blocksize = None
        """Size of the block for this index."""
        self.slicesize = None
        """Size of the slice for this index."""
        self.chunksize = None
        """Size of the chunk for this index."""
        self.tmpfilename = None
        """Filename for temporary bounds."""
        self.opt_search_types = opt_search_types
        """The types for which and optimized search has been implemented."""
        self.noverlaps = -1
        """The number of overlaps in an index.  0 means a completely
        sorted index. -1 means that this number is not computed yet."""
        self.tprof = 0
        """Time counter for benchmarking purposes."""

        from .file import open_file
        self._openFile = open_file
        """The `open_file()` function, to avoid a circular import."""

        super().__init__(parentnode, name, title, new, filters)

    def _g_post_init_hook(self):
        if self._v_new:
            # The version for newly created indexes
            self._v_version = obversion
        super()._g_post_init_hook()

        # Index arrays must only be created for new indexes
        if not self._v_new:
            idxversion = self._v_version
            # Set-up some variables from info on disk and return
            attrs = self._v_attrs
            # Coerce NumPy scalars to Python scalars in order
            # to avoid undesired upcasting operations.
            self.superblocksize = int(attrs.superblocksize)
            self.blocksize = int(attrs.blocksize)
            self.slicesize = int(attrs.slicesize)
            self.chunksize = int(attrs.chunksize)
            self.blocksizes = (self.superblocksize, self.blocksize,
                               self.slicesize, self.chunksize)
            self.optlevel = int(attrs.optlevel)
            sorted = self.sorted
            indices = self.indices
            self.dtype = sorted.atom.dtype
            self.type = sorted.atom.type
            self.indsize = indices.atom.itemsize
            # Some sanity checks for slicesize, chunksize and indsize
            assert self.slicesize == indices.shape[1], "Wrong slicesize"
            assert self.chunksize == indices._v_chunkshape[
                1], "Wrong chunksize"
            assert self.indsize in (1, 2, 4, 8), "Wrong indices itemsize"
            if idxversion > "2.0":
                self.reduction = int(attrs.reduction)
                nelementsSLR = int(self.sortedLR.attrs.nelements)
                nelementsILR = int(self.indicesLR.attrs.nelements)
            else:
                self.reduction = 1
                nelementsILR = self.indicesLR[-1]
                nelementsSLR = nelementsILR
            self.nrows = sorted.nrows
            self.nelements = self.nrows * self.slicesize + nelementsILR
            self.nelementsSLR = nelementsSLR
            self.nelementsILR = nelementsILR
            if nelementsILR > 0:
                self.nrows += 1
            # Get the bounds as a cache (this has to remain here!)
            rchunksize = self.chunksize // self.reduction
            nboundsLR = (nelementsSLR - 1) // rchunksize
            if nboundsLR < 0:
                nboundsLR = 0  # correction for -1 bounds
            nboundsLR += 2  # bounds + begin + end
            # All bounds values (+begin + end) are at the end of sortedLR
            self.bebounds = self.sortedLR[
                nelementsSLR:nelementsSLR + nboundsLR]
            return

        # The index is new. Initialize the values
        self.nrows = 0
        self.nelements = 0
        self.nelementsSLR = 0
        self.nelementsILR = 0

        # The atom
        atom = Atom.from_dtype(self.dtype)

        # The filters
        filters = self.filters

        # Compute the superblocksize, blocksize, slicesize and chunksize values
        # (in case these parameters haven't been passed to the constructor)
        if self.blocksizes is None:
            self.blocksizes = calc_chunksize(
                self.expectedrows, self.optlevel, self.indsize, node=self)
        (self.superblocksize, self.blocksize,
         self.slicesize, self.chunksize) = self.blocksizes
        if debug:
            print("blocksizes:", self.blocksizes)
        # Compute the reduction level
        self.reduction = get_reduction_level(
            self.indsize, self.optlevel, self.slicesize, self.chunksize)
        rchunksize = self.chunksize // self.reduction
        rslicesize = self.slicesize // self.reduction

        # Save them on disk as attributes
        self._v_attrs.superblocksize = np.uint64(self.superblocksize)
        self._v_attrs.blocksize = np.uint64(self.blocksize)
        self._v_attrs.slicesize = np.uint32(self.slicesize)
        self._v_attrs.chunksize = np.uint32(self.chunksize)
        # Save the optlevel as well
        self._v_attrs.optlevel = self.optlevel
        # Save the reduction level
        self._v_attrs.reduction = self.reduction

        # Create the IndexArray for sorted values
        sorted = IndexArray(self, 'sorted', atom, "Sorted Values",
                            filters, self.byteorder)

        # Create the IndexArray for index values
        IndexArray(self, 'indices', UIntAtom(itemsize=self.indsize),
                   "Number of chunk in table", filters, self.byteorder)

        # Create the cache for range values  (1st order cache)
        CacheArray(self, 'ranges', atom, (0, 2), "Range Values", filters,
                   self.expectedrows // self.slicesize,
                   byteorder=self.byteorder)
        # median ranges
        EArray(self, 'mranges', atom, (0,), "Median ranges", filters,
               byteorder=self.byteorder, _log=False)

        # Create the cache for boundary values (2nd order cache)
        nbounds_inslice = (rslicesize - 1) // rchunksize
        CacheArray(self, 'bounds', atom, (0, nbounds_inslice),
                   "Boundary Values", filters, self.nchunks,
                   (1, nbounds_inslice), byteorder=self.byteorder)

        # begin, end & median bounds (only for numerical types)
        EArray(self, 'abounds', atom, (0,), "Start bounds", filters,
               byteorder=self.byteorder, _log=False)
        EArray(self, 'zbounds', atom, (0,), "End bounds", filters,
               byteorder=self.byteorder, _log=False)
        EArray(self, 'mbounds', atom, (0,), "Median bounds", filters,
               byteorder=self.byteorder, _log=False)

        # Create the Array for last (sorted) row values + bounds
        shape = (rslicesize + 2 + nbounds_inslice,)
        sortedLR = LastRowArray(self, 'sortedLR', atom, shape,
                                "Last Row sorted values + bounds",
                                filters, (rchunksize,),
                                byteorder=self.byteorder)

        # Create the Array for the number of chunk in last row
        shape = (self.slicesize,)     # enough for indexes and length
        indicesLR = LastRowArray(self, 'indicesLR',
                                 UIntAtom(itemsize=self.indsize),
                                 shape, "Last Row indices",
                                 filters, (self.chunksize,),
                                 byteorder=self.byteorder)

        # The number of elements in LR will be initialized here
        sortedLR.attrs.nelements = 0
        indicesLR.attrs.nelements = 0

        # All bounds values (+begin + end) are uninitialized in creation time
        self.bebounds = None

        # The starts and lengths initialization
        self.starts = np.empty(shape=self.nrows, dtype=np.int32)
        """Where the values fulfiling conditions starts for every slice."""
        self.lengths = np.empty(shape=self.nrows, dtype=np.int32)
        """Lengths of the values fulfilling conditions for every slice."""

        # Finally, create a temporary file for indexes if needed
        if self.temp_required:
            self.create_temp()

    def initial_append(self, xarr, nrow, reduction):
        """Compute an initial indices arrays for data to be indexed."""

        if profile:
            tref = clock()
        if profile:
            show_stats("Entering initial_append", tref)
        arr = xarr.pop()
        indsize = self.indsize
        slicesize = self.slicesize
        nelementsILR = self.nelementsILR
        if profile:
            show_stats("Before creating idx", tref)
        if indsize == 8:
            idx = np.arange(0, len(arr), dtype="uint64") + nrow * slicesize
        elif indsize == 4:
            # For medium (32-bit) all the rows in tables should be
            # directly reachable.  But as len(arr) < 2**31, we can
            # choose uint32 for representing indices.  In this way, we
            # consume far less memory during the keysort process.  The
            # offset will be added in self.final_idx32() later on.
            #
            # This optimization also prevents the values in LR to
            # participate in the ``swap_chunks`` process, and this is
            # the main reason to not allow the medium indexes to create
            # completely sorted indexes.  However, I don't find this to
            # be a big limitation, as probably fully indexes are much
            # more suitable for producing completely sorted indexes
            # because in this case the indices part is usable for
            # getting the reverse indices of the index, and I forsee
            # this to be a common requirement in many operations (for
            # example, in table sorts).
            #
            # F. Alted 2008-09-15
            idx = np.arange(0, len(arr), dtype="uint32")
        else:
            idx = np.empty(len(arr), "uint%d" % (indsize * 8))
            lbucket = self.lbucket
            # Fill the idx with the bucket indices
            offset = int(lbucket - ((nrow * (slicesize % lbucket)) % lbucket))
            idx[0:offset] = 0
            for i in range(offset, slicesize, lbucket):
                idx[i:i + lbucket] = (i + lbucket - 1) // lbucket
            if indsize == 2:
                # Add a second offset in this case
                # First normalize the number of rows
                offset2 = (nrow % self.nslicesblock) * slicesize // lbucket
                idx += offset2
        # Add the last row at the beginning of arr & idx (if needed)
        if (indsize == 8 and nelementsILR > 0):
            # It is possible that the values in LR are already sorted.
            # Fetch them and override existing values in arr and idx.
            assert len(arr) > nelementsILR
            self.read_slice_lr(self.sortedLR, arr[:nelementsILR])
            self.read_slice_lr(self.indicesLR, idx[:nelementsILR])
        # In-place sorting
        if profile:
            show_stats("Before keysort", tref)
        indexesextension.keysort(arr, idx)
        larr = arr[-1]
        if reduction > 1:
            # It's important to do a copy() here in order to ensure that
            # sorted._append() will receive a contiguous array.
            if profile:
                show_stats("Before reduction", tref)
            reduc = arr[::reduction].copy()
            if profile:
                show_stats("After reduction", tref)
            arr = reduc
            if profile:
                show_stats("After arr <-- reduc", tref)
        # A completely sorted index is not longer possible after an
        # append of an index with already one slice.
        if nrow > 0:
            self._v_attrs.is_csi = False
        if profile:
            show_stats("Exiting initial_append", tref)
        return larr, arr, idx

    def final_idx32(self, idx, offset):
        """Perform final operations in 32-bit indices."""

        if profile:
            tref = clock()
        if profile:
            show_stats("Entering final_idx32", tref)
        # Do an upcast first in order to add the offset.
        idx = idx.astype('uint64')
        idx += offset
        # The next partition is valid up to table sizes of
        # 2**30 * 2**18 = 2**48 bytes, that is, 256 Tera-elements,
        # which should be a safe figure, at least for a while.
        idx //= self.lbucket
        # After the division, we can downsize the indexes to 'uint32'
        idx = idx.astype('uint32')
        if profile:
            show_stats("Exiting final_idx32", tref)
        return idx

    def append(self, xarr, update=False):
        """Append the array to the index objects."""

        if profile:
            tref = clock()
        if profile:
            show_stats("Entering append", tref)
        if not update and self.temp_required:
            where = self.tmp
            # The reduction will take place *after* the optimization process
            reduction = 1
        else:
            where = self
            reduction = self.reduction
        sorted = where.sorted
        indices = where.indices
        ranges = where.ranges
        mranges = where.mranges
        bounds = where.bounds
        mbounds = where.mbounds
        abounds = where.abounds
        zbounds = where.zbounds
        sortedLR = where.sortedLR
        indicesLR = where.indicesLR
        nrows = sorted.nrows  # before sorted.append()
        larr, arr, idx = self.initial_append(xarr, nrows, reduction)
        # Save the sorted array
        sorted.append(arr.reshape(1, arr.size))
        cs = self.chunksize // reduction
        ncs = self.nchunkslice
        # Save ranges & bounds
        ranges.append([[arr[0], larr]])
        bounds.append([arr[cs::cs]])
        abounds.append(arr[0::cs])
        zbounds.append(arr[cs - 1::cs])
        # Compute the medians
        smedian = arr[cs // 2::cs]
        mbounds.append(smedian)
        mranges.append([smedian[ncs // 2]])
        if profile:
            show_stats("Before deleting arr & smedian", tref)
        del arr, smedian   # delete references
        if profile:
            show_stats("After deleting arr & smedian", tref)
        # Now that arr is gone, we can upcast the indices and add the offset
        if self.indsize == 4:
            idx = self.final_idx32(idx, nrows * self.slicesize)
        indices.append(idx.reshape(1, idx.size))
        if profile:
            show_stats("Before deleting idx", tref)
        del idx
        # Update counters after a successful append
        self.nrows = nrows + 1
        self.nelements = self.nrows * self.slicesize
        self.nelementsSLR = 0  # reset the counter of the last row index to 0
        self.nelementsILR = 0  # reset the counter of the last row index to 0
        # The number of elements will be saved as an attribute.
        # This is necessary in case the LR arrays can remember its values
        # after a possible node preemtion/reload.
        sortedLR.attrs.nelements = self.nelementsSLR
        indicesLR.attrs.nelements = self.nelementsILR
        self.dirtycache = True   # the cache is dirty now
        if profile:
            show_stats("Exiting append", tref)

    def append_last_row(self, xarr, update=False):
        """Append the array to the last row index objects."""

        if profile:
            tref = clock()
        if profile:
            show_stats("Entering appendLR", tref)
        # compute the elements in the last row sorted & bounds array
        nrows = self.nslices
        if not update and self.temp_required:
            where = self.tmp
            # The reduction will take place *after* the optimization process
            reduction = 1
        else:
            where = self
            reduction = self.reduction
        indicesLR = where.indicesLR
        sortedLR = where.sortedLR
        larr, arr, idx = self.initial_append(xarr, nrows, reduction)
        nelementsSLR = len(arr)
        nelementsILR = len(idx)
        # Build the cache of bounds
        rchunksize = self.chunksize // reduction
        self.bebounds = np.concatenate((arr[::rchunksize], [larr]))
        # The number of elements will be saved as an attribute
        sortedLR.attrs.nelements = nelementsSLR
        indicesLR.attrs.nelements = nelementsILR
        # Save the number of elements, bounds and sorted values
        # at the end of the sorted array
        offset2 = len(self.bebounds)
        sortedLR[nelementsSLR:nelementsSLR + offset2] = self.bebounds
        sortedLR[:nelementsSLR] = arr
        del arr
        # Now that arr is gone, we can upcast the indices and add the offset
        if self.indsize == 4:
            idx = self.final_idx32(idx, nrows * self.slicesize)
        # Save the reverse index array
        indicesLR[:len(idx)] = idx
        del idx
        # Update counters after a successful append
        self.nrows = nrows + 1
        self.nelements = nrows * self.slicesize + nelementsILR
        self.nelementsILR = nelementsILR
        self.nelementsSLR = nelementsSLR
        self.dirtycache = True   # the cache is dirty now
        if profile:
            show_stats("Exiting appendLR", tref)

    def optimize(self, verbose=False):
        """Optimize an index so as to allow faster searches.

        verbose
            If True, messages about the progress of the
            optimization process are printed out.

        """

        if not self.temp_required:
            return

        if verbose:
            self.verbose = True
        else:
            self.verbose = debug

        # Initialize last_tover and last_nover
        self.last_tover = 0
        self.last_nover = 0

        # Compute the correct optimizations for current optim level
        opts = calcoptlevels(self.nblocks, self.optlevel, self.indsize)
        optmedian, optstarts, optstops, optfull = opts

        if debug:
            print("optvalues:", opts)

        self.create_temp2()
        # Start the optimization process
        while True:
            if optfull:
                for niter in range(optfull):
                    if self.swap('chunks', 'median'):
                        break
                    if self.nblocks > 1:
                        # Swap slices only in the case that we have
                        # several blocks
                        if self.swap('slices', 'median'):
                            break
                        if self.swap('chunks', 'median'):
                            break
                    if self.swap('chunks', 'start'):
                        break
                    if self.swap('chunks', 'stop'):
                        break
            else:
                if optmedian:
                    if self.swap('chunks', 'median'):
                        break
                if optstarts:
                    if self.swap('chunks', 'start'):
                        break
                if optstops:
                    if self.swap('chunks', 'stop'):
                        break
            break  # If we reach this, exit the loop

        # Check if we require a complete sort.  Important: this step
        # should be carried out *after* the optimization process has
        # been completed (this is to guarantee that the complete sort
        # does not take too much memory).
        if self.want_complete_sort:
            if self.noverlaps > 0:
                self.do_complete_sort()
            # Check that we have effectively achieved the complete sort
            if self.noverlaps > 0:
                warnings.warn(
                    "OPSI was not able to achieve a completely sorted index."
                    "  Please report this to the authors.", UserWarning)

        # Close and delete the temporal optimization index file
        self.cleanup_temp()
        return

    def do_complete_sort(self):
        """Bring an already optimized index into a complete sorted state."""

        if self.verbose:
            t1 = clock()
            c1 = cpuclock()
        ss = self.slicesize
        tmp = self.tmp
        ranges = tmp.ranges[:]
        nslices = self.nslices

        nelementsLR = self.nelementsILR
        if nelementsLR > 0:
            # Add the ranges corresponding to the last row
            rangeslr = np.array([self.bebounds[0], self.bebounds[-1]])
            ranges = np.concatenate((ranges, [rangeslr]))
            nslices += 1

        sorted = tmp.sorted
        indices = tmp.indices
        sortedLR = tmp.sortedLR
        indicesLR = tmp.indicesLR
        sremain = np.array([], dtype=self.dtype)
        iremain = np.array([], dtype='u%d' % self.indsize)
        starts = np.zeros(shape=nslices, dtype=np.int_)
        for i in range(nslices):
            # Find the overlapping elements for slice i
            sover = np.array([], dtype=self.dtype)
            iover = np.array([], dtype='u%d' % self.indsize)
            prev_end = ranges[i, 1]
            for j in range(i + 1, nslices):
                stj = starts[j]
                if ((j < self.nslices and stj == ss) or
                        (j == self.nslices and stj == nelementsLR)):
                    # This slice has been already dealt with
                    continue
                if j < self.nslices:
                    assert stj < ss, \
                        "Two slices cannot overlap completely at this stage!"
                    next_beg = sorted[j, stj]
                else:
                    assert stj < nelementsLR, \
                        "Two slices cannot overlap completely at this stage!"
                    next_beg = sortedLR[stj]
                next_end = ranges[j, 1]
                if prev_end > next_end:
                    # Complete overlapping case
                    if j < self.nslices:
                        sover = np.concatenate((sover, sorted[j, stj:]))
                        iover = np.concatenate((iover, indices[j, stj:]))
                        starts[j] = ss
                    else:
                        n = nelementsLR
                        sover = np.concatenate((sover, sortedLR[stj:n]))
                        iover = np.concatenate((iover, indicesLR[stj:n]))
                        starts[j] = nelementsLR
                elif prev_end > next_beg:
                    idx = self.search_item_lt(tmp, prev_end, j, ranges[j], stj)
                    if j < self.nslices:
                        sover = np.concatenate((sover, sorted[j, stj:idx]))
                        iover = np.concatenate((iover, indices[j, stj:idx]))
                    else:
                        sover = np.concatenate((sover, sortedLR[stj:idx]))
                        iover = np.concatenate((iover, indicesLR[stj:idx]))
                    starts[j] = idx
            # Build the extended slices to sort out
            if i < self.nslices:
                ssorted = np.concatenate(
                    (sremain, sorted[i, starts[i]:], sover))
                sindices = np.concatenate(
                    (iremain, indices[i, starts[i]:], iover))
            else:
                ssorted = np.concatenate(
                    (sremain, sortedLR[starts[i]:nelementsLR], sover))
                sindices = np.concatenate(
                    (iremain, indicesLR[starts[i]:nelementsLR], iover))
            # Sort the extended slices
            indexesextension.keysort(ssorted, sindices)
            # Save the first elements of extended slices in the slice i
            if i < self.nslices:
                sorted[i] = ssorted[:ss]
                indices[i] = sindices[:ss]
                # Update caches for this slice
                self.update_caches(i, ssorted[:ss])
                # Save the remaining values in a separate array
                send = len(sover) + len(sremain)
                sremain = ssorted[ss:ss + send]
                iremain = sindices[ss:ss + send]
            else:
                # Still some elements remain for the last row
                n = len(ssorted)
                assert n == nelementsLR
                send = 0
                sortedLR[:n] = ssorted
                indicesLR[:n] = sindices
                # Update the caches for last row
                sortedlr = sortedLR[:nelementsLR]
                bebounds = np.concatenate(
                    (sortedlr[::self.chunksize], [sortedlr[-1]]))
                sortedLR[nelementsLR:nelementsLR + len(bebounds)] = bebounds
                self.bebounds = bebounds

        # Verify that we have dealt with all the remaining values
        assert send == 0

        # Compute the overlaps in order to verify that we have achieved
        # a complete sort.  This has to be executed always (and not only
        # in verbose mode!).
        self.compute_overlaps(self.tmp, "do_complete_sort()", self.verbose)
        if self.verbose:
            print(f"time: {clock() - t1:.4f}. clock: {cpuclock() - c1:.4f}")

    def swap(self, what, mode=None):
        """Swap chunks or slices using a certain bounds reference."""

        # Thresholds for avoiding continuing the optimization
        # thnover = 4 * self.slicesize  # minimum number of overlapping
        #                               # elements
        thnover = 40
        thmult = 0.1      # minimum ratio of multiplicity (a 10%)
        thtover = 0.01    # minimum overlaping index for slices (a 1%)

        if self.verbose:
            t1 = clock()
            c1 = cpuclock()
        if what == "chunks":
            self.swap_chunks(mode)
        elif what == "slices":
            self.swap_slices(mode)
        if mode:
            message = f"swap_{what}({mode})"
        else:
            message = f"swap_{what}"
        (nover, mult, tover) = self.compute_overlaps(
            self.tmp, message, self.verbose)
        rmult = len(mult.nonzero()[0]) / len(mult)
        if self.verbose:
            print(f"time: {clock() - t1:.4f}. clock: {cpuclock() - c1:.4f}")
        # Check that entropy is actually decreasing
        if what == "chunks" and self.last_tover > 0 and self.last_nover > 0:
            tover_var = (self.last_tover - tover) / self.last_tover
            nover_var = (self.last_nover - nover) / self.last_nover
            if tover_var < 0.05 and nover_var < 0.05:
                # Less than a 5% of improvement is too few
                return True
        self.last_tover = tover
        self.last_nover = nover
        # Check if some threshold has met
        if nover < thnover:
            return True
        if rmult < thmult:
            return True
        # Additional check for the overlap ratio
        if 0 <= tover < thtover:
            return True
        return False

    def create_temp(self):
        """Create some temporary objects for slice sorting purposes."""

        # The index will be dirty during the index optimization process
        self.dirty = True
        # Build the name of the temporary file
        fd, self.tmpfilename = tempfile.mkstemp(
            ".tmp", "pytables-", self.tmp_dir)
        # Close the file descriptor so as to avoid leaks
        os.close(fd)
        # Create the proper PyTables file
        self.tmpfile = self._openFile(self.tmpfilename, "w")
        self.tmp = tmp = self.tmpfile.root
        cs = self.chunksize
        ss = self.slicesize
        filters = self.filters
        # temporary sorted & indices arrays
        shape = (0, ss)
        atom = Atom.from_dtype(self.dtype)
        EArray(tmp, 'sorted', atom, shape,
               "Temporary sorted", filters, chunkshape=(1, cs))
        EArray(tmp, 'indices', UIntAtom(itemsize=self.indsize), shape,
               "Temporary indices", filters, chunkshape=(1, cs))
        # temporary bounds
        nbounds_inslice = (ss - 1) // cs
        shape = (0, nbounds_inslice)
        EArray(tmp, 'bounds', atom, shape, "Temp chunk bounds",
               filters, chunkshape=(cs, nbounds_inslice))
        shape = (0,)
        EArray(tmp, 'abounds', atom, shape, "Temp start bounds",
               filters, chunkshape=(cs,))
        EArray(tmp, 'zbounds', atom, shape, "Temp end bounds",
               filters, chunkshape=(cs,))
        EArray(tmp, 'mbounds', atom, shape, "Median bounds",
               filters, chunkshape=(cs,))
        # temporary ranges
        EArray(tmp, 'ranges', atom, (0, 2),
               "Temporary range values", filters, chunkshape=(cs, 2))
        EArray(tmp, 'mranges', atom, (0,),
               "Median ranges", filters, chunkshape=(cs,))
        # temporary last row (sorted)
        shape = (ss + 2 + nbounds_inslice,)
        CArray(tmp, 'sortedLR', atom, shape,
               "Temp Last Row sorted values + bounds",
               filters, chunkshape=(cs,))
        # temporary last row (indices)
        shape = (ss,)
        CArray(tmp, 'indicesLR',
               UIntAtom(itemsize=self.indsize),
               shape, "Temp Last Row indices",
               filters, chunkshape=(cs,))

    def create_temp2(self):
        """Create some temporary objects for slice sorting purposes."""

        # The algorithms for doing the swap can be optimized so that
        # one should be necessary to create temporaries for keeping just
        # the contents of a single superblock.
        # F. Alted 2007-01-03
        cs = self.chunksize
        ss = self.slicesize
        filters = self.filters
        # temporary sorted & indices arrays
        shape = (self.nslices, ss)
        atom = Atom.from_dtype(self.dtype)
        tmp = self.tmp
        CArray(tmp, 'sorted2', atom, shape,
               "Temporary sorted 2", filters, chunkshape=(1, cs))
        CArray(tmp, 'indices2', UIntAtom(itemsize=self.indsize), shape,
               "Temporary indices 2", filters, chunkshape=(1, cs))
        # temporary bounds
        nbounds_inslice = (ss - 1) // cs
        shape = (self.nslices, nbounds_inslice)
        CArray(tmp, 'bounds2', atom, shape, "Temp chunk bounds 2",
               filters, chunkshape=(cs, nbounds_inslice))
        shape = (self.nchunks,)
        CArray(tmp, 'abounds2', atom, shape, "Temp start bounds 2",
               filters, chunkshape=(cs,))
        CArray(tmp, 'zbounds2', atom, shape, "Temp end bounds 2",
               filters, chunkshape=(cs,))
        CArray(tmp, 'mbounds2', atom, shape, "Median bounds 2",
               filters, chunkshape=(cs,))
        # temporary ranges
        CArray(tmp, 'ranges2', atom, (self.nslices, 2),
               "Temporary range values 2", filters, chunkshape=(cs, 2))
        CArray(tmp, 'mranges2', atom, (self.nslices,),
               "Median ranges 2", filters, chunkshape=(cs,))

    def cleanup_temp(self):
        """Copy the data and delete the temporaries for sorting purposes."""

        if self.verbose:
            print("Copying temporary data...")
        # tmp -> index
        reduction = self.reduction
        cs = self.chunksize // reduction
        ncs = self.nchunkslice
        tmp = self.tmp
        for i in range(self.nslices):
            # Copy sorted & indices slices
            sorted = tmp.sorted[i][::reduction].copy()
            self.sorted.append(sorted.reshape(1, sorted.size))
            # Compute ranges
            self.ranges.append([[sorted[0], sorted[-1]]])
            # Compute chunk bounds
            self.bounds.append([sorted[cs::cs]])
            # Compute start, stop & median bounds and ranges
            self.abounds.append(sorted[0::cs])
            self.zbounds.append(sorted[cs - 1::cs])
            smedian = sorted[cs // 2::cs]
            self.mbounds.append(smedian)
            self.mranges.append([smedian[ncs // 2]])
            del sorted, smedian   # delete references
            # Now that sorted is gone, we can copy the indices
            indices = tmp.indices[i]
            self.indices.append(indices.reshape(1, indices.size))

        # Now it is the last row turn (if needed)
        if self.nelementsSLR > 0:
            # First, the sorted values
            sortedLR = self.sortedLR
            indicesLR = self.indicesLR
            nelementsLR = self.nelementsILR
            sortedlr = tmp.sortedLR[:nelementsLR][::reduction].copy()
            nelementsSLR = len(sortedlr)
            sortedLR[:nelementsSLR] = sortedlr
            # Now, the bounds
            self.bebounds = np.concatenate((sortedlr[::cs], [sortedlr[-1]]))
            offset2 = len(self.bebounds)
            sortedLR[nelementsSLR:nelementsSLR + offset2] = self.bebounds
            # Finally, the indices
            indicesLR[:] = tmp.indicesLR[:]
            # Update the number of (reduced) sorted elements
            self.nelementsSLR = nelementsSLR
        # The number of elements will be saved as an attribute
        self.sortedLR.attrs.nelements = self.nelementsSLR
        self.indicesLR.attrs.nelements = self.nelementsILR

        if self.verbose:
            print("Deleting temporaries...")
        self.tmp = None
        self.tmpfile.close()
        Path(self.tmpfilename).unlink()
        self.tmpfilename = None

        # The optimization process has finished, and the index is ok now
        self.dirty = False
        # ...but the memory data cache is dirty now
        self.dirtycache = True

    def get_neworder(self, neworder, src_disk, tmp_disk,
                     lastrow, nslices, offset, dtype):
        """Get sorted & indices values in new order."""

        cs = self.chunksize
        ncs = ncs2 = self.nchunkslice
        self_nslices = self.nslices
        tmp = np.empty(shape=self.slicesize, dtype=dtype)
        for i in range(nslices):
            ns = offset + i
            if ns == self_nslices:
                # The number of complete chunks in the last row
                ncs2 = self.nelementsILR // cs
            # Get slices in new order
            for j in range(ncs2):
                idx = neworder[i * ncs + j]
                ins = idx // ncs
                inc = (idx - ins * ncs) * cs
                ins += offset
                nc = j * cs
                if ins == self_nslices:
                    tmp[nc:nc + cs] = lastrow[inc:inc + cs]
                else:
                    tmp[nc:nc + cs] = src_disk[ins, inc:inc + cs]
            if ns == self_nslices:
                # The number of complete chunks in the last row
                lastrow[:ncs2 * cs] = tmp[:ncs2 * cs]
                # The elements in the last chunk of the last row will
                # participate in the global reordering later on, during
                # the phase of sorting of *two* slices at a time
                # (including the last row slice, see
                # self.reorder_slices()).  The caches for last row will
                # be updated in self.reorder_slices() too.
                # F. Altet 2008-08-25
            else:
                tmp_disk[ns] = tmp

    def swap_chunks(self, mode="median"):
        """Swap & reorder the different chunks in a block."""

        boundsnames = {
            'start': 'abounds', 'stop': 'zbounds', 'median': 'mbounds'}
        tmp = self.tmp
        sorted = tmp.sorted
        indices = tmp.indices
        tmp_sorted = tmp.sorted2
        tmp_indices = tmp.indices2
        sortedLR = tmp.sortedLR
        indicesLR = tmp.indicesLR
        cs = self.chunksize
        ncs = self.nchunkslice
        nsb = self.nslicesblock
        ncb = ncs * nsb
        ncb2 = ncb
        boundsobj = tmp._f_get_child(boundsnames[mode])
        can_cross_bbounds = (self.indsize == 8 and self.nelementsILR > 0)
        for nblock in range(self.nblocks):
            # Protection for last block having less chunks than ncb
            remainingchunks = self.nchunks - nblock * ncb
            if remainingchunks < ncb:
                ncb2 = remainingchunks
            if ncb2 <= 1:
                # if only zero or one chunks remains we are done
                break
            nslices = ncb2 // ncs
            bounds = boundsobj[nblock * ncb:nblock * ncb + ncb2]
            # Do this only if lastrow elements can cross block boundaries
            if (nblock == self.nblocks - 1 and  # last block
                    can_cross_bbounds):
                nslices += 1
                ul = self.nelementsILR // cs
                bounds = np.concatenate((bounds, self.bebounds[:ul]))
            sbounds_idx = bounds.argsort(kind=defsort)
            offset = int(nblock * nsb)
            # Swap sorted and indices following the new order
            self.get_neworder(sbounds_idx, sorted, tmp_sorted, sortedLR,
                              nslices, offset, self.dtype)
            self.get_neworder(sbounds_idx, indices, tmp_indices, indicesLR,
                              nslices, offset, 'u%d' % self.indsize)
        # Reorder completely the index at slice level
        self.reorder_slices(tmp=True)

    def read_slice(self, where, nslice, buffer, start=0):
        """Read a slice from the `where` dataset and put it in `buffer`."""

        # Create the buffers for specifying the coordinates
        self.startl = np.array([nslice, start], np.uint64)
        self.stopl = np.array([nslice + 1, start + buffer.size], np.uint64)
        self.stepl = np.ones(shape=2, dtype=np.uint64)
        where._g_read_slice(self.startl, self.stopl, self.stepl, buffer)

    def write_slice(self, where, nslice, buffer, start=0):
        """Write a `slice` to the `where` dataset with the `buffer` data."""

        self.startl = np.array([nslice, start], np.uint64)
        self.stopl = np.array([nslice + 1, start + buffer.size], np.uint64)
        self.stepl = np.ones(shape=2, dtype=np.uint64)
        countl = self.stopl - self.startl   # (1, self.slicesize)
        where._g_write_slice(self.startl, self.stepl, countl, buffer)

    # Read version for LastRow
    def read_slice_lr(self, where, buffer, start=0):
        """Read a slice from the `where` dataset and put it in `buffer`."""

        startl = np.array([start], dtype=np.uint64)
        stopl = np.array([start + buffer.size], dtype=np.uint64)
        stepl = np.array([1], dtype=np.uint64)
        where._g_read_slice(startl, stopl, stepl, buffer)

    # Write version for LastRow
    def write_sliceLR(self, where, buffer, start=0):
        """Write a slice from the `where` dataset with the `buffer` data."""

        startl = np.array([start], dtype=np.uint64)
        countl = np.array([start + buffer.size], dtype=np.uint64)
        stepl = np.array([1], dtype=np.uint64)
        where._g_write_slice(startl, stepl, countl, buffer)

    def reorder_slice(self, nslice, sorted, indices, ssorted, sindices,
                      tmp_sorted, tmp_indices):
        """Copy & reorder the slice in source to final destination."""

        ss = self.slicesize
        # Load the second part in buffers
        self.read_slice(tmp_sorted, nslice, ssorted[ss:])
        self.read_slice(tmp_indices, nslice, sindices[ss:])
        indexesextension.keysort(ssorted, sindices)
        # Write the first part of the buffers to the regular leaves
        self.write_slice(sorted, nslice - 1, ssorted[:ss])
        self.write_slice(indices, nslice - 1, sindices[:ss])
        # Update caches
        self.update_caches(nslice - 1, ssorted[:ss])
        # Shift the slice in the end to the beginning
        ssorted[:ss] = ssorted[ss:]
        sindices[:ss] = sindices[ss:]

    def update_caches(self, nslice, ssorted):
        """Update the caches for faster lookups."""

        cs = self.chunksize
        ncs = self.nchunkslice
        tmp = self.tmp
        # update first & second cache bounds (ranges & bounds)
        tmp.ranges[nslice] = ssorted[[0, -1]]
        tmp.bounds[nslice] = ssorted[cs::cs]
        # update start & stop bounds
        tmp.abounds[nslice * ncs:(nslice + 1) * ncs] = ssorted[0::cs]
        tmp.zbounds[nslice * ncs:(nslice + 1) * ncs] = ssorted[cs - 1::cs]
        # update median bounds
        smedian = ssorted[cs // 2::cs]
        tmp.mbounds[nslice * ncs:(nslice + 1) * ncs] = smedian
        tmp.mranges[nslice] = smedian[ncs // 2]

    def reorder_slices(self, tmp):
        """Reorder completely the index at slice level.

        This method has to maintain the locality of elements in the
        ambit of ``blocks``, i.e. an element of a ``block`` cannot be
        sent to another ``block`` during this reordering.  This is
        *critical* for ``light`` indexes to be able to use this.

        This version of reorder_slices is optimized in that *two*
        complete slices are taken at a time (including the last row
        slice) so as to sort them.  Then, each new slice that is read is
        put at the end of this two-slice buffer, while the previous one
        is moved to the beginning of the buffer.  This is in order to
        better reduce the entropy of the regular part (i.e. all except
        the last row) of the index.

        A secondary effect of this is that it takes at least *twice* of
        memory than a previous version of reorder_slices() that only
        reorders on a slice-by-slice basis.  However, as this is more
        efficient than the old version, one can configure the slicesize
        to be smaller, so the memory consumption is barely similar.

        """

        tmp = self.tmp
        sorted = tmp.sorted
        indices = tmp.indices
        if tmp:
            tmp_sorted = tmp.sorted2
            tmp_indices = tmp.indices2
        else:
            tmp_sorted = tmp.sorted
            tmp_indices = tmp.indices
        cs = self.chunksize
        ss = self.slicesize
        nsb = self.blocksize // self.slicesize
        nslices = self.nslices
        nblocks = self.nblocks
        nelementsLR = self.nelementsILR
        # Create the buffer for reordering 2 slices at a time
        ssorted = np.empty(shape=ss * 2, dtype=self.dtype)
        sindices = np.empty(shape=ss * 2, dtype=np.dtype('u%d' % self.indsize))

        if self.indsize == 8:
            # Bootstrap the process for reordering
            # Read the first slice in buffers
            self.read_slice(tmp_sorted, 0, ssorted[:ss])
            self.read_slice(tmp_indices, 0, sindices[:ss])

            nslice = 0   # Just in case the loop behind executes nothing
            # Loop over the remainding slices in block
            for nslice in range(1, sorted.nrows):
                self.reorder_slice(nslice, sorted, indices,
                                   ssorted, sindices,
                                   tmp_sorted, tmp_indices)

            # End the process (enrolling the lastrow if necessary)
            if nelementsLR > 0:
                sortedLR = self.tmp.sortedLR
                indicesLR = self.tmp.indicesLR
                # Shrink the ssorted and sindices arrays to the minimum
                ssorted2 = ssorted[:ss + nelementsLR]
                sortedlr = ssorted2[ss:]
                sindices2 = sindices[:ss + nelementsLR]
                indiceslr = sindices2[ss:]
                # Read the last row info in the second part of the buffer
                self.read_slice_lr(sortedLR, sortedlr)
                self.read_slice_lr(indicesLR, indiceslr)
                indexesextension.keysort(ssorted2, sindices2)
                # Write the second part of the buffers to the lastrow indices
                self.write_sliceLR(sortedLR, sortedlr)
                self.write_sliceLR(indicesLR, indiceslr)
                # Update the caches for last row
                bebounds = np.concatenate((sortedlr[::cs], [sortedlr[-1]]))
                sortedLR[nelementsLR:nelementsLR + len(bebounds)] = bebounds
                self.bebounds = bebounds
            # Write the first part of the buffers to the regular leaves
            self.write_slice(sorted, nslice, ssorted[:ss])
            self.write_slice(indices, nslice, sindices[:ss])
            # Update caches for this slice
            self.update_caches(nslice, ssorted[:ss])
        else:
            # Iterate over each block.  No data should cross block
            # boundaries to avoid adressing problems with short indices.
            for nb in range(nblocks):
                # Bootstrap the process for reordering
                # Read the first slice in buffers
                nrow = nb * nsb
                self.read_slice(tmp_sorted, nrow, ssorted[:ss])
                self.read_slice(tmp_indices, nrow, sindices[:ss])

                # Loop over the remainding slices in block
                lrb = nrow + nsb
                if lrb > nslices:
                    lrb = nslices
                nslice = nrow   # Just in case the loop behind executes nothing
                for nslice in range(nrow + 1, lrb):
                    self.reorder_slice(nslice, sorted, indices,
                                       ssorted, sindices,
                                       tmp_sorted, tmp_indices)

                # Write the first part of the buffers to the regular leaves
                self.write_slice(sorted, nslice, ssorted[:ss])
                self.write_slice(indices, nslice, sindices[:ss])
                # Update caches for this slice
                self.update_caches(nslice, ssorted[:ss])

    def swap_slices(self, mode="median"):
        """Swap slices in a superblock."""

        tmp = self.tmp
        sorted = tmp.sorted
        indices = tmp.indices
        tmp_sorted = tmp.sorted2
        tmp_indices = tmp.indices2
        ncs = self.nchunkslice
        nss = self.superblocksize // self.slicesize
        nss2 = nss
        for sblock in range(self.nsuperblocks):
            # Protection for last superblock having less slices than nss
            remainingslices = self.nslices - sblock * nss
            if remainingslices < nss:
                nss2 = remainingslices
            if nss2 <= 1:
                break
            if mode == "start":
                ranges = tmp.ranges[sblock * nss:sblock * nss + nss2, 0]
            elif mode == "stop":
                ranges = tmp.ranges[sblock * nss:sblock * nss + nss2, 1]
            elif mode == "median":
                ranges = tmp.mranges[sblock * nss:sblock * nss + nss2]
            sranges_idx = ranges.argsort(kind=defsort)
            # Don't swap the superblock at all if one doesn't need to
            ndiff = (sranges_idx != np.arange(nss2)).sum() / 2
            if ndiff * 50 < nss2:
                # The number of slices to rearrange is less than 2.5%,
                # so skip the reordering of this superblock
                # (too expensive for such a little improvement)
                if self.verbose:
                    print("skipping reordering of superblock ->", sblock)
                continue
            ns = sblock * nss2
            # Swap sorted and indices slices following the new order
            for i in range(nss2):
                idx = sranges_idx[i]
                # Swap sorted & indices slices
                oi = ns + i
                oidx = ns + idx
                tmp_sorted[oi] = sorted[oidx]
                tmp_indices[oi] = indices[oidx]
                # Swap start, stop & median ranges
                tmp.ranges2[oi] = tmp.ranges[oidx]
                tmp.mranges2[oi] = tmp.mranges[oidx]
                # Swap chunk bounds
                tmp.bounds2[oi] = tmp.bounds[oidx]
                # Swap start, stop & median bounds
                j = oi * ncs
                jn = (oi + 1) * ncs
                xj = oidx * ncs
                xjn = (oidx + 1) * ncs
                tmp.abounds2[j:jn] = tmp.abounds[xj:xjn]
                tmp.zbounds2[j:jn] = tmp.zbounds[xj:xjn]
                tmp.mbounds2[j:jn] = tmp.mbounds[xj:xjn]
            # tmp -> originals
            for i in range(nss2):
                # Copy sorted & indices slices
                oi = ns + i
                sorted[oi] = tmp_sorted[oi]
                indices[oi] = tmp_indices[oi]
                # Copy start, stop & median ranges
                tmp.ranges[oi] = tmp.ranges2[oi]
                tmp.mranges[oi] = tmp.mranges2[oi]
                # Copy chunk bounds
                tmp.bounds[oi] = tmp.bounds2[oi]
                # Copy start, stop & median bounds
                j = oi * ncs
                jn = (oi + 1) * ncs
                tmp.abounds[j:jn] = tmp.abounds2[j:jn]
                tmp.zbounds[j:jn] = tmp.zbounds2[j:jn]
                tmp.mbounds[j:jn] = tmp.mbounds2[j:jn]

    def search_item_lt(self, where, item, nslice, limits, start=0):
        """Search a single item in a specific sorted slice."""

        # This method will only works under the assumtion that item
        # *is to be found* in the nslice.
        assert nan_aware_lt(limits[0], item) and nan_aware_le(item, limits[1])
        cs = self.chunksize
        ss = self.slicesize
        nelementsLR = self.nelementsILR
        bstart = start // cs

        # Find the chunk
        if nslice < self.nslices:
            nchunk = bisect_left(where.bounds[nslice], item, bstart)
        else:
            # We need to subtract 1 chunk here because bebounds
            # has a leading value
            nchunk = bisect_left(self.bebounds, item, bstart) - 1
        assert nchunk >= 0

        # Find the element in chunk
        pos = nchunk * cs
        if nslice < self.nslices:
            pos += bisect_left(where.sorted[nslice, pos:pos + cs], item)
            assert pos <= ss
        else:
            end = pos + cs
            if end > nelementsLR:
                end = nelementsLR
            pos += bisect_left(self.sortedLR[pos:end], item)
            assert pos <= nelementsLR
        assert pos > 0
        return pos

    def compute_overlaps_finegrain(self, where, message, verbose):
        """Compute some statistics about overlaping of slices in index.

        It returns the following info:

        noverlaps : int
            The total number of elements that overlaps in index.
        multiplicity : array of int
            The number of times that a concrete slice overlaps with any other.
        toverlap : float
            An ovelap index: the sum of the values in segment slices that
            overlaps divided by the entire range of values.  This index is only
            computed for numerical types.

        """

        ss = self.slicesize
        ranges = where.ranges[:]
        sorted = where.sorted
        sortedLR = where.sortedLR
        nslices = self.nslices
        nelementsLR = self.nelementsILR
        if nelementsLR > 0:
            # Add the ranges corresponding to the last row
            rangeslr = np.array([self.bebounds[0], self.bebounds[-1]])
            ranges = np.concatenate((ranges, [rangeslr]))
            nslices += 1
        soverlap = 0
        toverlap = -1
        multiplicity = np.zeros(shape=nslices, dtype="int_")
        overlaps = multiplicity.copy()
        starts = multiplicity.copy()
        for i in range(nslices):
            prev_end = ranges[i, 1]
            for j in range(i + 1, nslices):
                stj = starts[j]
                assert stj <= ss
                if stj == ss:
                    # This slice has already been counted
                    continue
                if j < self.nslices:
                    next_beg = sorted[j, stj]
                else:
                    next_beg = sortedLR[stj]
                next_end = ranges[j, 1]
                if prev_end > next_end:
                    # Complete overlapping case
                    multiplicity[j - i] += 1
                    if j < self.nslices:
                        overlaps[i] += ss - stj
                        starts[j] = ss   # a sentinel
                    else:
                        overlaps[i] += nelementsLR - stj
                        starts[j] = nelementsLR   # a sentinel
                elif prev_end > next_beg:
                    multiplicity[j - i] += 1
                    idx = self.search_item_lt(
                        where, prev_end, j, ranges[j], stj)
                    nelem = idx - stj
                    overlaps[i] += nelem
                    starts[j] = idx
                    if self.type != "string":
                        # Convert ranges into floats in order to allow
                        # doing operations with them without overflows
                        soverlap += float(ranges[i, 1]) - float(ranges[j, 0])

        # Return the overlap as the ratio between overlaps and entire range
        if self.type != "string":
            erange = float(ranges[-1, 1]) - float(ranges[0, 0])
            # Check that there is an effective range of values
            # Beware, erange can be negative in situations where
            # the values are suffering overflow. This can happen
            # specially on big signed integer values (on overflows,
            # the end value will become negative!).
            # Also, there is no way to compute overlap ratios for
            # non-numerical types. So, be careful and always check
            # that toverlap has a positive value (it must have been
            # initialized to -1. before) before using it.
            # F. Alted 2007-01-19
            if erange > 0:
                toverlap = soverlap / erange
        if verbose and message != "init":
            print("toverlap (%s):" % message, toverlap)
            print("multiplicity:\n", multiplicity, multiplicity.sum())
            print("overlaps:\n", overlaps, overlaps.sum())
        noverlaps = overlaps.sum()
        # For full indexes, set the 'is_csi' flag
        if self.indsize == 8 and self._v_file._iswritable():
            self._v_attrs.is_csi = (noverlaps == 0)
        # Save the number of overlaps for future references
        self.noverlaps = noverlaps
        return (noverlaps, multiplicity, toverlap)

    def compute_overlaps(self, where, message, verbose):
        """Compute some statistics about overlaping of slices in index.

        It returns the following info:

        noverlaps : int
            The total number of slices that overlaps in index.
        multiplicity : array of int
            The number of times that a concrete slice overlaps with any other.
        toverlap : float
            An ovelap index: the sum of the values in segment slices that
            overlaps divided by the entire range of values.  This index is only
            computed for numerical types.

        """

        ranges = where.ranges[:]
        nslices = self.nslices
        if self.nelementsILR > 0:
            # Add the ranges corresponding to the last row
            rangeslr = np.array([self.bebounds[0], self.bebounds[-1]])
            ranges = np.concatenate((ranges, [rangeslr]))
            nslices += 1
        noverlaps = 0
        soverlap = 0
        toverlap = -1
        multiplicity = np.zeros(shape=nslices, dtype="int_")
        for i in range(nslices):
            for j in range(i + 1, nslices):
                if ranges[i, 1] > ranges[j, 0]:
                    noverlaps += 1
                    multiplicity[j - i] += 1
                    if self.type != "string":
                        # Convert ranges into floats in order to allow
                        # doing operations with them without overflows
                        soverlap += float(ranges[i, 1]) - float(ranges[j, 0])

        # Return the overlap as the ratio between overlaps and entire range
        if self.type != "string":
            erange = float(ranges[-1, 1]) - float(ranges[0, 0])
            # Check that there is an effective range of values
            # Beware, erange can be negative in situations where
            # the values are suffering overflow. This can happen
            # specially on big signed integer values (on overflows,
            # the end value will become negative!).
            # Also, there is no way to compute overlap ratios for
            # non-numerical types. So, be careful and always check
            # that toverlap has a positive value (it must have been
            # initialized to -1. before) before using it.
            # F. Altet 2007-01-19
            if erange > 0:
                toverlap = soverlap / erange
        if verbose:
            print("overlaps (%s):" % message, noverlaps, toverlap)
            print(multiplicity)
        # For full indexes, set the 'is_csi' flag
        if self.indsize == 8 and self._v_file._iswritable():
            self._v_attrs.is_csi = (noverlaps == 0)
        # Save the number of overlaps for future references
        self.noverlaps = noverlaps
        return (noverlaps, multiplicity, toverlap)

    def read_sorted_indices(self, what, start, stop, step):
        """Return the sorted or indices values in the specified range."""
        (start, stop, step) = self._process_range(start, stop, step)
        if start >= stop:
            return np.empty(0, self.dtype)
        # Correction for negative values of step (reverse indices)
        if step < 0:
            tmp = start
            start = self.nelements - stop
            stop = self.nelements - tmp
        if what == "sorted":
            values = self.sorted
            valuesLR = self.sortedLR
            buffer_ = np.empty(stop - start, dtype=self.dtype)
        else:
            values = self.indices
            valuesLR = self.indicesLR
            buffer_ = np.empty(stop - start, dtype="u%d" % self.indsize)
        ss = self.slicesize
        nrow_start = start // ss
        istart = start % ss
        nrow_stop = stop // ss
        tlen = stop - start
        bstart = 0
        ilen = 0
        for nrow in range(nrow_start, nrow_stop + 1):
            blen = ss - istart
            if ilen + blen > tlen:
                blen = tlen - ilen
            if blen <= 0:
                break
            if nrow < self.nslices:
                self.read_slice(
                    values, nrow, buffer_[bstart:bstart + blen], istart)
            else:
                self.read_slice_lr(
                    valuesLR, buffer_[bstart:bstart + blen], istart)
            istart = 0
            bstart += blen
            ilen += blen
        return buffer_[::step]

    def read_sorted(self, start=None, stop=None, step=None):
        """Return the sorted values of index in the specified range.

        The meaning of the start, stop and step arguments is the same as in
        :meth:`Table.read_sorted`.

        """

        return self.read_sorted_indices('sorted', start, stop, step)

    def read_indices(self, start=None, stop=None, step=None):
        """Return the indices values of index in the specified range.

        The meaning of the start, stop and step arguments is the same as in
        :meth:`Table.read_sorted`.

        """

        return self.read_sorted_indices('indices', start, stop, step)

    def _process_range(self, start, stop, step):
        """Get a range specifc for the index usage."""

        if start is not None and stop is None:
            # Special case for the behaviour of PyTables iterators
            stop = idx2long(start + 1)
        if start is None:
            start = 0
        else:
            start = idx2long(start)
        if stop is None:
            stop = idx2long(self.nelements)
        else:
            stop = idx2long(stop)
        if step is None:
            step = 1
        else:
            step = idx2long(step)
        return (start, stop, step)

    def __getitem__(self, key):
        """Return the indices values of index in the specified range.

        If key argument is an integer, the corresponding index is returned.  If
        key is a slice, the range of indices determined by it is returned.  A
        negative value of step in slice is supported, meaning that the results
        will be returned in reverse order.

        This method is equivalent to :meth:`Index.read_indices`.

        """

        if is_idx(key):
            key = operator.index(key)

            if key < 0:
                # To support negative values
                key += self.nelements
            return self.read_indices(key, key + 1, 1)[0]
        elif isinstance(key, slice):
            return self.read_indices(key.start, key.stop, key.step)

    def __len__(self):
        return self.nelements

    def restorecache(self):
        """Clean the limits cache and resize starts and lengths arrays"""

        params = self._v_file.params
        # The sorted IndexArray is absolutely required to be in memory
        # at the same time than the Index instance, so create a strong
        # reference to it.  We are not introducing leaks because the
        # strong reference will disappear when this Index instance is
        # to be closed.
        self._sorted = self.sorted
        self._sorted.boundscache = ObjectCache(params['BOUNDS_MAX_SLOTS'],
                                               params['BOUNDS_MAX_SIZE'],
                                               'non-opt types bounds')
        self.sorted.boundscache = ObjectCache(params['BOUNDS_MAX_SLOTS'],
                                              params['BOUNDS_MAX_SIZE'],
                                              'non-opt types bounds')
        """A cache for the bounds (2nd hash) data. Only used for
        non-optimized types searches."""
        self.limboundscache = ObjectCache(params['LIMBOUNDS_MAX_SLOTS'],
                                          params['LIMBOUNDS_MAX_SIZE'],
                                          'bounding limits')
        """A cache for bounding limits."""
        self.sortedLRcache = ObjectCache(params['SORTEDLR_MAX_SLOTS'],
                                         params['SORTEDLR_MAX_SIZE'],
                                         'last row chunks')
        """A cache for the last row chunks. Only used for searches in
        the last row, and mainly useful for small indexes."""
        self.starts = np.empty(shape=self.nrows, dtype=np.int32)
        self.lengths = np.empty(shape=self.nrows, dtype=np.int32)
        self.sorted._init_sorted_slice(self)
        self.dirtycache = False

    def search(self, item):
        """Do a binary search in this index for an item."""

        if profile:
            tref = clock()
        if profile:
            show_stats("Entering search", tref)

        if self.dirtycache:
            self.restorecache()

        # An empty item or if left limit is larger than the right one
        # means that the number of records is always going to be empty,
        # so we avoid further computation (including looking up the
        # limits cache).
        if not item or item[0] > item[1]:
            self.starts[:] = 0
            self.lengths[:] = 0
            return 0

        tlen = 0
        # Check whether the item tuple is in the limits cache or not
        nslot = self.limboundscache.getslot(item)
        if nslot >= 0:
            startlengths = self.limboundscache.getitem(nslot)
            # Reset the lengths array (not necessary for starts)
            self.lengths[:] = 0
            # Now, set the interesting rows
            for nrow2, start, length in startlengths:
                self.starts[nrow2] = start
                self.lengths[nrow2] = length
                tlen = tlen + length
            return tlen
        # The item is not in cache. Do the real lookup.
        sorted = self.sorted
        if self.nslices > 0:
            if self.type in self.opt_search_types:
                # The next are optimizations. However, they hide the
                # CPU functions consumptions from python profiles.
                # You may want to de-activate them during profiling.
                if self.type == "int32":
                    tlen = sorted._search_bin_na_i(*item)
                elif self.type == "int64":
                    tlen = sorted._search_bin_na_ll(*item)
                elif self.type == "float16":
                    tlen = sorted._search_bin_na_e(*item)
                elif self.type == "float32":
                    tlen = sorted._search_bin_na_f(*item)
                elif self.type == "float64":
                    tlen = sorted._search_bin_na_d(*item)
                elif self.type == "float96":
                    tlen = sorted._search_bin_na_g(*item)
                elif self.type == "float128":
                    tlen = sorted._search_bin_na_g(*item)
                elif self.type == "uint32":
                    tlen = sorted._search_bin_na_ui(*item)
                elif self.type == "uint64":
                    tlen = sorted._search_bin_na_ull(*item)
                elif self.type == "int8":
                    tlen = sorted._search_bin_na_b(*item)
                elif self.type == "int16":
                    tlen = sorted._search_bin_na_s(*item)
                elif self.type == "uint8":
                    tlen = sorted._search_bin_na_ub(*item)
                elif self.type == "uint16":
                    tlen = sorted._search_bin_na_us(*item)
                else:
                    assert False, "This can't happen!"
            else:
                tlen = self.search_scalar(item, sorted)
        # Get possible remaining values in last row
        if self.nelementsSLR > 0:
            # Look for more indexes in the last row
            (start, stop) = self.search_last_row(item)
            self.starts[-1] = start
            self.lengths[-1] = stop - start
            tlen += stop - start

        if self.limboundscache.couldenablecache():
            # Get a startlengths tuple and save it in cache.
            # This is quite slow, but it is a good way to compress
            # the bounds info. Moreover, the .couldenablecache()
            # is doing a good work so as to avoid computing this
            # when it is not necessary to do it.
            startlengths = []
            for nrow, length in enumerate(self.lengths):
                if length > 0:
                    startlengths.append((nrow, self.starts[nrow], length))
            # Compute the size of the recarray (aproximately)
            # The +1 at the end is important to avoid 0 lengths
            # (remember, the object headers take some space)
            size = len(startlengths) * 8 * 2 + 1
            # Put this startlengths list in cache
            self.limboundscache.setitem(item, startlengths, size)

        if profile:
            show_stats("Exiting search", tref)
        return tlen

    # This is an scalar version of search. It works with strings as well.
    def search_scalar(self, item, sorted):
        """Do a binary search in this index for an item."""

        tlen = 0
        # Do the lookup for values fullfilling the conditions
        for i in range(self.nslices):
            (start, stop) = sorted._search_bin(i, item)
            self.starts[i] = start
            self.lengths[i] = stop - start
            tlen += stop - start
        return tlen

    def search_last_row(self, item):
        # Variable initialization
        item1, item2 = item
        bebounds = self.bebounds
        b0, b1 = bebounds[0], bebounds[-1]
        bounds = bebounds[1:-1]
        itemsize = self.dtype.itemsize
        sortedLRcache = self.sortedLRcache
        hi = self.nelementsSLR               # maximum number of elements
        rchunksize = self.chunksize // self.reduction

        nchunk = -1
        # Lookup for item1
        if nan_aware_gt(item1, b0):
            if nan_aware_le(item1, b1):
                # Search the appropriate chunk in bounds cache
                nchunk = bisect_left(bounds, item1)
                # Lookup for this chunk in cache
                nslot = sortedLRcache.getslot(nchunk)
                if nslot >= 0:
                    chunk = sortedLRcache.getitem(nslot)
                else:
                    begin = rchunksize * nchunk
                    end = rchunksize * (nchunk + 1)
                    if end > hi:
                        end = hi
                    # Read the chunk from disk
                    chunk = self.sortedLR._read_sorted_slice(
                        self.sorted, begin, end)
                    # Put it in cache.  It's important to *copy*
                    # the buffer, as it is reused in future reads!
                    sortedLRcache.setitem(nchunk, chunk.copy(),
                                          (end - begin) * itemsize)
                start = bisect_left(chunk, item1)
                start += rchunksize * nchunk
            else:
                start = hi
        else:
            start = 0
        # Lookup for item2
        if nan_aware_ge(item2, b0):
            if nan_aware_lt(item2, b1):
                # Search the appropriate chunk in bounds cache
                nchunk2 = bisect_right(bounds, item2)
                if nchunk2 != nchunk:
                    # Lookup for this chunk in cache
                    nslot = sortedLRcache.getslot(nchunk2)
                    if nslot >= 0:
                        chunk = sortedLRcache.getitem(nslot)
                    else:
                        begin = rchunksize * nchunk2
                        end = rchunksize * (nchunk2 + 1)
                        if end > hi:
                            end = hi
                        # Read the chunk from disk
                        chunk = self.sortedLR._read_sorted_slice(
                            self.sorted, begin, end)
                        # Put it in cache.  It's important to *copy*
                        # the buffer, as it is reused in future reads!
                        # See bug #60 in xot.carabos.com
                        sortedLRcache.setitem(nchunk2, chunk.copy(),
                                              (end - begin) * itemsize)
                stop = bisect_right(chunk, item2)
                stop += rchunksize * nchunk2
            else:
                stop = hi
        else:
            stop = 0
        return (start, stop)

    def get_chunkmap(self):
        """Compute a map with the interesting chunks in index."""

        if profile:
            tref = clock()
        if profile:
            show_stats("Entering get_chunkmap", tref)
        ss = self.slicesize
        nsb = self.nslicesblock
        nslices = self.nslices
        lbucket = self.lbucket
        indsize = self.indsize
        bucketsinblock = self.blocksize / lbucket
        nchunks = math.ceil(self.nelements / lbucket)
        chunkmap = np.zeros(shape=nchunks, dtype="bool")
        reduction = self.reduction
        starts = (self.starts - 1) * reduction + 1
        stops = (self.starts + self.lengths) * reduction
        starts[starts < 0] = 0    # All negative values set to zero
        indices = self.indices
        for nslice in range(self.nrows):
            start = starts[nslice]
            stop = stops[nslice]
            if stop > start:
                idx = np.empty(shape=stop - start, dtype='u%d' % indsize)
                if nslice < nslices:
                    indices._read_index_slice(nslice, start, stop, idx)
                else:
                    self.indicesLR._read_index_slice(start, stop, idx)
                if indsize == 8:
                    idx //= lbucket
                elif indsize == 2:
                    # The chunkmap size cannot be never larger than 'int_'
                    idx = idx.astype("int_")
                    offset = int((nslice // nsb) * bucketsinblock)
                    idx += offset
                elif indsize == 1:
                    # The chunkmap size cannot be never larger than 'int_'
                    idx = idx.astype("int_")
                    offset = (nslice * ss) // lbucket
                    idx += offset
                chunkmap[idx] = True
        # The case lbucket < nrowsinchunk should only happen in tests
        nrowsinchunk = self.nrowsinchunk
        if lbucket != nrowsinchunk:
            # Map the 'coarse grain' chunkmap into the 'true' chunkmap
            nelements = self.nelements
            tnchunks = math.ceil(nelements / nrowsinchunk)
            tchunkmap = np.zeros(shape=tnchunks, dtype="bool")
            ratio = lbucket / nrowsinchunk
            idx = chunkmap.nonzero()[0]
            starts = (idx * ratio).astype('int_')
            stops = np.ceil((idx + 1) * ratio).astype('int_')
            for start, stop in zip(starts, stops):
                tchunkmap[start:stop] = True
            chunkmap = tchunkmap
        if profile:
            show_stats("Exiting get_chunkmap", tref)
        return chunkmap

    def get_lookup_range(self, ops, limits):
        assert len(ops) in [1, 2]
        assert len(limits) in [1, 2]
        assert len(ops) == len(limits)

        column = self.column
        coldtype = column.dtype.base
        itemsize = coldtype.itemsize

        if len(limits) == 1:
            assert ops[0] in ['lt', 'le', 'eq', 'ge', 'gt']
            limit = limits[0]
            op = ops[0]
            if op == 'lt':
                range_ = (inftype(coldtype, itemsize, sign=-1),
                          nextafter(limit, -1, coldtype, itemsize))
            elif op == 'le':
                range_ = (inftype(coldtype, itemsize, sign=-1),
                          limit)
            elif op == 'gt':
                range_ = (nextafter(limit, +1, coldtype, itemsize),
                          inftype(coldtype, itemsize, sign=+1))
            elif op == 'ge':
                range_ = (limit,
                          inftype(coldtype, itemsize, sign=+1))
            elif op == 'eq':
                range_ = (limit, limit)

        elif len(limits) == 2:
            assert ops[0] in ('gt', 'ge') and ops[1] in ('lt', 'le')

            lower, upper = limits
            if lower > upper:
                # ``a <[=] x <[=] b`` is always false if ``a > b``.
                return ()

            if ops == ('gt', 'lt'):  # lower < col < upper
                range_ = (nextafter(lower, +1, coldtype, itemsize),
                          nextafter(upper, -1, coldtype, itemsize))
            elif ops == ('ge', 'lt'):  # lower <= col < upper
                range_ = (lower, nextafter(upper, -1, coldtype, itemsize))
            elif ops == ('gt', 'le'):  # lower < col <= upper
                range_ = (nextafter(lower, +1, coldtype, itemsize), upper)
            elif ops == ('ge', 'le'):  # lower <= col <= upper
                range_ = (lower, upper)

        return range_

    def _f_remove(self, recursive=False):
        """Remove this Index object."""

        # Index removal is always recursive,
        # no matter what `recursive` says.
        super()._f_remove(True)

    def __str__(self):
        """This provides a more compact representation than __repr__"""

        # The filters
        filters = []
        if self.filters.complevel:
            if self.filters.shuffle:
                filters.append('shuffle')
            if self.filters.bitshuffle:
                filters.append('bitshuffle')
            filters.append(f'{self.filters.complib}({self.filters.complevel})')
        return (f"Index({self.optlevel}, "
                f"{self.kind}{', '.join(filters)}).is_csi={self.is_csi}")

    def __repr__(self):
        """This provides more metainfo than standard __repr__"""

        cpathname = f"{self.table._v_pathname}.cols.{self.column.pathname}"
        retstr = f"""{self._v_pathname} (Index for column {cpathname})
  optlevel := {self.optlevel}
  kind := {self.kind}
  filters := {self.filters}
  is_csi := {self.is_csi}
  nelements := {self.nelements}
  chunksize := {self.chunksize}
  slicesize := {self.slicesize}
  blocksize := {self.blocksize}
  superblocksize := {self.superblocksize}
  dirty := {self.dirty}
  byteorder := {self.byteorder!r}
    sorted := {self.sorted}
    indices := {self.indices}
    ranges := {self.ranges}
    bounds := {self.bounds}
    sortedLR := {self.sortedLR}
    indicesLR := {self.indicesLR}"""
        return retstr


class IndexesDescG(NotLoggedMixin, Group):
    _c_classid = 'DINDEX'

    def _g_width_warning(self):
        warnings.warn(
            "the number of indexed columns on a single description group "
            "is exceeding the recommended maximum (%d); "
            "be ready to see PyTables asking for *lots* of memory "
            "and possibly slow I/O" % self._v_max_group_width,
            PerformanceWarning)


class IndexesTableG(NotLoggedMixin, Group):
    _c_classid = 'TINDEX'

    @property
    def auto(self):
        if 'AUTO_INDEX' not in self._v_attrs:
            return default_auto_index
        return self._v_attrs.AUTO_INDEX

    @auto.setter
    def auto(self, auto):
        self._v_attrs.AUTO_INDEX = bool(auto)

    @auto.deleter
    def auto(self):
        del self._v_attrs.AUTO_INDEX

    def _g_width_warning(self):
        warnings.warn(
            "the number of indexed columns on a single table "
            "is exceeding the recommended maximum (%d); "
            "be ready to see PyTables asking for *lots* of memory "
            "and possibly slow I/O" % self._v_max_group_width,
            PerformanceWarning)

    def _g_check_name(self, name):
        if not name.startswith('_i_'):
            raise ValueError(
                "names of index groups must start with ``_i_``: %s" % name)

    @property
    def table(self):
        """Accessor for the `Table` object of this `IndexesTableG`
        container."""
        names = self._v_pathname.split("/")
        tablename = names.pop()[3:]   # "_i_" is at the beginning
        parentpathname = "/".join(names)
        tablepathname = join_path(parentpathname, tablename)
        table = self._v_file._get_node(tablepathname)
        return table


class OldIndex(NotLoggedMixin, Group):
    """This is meant to hide indexes of PyTables 1.x files."""

    _c_classid = 'CINDEX'
