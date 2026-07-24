########################################################################
#
# License: BSD
# Created: May 18, 2006
# Author:  Francesc Alted - faltet@pytables.com
#
# $Id$
#
########################################################################

"""cython interface for keeping indexes classes.

Classes (type extensions):

    IndexArray
    CacheArray
    LastRowArray

Functions:

    keysort

Misc variables:

"""

import numpy as np
import cython

cimport numpy as cnp

from .exceptions import HDF5ExtError

# Types, constants, functions, classes & other objects from everywhere
from numpy cimport (
    import_array,
    ndarray,
    npy_int8,
    npy_int16,
    npy_int32,
    npy_int64,
    npy_uint8,
    npy_uint16,
    npy_uint32,
    npy_uint64,
    npy_float32,
    npy_float64,
    npy_float,
    npy_double,
    npy_longdouble,
    PyArray_BYTES,
    PyArray_DATA,
)

from .hdf5extension cimport Array

# These two types are defined in npy_common.h but not in cython's numpy.pxd
ctypedef unsigned char npy_bool
ctypedef npy_uint16 npy_float16

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, strncmp

from .definitions cimport hid_t, herr_t, hsize_t, H5Screate_simple, H5Sclose
from .lrucacheextension cimport NumCache

#-------------------------------------------------------------------

# External C functions

# Functions for optimized operations with ARRAY for indexing purposes
cdef extern from "H5ARRAY-opt.h" nogil:
  herr_t H5ARRAYOinit_readSlice(
    hid_t dataset_id, hid_t *mem_space_id, hsize_t count)
  herr_t H5ARRAYOread_readSlice(
    hid_t dataset_id, hid_t type_id,
    hsize_t irow, hsize_t start, hsize_t stop, void *data)
  herr_t H5ARRAYOread_readSortedSlice(
    hid_t dataset_id, hid_t mem_space_id, hid_t type_id,
    hsize_t irow, hsize_t start, hsize_t stop, void *data)
  herr_t H5ARRAYOread_readBoundsSlice(
    hid_t dataset_id, hid_t mem_space_id, hid_t type_id,
    hsize_t irow, hsize_t start, hsize_t stop, void *data)
  herr_t H5ARRAYOreadSliceLR(
    hid_t dataset_id, hid_t type_id, hsize_t start, hsize_t stop, void *data)


# Functions for optimized operations for dealing with indexes
cdef extern from "idx-opt.h" nogil:
  int bisect_left_b(npy_int8 *a, long x, int hi, int offset)
  int bisect_left_ub(npy_uint8 *a, long x, int hi, int offset)
  int bisect_right_b(npy_int8 *a, long x, int hi, int offset)
  int bisect_right_ub(npy_uint8 *a, long x, int hi, int offset)
  int bisect_left_s(npy_int16 *a, long x, int hi, int offset)
  int bisect_left_us(npy_uint16 *a, long x, int hi, int offset)
  int bisect_right_s(npy_int16 *a, long x, int hi, int offset)
  int bisect_right_us(npy_uint16 *a, long x, int hi, int offset)
  int bisect_left_i(npy_int32 *a, long x, int hi, int offset)
  int bisect_left_ui(npy_uint32 *a, npy_uint32 x, int hi, int offset)
  int bisect_right_i(npy_int32 *a, long x, int hi, int offset)
  int bisect_right_ui(npy_uint32 *a, npy_uint32 x, int hi, int offset)
  int bisect_left_ll(npy_int64 *a, npy_int64 x, int hi, int offset)
  int bisect_left_ull(npy_uint64 *a, npy_uint64 x, int hi, int offset)
  int bisect_right_ll(npy_int64 *a, npy_int64 x, int hi, int offset)
  int bisect_right_ull(npy_uint64 *a, npy_uint64 x, int hi, int offset)
  int bisect_left_e(npy_float16 *a, npy_float64 x, int hi, int offset)
  int bisect_right_e(npy_float16 *a, npy_float64 x, int hi, int offset)
  int bisect_left_f(npy_float32 *a, npy_float64 x, int hi, int offset)
  int bisect_right_f(npy_float32 *a, npy_float64 x, int hi, int offset)
  int bisect_left_d(npy_float64 *a, npy_float64 x, int hi, int offset)
  int bisect_right_d(npy_float64 *a, npy_float64 x, int hi, int offset)
  int bisect_left_g(npy_longdouble *a, npy_longdouble x, int hi, int offset)
  int bisect_right_g(npy_longdouble *a, npy_longdouble x, int hi, int offset)


#----------------------------------------------------------------------------

# Initialization code

# The numpy API requires this function to be called before
# using any numpy facilities in an extension module.
import_array()

#---------------------------------------------------------------------------

ctypedef fused floating_type:
    npy_float32
    npy_float64
    npy_longdouble


ctypedef fused number_type:
    npy_int8
    npy_int16
    npy_int32
    npy_int64

    npy_uint8
    npy_uint16
    npy_uint32
    npy_uint64

    npy_float32
    npy_float64
    npy_longdouble

#===========================================================================
# Functions
#===========================================================================

#---------------------------------------------------------------------------
# keysort
#---------------------------------------------------------------------------

DEF PYA_QS_STACK = 100
DEF SMALL_QUICKSORT = 15

def keysort(ndarray array1, ndarray array2):
    """Sort array1 in-place. array2 is also sorted following the array1 order.

    array1 can be of any type, except complex or string.  array2 may be made of
    elements on any size.

    """
    cdef size_t size = cnp.PyArray_SIZE(array1)
    cdef size_t elsize1 = cnp.PyArray_ITEMSIZE(array1)
    cdef size_t elsize2 = cnp.PyArray_ITEMSIZE(array2)
    cdef int type_num = cnp.PyArray_TYPE(array1)

    # floating types
    if type_num == cnp.NPY_FLOAT16:
        _keysort[npy_float16](<npy_float16*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_FLOAT32:
        _keysort[npy_float32](<npy_float32*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_FLOAT64:
        _keysort[npy_float64](<npy_float64*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_LONGDOUBLE:
        _keysort[npy_longdouble](<npy_longdouble*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    # signed integer types
    elif type_num == cnp.NPY_INT8:
        _keysort[npy_int8](<npy_int8*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_INT16:
        _keysort[npy_int16](<npy_int16*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_INT32:
        _keysort[npy_int32](<npy_int32*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_INT64:
        _keysort[npy_int64](<npy_int64*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    # unsigned integer types
    elif type_num == cnp.NPY_UINT8:
        _keysort[npy_uint8](<npy_uint8*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_UINT16:
        _keysort[npy_uint16](<npy_uint16*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_UINT32:
        _keysort[npy_uint32](<npy_uint32*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_UINT64:
        _keysort[npy_uint64](<npy_uint64*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    # other
    elif type_num == cnp.NPY_BOOL:
        _keysort[npy_bool](<npy_bool*>PyArray_DATA(array1), PyArray_BYTES(array2), elsize2, size)
    elif type_num == cnp.NPY_STRING:
        _keysort_string(PyArray_BYTES(array1), elsize1, PyArray_BYTES(array2), elsize2, size)
    else:
        raise ValueError("Unknown array datatype")


cdef inline void swap_bytes(char *x, char *y, size_t n) noexcept nogil:
    if n == 8:
        (<npy_int64*>x)[0], (<npy_int64*>y)[0] = (<npy_int64*>y)[0], (<npy_int64*>x)[0]
    elif n == 4:
        (<npy_int32*>x)[0], (<npy_int32*>y)[0] = (<npy_int32*>y)[0], (<npy_int32*>x)[0]
    elif n == 2:
        (<npy_int16*>x)[0], (<npy_int16*>y)[0] = (<npy_int16*>y)[0], (<npy_int16*>x)[0]
    else:
        for i in range(n):
            x[i], y[i] = y[i], x[i]


cdef inline int less_than(number_type* a, number_type* b) nogil:
    if number_type in floating_type:
        return a[0] < b[0] or (b[0] != b[0] and a[0] == a[0])
    else:
        return a[0] < b[0]


@cython.cdivision(True)
cdef void _keysort(number_type* start1, char* start2, size_t elsize2, size_t n) noexcept nogil:
    cdef number_type *pl = start1
    cdef number_type *pr = start1 + (n - 1)

    cdef char *ipl = start2
    cdef char *ipr = start2 + (n - 1) * elsize2

    cdef number_type vp
    cdef char *ivp = <char *> malloc(elsize2)

    cdef number_type *stack[PYA_QS_STACK]
    cdef number_type **sptr = stack

    cdef char *istack[PYA_QS_STACK]
    cdef char **isptr = istack

    cdef size_t stack_index = 0

    cdef number_type *pm
    cdef number_type *pi
    cdef number_type *pj
    cdef number_type *pt
    cdef char *ipm
    cdef char *ipi
    cdef char *ipj
    cdef char *ipt

    while True:
        while pr - pl > SMALL_QUICKSORT:
            pm  = pl + ((pr - pl) >> 1)
            ipm  = ipl + ((ipr - ipl)//elsize2 >> 1)*elsize2

            if less_than(pm, pl):
                pm[0], pl[0] =  pl[0], pm[0]
                swap_bytes(ipm, ipl, elsize2)

            if less_than(pr, pm):
                pr[0], pm[0] =  pm[0], pr[0]
                swap_bytes(ipr, ipm, elsize2)

            if less_than(pm, pl):
                pm[0], pl[0] =  pl[0], pm[0]
                swap_bytes(ipm, ipl, elsize2)

            vp = pm[0]

            pi = pl
            ipi = ipl

            pj = pr - 1
            ipj = ipr - elsize2

            pm[0], pj[0] = pj[0], pm[0]
            swap_bytes(ipm, ipj, elsize2)

            while True:
                pi += 1
                ipi += elsize2
                while less_than(pi, &vp):
                    pi += 1
                    ipi += elsize2

                pj -= 1
                ipj -= elsize2
                while less_than(&vp, pj):
                    pj -= 1
                    ipj -= elsize2

                if pi >= pj:
                    break

                pi[0], pj[0] = pj[0], pi[0]
                swap_bytes(ipi, ipj, elsize2)

            pi[0], (pr-1)[0] = (pr-1)[0], pi[0]
            swap_bytes(ipi, ipr-elsize2, elsize2)

            # push largest partition on stack and proceed with the other
            if (pi - pl) < (pr - pi):
                sptr[0] = pi + 1
                sptr[1] = pr
                sptr += 2

                isptr[0] = ipi + elsize2
                isptr[1] = ipr
                isptr += 2

                pr = pi - 1
                ipr = ipi - elsize2
            else:
                sptr[0] = pl
                sptr[1] = pi - 1
                sptr += 2

                isptr[0] = ipl
                isptr[1] = ipi - elsize2
                isptr += 2

                pl = pi + 1
                ipl = ipi + elsize2

        pi = pl + 1
        ipi = ipl + elsize2
        while pi <= pr:
            vp = pi[0]
            memcpy(ivp, ipi, elsize2)

            pj = pi
            pt = pi - 1

            ipj = ipi
            ipt = ipi - elsize2

            while pj > pl and less_than(&vp, pt):
                pj[0] = pt[0]
                pj -= 1
                pt -= 1

                memcpy(ipj, ipt, elsize2)
                ipj -= elsize2
                ipt -= elsize2

            pj[0] = vp
            memcpy(ipj, ivp, elsize2)

            pi += 1
            ipi += elsize2

        if sptr == stack:
            break

        sptr -= 2
        pl = sptr[0]
        pr = sptr[1]

        isptr -= 2
        ipl = isptr[0]
        ipr = isptr[1]

    free(ivp)


@cython.cdivision(True)
cdef void _keysort_string(char* start1, size_t ss, char* start2, size_t ts, size_t n) noexcept nogil:
    cdef char *pl = start1
    cdef char *pr = start1 + (n - 1) * ss

    cdef char *ipl = start2
    cdef char *ipr = start2 + (n - 1) * ts

    cdef char *vp = <char *>malloc(ss)
    cdef char *ivp = <char *>malloc(ts)

    cdef char *stack[PYA_QS_STACK]
    cdef char **sptr = stack

    cdef char *istack[PYA_QS_STACK]
    cdef char **isptr = istack

    cdef size_t stack_index = 0

    cdef char *pm
    cdef char *pi
    cdef char *pj
    cdef char *pt

    cdef char *ipm
    cdef char *ipi
    cdef char *ipj
    cdef char *ipt

    while True:
        while pr - pl > <long>(SMALL_QUICKSORT * ss):
            pm  = pl + ((pr - pl)//ss >> 1)*ss
            ipm  = ipl + ((ipr - ipl)//ts >> 1)*ts

            if strncmp(pm, pl, ss) < 0:
                swap_bytes(pm, pl, ss)
                swap_bytes(ipm, ipl, ts)

            if strncmp(pr, pm, ss) < 0:
                swap_bytes(pr, pm, ss)
                swap_bytes(ipr, ipm, ts)

            if strncmp(pm, pl, ss) < 0:
                swap_bytes(pm, pl, ss)
                swap_bytes(ipm, ipl, ts)

            memcpy(vp, pm, ss)

            pi = pl
            ipi = ipl

            pj = pr - ss
            ipj = ipr - ts

            swap_bytes(pm, pj, ss)
            swap_bytes(ipm, ipj, ts)

            while True:
                pi += ss
                ipi += ts
                while strncmp(pi, vp, ss) < 0:
                    pi += ss
                    ipi += ts

                pj -= ss
                ipj -= ts
                while strncmp(vp, pj, ss) < 0:
                    pj -= ss
                    ipj -= ts

                if pi >= pj:
                    break

                swap_bytes(pi, pj, ss)
                swap_bytes(ipi, ipj, ts)

            swap_bytes(pi, pr-ss, ss)
            swap_bytes(ipi, ipr-ts, ts)

            # push largest partition on stack and proceed with the other
            if (pi - pl) < (pr - pi):
                sptr[0] = pi + ss
                sptr[1] = pr
                sptr += 2

                isptr[0] = ipi + ts
                isptr[1] = ipr
                isptr += 2

                pr = pi - ss
                ipr = ipi - ts
            else:
                sptr[0] = pl
                sptr[1] = pi - ss
                sptr += 2

                isptr[0] = ipl
                isptr[1] = ipi - ts
                isptr += 2

                pl = pi + ss
                ipl = ipi + ts

        pi = pl + ss
        ipi = ipl + ts

        while pi <= pr:
            memcpy(vp, pi, ss)
            memcpy(ivp, ipi, ts)

            pj = pi
            pt = pi - ss

            ipj = ipi
            ipt = ipi - ts

            while pj > pl and strncmp(vp, pt, ss) < 0:
                memcpy(pj, pt, ss)
                pj -= ss
                pt -= ss

                memcpy(ipj, ipt, ts)
                ipj -= ts
                ipt -= ts

            memcpy(pj, vp, ss)
            memcpy(ipj, ivp, ts)

            pi += ss
            ipi += ts

        if sptr == stack:
            break

        sptr -= 2
        pl = sptr[0]
        pr = sptr[1]

        isptr -= 2
        ipl = isptr[0]
        ipr = isptr[1]

    free(vp)
    free(ivp)

#---------------------------------------------------------------------------
# bisect
#---------------------------------------------------------------------------

# This has been copied from the standard module bisect.
# Checks for the values out of limits has been added at the beginning
# because I forsee that this should be a very common case.
# 2004-05-20
def _bisect_left(a, x, int hi):
  """Return the index where to insert item x in list a, assuming a is sorted.

  The return value i is such that all e in a[:i] have e < x, and all e in
  a[i:] have e >= x.  So if x already appears in the list, i points just
  before the leftmost x already there.

  """

  cdef int lo, mid

  lo = 0
  if x <= a[0]: return 0
  if a[-1] < x: return hi
  while lo < hi:
      mid = (lo+hi)//2
      if a[mid] < x: lo = mid+1
      else: hi = mid
  return lo


def _bisect_right(a, x, int hi):
  """Return the index where to insert item x in list a, assuming a is sorted.

  The return value i is such that all e in a[:i] have e <= x, and all e in
  a[i:] have e > x.  So if x already appears in the list, i points just
  beyond the rightmost x already there.

  """

  cdef int lo, mid

  lo = 0
  if x < a[0]: return 0
  if a[-1] <= x: return hi
  while lo < hi:
    mid = (lo+hi)//2
    if x < a[mid]: hi = mid
    else: lo = mid+1
  return lo


#===========================================================================
# Classes
#===========================================================================



cdef class Index:
  pass


cdef class CacheArray(Array):
  """Container for keeping index caches of 1st and 2nd level."""

  cdef hid_t mem_space_id

  cdef initread(self, int nbounds):
    # "Actions to accelerate the reads afterwards."

    # Precompute the mem_space_id
    if (H5ARRAYOinit_readSlice(self.dataset_id, &self.mem_space_id,
                               nbounds) < 0):
      raise HDF5ExtError("Problems initializing the bounds array data.")
    return

  cdef read_slice(self, hsize_t nrow, hsize_t start, hsize_t stop, void *rbuf):
    # "Read an slice of bounds."

    if (H5ARRAYOread_readBoundsSlice(
      self.dataset_id, self.mem_space_id, self.type_id,
      nrow, start, stop, rbuf) < 0):
      raise HDF5ExtError("Problems reading the bounds array data.")
    return

  def _g_close(self):
    super()._g_close()
    # Release specific resources of this class
    if self.mem_space_id > 0:
      H5Sclose(self.mem_space_id)


cdef class IndexArray(Array):
  """Container for keeping sorted and indices values."""

  cdef void    *rbufst
  cdef void    *rbufln
  cdef void    *rbufrv
  cdef void    *rbufbc
  cdef void    *rbuflb
  cdef hid_t   mem_space_id
  cdef int     l_chunksize, l_slicesize, nbounds, indsize
  cdef CacheArray bounds_ext
  cdef NumCache boundscache, sortedcache
  cdef ndarray bufferbc, bufferlb

  def _read_index_slice(self, hsize_t irow, hsize_t start, hsize_t stop,
                      ndarray idx):
    cdef herr_t ret
    cdef void *buf = PyArray_DATA(idx)

    # Do the physical read
    with nogil:
        ret = H5ARRAYOread_readSlice(self.dataset_id, self.type_id,
                                     irow, start, stop, buf)

    if ret < 0:
      raise HDF5ExtError("Problems reading the index indices.")


  def _init_sorted_slice(self, index):
    """Initialize the structures for doing a binary search."""

    cdef long ndims
    cdef int  rank, buflen, cachesize
    cdef char *bname
    cdef hsize_t count[2]
    cdef ndarray starts, lengths, rvcache
    cdef object maxslots, rowsize

    dtype = self.atom.dtype
    # Create the buffer for reading sorted data chunks if not created yet
    if <object>self.bufferlb is None:
      # Internal buffers
      self.bufferlb = np.empty(dtype=dtype, shape=self.chunksize)
      # Get the pointers to the different buffer data areas
      self.rbuflb = PyArray_DATA(self.bufferlb)
      # Init structures for accelerating sorted array reads
      rank = 2
      count[0] = 1
      count[1] = self.chunksize
      self.mem_space_id = H5Screate_simple(rank, count, NULL)
      # Cache some counters in local extension variables
      self.l_chunksize = self.chunksize
      self.l_slicesize = self.slicesize

    # Get the addresses of buffer data
    starts = index.starts
    lengths = index.lengths
    self.rbufst = PyArray_DATA(starts)
    self.rbufln = PyArray_DATA(lengths)
    # The 1st cache is loaded completely in memory and needs to be reloaded
    rvcache = index.ranges[:]
    self.rbufrv = PyArray_DATA(rvcache)
    index.rvcache = <object>rvcache
    # Init the bounds array for reading
    self.nbounds = index.bounds.shape[1]
    self.bounds_ext = <CacheArray>index.bounds
    self.bounds_ext.initread(self.nbounds)
    if str(dtype) in self._v_parent.opt_search_types:
      # The next caches should be defined only for optimized search types.
      # The 2nd level cache will replace the already existing ObjectCache and
      # already bound to the boundscache attribute. This way, the cache will
      # not be duplicated (I know, this smells badly, but anyway).
      params = self._v_file.params
      rowsize = (self.bounds_ext._v_chunkshape[1] * dtype.itemsize)
      maxslots = params['BOUNDS_MAX_SIZE'] // rowsize
      self.boundscache = <NumCache>NumCache(
        (maxslots, self.nbounds), dtype, 'non-opt types bounds')
      self.bufferbc = np.empty(dtype=dtype, shape=self.nbounds)
      # Get the pointer for the internal buffer for 2nd level cache
      self.rbufbc = PyArray_DATA(self.bufferbc)
      # Another NumCache for the sorted values
      rowsize = (self.chunksize*dtype.itemsize)
      maxslots = params['SORTED_MAX_SIZE'] // (self.chunksize*dtype.itemsize)
      self.sortedcache = <NumCache>NumCache(
        (maxslots, self.chunksize), dtype, 'sorted')



  cdef void *_g_read_sorted_slice(self, hsize_t irow, hsize_t start,
                                hsize_t stop):
    """Read the sorted part of an index."""

    with nogil:
        ret = H5ARRAYOread_readSortedSlice(
          self.dataset_id, self.mem_space_id, self.type_id,
          irow, start, stop, self.rbuflb)

    if ret < 0:
      raise HDF5ExtError("Problems reading the array data.")

    return self.rbuflb

  # can't time machine since this function is cdef'd
  #_g_read_sorted_slice = prveious_api(_g_read_sorted_slice)

  # This is callable from python
  def _read_sorted_slice(self, hsize_t irow, hsize_t start, hsize_t stop):
    """Read the sorted part of an index."""

    self._g_read_sorted_slice(irow, start, stop)
    return self.bufferlb


  cdef void *get_lru_bounds(self, int nrow, int nbounds):
    """Get the bounds from the cache, or read them."""

    cdef void *vpointer
    cdef long nslot

    nslot = self.boundscache.getslot_(nrow)
    if nslot >= 0:
      vpointer = self.boundscache.getitem1_(nslot)
    else:
      # Bounds row is not in cache. Read it and put it in the LRU cache.
      self.bounds_ext.read_slice(nrow, 0, nbounds, self.rbufbc)
      self.boundscache.setitem_(nrow, self.rbufbc, 0)
      vpointer = self.rbufbc
    return vpointer

  # can't time machine since get_lru_bounds() function is cdef'd

  cdef void *get_lru_sorted(self, int nrow, int ncs, int nchunk, int cs):
    """Get the sorted row from the cache or read it."""

    cdef void *vpointer
    cdef npy_int64 nckey
    cdef long nslot
    cdef hsize_t start, stop

    # Compute the number of chunk read and use it as the key for the cache.
    nckey = nrow*ncs+nchunk
    nslot = self.sortedcache.getslot_(nckey)
    if nslot >= 0:
      vpointer = self.sortedcache.getitem1_(nslot)
    else:
      # The sorted chunk is not in cache. Read it and put it in the LRU cache.
      start = cs*nchunk
      stop = cs*(nchunk+1)
      vpointer = self._g_read_sorted_slice(nrow, start, stop)
      self.sortedcache.setitem_(nckey, vpointer, 0)
    return vpointer

  # can't time machine since get_lru_sorted() function is cdef'd

  # Optimized version for int8
  def _search_bin_na_b(self, long item1, long item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_int8 *rbufrv
    cdef npy_int8 *rbufbc = NULL
    cdef npy_int8 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_int8 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1

      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_int8 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_b(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_int8 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_b(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_int8 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_b(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_int8 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_b(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for uint8
  def _search_bin_na_ub(self, long item1, long item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_uint8 *rbufrv
    cdef npy_uint8 *rbufbc = NULL
    cdef npy_uint8 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_uint8 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1

      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_uint8 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_ub(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_uint8 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_ub(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_uint8 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_ub(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_uint8 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_ub(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for int16
  def _search_bin_na_s(self, long item1, long item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_int16 *rbufrv
    cdef npy_int16 *rbufbc = NULL
    cdef npy_int16 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_int16 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_int16 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_s(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_int16 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_s(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_int16 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_s(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_int16 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_s(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for uint16
  def _search_bin_na_us(self, long item1, long item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_uint16 *rbufrv
    cdef npy_uint16 *rbufbc = NULL
    cdef npy_uint16 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_uint16 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_uint16 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_us(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_uint16 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_us(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_uint16 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_us(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_uint16 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_us(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for int32
  def _search_bin_na_i(self, long item1, long item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_int32 *rbufrv
    cdef npy_int32 *rbufbc = NULL
    cdef npy_int32 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_int32 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_int32 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_i(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_int32 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_i(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_int32 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_i(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_int32 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_i(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for uint32
  def _search_bin_na_ui(self, npy_uint32 item1, npy_uint32 item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_uint32 *rbufrv
    cdef npy_uint32 *rbufbc = NULL
    cdef npy_uint32 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_uint32 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_uint32 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_ui(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_uint32 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_ui(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_uint32 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_ui(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_uint32 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_ui(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for int64
  def _search_bin_na_ll(self, npy_int64 item1, npy_int64 item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_int64 *rbufrv
    cdef npy_int64 *rbufbc = NULL
    cdef npy_int64 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_int64 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_int64 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_ll(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_int64 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_ll(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_int64 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_ll(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_int64 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_ll(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for uint64
  def _search_bin_na_ull(self, npy_uint64 item1, npy_uint64 item2):
    cdef int cs, ss, ncs, nrow, nrows, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_uint64 *rbufrv
    cdef npy_uint64 *rbufbc = NULL
    cdef npy_uint64 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    rbufrv = <npy_uint64 *>self.rbufrv
    tlength = 0
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_uint64 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_ull(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_uint64 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_ull(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_uint64 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_ull(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_uint64 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_ull(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for float16
  def _search_bin_na_e(self, npy_float64 item1, npy_float64 item2):
    cdef int cs, ss, ncs, nrow, nrows, nrow2, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_float16 *rbufrv
    cdef npy_float16 *rbufbc = NULL
    cdef npy_float16 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    tlength = 0
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln
    # Limits not in cache, do a lookup
    rbufrv = <npy_float16 *>self.rbufrv
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1

      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_float16 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_e(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_float16 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_e(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_float16 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_e(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_float16 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_e(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for float32
  def _search_bin_na_f(self, npy_float64 item1, npy_float64 item2):
    cdef int cs, ss, ncs, nrow, nrows, nrow2, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln
    # Variables with specific type
    cdef npy_float32 *rbufrv
    cdef npy_float32 *rbufbc = NULL
    cdef npy_float32 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    tlength = 0
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln

    # Limits not in cache, do a lookup
    rbufrv = <npy_float32 *>self.rbufrv
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1
      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_float32 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_f(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_float32 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_f(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_float32 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_f(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_float32 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_f(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for float64
  def _search_bin_na_d(self, npy_float64 item1, npy_float64 item2):
    cdef int cs, ss, ncs, nrow, nrows, nrow2, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_float64 *rbufrv
    cdef npy_float64 *rbufbc = NULL
    cdef npy_float64 *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    tlength = 0
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln

    # Limits not in cache, do a lookup
    rbufrv = <npy_float64 *>self.rbufrv
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1

      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_float64 *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_d(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_float64 *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_d(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_float64 *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_d(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_float64 *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_d(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  # Optimized version for npy_longdouble/float96/float128
  def _search_bin_na_g(self, npy_longdouble item1, npy_longdouble item2):
    cdef int cs, ss, ncs, nrow, nrows, nrow2, nbounds, rvrow
    cdef int start, stop, tlength, length, bread, nchunk, nchunk2
    cdef int *rbufst
    cdef int *rbufln

    # Variables with specific type
    cdef npy_longdouble *rbufrv
    cdef npy_longdouble *rbufbc = NULL
    cdef npy_longdouble *rbuflb = NULL

    cs = self.l_chunksize
    ss = self.l_slicesize
    ncs = ss // cs
    nbounds = self.nbounds
    nrows = self.nrows
    tlength = 0
    rbufst = <int *>self.rbufst
    rbufln = <int *>self.rbufln

    # Limits not in cache, do a lookup
    rbufrv = <npy_longdouble *>self.rbufrv
    for nrow from 0 <= nrow < nrows:
      rvrow = nrow*2
      bread = 0
      nchunk = -1

      # Look if item1 is in this row
      if item1 > rbufrv[rvrow]:
        if item1 <= rbufrv[rvrow+1]:
          # Get the bounds row from the LRU cache or read them.
          rbufbc = <npy_longdouble *>self.get_lru_bounds(nrow, nbounds)
          bread = 1
          nchunk = bisect_left_g(rbufbc, item1, nbounds, 0)
          # Get the sorted row from the LRU cache or read it.
          rbuflb = <npy_longdouble *>self.get_lru_sorted(nrow, ncs, nchunk, cs)
          start = bisect_left_g(rbuflb, item1, cs, 0) + cs*nchunk
        else:
          start = ss
      else:
        start = 0
      # Now, for item2
      if item2 >= rbufrv[rvrow]:
        if item2 < rbufrv[rvrow+1]:
          if not bread:
            # Get the bounds row from the LRU cache or read them.
            rbufbc = <npy_longdouble *>self.get_lru_bounds(nrow, nbounds)
          nchunk2 = bisect_right_g(rbufbc, item2, nbounds, 0)
          if nchunk2 <> nchunk:
            # Get the sorted row from the LRU cache or read it.
            rbuflb = <npy_longdouble *>self.get_lru_sorted(nrow, ncs, nchunk2, cs)
          stop = bisect_right_g(rbuflb, item2, cs, 0) + cs*nchunk2
        else:
          stop = ss
      else:
        stop = 0
      length = stop - start
      tlength = tlength + length
      rbufst[nrow] = start
      rbufln[nrow] = length
    return tlength


  def _g_close(self):
    super()._g_close()
    # Release specific resources of this class
    if self.mem_space_id > 0:
      H5Sclose(self.mem_space_id)


cdef class LastRowArray(Array):
  """
  Container for keeping sorted and indices values of last rows of an index.
  """

  def _read_index_slice(self, hsize_t start, hsize_t stop, ndarray idx):
    """Read the reverse index part of an LR index."""

    cdef void *buf = PyArray_DATA(idx)
    with nogil:
        ret = H5ARRAYOreadSliceLR(self.dataset_id, self.type_id,
                                  start, stop, buf)

    if ret < 0:
      raise HDF5ExtError("Problems reading the index data in Last Row.")


  def _read_sorted_slice(self, IndexArray sorted, hsize_t start, hsize_t stop):
    """Read the sorted part of an LR index."""

    cdef void  *rbuflb

    rbuflb = sorted.rbuflb  # direct access to rbuflb: very fast.
    with nogil:
        ret = H5ARRAYOreadSliceLR(self.dataset_id, self.type_id,
                                  start, stop, rbuflb)

    if ret < 0:
      raise HDF5ExtError("Problems reading the index data.")
    return sorted.bufferlb[:stop-start]



## Local Variables:
## mode: python
## py-indent-offset: 2
## tab-width: 2
## fill-column: 78
## End:
