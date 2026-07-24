########################################################################
#
# License: BSD
# Created:
# Author:  Francesc Alted - faltet@pytables.com
#
# $Id$
#
########################################################################

from numpy cimport ndarray


# Declaration of instance variables for shared classes
# The NodeCache class is useful for caching general objects (like Nodes).
cdef class NodeCache:
  cdef readonly long nslots
  cdef long nextslot
  cdef object nodes, paths
  cdef object setitem(self, object path, object node)
  cdef long getslot(self, object path)
  cdef object cpop(self, object path)


# Base class for other caches
cdef class BaseCache:
  cdef int iscachedisabled, incsetcount
  cdef long setcount, getcount, containscount
  cdef long disablecyclecount, disableeverycycles
  cdef long enablecyclecount, enableeverycycles
  cdef double nprobes, hitratio
  cdef long seqn_, nextslot, nslots
  cdef long *ratimes
  cdef double lowesthr
  cdef ndarray atimes
  cdef object name
  cdef int checkhitratio(self)
  cdef int couldenablecache_(self)
  cdef long incseqn(self)


#  Helper class for ObjectCache
cdef class ObjectNode:
  cdef object key, obj
  cdef long nslot


# The ObjectCache class is useful for general python objects
cdef class ObjectCache(BaseCache):
  cdef long maxcachesize, cachesize, maxobjsize
  cdef long *rsizes
  cdef ndarray sizes
  cdef object __list, __dict
  cdef ObjectNode mrunode
  cdef removeslot_(self, long nslot)
  cdef clearcache_(self)
  cdef updateslot_(self, long nslot, long size, object key, object value)
  cdef long setitem_(self, object key, object value, long size)
  cdef long getslot_(self, object key)
  cdef object getitem_(self, long nslot)


# The NumCache class is useful for caching numerical data in an efficient way
cdef class NumCache(BaseCache):
  cdef long itemsize, slotsize
  cdef ndarray cacheobj, keys
  cdef void *rcache
  cdef long long *rkeys
  cdef object __dict
  cdef void *getaddrslot_(self, long nslot)
  cdef long setitem_(self, long long key, void *data, long start)
  cdef long setitem1_(self, long long key)
  cdef long getslot_(self, long long key)
  cdef getitem_(self, long nslot, void *data, long start)
  cdef void *getitem1_(self, long nslot)


## Local Variables:
## mode: python
## py-indent-offset: 2
## tab-width: 2
## fill-column: 78
## End:
