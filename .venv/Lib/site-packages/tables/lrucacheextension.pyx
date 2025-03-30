########################################################################
#
# License: BSD
# Created: Aug 13, 2006
# Author:  Francesc Alted - faltet@pytables.com
#
# $Id: $
#
########################################################################

"""Cython interface for several LRU cache systems.

Classes (type extensions):

    NodeCache
    ObjectCache
    NumCache

Functions:

Misc variables:

"""

cdef extern from "Python.h":
    int PyUnicode_Compare(object, object)

import sys

import numpy as np

from numpy cimport import_array, ndarray, PyArray_DATA
from libc.string cimport memcpy, strcmp
from cpython.unicode cimport PyUnicode_Check

from .parameters import (
    DISABLE_EVERY_CYCLES,
    ENABLE_EVERY_CYCLES,
    LOWEST_HIT_RATIO,
)

#----------------------------------------------------------------------------
# Initialization code.
# The numpy API requires this function to be called before
# using any numpy facilities in an extension module.
import_array()
#----------------------------------------------------------------------------


# ------- Minimalist NodeCache for nodes in PyTables ---------

# The next NodeCache code relies on the fact that a node that is
# fetched from the cache will be removed from it. Said in other words:
# "A node cannot be alive and dead at the same time."

# Thanks to the above behaviour, the next code has been stripped down
# to a bare minimum (the info in cache is kept in just 2 lists).

#*********************** Important note! *****************************
# The code behind has been carefully tuned to serve the needs of
# PyTables cache for nodes. As a consequence, it is no longer
# appropriate as a general LRU cache implementation. You have been
# warned!.  F. Alted 2006-08-08
#*********************************************************************


cdef class NodeCache:
  """Least-Recently-Used (LRU) cache for PyTables nodes."""

  def __init__(self, nslots):
    """Maximum nslots of the cache.

    If more than 'nslots' elements are added to the cache,
    the least-recently-used ones will be discarded.

    """

    if nslots < 0:
      raise ValueError("Negative number (%s) of slots!" % nslots)
    self.nslots = nslots
    self.nextslot = 0
    self.nodes = []
    self.paths = []

  def __len__(self):
    return len(self.nodes)

  def __setitem__(self, path, node):
    self.setitem(path, node)

  cdef setitem(self, object path, object node):
    """Puts a new node in the node list."""

    if self.nslots == 0:   # Oops, the cache is set to empty
      return
    # Check if we are growing out of space
    if self.nextslot == self.nslots:
      # It is critical to reduce nextslot *before* the preemption of
      # the LRU node.  If not, this can lead with problems in situations
      # with very small caches (length 1 or so).
      # F. Alted 2008-10-22
      self.nextslot = self.nextslot - 1
      # Remove the LRU node and path (the start of the lists)
      del self.nodes[0]
      del self.paths[0]
    # The equality protection has been put for situations in which a
    # node is being preempted and added simultaneously (with very small
    # caches).
    if len(self.nodes) == len(self.paths):
      # Add the node and path to the end of its lists
      self.nodes.append(node)
      self.paths.append(path)
      self.nextslot = self.nextslot + 1

  def __contains__(self, path):
    if self.getslot(path) == -1:
      return 0
    else:
      return 1

  cdef long getslot(self, object path):
    """Checks whether path is in this cache or not."""

    cdef long i, nslot, compare

    nslot = -1  # -1 means not found
    if PyUnicode_Check(path):
        # Start looking from the trailing values (most recently used)
        for i from self.nextslot > i >= 0:
          #if strcmp(<char *>encoded_path, <char *>self.paths[i]) == 0:
          if PyUnicode_Compare(path, self.paths[i]) == 0:
            nslot = i
            break
    else:
        # Start looking from the trailing values (most recently used)
        for i from self.nextslot > i >= 0:
          #if strcmp(<char *>path, <char *>self.paths[i]) == 0:
          if PyUnicode_Check(self.paths[i]):
            compare = PyUnicode_Compare(path, self.paths[i])
          else:
            compare = strcmp(<char *>path, <char *>self.paths[i])
          if compare == 0:
            nslot = i
            break

    return nslot

  __marker = object()

  def pop(self, path, d=__marker):
    try:
      node = self.cpop(path)
    except KeyError:
      if d is not self.__marker:
        return d
      else:
        raise
    else:
      return node

  cdef object cpop(self, object path):
    cdef long nslot

    nslot = self.getslot(path)
    if nslot == -1:
        raise KeyError(path)
    else:
        node = self.nodes[nslot]
        del self.nodes[nslot]
        del self.paths[nslot]
        self.nextslot = self.nextslot - 1
    return node

  def __iter__(self):
    # Do a copy of the paths list because it can be modified in the middle of
    # the iterator!
    copy = self.paths[:]
    return iter(copy)

  def __repr__(self):
    return "<%s (%d elements)>" % (str(self.__class__), len(self.paths))


########################################################################
# Common code for other LRU cache classes
########################################################################

cdef class BaseCache:
  """Base class that implements automatic probing/disabling of the cache."""

  def __init__(self, long nslots, object name):

    if nslots < 0:
      raise ValueError("Negative number (%s) of slots!" % nslots)
    self.setcount = 0;  self.getcount = 0;  self.containscount = 0
    self.enablecyclecount = 0;  self.disablecyclecount = 0
    self.iscachedisabled = False  # Cache is enabled by default
    self.disableeverycycles = DISABLE_EVERY_CYCLES
    self.enableeverycycles = ENABLE_EVERY_CYCLES
    self.lowesthr = LOWEST_HIT_RATIO
    self.nprobes = 0.0;  self.hitratio = 0.0
    self.nslots = nslots
    self.seqn_ = 0;  self.nextslot = 0
    self.name = name
    self.incsetcount = False
    # The array for keeping the access times (using long ints here)
    self.atimes = <ndarray>np.zeros(shape=nslots, dtype=np.int_)
    self.ratimes = <long *>PyArray_DATA(self.atimes)

  def __len__(self):
    return self.nslots

  # Machinery for determining whether the hit ratio is being effective
  # or not.  If not, the cache will be disabled. The efficency will be
  # checked every cycle (the time that the cache would be refilled
  # completely).  In situations where the cache is not being re-filled
  # (i.e. it is not enabled) for a long time, it is forced to be
  # re-enabled when a certain number of cycles has passed so as to
  # check whether a new scenario where the cache can be useful again
  # has come.
  # F. Alted 2006-08-09
  cdef int checkhitratio(self):
    cdef double hitratio
    cdef long nslot

    if self.setcount > self.nslots:
      self.disablecyclecount = self.disablecyclecount + 1
      self.enablecyclecount = self.enablecyclecount + 1
      self.nprobes = self.nprobes + 1
      hitratio = <double>self.getcount / self.containscount
      self.hitratio = self.hitratio + hitratio
      # Reset the hit counters
      self.setcount = 0;  self.getcount = 0;  self.containscount = 0
      if (not self.iscachedisabled and
          self.disablecyclecount >= self.disableeverycycles):
        # Check whether the cache is being effective or not
        if hitratio < self.lowesthr:
          # Hit ratio is low. Disable the cache.
          self.iscachedisabled = True
        else:
          # Hit ratio is acceptable. (Re-)Enable the cache.
          self.iscachedisabled = False
        self.disablecyclecount = 0
      if self.enablecyclecount >= self.enableeverycycles:
        # We have reached the time for forcing the cache to act again
        self.iscachedisabled = False
        self.enablecyclecount = 0
    return not self.iscachedisabled

  def couldenablecache(self):
    return self.couldenablecache_()

  # Check whether the cache is enabled or *could* be enabled in the next
  # setitem operation. This method can be used in order to probe whether
  # an (expensive) operation to be done before a .setitem() is worth the
  # effort or not.
  cdef int couldenablecache_(self):

    if self.nslots == 0:
      return False
    # Increment setitem because it can be that .setitem() doesn't
    # get called after calling this.
    self.setcount = self.setcount + 1;  self.incsetcount = True
    if self.iscachedisabled:
      if self.setcount == self.nslots:
        # The cache *could* be enabled in the next setitem operation
        return True
      else:
        return False
    else:
      return True

  # Increase the access time (implemented as a C long sequence)
  cdef long incseqn(self):

    self.seqn_ = self.seqn_ + 1
    if self.seqn_ < 0:
      # Ooops, the counter has run out of range! Reset all the access times.
      self.atimes[:] = sys.maxsize
      # Set the counter to 1 (to indicate that it is newer than existing ones)
      self.seqn_ = 1
    return self.seqn_

  def __repr__(self):
    return "<%s(%s) (%d elements)>" % (self.name, str(self.__class__),
                                       self.nslots)


########################################################################
#  Helper class for ObjectCache
########################################################################

cdef class ObjectNode:
  """Record of a cached value. Not for public consumption."""

  def __init__(self, object key, object obj, long nslot):
    object.__init__(self)
    self.key = key
    self.obj = obj
    self.nslot = nslot

  def __repr__(self):
    return "<%s %s (slot #%s) => %s>" % (self.__class__, self.key, self.nslot,
                                         self.object)


########################################################################
#  Minimalistic LRU cache implementation for general python objects
#        This is a *true* general lru cache for python objects
########################################################################

cdef class ObjectCache(BaseCache):
  """Least-Recently-Used (LRU) cache specific for python objects."""

  def __init__(self, long nslots, long maxcachesize, object name):
    """Maximum size of the cache.

    If more than 'nslots' elements are added to the cache,
    the least-recently-used ones will be discarded.

    Parameters:
    nslots - The number of slots in cache
    name - A descriptive name for this cache

    """

    super().__init__(nslots, name)
    self.cachesize = 0
    self.maxcachesize = maxcachesize
    # maxobjsize will be the same as the maximum cache size
    self.maxobjsize = maxcachesize
    self.__list = [None]*nslots
    self.__dict = {}
    self.mrunode = <ObjectNode>None   # Most Recent Used node
    # The array for keeping the object size (using long ints here)
    self.sizes = <ndarray>np.zeros(shape=nslots, dtype=np.int_)
    self.rsizes = <long *>PyArray_DATA(self.sizes)

  # Clear cache
  cdef clearcache_(self):
    self.__list = [None]*self.nslots
    self.__dict = {}
    self.mrunode = <ObjectNode>None
    self.cachesize = 0
    self.nextslot = 0
    self.seqn_ = 0

  # Remove a slot (if it exists in cache)
  cdef removeslot_(self, long nslot):
    cdef ObjectNode node

    assert nslot < self.nslots, "Attempting to remove beyond cache capacity."
    node = self.__list[nslot]
    if node is not None:
      self.__list[nslot] = None
      del self.__dict[node.key]
      self.cachesize = self.cachesize - self.rsizes[nslot]
      self.rsizes[nslot] = 0
      if self.mrunode and self.mrunode.nslot == nslot:
        self.mrunode = <ObjectNode>None
    # The next slot to be updated will be this one
    self.nextslot = nslot

  # Update a slot
  cdef updateslot_(self, long nslot, long size, object key, object value):
    cdef ObjectNode node, oldnode
    cdef long nslot1, nslot2
    cdef object lruidx

    assert nslot < self.nslots, "Number of nodes exceeding cache capacity."
    # Remove the previous nslot
    self.removeslot_(nslot)
    # Protection against too large data cache size
    while size + self.cachesize > self.maxcachesize:
      # Remove the LRU node among the 10 largest ones
      largidx = self.sizes.argsort()[-10:]
      nslot1 = self.atimes[largidx].argmin()
      nslot2 = largidx[nslot1]
      self.removeslot_(nslot2)
    # Insert the new one
    node = ObjectNode(key, value, nslot)
    self.ratimes[nslot] = self.incseqn()
    self.rsizes[nslot] = size
    self.__list[nslot] = node
    self.__dict[key] = node
    self.mrunode = node
    self.cachesize = self.cachesize + size
    # The next slot to update will be the LRU
    self.nextslot = self.atimes.argmin()

  # Put the object to the data in cache (for Python calls)
  def setitem(self, object key, object value, object size):
    return self.setitem_(key, value, size)

  # Put the object in cache (for cython calls)
  # size can be the exact size of the value object or an estimation.
  cdef long setitem_(self, object key, object value, long size):
    cdef long nslot

    if self.nslots == 0:   # The cache has been set to empty
      return -1
    nslot = -1
    # Perhaps setcount has been already incremented in couldenablecache()
    if not self.incsetcount:
      self.setcount = self.setcount + 1
    else:
      self.incsetcount = False
    if size > self.maxobjsize:  # Check if the object is too large
      return -1
    if self.checkhitratio():
      nslot = self.nextslot
      self.updateslot_(nslot, size, key, value)
    else:
      # Empty the cache because it is not effective and it is taking space
      self.clearcache_()
    return nslot

  # Tells whether the key is in cache or not
  def __contains__(self, object key):
    return self.__dict.has_key(key)

  # Tells in which slot the key is. If not found, -1 is returned.
  def getslot(self, object key):
    return self.getslot_(key)

  # Tells in which slot the key is. If not found, -1 is returned.
  cdef long getslot_(self, object key):
    cdef ObjectNode node

    if self.nslots == 0:   # The cache has been set to empty
      return -1
    self.containscount = self.containscount + 1
    # Give a chance to the MRU node
    node = self.mrunode
    if node and node.key == key:
      return node.nslot
    # No luck. Look in the dictionary.
    node = self.__dict.get(key)
    if node is <ObjectNode>None:
      return -1
    return node.nslot

  # Return the object to the data in cache (for Python calls)
  def getitem(self, object nslot):
    return self.getitem_(nslot)

  # Return the object to the data in cache (for cython calls)
  cdef object getitem_(self, long nslot):
    cdef ObjectNode node

    self.getcount = self.getcount + 1
    node = self.__list[nslot]
    self.ratimes[nslot] = self.incseqn()
    self.mrunode = node
    return node.obj

  def __repr__(self):
    if self.nprobes > 0:
      hitratio = self.hitratio / self.nprobes
    else:
      hitratio = <double>self.getcount / self.containscount
    return """<%s(%s)
  (%d maxslots, %d slots used, %.3f KB cachesize,
  hit ratio: %.3f, disabled? %s)>
  """ % (self.name, str(self.__class__), self.nslots, self.nextslot,
         self.cachesize / 1024., hitratio, self.iscachedisabled)


###################################################################
#  Minimalistic LRU cache implementation for numerical data
###################################################################
# The next code is more efficient in situations where efficiency is low.
###################################################################

#*********************** Important note! ****************************
# The code behind has been carefully tuned to serve the needs of
# caching numerical data. As a consequence, it is no longer appropriate
# as a general LRU cache implementation. You have been warned!.
# F. Alted 2006-08-09
#********************************************************************

cdef class NumCache(BaseCache):
  """Least-Recently-Used (LRU) cache specific for Numerical data."""

  def __init__(self, object shape, object dtype, object name):
    """Maximum size of the cache.

    If more than 'nslots' elements are added to the cache,
    the least-recently-used ones will be discarded.

    Parameters:
    shape - The rectangular shape of the cache (nslots, nelemsperslot)
    itemsize - The size of the element base in cache
    name - A descriptive name for this cache

    """

    cdef long nslots

    nslots = shape[0];  self.slotsize = shape[1]
    if nslots >= 1<<16:
      # nslots can't be higher than 2**16. Will silently trunk the number.
      nslots = <long>((1<<16)-1)  # Cast makes cython happy here
    super().__init__(nslots, name)
    self.itemsize = dtype.itemsize
    self.__dict = {}
    # The cache object where all data will go
    # The last slot is to allow the setitem1_ method to still return
    # a valid scratch area for writing purposes
    self.cacheobj = <ndarray>np.empty(shape=(nslots+1, self.slotsize),
                                         dtype=dtype)
    self.rcache = PyArray_DATA(self.cacheobj)
    # The array for keeping the keys of slots
    self.keys = <ndarray>(-np.ones(shape=nslots, dtype=np.int64))
    self.rkeys = <long long *>PyArray_DATA(self.keys)

  # Returns the address of nslot
  cdef void *getaddrslot_(self, long nslot):
    if nslot >= 0:
      return <char *>self.rcache + nslot * self.slotsize * self.itemsize
    else:
      return <char *>self.rcache + self.nslots * self.slotsize * self.itemsize

  def setitem(self, long long key, ndarray nparr, long start):
    return self.setitem_(key, PyArray_DATA(nparr), start)

  # Copy the new data into a cache slot
  cdef long setitem_(self, long long key, void *data, long start):
    cdef long nslot

    nslot = self.setitem1_(key)
    if nslot >= 0:
      # Copy the data to cache
      memcpy(<char *>self.rcache + nslot * self.slotsize * self.itemsize,
             <char *>data + start * self.itemsize,
             self.slotsize * self.itemsize)
    return nslot

  # Return a cache data pointer appropriate to save data.
  # Even if the cache is disabled, this will return a -1, which is
  # the last element in the cache.
  # This version avoids a memcpy of data, but the user should be
  # aware that data in nslot cannot be overwritten!
  cdef long setitem1_(self, long long key):
    cdef long nslot
    cdef object key2

    if self.nslots == 0:   # Oops, the cache is set to empty
      return -1
    # Perhaps setcount has been already incremented in couldenablecache()
    if not self.incsetcount:
      self.setcount = self.setcount + 1
    else:
      self.incsetcount = False
    nslot = -1
    if self.checkhitratio():
      # Check if we are growing out of space
      if self.nextslot == self.nslots:
        # Get the least recently used slot
        nslot = self.atimes.argmin()
        # Remove the slot from the dict
        key2 = self.keys[nslot]
        del self.__dict[key2]
        self.nextslot = self.nextslot - 1
      else:
        # Get the next slot available
        nslot = self.nextslot
      # Insert the slot in the dictionary
      self.__dict[key] = nslot
      self.keys[nslot] = key
      self.ratimes[nslot] = self.incseqn()
      self.nextslot = self.nextslot + 1
      # The next reduces the performance of the cache in scenarios where
      # the efficicency is near to zero.  I don't understand exactly why.
      # F. Alted 24-03-2008
    elif self.nextslot > 0:
      # Empty the cache if needed
      self.__dict.clear()
      self.nextslot = 0
    return nslot

  def getslot(self, long long key):
    return self.getslot_(key)

  # Tells in which slot key is. If not found, -1 is returned.
  cdef long getslot_(self, long long key):
    cdef object nslot

    self.containscount = self.containscount + 1
    if self.nextslot == 0:   # No chances for finding a slot
      return -1
    try:
      nslot = self.__dict[key]
    except KeyError:
      return -1
    return nslot

  def getitem(self, long nslot, ndarray nparr, long start):
    self.getitem_(nslot, PyArray_DATA(nparr), start)

  # This version copies data in cache to data+start.
  # The user should be responsible to provide a large enough data buffer
  # to keep all the data.
  cdef getitem_(self, long nslot, void *data, long start):
    cdef void *cachedata

    cachedata = self.getitem1_(nslot)
    # Copy the data in cache to destination
    memcpy(<char *>data + start * self.itemsize, cachedata,
           self.slotsize * self.itemsize)

  # Return the pointer to the data in cache
  # This version avoids a memcpy of data, but the user should be
  # aware that data in nslot cannot be overwritten!
  cdef void *getitem1_(self, long nslot):

    self.getcount = self.getcount + 1
    self.ratimes[nslot] = self.incseqn()
    return <char *>self.rcache + nslot * self.slotsize * self.itemsize

  def __repr__(self):
    cachesize = (self.nslots * self.slotsize * self.itemsize) / 1024.
    if self.nprobes > 0:
      hitratio = self.hitratio / self.nprobes
    elif self.containscount > 0:
      hitratio = <double>self.getcount / self.containscount
    else:
      hitratio = np.nan
    return """<%s(%s)
  (%d maxslots, %d slots used, %.3f KB cachesize,
  hit ratio: %.3f, disabled? %s)>
  """ % (self.name, str(self.__class__), self.nslots, self.nextslot,
         cachesize, hitratio, self.iscachedisabled)


## Local Variables:
## mode: python
## py-indent-offset: 2
## tab-width: 2
## fill-column: 78
## End:
