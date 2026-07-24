########################################################################
#
# License: BSD
# Created: June 17, 2005
# Author:  Francesc Alted - faltet@pytables.com
#
# $Id$
#
########################################################################

"""Here is where Table and Row extension types live.

Classes (type extensions):

    Table
    Row

Functions:

Misc variables:

"""
import os
import sys
import math
import platform
from time import time

import numpy as np

from .utils import SizeType
from .conditions import call_on_recarr
from .exceptions import HDF5ExtError
from .description import Col
from .utilsextension import (
    get_nested_field,
    atom_from_hdf5_type,
    create_nested_type,
    hdf5_to_np_ext_type,
    platform_byteorder,
    pttype_to_hdf5,
    pt_special_kinds,
    npext_prefixes_to_ptkinds,
    hdf5_class_to_string,
    H5T_STD_I64,
)

from numpy cimport (
    import_array,
    ndarray,
    npy_intp,
    PyArray_GETITEM,
    PyArray_SETITEM,
    PyArray_BYTES,
    PyArray_DATA,
    PyArray_NDIM,
    PyArray_STRIDE,
)
from cpython cimport PyErr_Clear
from libc.stdio cimport snprintf
from libc.stdint cimport int32_t
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, strdup, strcmp, strlen

from .definitions cimport (
    hid_t,
    herr_t,
    hsize_t,
    haddr_t,
    htri_t,
    hbool_t,
    H5F_ACC_RDONLY,
    H5P_DEFAULT,
    H5D_CHUNKED,
    H5T_DIR_DEFAULT,
    H5F_SCOPE_LOCAL,
    H5F_SCOPE_GLOBAL,
    H5T_COMPOUND,
    H5Tget_order,
    H5Fflush,
    H5Dget_create_plist,
    H5T_ORDER_LE,
    H5D_layout_t,
    H5Dopen,
    H5Dclose,
    H5Dread,
    H5Dget_type,
    H5Dget_space,
    H5Pget_layout,
    H5Pget_chunk,
    H5Pclose,
    H5Sget_simple_extent_ndims,
    H5Sget_simple_extent_dims,
    H5Sclose,
    H5T_class_t,
    H5Tget_size,
    H5Tset_size,
    H5Tcreate,
    H5Tcopy,
    H5Tclose,
    H5Tget_nmembers,
    H5Tget_member_name,
    H5Tget_member_type,
    H5Tget_native_type,
    H5Tget_member_offset,
    H5Tinsert,
    H5Tget_class,
    H5Tget_super,
    H5Tget_offset,
    H5T_cset_t,
    H5T_CSET_ASCII,
    H5T_CSET_UTF8,
    H5ATTRset_attribute_string,
    H5ATTRset_attribute,
    get_len_of_range,
    get_order,
    set_order,
    is_complex,
    conv_float64_timeval32,
    truncate_dset,
    H5free_memory,
)

# numpy functions & objects
from .hdf5extension cimport Leaf
from .utilsextension cimport get_native_type, cstr_to_pystr
from .lrucacheextension cimport ObjectCache, NumCache

#-----------------------------------------------------------------

# Optimized HDF5 API for PyTables
cdef extern from "H5TB-opt.h" nogil:

  ctypedef struct chunk_iter_op:
    size_t itemsize
    size_t chunkshape
    haddr_t *addrs

  int fill_chunk_addrs(hid_t dataset_id, hsize_t nchunks, chunk_iter_op *chunk_op)
  int clean_chunk_addrs(chunk_iter_op *chunk_op)

  herr_t H5TBOmake_table( char *table_title, hid_t loc_id, char *dset_name,
                          char *version, char *class_,
                          hid_t mem_type_id, hsize_t nrecords,
                          hsize_t chunk_size, hsize_t block_size,
                          void *fill_data, int compress,
                          char *complib, int shuffle, int fletcher32,
                          hbool_t track_times, hbool_t blosc2_support,
                          void *data )

  herr_t H5TBOread_records( char* filename, hbool_t blosc2_support,
                            chunk_iter_op chunk_op,
                            hid_t dataset_id, hid_t mem_type_id,
                            hsize_t start, hsize_t nrecords, void *data )

  herr_t H5TBOread_elements( hid_t dataset_id, hid_t mem_type_id,
                             hsize_t nrecords, void *coords, void *data )

  herr_t H5TBOappend_records( hbool_t blosc2_support, hid_t dataset_id,
                              hid_t mem_type_id, hsize_t start,
                              hsize_t nrecords, void *data )

  herr_t H5TBOwrite_records ( hbool_t blosc2_support, hid_t dataset_id,
                              hid_t mem_type_id, hsize_t start,
                              hsize_t nrecords, hsize_t step, void *data )

  herr_t write_records_blosc2( hid_t dataset_id, hid_t mem_type_id,
                               hsize_t start, hsize_t nrecords,
                               const void *data )

  herr_t H5TBOwrite_elements( hid_t dataset_id, hid_t mem_type_id,
                              hsize_t nrecords, void *coords, void *data )

  herr_t H5TBOdelete_records( char* filename, hbool_t blosc2_support,
                              chunk_iter_op chunk_op,
                              hid_t dataset_id, hid_t mem_type_id,
                              hsize_t ntotal_records, size_t src_size,
                              hsize_t start, hsize_t nrecords,
                              hsize_t maxtuples )


#----------------------------------------------------------------------------

# Initialization code

# The numpy API requires this function to be called before
# using any numpy facilities in an extension module.
import_array()

#-------------------------------------------------------------


# Private functions
cdef get_nested_field_cache(recarray, fieldname, fieldcache):
  """Get the maybe nested field named `fieldname` from the `recarray`.

  The `fieldname` may be a simple field name or a nested field name with
  slah-separated components. It can also be an integer specifying the position
  of the field.

  """

  try:
    field = fieldcache[fieldname]
  except KeyError:
    # Check whether fieldname is an integer and if so, get the field
    # straight from the recarray dictionary (it can't be anywhere else)
    if isinstance(fieldname, int):
      field = recarray[fieldname]
    else:
      field = get_nested_field(recarray, fieldname)
    fieldcache[fieldname] = field
  return field


cdef join_path(object parent, object name):
  if parent == "":
    return name
  else:
    return parent + '/' + name


# Public classes

cdef class Table(Leaf):
  # instance variables
  cdef void *wbuf
  cdef chunk_iter_op chunk_op
  cdef hbool_t blosc2_support_read
  cdef hbool_t blosc2_support_write

  def _create_table(self, title, complib, obversion):
    cdef int     offset
    cdef int     ret
    cdef long    buflen
    cdef hid_t   oid
    cdef void    *data
    cdef hsize_t nrows
    cdef bytes   class_
    cdef ndarray wdflts
    cdef void    *fill_data
    cdef ndarray recarr
    cdef object  name
    cdef bytes encoded_title, encoded_complib, encoded_obversion
    cdef char *ctitle = NULL
    cdef char *cobversion = NULL
    cdef bytes encoded_name
    cdef char fieldname[128]
    cdef int i
    cdef H5T_cset_t cset = H5T_CSET_ASCII

    encoded_title = title.encode('utf-8')
    encoded_complib = complib.encode('utf-8')
    encoded_obversion = obversion.encode('utf-8')
    encoded_name = self.name.encode('utf-8')

    # Get the C pointer
    ctitle = encoded_title
    cobversion = encoded_obversion

    # Compute the complete compound datatype based on the table description
    self.disk_type_id = create_nested_type(self.description, self.byteorder)
    #self.type_id = H5Tcopy(self.disk_type_id)
    # A H5Tcopy only is not enough, as we want the in-memory type to be
    # in the byteorder of the machine (sys.byteorder).
    self.type_id = create_nested_type(self.description, sys.byteorder)

    # The fill values area
    wdflts = self._v_wdflts
    if wdflts is None:
      fill_data = NULL
    else:
      fill_data = PyArray_DATA(wdflts)

    # test if there is data to be saved initially
    if self._v_recarray is not None:
      recarr = self._v_recarray
      data = PyArray_DATA(recarr)
    else:
      data = NULL

    # Decide whether Blosc2 optimized operations can be used.
    self.blosc2_support_write = (
        (self.byteorder == sys.byteorder) and
        (not self.filters.fletcher32) and
        (self.filters.complib is not None) and
        (self.filters.complib.startswith("blosc2")))
    # For reading, Windows does not support re-opening a file twice
    # in not read-only mode (for good reason), so we cannot use the
    # blosc2 opt
    self.blosc2_support_read = (
        self.blosc2_support_write and
        ((platform.system().lower() != 'windows') or
         (self._v_file.mode == 'r')))

    class_ = self._c_classid.encode('utf-8')
    cdef hsize_t blocksize = int(os.environ.get("PT_DEFAULT_B2_BLOCKSIZE", "0"))
    self.dataset_id = H5TBOmake_table(ctitle, self.parent_id, encoded_name,
                                      cobversion, class_, self.disk_type_id,
                                      self.nrows, self.chunkshape[0],
                                      blocksize, fill_data,
                                      self.filters.complevel, encoded_complib,
                                      self.filters.shuffle_bitshuffle,
                                      self.filters.fletcher32,
                                      self._want_track_times,
                                      self.blosc2_support_write, data)
    if self.dataset_id < 0:
      raise HDF5ExtError("Problems creating the table")

    if self._v_file.params['PYTABLES_SYS_ATTRS']:
      cset = H5T_CSET_UTF8
      # Set the conforming table attributes
      # Attach the CLASS attribute
      ret = H5ATTRset_attribute_string(self.dataset_id, "CLASS", class_,
                                       len(class_), cset)
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in table:\n %s." %
                           ("CLASS", self.name))
      # Attach the VERSION attribute
      ret = H5ATTRset_attribute_string(self.dataset_id, "VERSION", cobversion,
                                       len(encoded_obversion), cset)
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in table:\n %s." %
                           ("VERSION", self.name))
      # Attach the TITLE attribute
      ret = H5ATTRset_attribute_string(self.dataset_id, "TITLE", ctitle,
                                       len(encoded_title), cset)
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in table:\n %s." %
                           ("TITLE", self.name))
      # Attach the NROWS attribute
      nrows = self.nrows
      ret = H5ATTRset_attribute(self.dataset_id, "NROWS", H5T_STD_I64,
                                0, NULL, <char *>&nrows)
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in table:\n %s." %
                           ("NROWS", self.name))

      # Attach the FIELD_N_NAME attributes
      # We write only the first level names
      for i, name in enumerate(self.description._v_names):
        snprintf(fieldname, 128, "FIELD_%d_NAME", i)
        encoded_name = name.encode('utf-8')
        ret = H5ATTRset_attribute_string(self.dataset_id, fieldname,
                                         encoded_name, len(encoded_name),
                                         cset)
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in table:\n %s." %
                           (fieldname, self.name))

    # If created in PyTables, the table is always chunked
    self._chunked = True  # Accessible from python

    # Initialize blosc2 struct for chunk addresses
    self.chunk_op = chunk_iter_op(self.description._v_itemsize, self.chunkshape[0], NULL)

    # Finally, return the object identifier.
    return self.dataset_id


  cdef get_nested_type(self, hid_t type_id, hid_t native_type_id,
                       object colpath, object field_byteorders):
    """Open a nested type and return a nested dictionary as description."""

    cdef hid_t   member_type_id, native_member_type_id, member_offset
    cdef hsize_t nfields, i
    cdef hsize_t dims[1]
    cdef size_t  itemsize
    cdef char    *c_colname
    cdef H5T_class_t class_id
    cdef char    c_byteorder2[11]  # "irrelevant" fits easily here
    cdef char    *sys_byteorder
    cdef object  desc, colobj, colpath2, typeclassname, typeclass
    cdef object  byteorder
    cdef str     colname, byteorder2

    offset = 0
    desc = {}
    # Get the number of members
    nfields = H5Tget_nmembers(type_id)

    # Iterate through fields to get the correct order that elements may appear in
    # The object type can be stored not in order, so order based on the offset in the data
    position_order = []
    for i in range(nfields):
      member_offset = H5Tget_member_offset(type_id, i)
      position_order.append((member_offset, i))

    position_order.sort()

    # Iterate thru the members
    for pos, i in enumerate([x[1] for x in position_order]):
      # Get the member name
      c_colname = H5Tget_member_name(type_id, i)
      colname = cstr_to_pystr(c_colname)

      # Get the member type
      member_type_id = H5Tget_member_type(type_id, i)
      # Get the member offset
      member_offset = H5Tget_member_offset(type_id, i)
      # Get the HDF5 class
      class_id = H5Tget_class(member_type_id)
      if class_id == H5T_COMPOUND and not is_complex(member_type_id):
        colpath2 = join_path(colpath, colname)
        # Create the native data in-memory
        itemsize = H5Tget_size(member_type_id)
        native_member_type_id = H5Tcreate(H5T_COMPOUND, itemsize)
        desc[colname], itemsize = self.get_nested_type(
          member_type_id, native_member_type_id, colpath2, field_byteorders)
        desc[colname]["_v_pos"] = pos
        desc[colname]["_v_offset"] = member_offset
      else:
        # Get the member format and the corresponding Col object
        try:
          native_member_type_id = get_native_type(member_type_id)
          atom = atom_from_hdf5_type(native_member_type_id)
          colobj = Col.from_atom(atom, pos=pos, _offset=member_offset)
          itemsize = H5Tget_size(native_member_type_id)
        except TypeError, te:
          # Re-raise TypeError again with more info
          raise TypeError(
            ("table ``%s``, column ``%s``: %%s" % (self.name, colname))
            % te.args[0])
        desc[colname] = colobj
        # For time kinds, save the byteorder of the column
        # (useful for conversion of time datatypes later on)
        if colobj.kind == "time":
          colobj._byteorder = H5Tget_order(member_type_id)
          if colobj._byteorder == H5T_ORDER_LE:
            field_byteorders.append("little")
          else:
            field_byteorders.append("big")
        elif colobj.kind in ['int', 'uint', 'float', 'complex', 'enum']:
          # Keep track of the byteorder for this column
          get_order(member_type_id, c_byteorder2)
          byteorder2 = cstr_to_pystr(c_byteorder2)
          if byteorder2 in ["little", "big"]:
            field_byteorders.append(byteorder2)

      # Insert the native member
      H5Tinsert(native_type_id, c_colname, member_offset, native_member_type_id)
      # Update the offset
      offset = offset + itemsize
      # Release resources
      H5Tclose(native_member_type_id)
      H5Tclose(member_type_id)
      H5free_memory(c_colname)

    # set the byteorder and other things (just in top level)
    if colpath == "":
      # Compute a byteorder for the entire table
      if len(field_byteorders) > 0:
        field_byteorders = np.array(field_byteorders)
        # Cython doesn't interpret well the extended comparison
        # operators so this: field_byteorders == "little" doesn't work
        # as expected
        if np.all(field_byteorders.__eq__("little")):
          byteorder = "little"
        elif np.all(field_byteorders.__eq__("big")):
          byteorder = "big"
        else:  # Yes! someone has done it!
          byteorder = "mixed"
      else:
        byteorder = "irrelevant"
      self.byteorder = byteorder

    return desc, offset

  def _get_info(self):
    """Get info from a table on disk."""

    cdef hid_t   space_id, plist
    cdef size_t  type_size, size2
    cdef hsize_t dims[1]        # enough for unidimensional tables
    cdef hsize_t chunksize[1]
    cdef H5D_layout_t layout
    cdef bytes encoded_name

    # Open the dataset
    encoded_name = self.name.encode('utf-8')
    self.dataset_id = H5Dopen(self.parent_id, encoded_name, H5P_DEFAULT)
    if self.dataset_id < 0:
      raise HDF5ExtError("Non-existing node ``%s`` under ``%s``" %
                         (self.name, self._v_parent._v_pathname))

    # Get the datatype on disk
    self.disk_type_id = H5Dget_type(self.dataset_id)
    if H5Tget_class(self.disk_type_id) != H5T_COMPOUND:
        raise ValueError("Node ``%s`` is not a Table object" %
                         (self._v_parent._v_leaves[self.name]._v_pathname))
    # Get the number of rows
    space_id = H5Dget_space(self.dataset_id)
    H5Sget_simple_extent_dims(space_id, dims, NULL)
    self.nrows = SizeType(dims[0])
    # Free resources
    H5Sclose(space_id)

    # Get the layout of the datatype
    plist = H5Dget_create_plist(self.dataset_id)
    layout = H5Pget_layout(plist)
    if layout == H5D_CHUNKED:
      self._chunked = 1
      # Get the chunksize
      H5Pget_chunk(plist, 1, chunksize)
    else:
      self._chunked = 0
      chunksize[0] = 0
    H5Pclose(plist)

    # Get the type size
    type_size = H5Tget_size(self.disk_type_id)
    # Create the native data in-memory
    self.type_id = H5Tcreate(H5T_COMPOUND, type_size)
    # Fill-up the (nested) native type and description
    desc, offset = self.get_nested_type(self.disk_type_id, self.type_id, "", [])

    if desc == {}:
      raise HDF5ExtError("Problems getting desciption for table %s", self.name)

    if offset < type_size:
      # Trailing padding, set the itemsize to the correct type_size (see #765)
      desc['_v_itemsize'] = type_size

    # Initialize blosc2 struct for chunk addresses
    self.chunk_op = chunk_iter_op(type_size, chunksize[0], NULL)

    # Return the object ID and the description
    return (self.dataset_id, desc, SizeType(chunksize[0]))

  cdef _convert_time64_(self, ndarray nparr, hsize_t nrecords, int sense):
    """Converts a NumPy of Time64 elements between NumPy and HDF5 formats.

    NumPy to HDF5 conversion is performed when 'sense' is 0.  Otherwise, HDF5
    to NumPy conversion is performed.  The conversion is done in place,
    i.e. 'nparr' is modified.

    """

    cdef void *t64buf
    cdef long byteoffset
    cdef npy_intp bytestride, nelements

    byteoffset = 0   # NumPy objects doesn't have an offset
    bytestride = PyArray_STRIDE(nparr, 0)  # supports multi-dimensional recarray
    # Compute the number of elements in the multidimensional cell
    nelements = nparr.size // len(nparr)
    t64buf = PyArray_DATA(nparr)

    conv_float64_timeval32(
      t64buf, byteoffset, bytestride, nrecords, nelements, sense)

  cpdef _convert_types(self, ndarray recarr, hsize_t nrecords, int sense):
    """Converts columns in 'recarr' between NumPy and HDF5 formats.

    NumPy to HDF5 conversion is performed when 'sense' is 0.  Otherwise, HDF5
    to NumPy conversion is performed.  The conversion is done in place,
    i.e. 'recarr' is modified.

    """

    # For reading, first swap the byteorder by hand
    # (this is not currently supported by HDF5)
    if sense == 1:
      for colpathname in self.colpathnames:
        if self.coltypes[colpathname] in ["time32", "time64"]:
          colobj = self.coldescrs[colpathname]
          if hasattr(colobj, "_byteorder"):
            if colobj._byteorder != platform_byteorder:
              column = get_nested_field(recarr, colpathname)
              # Do an *inplace* byteswapping
              column.byteswap(True)

    # This should be generalised to support other type conversions.
    for t64cname in self._time64colnames:
      column = get_nested_field(recarr, t64cname)
      self._convert_time64_(column, nrecords, sense)

  def _open_append(self, ndarray recarr):
    self._v_recarray = <object>recarr
    # Get the pointer to the buffer data area
    self.wbuf = PyArray_DATA(recarr)

  def _append_records(self, hsize_t nrecords):
    cdef int ret
    cdef hsize_t nrows

    # Clean address cache
    self._clean_chunk_addrs()

    # Convert some NumPy types to HDF5 before storing.
    self._convert_types(self._v_recarray, nrecords, 0)

    nrows = self.nrows
    # release GIL (allow other threads to use the Python interpreter)
    with nogil:
        # Append the records:
        ret = H5TBOappend_records(self.blosc2_support_write, self.dataset_id,
                                  self.type_id, nrows, nrecords, self.wbuf)

    if ret < 0:
      raise HDF5ExtError("Problems appending the records.")

    self.nrows = self.nrows + nrecords

  def _close_append(self):
    cdef hsize_t nrows

    if self._v_file.params['PYTABLES_SYS_ATTRS']:
      # Update the NROWS attribute
      nrows = self.nrows
      if (H5ATTRset_attribute(self.dataset_id, "NROWS", H5T_STD_I64,
                              0, NULL, <char *>&nrows) < 0):
        raise HDF5ExtError("Problems setting the NROWS attribute.")

    # Set the caches to dirty (in fact, and for the append case,
    # it should be only the caches based on limits, but anyway)
    self._dirtycache = True
    self._clean_chunk_addrs()
    # Delete the reference to recarray as we doesn't need it anymore
    self._v_recarray = None

  def _update_records(self, hsize_t start, hsize_t stop,
                      hsize_t step, ndarray recarr):
    cdef herr_t ret
    cdef void *rbuf
    cdef hsize_t nrecords, nrows

    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(recarr)

    # Compute the number of records to update
    nrecords = len(recarr)
    nrows = get_len_of_range(start, stop, step)
    if nrecords > nrows:
      nrecords = nrows

    # Convert some NumPy types to HDF5 before storing.
    self._convert_types(recarr, nrecords, 0)
    # Update the records:
    with nogil:
        ret = H5TBOwrite_records(self.blosc2_support_write and (step == 1), self.dataset_id,
                                 self.type_id, start, nrecords, step, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems updating the records.")

    # Set the caches to dirty
    self._dirtycache = True
    self._clean_chunk_addrs()


  def _update_elements(self, hsize_t nrecords, ndarray coords,
                       ndarray recarr):
    cdef herr_t ret
    cdef void *rbuf
    cdef void *rcoords

    # Get the chunk of the coords that correspond to a buffer
    rcoords = PyArray_DATA(coords)

    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(recarr)

    # Convert some NumPy types to HDF5 before storing.
    self._convert_types(recarr, nrecords, 0)

    # Update the records:
    with nogil:
        ret = H5TBOwrite_elements(self.dataset_id, self.type_id,
                                  nrecords, rcoords, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems updating the records.")

    # Set the caches to dirty
    self._dirtycache = True
    self._clean_chunk_addrs()


  def _read_records(self, hsize_t start, hsize_t nrecords, ndarray recarr):
    cdef void *rbuf
    cdef int ret
    cdef bytes fname = self._v_file.filename.encode('utf8')
    cdef char* filename = fname

    if self.blosc2_support_read:
      # Grab the addresses for the blosc2 frames (HDF5 chunks)
      nchunks = math.ceil(self.nrows / self.chunkshape[0])
      fill_chunk_addrs(self.dataset_id, nchunks, &self.chunk_op)

    # Correct the number of records to read, if needed
    if (start + nrecords) > self.nrows:
      nrecords = self.nrows - start

    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(recarr)

    # Read the records from disk
    with nogil:
        ret = H5TBOread_records(filename, self.blosc2_support_read, self.chunk_op,
                                self.dataset_id, self.type_id, start,
                                nrecords, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems reading records.")

    # Convert some HDF5 types to NumPy after reading.
    self._convert_types(recarr, nrecords, 1)

    return nrecords

  cdef hsize_t _read_chunk(self, hsize_t nchunk, ndarray iobuf, long cstart):
    cdef long nslot
    cdef hsize_t start, nrecords, chunkshape
    cdef int ret
    cdef void *rbuf
    cdef NumCache chunkcache
    cdef bytes fname = self._v_file.filename.encode('utf8')
    cdef char* filename = fname

    if self.blosc2_support_read:
      # Grab the addresses for the blosc2 frames (HDF5 chunks)
      nchunks = math.ceil(self.nrows / self.chunkshape[0])
      fill_chunk_addrs(self.dataset_id, nchunks, &self.chunk_op)

    chunkcache = self._chunkcache
    chunkshape = chunkcache.slotsize
    # Correct the number of records to read, if needed
    start = nchunk*chunkshape
    nrecords = chunkshape
    if (start + nrecords) > self.nrows:
      nrecords = self.nrows - start
    rbuf = PyArray_BYTES(iobuf) + cstart * chunkcache.itemsize
    # Try to see if the chunk is in cache
    nslot = chunkcache.getslot_(nchunk)
    if nslot >= 0:
      chunkcache.getitem_(nslot, rbuf, 0)
    else:
      # Chunk is not in cache. Read it and put it in the LRU cache.
      with nogil:
          ret = H5TBOread_records(filename, self.blosc2_support_read, self.chunk_op,
                                  self.dataset_id, self.type_id, start,
                                  nrecords, rbuf)

      if ret < 0:
        raise HDF5ExtError("Problems reading chunk records.")
      nslot = chunkcache.setitem_(nchunk, rbuf, 0)
    return nrecords

  def _read_elements(self, ndarray coords, ndarray recarr):
    cdef long nrecords
    cdef void *rbuf
    cdef void *rbuf2
    cdef int ret

    # Get the chunk of the coords that correspond to a buffer
    nrecords = coords.size
    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(recarr)
    # Get the pointer to the buffer coords area
    rbuf2 = PyArray_DATA(coords)

    with nogil:
        ret = H5TBOread_elements(self.dataset_id, self.type_id,
                                 nrecords, rbuf2, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems reading records.")

    # Convert some HDF5 types to NumPy after reading.
    self._convert_types(recarr, nrecords, 1)

    return nrecords

  def _remove_rows(self, hsize_t start, hsize_t stop, long step):
    cdef size_t rowsize
    cdef hsize_t nrecords=0, nrecords2
    cdef hsize_t i
    cdef bytes fname = self._v_file.filename.encode('utf8')
    cdef char* filename = fname

    if step == 1:
      nrecords = stop - start
      rowsize = self.rowsize
      # Using self.disk_type_id should be faster (i.e. less conversions)
      if (H5TBOdelete_records(filename, self.blosc2_support_read, self.chunk_op,
                              self.dataset_id,
                              self.disk_type_id, self.nrows, rowsize,
                              start, nrecords, self.nrowsinbuf) < 0):
        raise HDF5ExtError("Problems deleting records.")

      self.nrows = self.nrows - nrecords
      if self._v_file.params['PYTABLES_SYS_ATTRS']:
        # Attach the NROWS attribute
        nrecords2 = self.nrows
        H5ATTRset_attribute(self.dataset_id, "NROWS", H5T_STD_I64,
                            0, NULL, <char *>&nrecords2)
      # Set the caches to dirty
      self._dirtycache = True
      self._clean_chunk_addrs()
    elif step == -1:
      nrecords = self._remove_rows(stop+1, start+1, 1)
    elif step >= 1:
      # always want to go through the space backwards
      for i in range(stop - step, <ssize_t>start - step, -step):
        nrecords += self._remove_rows(i, i+1, 1)
    elif step <= -1:
      # always want to go through the space backwards
      for i in range(start, stop, step):
        nrecords += self._remove_rows(i, i+1, 1)
    else:
      raise ValueError("step size may not be 0.")

    # Return the number of records removed
    return nrecords

  # Clean address cache
  def _clean_chunk_addrs(self):
    clean_chunk_addrs(&self.chunk_op)


cdef class Row:
  """Table row iterator and field accessor.

  Instances of this class are used to fetch and set the values
  of individual table fields.  It works very much like a dictionary,
  where keys are the pathnames or positions (extended slicing is
  supported) of the fields in the associated table in a specific row.

  This class provides an *iterator interface*
  so that you can use the same Row instance to
  access successive table rows one after the other.  There are also
  some important methods that are useful for accessing, adding and
  modifying values in tables.

  .. rubric:: Row attributes

  .. attribute:: nrow

      The current row number.

      This property is useful for knowing which row is being dealt with in the
      middle of a loop or iterator.

  """

  cdef npy_intp _stride
  cdef long _row, _unsaved_nrows, _mod_nrows
  cdef long long start, absstep
  cdef long long stop, step, nextelement, _nrow, stopb  # has to be long long, not hsize_t, for negative step sizes
  cdef long long nrowsinbuf, nrows, nrowsread
  cdef long long chunksize, nchunksinbuf, totalchunks
  cdef long long startb, lenbuf
  cdef long long indexchunk
  cdef int     bufcounter, counter
  cdef int     exist_enum_cols
  cdef int     _riterator, _rowsize, _write_to_seqcache
  cdef int     wherecond, indexed
  cdef int     ro_filemode, chunked
  cdef int     _bufferinfo_done, sss_on
  cdef long long iterseq_max_elements
  cdef ndarray bufcoords, indexvalid, indexvalues, chunkmap
  cdef hsize_t *bufcoords_data
  cdef hsize_t *index_values_data
  cdef char    *chunkmap_data
  cdef char    *index_valid_data
  cdef object  dtype
  cdef object  iobuf, iobufcpy
  cdef object  wrec, wreccpy
  cdef object  wfields, rfields
  cdef object  coords
  cdef object  condfunc, condargs, condkwargs
  cdef object  mod_elements, colenums
  cdef object  rfieldscache, wfieldscache
  cdef object  iterseq
  cdef object  _table_file, _table_path
  cdef object  modified_fields
  cdef object  seqcache_key

  # The nrow() method has been converted into a property, which is handier
  property nrow:
    """The current row number.

    This property is useful for knowing which row is being dealt with in the
    middle of a loop or iterator.

    """

    def __get__(self):
      return SizeType(self._nrow)

  property table:
    def __get__(self):
        self._table_file._check_open()
        return self._table_file._get_node(self._table_path)

  def __cinit__(self, table):
    cdef int nfields, i
    # Location-dependent information.
    self._table_file = table._v_file
    self._table_path = table._v_pathname
    self._unsaved_nrows = 0
    self._mod_nrows = 0
    self._row = 0
    self._nrow = 0   # Useful in mod_append read iterators
    self._riterator = 0
    self._bufferinfo_done = 0
    # Some variables from table will be cached here
    if table._v_file.mode == 'r':
      self.ro_filemode = 1
    else:
      self.ro_filemode = 0
    self.chunked = table._chunked
    self.colenums = table._colenums
    self.exist_enum_cols = len(self.colenums)
    self.nrowsinbuf = table.nrowsinbuf
    self.chunksize = table.chunkshape[0]
    self.nchunksinbuf = self.nrowsinbuf // self.chunksize
    self.dtype = table._v_dtype
    self._new_buffer(table)
    self.mod_elements = None
    self.rfieldscache = {}
    self.wfieldscache = {}
    self.modified_fields = set()

  def _iter(self, start=0, stop=0, step=1, coords=None, chunkmap=None):
    """Return an iterator for traversiong the data in table."""
    self._init_loop(start, stop, step, coords, chunkmap)
    return iter(self)

  def __iter__(self):
    """Iterator that traverses all the data in the Table"""
    return self

  cdef _new_buffer(self, table):
    """Create the recarrays for I/O buffering"""

    wdflts = table._v_wdflts
    if wdflts is None:
      self.wrec = np.zeros(1, dtype=self.dtype)  # Defaults are zero
    else:
      self.wrec = table._v_wdflts.copy()
    self.wreccpy = self.wrec.copy()  # A copy of the defaults
    # Build the wfields dictionary for faster access to columns
    self.wfields = {}
    for name in self.dtype.names:
      self.wfields[name] = self.wrec[name]

    # Get the read buffer for this instance (it is private, remember!)
    buff = self.iobuf = table._get_container(self.nrowsinbuf)
    # Build the rfields dictionary for faster access to columns
    # This is quite fast, as it only takes around 5 us per column
    # in my laptop (Pentium 4 @ 2 GHz).
    # F. Alted 2006-08-18
    self.rfields = {}
    for i, name in enumerate(self.dtype.names):
      self.rfields[i] = buff[name]
      self.rfields[name] = buff[name]

    # Get the stride of these buffers
    self._stride = PyArray_STRIDE(buff, 0)
    # The rowsize
    self._rowsize = self.dtype.itemsize
    self.nrows = table.nrows  # This value may change

  cdef _init_loop(self, long long start, long long stop, long long step,
                 object coords, object chunkmap):
    """Initialization for the __iter__ iterator"""
    cdef Table table = self.table
    self._riterator = 1   # We are inside a read iterator
    self.start = start
    self.stop = stop
    self.step = step
    self.coords = coords
    self.startb = 0
    if step > 0:
      self._row = -1  # a sentinel
      self.nrowsread = start
    elif step < 0:
      self._row = 0
      self.nrowsread = 0
      self.nextelement = start
    self._nrow = start - self.step
    self.wherecond = 0
    self.indexed = 0
    self.nrows = table.nrows   # Update the row counter
    if table.blosc2_support_read:
      # Grab the addresses for the blosc2 frames (HDF5 chunks)
      nchunks = math.ceil(self.nrows / self.table.chunkshape[0])
      fill_chunk_addrs(table.dataset_id, nchunks, &table.chunk_op)

    if coords is not None and 0 < step:
      self.nrowsread = start
      self.nextelement = start
      self.stop = min(stop, len(coords))
      self.absstep = abs(step)
      return
    elif coords is not None and 0 > step:
      #self.nrowsread = 0
      #self.nextelement = start
      #self.stop = min(stop, len(coords))
      #self.stop = max(stop, start - len(coords))
      self.absstep = abs(step)
      return

    if table._where_condition:
      self.wherecond = 1
      #self.condkwargs = {'ex_uses_vml': True}
      self.condfunc, self.condargs, self.condkwargs = table._where_condition
      table._where_condition = None

    if table._use_index:
      # Indexing code depends on this condition (see #319)
      assert self.nrowsinbuf % self.chunksize == 0
      self.indexed = 1
      # Compute totalchunks here because self.nrows can change during the
      # life of a Row instance.
      self.totalchunks = self.nrows // self.chunksize
      if self.nrows % self.chunksize:
        self.totalchunks = self.totalchunks + 1
      self.nrowsread = 0
      self.nextelement = 0
      self.chunkmap = chunkmap
      self.chunkmap_data = PyArray_BYTES(self.chunkmap)
      table._use_index = False
      self.lenbuf = self.nrowsinbuf
      # Check if we have limitations on start, stop, step
      self.sss_on = (self.start > 0 or self.stop < self.nrows or self.step > 1)

    self.seqcache_key = table._seqcache_key
    table._seqcache_key = None
    if self.seqcache_key is not None:
      self._write_to_seqcache = 1
      self.iterseq_max_elements = table._v_file.params['ITERSEQ_MAX_ELEMENTS']
      self.iterseq = [] # all the row indexes, unless it would be longer than ITERSEQ_MAX_ELEMENTS
    else:
      self._write_to_seqcache = 0
      self.iterseq = None

  def __next__(self):
    """next() method for __iter__() that is called on each iteration"""

    if not self._riterator:
      # The iterator is already exhausted!
      raise StopIteration
    if self.indexed:
      return self.__next__indexed()
    elif self.coords is not None:
      return self.__next__coords()
    elif self.wherecond:
      return self.__next__inkernel()
    else:
      return self.__next__general()

  cdef __next__indexed(self):
    """The version of next() for indexed columns and a chunkmap."""

    cdef long recout, j, cs, vlen, rowsize
    cdef long long nchunksread
    cdef object tmp_range
    cdef Table table
    cdef ndarray iobuf
    cdef void *IObufData
    cdef long nslot
    cdef object seq
    cdef object seqcache

    assert self.nrowsinbuf >= self.chunksize
    while self.nextelement < self.stop:
      if self.nextelement >= self.nrowsread:
        # Skip until there is interesting information
        while self.start > self.nrowsread + self.nrowsinbuf:
          self.nrowsread = self.nrowsread + self.nrowsinbuf
          self.nextelement = self.nextelement + self.nrowsinbuf

        table = self.table
        iobuf = self.iobuf
        j = 0;  recout = 0;  cs = self.chunksize
        nchunksread = self.nrowsread // cs
        tmp_range = np.arange(0, cs, dtype='int64')
        self.bufcoords = np.empty(self.nrowsinbuf, dtype='int64')
        # Fetch valid chunks until the I/O buffer is full
        while nchunksread < self.totalchunks:
          if self.chunkmap_data[nchunksread]:
            self.bufcoords[j*cs:(j+1)*cs] = tmp_range + self.nrowsread
            # Not optimized read
            #  recout = recout + table._read_records(
            #    nchunksread*cs, cs, iobuf[j*cs:])
            #
            # Optimized read through the use of a chunk cache.  This cache has
            # more or less the same speed than the integrated HDF5 chunk
            # cache, but using the PyTables one has the advantage that the
            # user can easily change this parameter.
            recout = recout + table._read_chunk(nchunksread, iobuf, j*cs)
            j = j + 1
          self.nrowsread = (nchunksread+1)*cs
          if self.nrowsread > self.stop:
            self.nrowsread = self.stop
            break
          elif j == self.nchunksinbuf:
            break
          nchunksread = nchunksread + 1

        # Evaluate the condition on this table fragment.
        iobuf = iobuf[:recout]

        if len(iobuf) > 0:
          self.table._convert_types(iobuf, len(iobuf), 1)
        self.indexvalid = call_on_recarr(
          self.condfunc, self.condargs, iobuf, **self.condkwargs)
        self.index_valid_data = PyArray_BYTES(self.indexvalid)
        # Get the valid coordinates
        self.indexvalues = self.bufcoords[:recout][self.indexvalid]
        self.index_values_data = <hsize_t *>PyArray_DATA(self.indexvalues)
        self.lenbuf = self.indexvalues.size
        # Place the valid results at the beginning of the buffer
        iobuf[:self.lenbuf] = iobuf[self.indexvalid]

        # Initialize the internal buffer row counter
        self._row = -1

        if self._write_to_seqcache:
          # Feed the indexvalues into the seqcache
          seqcache = self.iterseq
          if self.lenbuf + len(seqcache) < self.iterseq_max_elements:
            seqcache.extend(self.indexvalues)
          else:
            self.iterseq = None
            self._write_to_seqcache = 0

      self._row = self._row + 1
      # Check whether we have read all the rows in buf
      if self._row == self.lenbuf:
        self.nextelement = self.nrowsread
        # Make _row to point to the last valid entry in buffer
        # (this is useful for accessing the last row after an iterator loop)
        self._row = self._row - 1
        continue
      self._nrow = self.index_values_data[self._row]
      # Check additional conditions on start, stop, step params
      if self.sss_on:
        if (self._nrow < self.start or self._nrow >= self.stop):
          self.nextelement = self.nextelement + 1
          continue
        if (self.step > 1 and
            ((self._nrow - self.start) % self.step > 0)):
          self.nextelement = self.nextelement + 1
          continue
      # Return this row
      self.nextelement = self._nrow + 1
      return self
    else:
      # All the elements have been read for this mode
      self._finish_riterator()

  cdef __next__coords(self):
    """The version of next() for user-required coordinates"""
    cdef int recout
    cdef long long lenbuf, nextelement
    cdef object tmp
    if 0 < self.step:
      while self.nextelement < self.stop:
        if self.nextelement >= self.nrowsread:
          # Correction for avoiding reading past self.stop
          if self.nrowsread+self.nrowsinbuf > self.stop:
            lenbuf = self.stop-self.nrowsread
          else:
            lenbuf = self.nrowsinbuf
          tmp = self.coords[self.nrowsread:self.nrowsread+lenbuf:self.step]
          # We have to get a contiguous buffer, so numpy.array is the way to go
          self.bufcoords = np.array(tmp, dtype="uint64")
          self._row = -1
          if self.bufcoords.size > 0:
            recout = self.table._read_elements(self.bufcoords, self.iobuf)
          else:
            recout = 0
          self.bufcoords_data = <hsize_t*>PyArray_DATA(self.bufcoords)
          self.nrowsread = self.nrowsread + lenbuf
          if recout == 0:
            # no items were read, skip out
            continue
        self._row = self._row + 1
        self._nrow = self.bufcoords_data[self._row]
        self.nextelement = self.nextelement + self.absstep
        return self
      else:
        # All the elements have been read for this mode
        self._finish_riterator()
    elif 0 > self.step:
      #print("self.nextelement = ", self.nextelement, self.start, self.nrowsread, self.nextelement <  self.start - self.nrowsread + 1)
      while self.nextelement > self.stop:
        if self.nextelement < self.start - (<long long> self.nrowsread) + 1:
          if 0 > self.nextelement - (<long long> self.nrowsinbuf) + 1:
            tmp = self.coords[0:self.nextelement + 1]
          else:
            tmp = self.coords[self.nextelement - (<long long> self.nrowsinbuf) + 1:self.nextelement + 1]
          self.bufcoords = np.array(tmp, dtype="uint64")
          recout = self.table._read_elements(self.bufcoords, self.iobuf)
          self.bufcoords_data = <hsize_t*>PyArray_DATA(self.bufcoords)
          self.nrowsread = self.nrowsread + self.nrowsinbuf
          self._row = len(self.bufcoords) - 1
        else:
          self._row = (self._row + self.step) % len(self.bufcoords)

        self._nrow = self.nextelement - self.step
        self.nextelement = self.nextelement + self.step
        # Return this value
        return self
      else:
        # All the elements have been read for this mode
        self._finish_riterator()
    else:
      self._finish_riterator()

  cdef __next__inkernel(self):
    """The version of next() in case of in-kernel conditions"""

    cdef hsize_t recout, correct
    cdef object numexpr_locals, colvar, col
    self.nextelement = self._nrow + self.step
    while self.nextelement < self.stop:
      if self.nextelement >= self.nrowsread:
        # Skip until there is interesting information
        while self.nextelement >= self.nrowsread + self.nrowsinbuf:
          self.nrowsread = self.nrowsread + self.nrowsinbuf
        # Compute the end for this iteration
        self.stopb = self.stop - self.nrowsread
        if self.stopb > self.nrowsinbuf:
          self.stopb = self.nrowsinbuf
        self._row = self.startb - self.step
        # Read a chunk
        recout = self.table._read_records(self.nextelement, self.nrowsinbuf,
                                          self.iobuf)
        self.nrowsread = self.nrowsread + recout
        self.indexchunk = -self.step

        # Evaluate the condition on this table fragment.
        self.indexvalid = call_on_recarr(
          self.condfunc, self.condargs, self.iobuf[:recout], **self.condkwargs)
        self.index_valid_data = PyArray_BYTES(self.indexvalid)

        # Is there any interesting information in this buffer?
        if not np.any(self.indexvalid):
          # No, so take the next one
          if self.step >= self.nrowsinbuf:
            self.nextelement = self.nextelement + self.step
          else:
            self.nextelement = self.nextelement + self.nrowsinbuf
            # Correction for step size > 1
            if self.step > 1:
              correct = (self.nextelement - self.start) % self.step
              self.nextelement = self.nextelement - correct
          continue

      self._row = self._row + self.step
      self._nrow = self.nextelement
      if self._row + self.step >= self.stopb:
        # Compute the start row for the next buffer
        self.startb = 0

      self.nextelement = self._nrow + self.step
      # Return only if this value is interesting
      self.indexchunk = self.indexchunk + self.step
      if self.index_valid_data[self.indexchunk]:
        return self
    else:
      self._finish_riterator()

  cdef __next__general(self):
    """The version of next() for the general cases"""
    cdef int recout
    if 0 < self.step:
      self.nextelement = self._nrow + self.step
      while self.nextelement < self.stop:
        if self.nextelement >= self.nrowsread:
          # Skip until there is interesting information
          while self.nextelement >= self.nrowsread + self.nrowsinbuf:
            self.nrowsread = self.nrowsread + self.nrowsinbuf
          # Compute the end for this iteration
          self.stopb = self.stop - self.nrowsread
          if self.stopb > self.nrowsinbuf:
            self.stopb = self.nrowsinbuf
          self._row = self.startb - self.step
          # Read a chunk
          recout = self.table._read_records(self.nrowsread, self.nrowsinbuf,
                                            self.iobuf)
          self.nrowsread = self.nrowsread + recout

        self._row = self._row + self.step
        self._nrow = self.nextelement
        if self._row + self.step >= self.stopb:
          # Compute the start row for the next buffer
          self.startb = (self._row + self.step) % self.nrowsinbuf

        self.nextelement = self._nrow + self.step
        # Return this value
        return self
      else:
        self._finish_riterator()
    elif 0 > self.step:
      self.stopb = -1
      while self.nextelement - 1 > self.stop:
        if self.nextelement < self.start - self.nrowsread + 1:
          # Read a chunk
          recout = self.table._read_records(self.nextelement - self.nrowsinbuf + 1,
                                            self.nrowsinbuf, self.iobuf)
          self.nrowsread = self.nrowsread + self.nrowsinbuf
          self._row = self.nrowsinbuf - 1
        else:
          self._row = (self._row + self.step) % self.nrowsinbuf

        self._nrow = self.nextelement - self.step
        self.nextelement = self.nextelement + self.step
        # Return this value
        return self
      else:
        self._finish_riterator()

  cdef _finish_riterator(self):
    """Clean-up things after iterator has been done"""
    cdef ObjectCache seqcache
    cdef Table table = self.table

    self.rfieldscache = {}     # empty rfields cache
    self.wfieldscache = {}     # empty wfields cache
    # Make a copy of the last read row in the private record
    # (this is useful for accessing the last row after an iterator loop)
    if self._row >= 0:
      self.wrec[:] = self.iobuf[self._row]
    if self._write_to_seqcache:
      seqcache = self.table._seqcache
      # Guessing iterseq size: Each element in self.iterseq should take at least 8 bytes
      seqcache.setitem_(self.seqcache_key, self.iterseq, len(self.iterseq) * 8)
    self._riterator = 0        # out of iterator
    self.iterseq = None        # empty seqcache-related things
    self.seqcache_key = None
    if self._mod_nrows > 0:    # Check if there is some modified row
      self._flush_mod_rows()     # Flush any possible modified row
    self.modified_fields = set()  # Empty the set of modified fields
    raise StopIteration        # end of iteration

  def _fill_col(self, result, start, stop, step, field):
    """Read a field from a table on disk and put the result in result"""

    cdef hsize_t startr, istartb
    cdef long long istart, inrowsinbuf, inextelement
    cdef long long stopr, istopb, i, j, inrowsread
    cdef long long istop, istep
    cdef object fields

    # We can't reuse existing buffers in this context
    self._init_loop(start, stop, step, None, None)
    istart, istop, istep = self.start, self.stop, self.step
    inrowsinbuf, inextelement, inrowsread = self.nrowsinbuf, istart, istart
    istartb, startr = self.startb, 0
    i = istart
    if 0 < istep:
      while i < istop:
        if (inextelement >= inrowsread + inrowsinbuf):
          inrowsread = inrowsread + inrowsinbuf
          i = i + inrowsinbuf
          continue
        # Compute the end for this iteration
        istopb = istop - inrowsread
        if istopb > inrowsinbuf:
          istopb = inrowsinbuf
        stopr = startr + ((istopb - istartb - 1) // istep) + 1
        # Read a chunk
        inrowsread = inrowsread + self.table._read_records(i, inrowsinbuf,
                                                           self.iobuf)
        # Assign the correct part to result
        fields = self.iobuf
        if field:
          fields = get_nested_field(fields, field)
        result[startr:stopr] = fields[istartb:istopb:istep]

        # Compute some indexes for the next iteration
        startr = stopr
        j = istartb + ((istopb - istartb - 1) // istep) * istep
        istartb = (j+istep) % inrowsinbuf
        inextelement = inextelement + istep
        i = i + inrowsinbuf
    elif istep < 0:
      inrowsinbuf = self.nrowsinbuf
      #istartb = self.startb
      istartb = self.nrowsinbuf - 1
      #istopb = self.stopb - 1
      istopb = -1
      startr = 0
      i = istart
      inextelement = istart
      inrowsread = 0
      while i-1 > istop:
        #if (inextelement <= inrowsread + inrowsinbuf):
        if (inextelement < i - inrowsinbuf):
          inrowsread = inrowsread + inrowsinbuf
          i = i - inrowsinbuf
          continue
        # Compute the end for this iteration
        # (we know we are going backward so try to keep indices positive)
        stopr = startr + (1 - istopb + istartb) // (-istep)
        # Read a chunk
        inrowsread = inrowsread + self.table._read_records(i - inrowsinbuf + 1,
                                                           inrowsinbuf, self.iobuf)
        # Assign the correct part to result
        fields = self.iobuf
        if field:
          fields = get_nested_field(fields, field)
        if istopb >= 0:
            result[startr:stopr] = fields[istartb:istopb:istep]
        else:
            result[startr:stopr] = fields[istartb::istep]

        # Compute some indexes for the next iteration
        startr = stopr
        istartb = (i - istartb)%inrowsinbuf
        inextelement = inextelement + istep
        i = i - inrowsinbuf
    self._riterator = 0  # out of iterator
    return


  def append(self):
    """Add a new row of data to the end of the dataset.

    Once you have filled the proper fields for the current
    row, calling this method actually appends the new data to the
    *output buffer* (which will eventually be
    dumped to disk).  If you have not set the value of a field, the
    default value of the column will be used.

    .. warning::

        After completion of the loop in which :meth:`Row.append` has
        been called, it is always convenient to make a call to
        :meth:`Table.flush` in order to avoid losing the last rows that
        may still remain in internal buffers.

    Examples
    --------

    ::

        row = table.row
        for i in xrange(nrows):
            row['col1'] = i-1
            row['col2'] = 'a'
            row['col3'] = -1.0
            row.append()
        table.flush()

    """
    cdef ndarray iobuf, wrec, wreccpy

    if self.ro_filemode:
      raise IOError("Attempt to write over a file opened in read-only mode")

    if not self.chunked:
      raise HDF5ExtError("You cannot append rows to a non-chunked table.",
                         h5tb=False)

    if self._riterator:
      raise NotImplementedError("You cannot append rows when in middle of a table iterator. If what you want is to update records, use Row.update() instead.")

    # Commit the private record into the write buffer
    # self.iobuf[self._unsaved_nrows] = self.wrec
    # The next is faster
    iobuf = <ndarray>self.iobuf; wrec = <ndarray>self.wrec
    memcpy(PyArray_BYTES(iobuf) + self._unsaved_nrows * self._stride,
           PyArray_BYTES(wrec), self._rowsize)
    # Restore the defaults for the private record
    # self.wrec[:] = self.wreccpy
    # The next is faster
    wreccpy = <ndarray>self.wreccpy
    memcpy(PyArray_BYTES(wrec), PyArray_BYTES(wreccpy), self._rowsize)
    self._unsaved_nrows = self._unsaved_nrows + 1
    # When the buffer is full, flush it
    if self._unsaved_nrows == self.nrowsinbuf:
      self._flush_buffered_rows()

  def _flush_buffered_rows(self):
    if self._unsaved_nrows > 0:
      self.table._save_buffered_rows(self.iobuf, self._unsaved_nrows)
      # Reset the buffer unsaved counter
      self._unsaved_nrows = 0


  def _get_unsaved_nrows(self):
    return self._unsaved_nrows


  def update(self):
    """Change the data of the current row in the dataset.

    This method allows you to modify values in a table when you are in the
    middle of a table iterator like :meth:`Table.iterrows` or
    :meth:`Table.where`.

    Once you have filled the proper fields for the current row, calling
    this method actually changes data in the *output buffer* (which will
    eventually be dumped to disk).  If you have not set the value of a
    field, its original value will be used.

    .. warning::

        After completion of the loop in which :meth:`Row.update` has
        been called, it is always convenient to make a call to
        :meth:`Table.flush` in order to avoid losing changed rows that
        may still remain in internal buffers.

    Examples
    --------

    ::

        for row in table.iterrows(step=10):
            row['col1'] = row.nrow
            row['col2'] = 'b'
            row['col3'] = 0.0
            row.update()
        table.flush()

    which modifies every tenth row in table.  Or::

        for row in table.where('col1 > 3'):
            row['col1'] = row.nrow
            row['col2'] = 'b'
            row['col3'] = 0.0
            row.update()
        table.flush()

    which just updates the rows with values bigger than 3 in the first
    column.

    """

    cdef ndarray iobufcpy, iobuf

    if self.ro_filemode:
      raise IOError("Attempt to write over a file opened in read-only mode")

    if not self._riterator:
      raise NotImplementedError("You are only allowed to update rows through the Row.update() method if you are in the middle of a table iterator.")

    if self.mod_elements is None:
      # Initialize an array for keeping the modified elements
      # (just in case Row.update() would be used)
      self.mod_elements = np.empty(shape=self.nrowsinbuf, dtype=SizeType)
      # We need a different copy for self.iobuf here
      self.iobufcpy = self.iobuf.copy()

    # Add this row to the list of elements to be modified
    self.mod_elements[self._mod_nrows] = self._nrow
    # Copy the current buffer row in input to the output buffer
    # self.iobufcpy[self._mod_nrows] = self.iobuf[self._row]
    # The next is faster
    iobufcpy = <ndarray>self.iobufcpy; iobuf = <ndarray>self.iobuf
    memcpy(PyArray_BYTES(iobufcpy) + self._mod_nrows * self._stride,
           PyArray_BYTES(iobuf) + self._row * self._stride, self._rowsize)
    # Increase the modified buffer count by one
    self._mod_nrows = self._mod_nrows + 1
    # No point writing seqcache -- Table.flush will invalidate it
    # since we no longer know whether this row will meet _where_condition
    self._write_to_seqcache = 0
    # When the buffer is full, flush it
    if self._mod_nrows == self.nrowsinbuf:
      self._flush_mod_rows()

  def _flush_mod_rows(self):
    """Flush any possible modified row using Row.update()"""

    table = self.table
    # Save the records on disk
    table._update_elements(self._mod_nrows, self.mod_elements, self.iobufcpy)
    # Reset the counter of modified rows to 0
    self._mod_nrows = 0
    # Mark the modified fields' indexes as dirty.
    table._mark_columns_as_dirty(self.modified_fields)


  def __contains__(self, item):
    """__contains__(item)

    A true value is returned if item is found in current row, false
    otherwise.

    """

    return item in self.fetch_all_fields()

  # This method is twice as faster than __getattr__ because there is
  # not a lookup in the local dictionary
  def __getitem__(self, key):
    """__getitem__(key)

    Get the row field specified by the `key`.

    The key can be a string (the name of the field), an integer (the
    position of the field) or a slice (the range of field positions). When
    key is a slice, the returned value is a *tuple* containing the values
    of the specified fields.

    Examples
    --------

    ::

        res = [row['var3'] for row in table.where('var2 < 20')]

    which selects the var3 field for all the rows that fulfil the
    condition. Or::

        res = [row[4] for row in table if row[1] < 20]

    which selects the field in the *4th* position for all the rows that
    fulfil the condition. Or::

        res = [row[:] for row in table if row['var2'] < 20]

    which selects the all the fields (in the form of a *tuple*) for all the
    rows that fulfil the condition. Or::

        res = [row[1::2] for row in table.iterrows(2, 3000, 3)]

    which selects all the fields in even positions (in the form of a
    *tuple*) for all the rows in the slice [2:3000:3].

    """

    cdef long offset
    cdef ndarray field
    cdef object row, fields, fieldscache

    if self._riterator:
      # If in the middle of an iterator loop, the user probably wants to
      # access the read buffer
      fieldscache = self.rfieldscache; fields = self.rfields
      offset = <long>self._row
    else:
      # We are not in an iterator loop, so the user probably wants to access
      # the write buffer
      fieldscache = self.wfieldscache; fields = self.wfields
      offset = 0

    try:
      # Check whether this object is in the cache dictionary
      field = fieldscache[key]
    except (KeyError, TypeError):
      try:
        # Try to get it from fields (str or int keys)
        field = get_nested_field_cache(fields, key, fieldscache)
      except TypeError:
        # No luck yet. Still, the key can be a slice.
        # Fetch the complete row and convert it into a tuple
        if self._riterator:
          row = self.iobuf[self._row].copy().item()
        else:
          row = self.wrec[0].copy().item()
        # Try with __getitem__()
        return row[key]

    if PyArray_NDIM(field) == 1:
      # For an scalar it is not needed a copy (immutable object)
      return PyArray_GETITEM(field, PyArray_BYTES(field) + offset * self._stride)
    else:
      # Do a copy of the array, so that it can be overwritten by the user
      # without damaging the internal self.rfields buffer
      return field[offset].copy()

  # This is slightly faster (around 3%) than __setattr__
  def __setitem__(self, object key, object value):
    """__setitem__(key, value)

    Set the key row field to the specified value.

    Differently from its __getitem__() counterpart, in this case key can
    only be a string (the name of the field). The changes done via
    __setitem__() will not take effect on the data on disk until any of the
    :meth:`Row.append` or :meth:`Row.update` methods are called.

    Examples
    --------

    ::

        for row in table.iterrows(step=10):
            row['col1'] = row.nrow
            row['col2'] = 'b'
            row['col3'] = 0.0
            row.update()
        table.flush()

    which modifies every tenth row in the table.

    """

    cdef int ret
    cdef long offset
    cdef ndarray field
    cdef object fields, fieldscache

    if self.ro_filemode:
      raise IOError("attempt to write over a file opened in read-only mode")

    if self._riterator:
      # If in the middle of an iterator loop, or *after*, the user
      # probably wants to access the read buffer
      fieldscache = self.rfieldscache; fields = self.rfields
      offset = <long>self._row
    else:
      # We are not in an iterator loop, so the user probably wants to access
      # the write buffer
      fieldscache = self.wfieldscache; fields = self.wfields
      offset = 0

    # Check validity of enumerated value.
    if self.exist_enum_cols:
      if key in self.colenums:
        enum = self.colenums[key]
        for cenval in np.asarray(value).flat:
          enum(cenval)  # raises ``ValueError`` on invalid values

    # Get the field to be modified
    field = get_nested_field_cache(fields, key, fieldscache)
    if key not in self.modified_fields:
      self.modified_fields.add(key)

    # Finally, try to set it to the value
    try:
      # Optimization for scalar values. This can optimize the writes
      # between a 10% and 100%, depending on the number of columns modified
      if PyArray_NDIM(field) == 1:
        ret = PyArray_SETITEM(field, PyArray_BYTES(field) + offset * self._stride, value)
        if ret < 0:
          PyErr_Clear()
          raise TypeError
      ##### End of optimization for scalar values
      else:
        field[offset] = value
    except TypeError:
      raise TypeError("invalid type (%s) for column ``%s``" % (type(value),
                                                               key))

  def fetch_all_fields(self):
    """Retrieve all the fields in the current row.

    Contrarily to row[:] (see :ref:`RowSpecialMethods`), this returns row
    data as a NumPy void scalar.  For instance::

        [row.fetch_all_fields() for row in table.where('col1 < 3')]

    will select all the rows that fulfill the given condition
    as a list of NumPy records.

    """

    # We need to do a cast for recognizing negative row numbers!
    if <signed long long>self._nrow < 0:
      return ("Warning: Row iterator has not been initialized for table:\n"
              "  %s\n"
              " You will normally want to use this method in iterator "
              "contexts." % self.table)

    # Always return a copy of the row so that new data that is written
    # in self.iobuf doesn't overwrite the original returned data.
    return self.iobuf[self._row].copy()

  def __str__(self):
    """Represent the record as an string"""

    # We need to do a cast for recognizing negative row numbers!
    if <signed long long>self._nrow < 0:
      return ("Warning: Row iterator has not been initialized for table:\n"
              "  %s\n"
              " You will normally want to use this object in iterator "
              "contexts." % self.table)

    tablepathname = self.table._v_pathname
    classname = self.__class__.__name__
    return "%s.row (%s), pointing to row #%d" %  (tablepathname, classname,
                                                  self._nrow)

  def __repr__(self):
    """Represent the record as an string"""

    return str(self)

## Local Variables:
## mode: python
## py-indent-offset: 2
## tab-width: 2
## fill-column: 78
## End:
