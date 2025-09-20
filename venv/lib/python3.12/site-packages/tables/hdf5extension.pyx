########################################################################
#
# License: BSD
# Created: September 21, 2002
# Author:  Francesc Alted - faltet@pytables.com
#
# $Id$
#
########################################################################

"""Cython interface between several PyTables classes and HDF5 library.

Classes (type extensions):

    File
    AttributeSet
    Node
    Leaf
    Group
    Array
    VLArray
    UnImplemented

Functions:

Misc variables:

"""

import os
import sys
import platform
import warnings
from collections import namedtuple

ObjInfo = namedtuple('ObjInfo', ['addr', 'rc'])
ObjTimestamps = namedtuple('ObjTimestamps', ['atime', 'mtime',
                                             'ctime', 'btime'])

import pickle

import numpy as np

from .atom import Atom
from .utils import check_file_access, byteorders, correct_byteorder, SizeType
from .exceptions import HDF5ExtError, DataTypeWarning
from .description import descr_from_dtype
from .utilsextension import (
    encode_filename,
    set_blosc_max_threads,
    set_blosc2_max_threads,
    atom_to_hdf5_type,
    atom_from_hdf5_type,
    hdf5_to_np_ext_type,
    create_nested_type,
    pttype_to_hdf5,
    pt_special_kinds,
    npext_prefixes_to_ptkinds,
    hdf5_class_to_string,
    platform_byteorder,
    get_filters,
)

# Types, constants, functions, classes & other objects from everywhere

from numpy cimport (
    import_array,
    ndarray,
    npy_intp,
    PyArray_BYTES,
    PyArray_DATA,
    PyArray_DIMS,
    PyArray_NDIM,
    PyArray_STRIDE,
)
from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen
from cpython.bytes cimport (
    PyBytes_AsString,
    PyBytes_FromStringAndSize,
    PyBytes_Check,
)
from cpython.unicode cimport PyUnicode_DecodeUTF8

from .definitions cimport (
    uintptr_t,
    hid_t,
    herr_t,
    hsize_t,
    hvl_t,
    uint32_t,
    H5S_seloper_t,
    H5D_FILL_VALUE_UNDEFINED,
    H5O_TYPE_UNKNOWN,
    H5O_TYPE_GROUP,
    H5O_TYPE_DATASET,
    H5O_TYPE_NAMED_DATATYPE,
    H5L_TYPE_ERROR,
    H5L_TYPE_HARD,
    H5L_TYPE_SOFT,
    H5L_TYPE_EXTERNAL,
    H5T_class_t,
    H5T_sign_t,
    H5T_NATIVE_INT,
    H5T_cset_t,
    H5T_CSET_ASCII,
    H5T_CSET_UTF8,
    H5F_SCOPE_GLOBAL,
    H5F_ACC_TRUNC,
    H5F_ACC_RDONLY,
    H5F_ACC_RDWR,
    H5P_DEFAULT,
    H5P_FILE_ACCESS,
    H5P_FILE_CREATE,
    H5T_DIR_DEFAULT,
    H5S_SELECT_SET,
    H5S_SELECT_AND,
    H5S_SELECT_NOTB,
    H5Fcreate,
    H5Fopen,
    H5Fclose,
    H5Fflush,
    H5Fget_vfd_handle,
    H5Fget_filesize,
    H5Fget_create_plist,
    H5Gcreate,
    H5Gopen,
    H5Gclose,
    H5Ldelete,
    H5Lmove,
    H5Dopen,
    H5Dclose,
    H5Dread,
    H5Dwrite,
    H5Dget_type,
    H5Dget_create_plist,
    H5Dget_space,
    H5Dvlen_reclaim,
    H5Dget_storage_size,
    H5Dvlen_get_buf_size,
    H5Dget_chunk_info_by_coord,
    haddr_t,
    HADDR_UNDEF,
    H5Dread_chunk,
    H5Dwrite_chunk,
    H5Tget_native_type,
    H5Tclose,
    H5Tis_variable_str,
    H5Tget_sign,
    H5Adelete,
    H5T_BITFIELD,
    H5T_INTEGER,
    H5T_FLOAT,
    H5T_STRING,
    H5Tget_order,
    H5Pcreate,
    H5Pset_cache,
    H5Pclose,
    H5Pget_userblock,
    H5Pset_userblock,
    H5Pset_fapl_sec2,
    H5Pset_fapl_log,
    H5Pset_fapl_stdio,
    H5Pset_fapl_core,
    H5Pset_fapl_split,
    H5Pget_obj_track_times,
    H5Sselect_all,
    H5Sselect_elements,
    H5Sselect_hyperslab,
    H5Screate_simple,
    H5Sclose,
    H5Oget_info,
    H5O_info_t,
    H5ATTRset_attribute,
    H5ATTRset_attribute_string,
    H5ATTRget_attribute,
    H5ATTRget_attribute_string,
    H5ATTRget_attribute_vlen_string_array,
    H5ATTRfind_attribute,
    H5ATTRget_type_ndims,
    H5ATTRget_dims,
    H5ARRAYget_ndims,
    H5ARRAYget_info,
    set_cache_size,
    get_objinfo,
    get_linkinfo,
    Giterate,
    Aiterate,
    H5UIget_info,
    get_len_of_range,
    conv_float64_timeval32,
    truncate_dset,
    H5_HAVE_DIRECT_DRIVER,
    pt_H5Pset_fapl_direct,
    H5_HAVE_WINDOWS_DRIVER,
    pt_H5Pset_fapl_windows,
    H5_HAVE_IMAGE_FILE,
    H5Pset_file_image,
    H5Fget_file_image,
    H5Tget_size,
    hobj_ref_t,
)


cdef int H5T_CSET_DEFAULT = 16

from .utilsextension cimport (
    malloc_dims,
    get_native_type,
    cstr_to_pystr,
    load_reference,
)

#-------------------------------------------------------------------

cdef extern from "Python.h":

    object PyByteArray_FromStringAndSize(char *s, Py_ssize_t len)

cdef extern from "H5ARRAY-opt.h" nogil:
  hid_t H5ARRAYOmake( hid_t loc_id,
                      const char *dset_name,
                      const char *obversion,
                      const int rank,
                      const hsize_t *dims,
                      int   extdim,
                      hid_t type_id,
                      hsize_t *dims_chunk,
                      hsize_t block_size,
                      void  *fill_data,
                      int   compress,
                      char  *complib,
                      int   shuffle,
                      int   fletcher32,
                      hbool_t track_times,
                      const void *data);


  herr_t H5ARRAYOreadSlice(char* filename,
                           hbool_t blosc2_support,
                           hid_t dataset_id,
                           hid_t type_id,
                           hsize_t *slice_start,
                           hsize_t *slice_stop,
                           hsize_t *slice_step,
                           void *slice_data);


# Functions from HDF5 ARRAY (this is not part of HDF5 HL; it's private)
cdef extern from "H5ARRAY.h" nogil:

  herr_t H5ARRAYmake(hid_t loc_id, char *dset_name, char *obversion,
                     int rank, hsize_t *dims, int extdim,
                     hid_t type_id, hsize_t *dims_chunk, void *fill_data,
                     int complevel, char  *complib, int shuffle,
                     int fletcher32, hbool_t track_times, void *data)

  herr_t H5ARRAYappend_records(hid_t dataset_id, hid_t type_id,
                               int rank, hsize_t *dims_orig,
                               hsize_t *dims_new, int extdim, void *data )

  herr_t H5ARRAYwrite_records(hid_t dataset_id, hid_t type_id,
                              int rank, hsize_t *start, hsize_t *step,
                              hsize_t *count, void *data)

  herr_t H5ARRAYread(hid_t dataset_id, hid_t type_id,
                     hsize_t start, hsize_t nrows, hsize_t step,
                     int extdim, void *data)

  herr_t H5ARRAYreadSlice(hid_t dataset_id, hid_t type_id,
                          hsize_t *start, hsize_t *stop,
                          hsize_t *step, void *data)

  herr_t H5ARRAYreadIndex(hid_t dataset_id, hid_t type_id, int notequal,
                          hsize_t *start, hsize_t *stop, hsize_t *step,
                          void *data)

  herr_t H5ARRAYget_chunkshape(hid_t dataset_id, int rank, hsize_t *dims_chunk)

  herr_t H5ARRAYget_fill_value( hid_t dataset_id, hid_t type_id,
                                int *status, void *value)


# Functions for dealing with VLArray objects
cdef extern from "H5VLARRAY.h" nogil:

  herr_t H5VLARRAYmake( hid_t loc_id, char *dset_name, char *obversion,
                        int rank, hsize_t *dims, hid_t type_id,
                        hsize_t chunk_size, void *fill_data, int complevel,
                        char *complib, int shuffle, int fletcher32,
                        hbool_t track_times, void *data)

  herr_t H5VLARRAYappend_records( hid_t dataset_id, hid_t type_id,
                                  int nobjects, hsize_t nrecords,
                                  void *data )

  herr_t H5VLARRAYmodify_records( hid_t dataset_id, hid_t type_id,
                                  hsize_t nrow, int nobjects,
                                  void *data )

  herr_t H5VLARRAYget_info( hid_t dataset_id, hid_t type_id,
                            hsize_t *nrecords, char *base_byteorder)


#----------------------------------------------------------------------------

# Initialization code

# The numpy API requires this function to be called before
# using any numpy facilities in an extension module.
import_array()

#---------------------------------------------------------------------------

# Helper functions

cdef hsize_t *npy_malloc_dims(int rank, npy_intp *pdims):
  """Returns a malloced hsize_t dims from a npy_intp *pdims."""

  cdef int i
  cdef hsize_t *dims

  dims = NULL
  if rank > 0:
    dims = <hsize_t *>malloc(rank * sizeof(hsize_t))
    for i from 0 <= i < rank:
      dims[i] = pdims[i]
  return dims


cdef object getshape(int rank, hsize_t *dims):
  """Return a shape (tuple) from a dims C array of rank dimensions."""

  cdef int i
  cdef object shape

  shape = []
  for i from 0 <= i < rank:
    shape.append(SizeType(dims[i]))

  return tuple(shape)


# Helper function for quickly fetch an attribute string
cdef object get_attribute_string_or_none(hid_t node_id, char* attr_name):
  """Returns a string/unicode attribute if it exists in node_id.

  It returns ``None`` in case it don't exists (or there have been problems
  reading it).

  """

  cdef char *attr_value
  cdef int cset = H5T_CSET_DEFAULT
  cdef object retvalue
  cdef hsize_t size

  attr_value = NULL
  retvalue = None   # Default value
  if H5ATTRfind_attribute(node_id, attr_name):
    size = H5ATTRget_attribute_string(node_id, attr_name, &attr_value, &cset)
    if size == 0:
      if cset == H5T_CSET_UTF8:
        retvalue = np.str_('')
      else:
        retvalue = np.bytes_(b'')
    elif cset == H5T_CSET_UTF8:
      retvalue = PyUnicode_DecodeUTF8(attr_value, size, NULL)
      retvalue = np.str_(retvalue)
    else:
      retvalue = PyBytes_FromStringAndSize(attr_value, size)
      # AV: oct 2012
      # since now we use the string size got form HDF5 we have to strip
      # trailing zeros used for padding.
      # The entire process is quite odd but due to a bug (??) in the way
      # numpy arrays are pickled in python 3 we can't assume that
      # strlen(attr_value) is the actual length of the attribute
      # and np.bytes_(attr_value) can give a truncated pickle string
      retvalue = retvalue.rstrip(b'\x00')
      retvalue = np.bytes_(retvalue)

    # Important to release attr_value, because it has been malloc'ed!
    if attr_value:
      free(<void *>attr_value)

  return retvalue


# Get the numpy dtype scalar attribute from an HDF5 type as fast as possible
cdef object get_dtype_scalar(hid_t type_id, H5T_class_t class_id,
                             size_t itemsize):
  cdef H5T_sign_t sign
  cdef object stype

  if class_id == H5T_BITFIELD:
    stype = "b1"
  elif class_id == H5T_INTEGER:
    # Get the sign
    sign = H5Tget_sign(type_id)
    if (sign > 0):
      stype = "i%s" % (itemsize)
    else:
      stype = "u%s" % (itemsize)
  elif class_id ==  H5T_FLOAT:
    stype = "f%s" % (itemsize)
  elif class_id ==  H5T_STRING:
    if H5Tis_variable_str(type_id):
      raise TypeError("variable length strings are not supported yet")
    stype = "S%s" % (itemsize)

  # Try to get a NumPy type.  If this can't be done, return None.
  try:
    ntype = np.dtype(stype)
  except TypeError:
    ntype = None
  return ntype


_supported_drivers = (
    "H5FD_SEC2",
    "H5FD_DIRECT",
    #"H5FD_LOG",
    "H5FD_WINDOWS",
    "H5FD_STDIO",
    "H5FD_CORE",
    #"H5FD_FAMILY",
    #"H5FD_MULTI",
    "H5FD_SPLIT",
    #"H5FD_MPIO",
    #"H5FD_MPIPOSIX",
    #"H5FD_STREAM",
)

HAVE_DIRECT_DRIVER = bool(H5_HAVE_DIRECT_DRIVER)
HAVE_WINDOWS_DRIVER = bool(H5_HAVE_WINDOWS_DRIVER)

# Type extensions declarations (these are subclassed by PyTables
# Python classes)

cdef class File:
  cdef hid_t   file_id
  cdef hid_t   access_plist
  cdef object  name

  def _g_new(self, name, pymode, **params):
    cdef herr_t err = 0
    cdef hid_t access_plist, create_plist = H5P_DEFAULT
    cdef hid_t meta_plist_id = H5P_DEFAULT, raw_plist_id = H5P_DEFAULT
    cdef size_t img_buf_len = 0, user_block_size = 0
    cdef void *img_buf_p = NULL
    cdef bytes encname
    #cdef bytes logfile_name

    # Check if we can handle the driver
    driver = params["DRIVER"]
    if driver is not None and driver not in _supported_drivers:
      raise ValueError("Invalid or not supported driver: '%s'" % driver)
    if driver == "H5FD_SPLIT":
      meta_ext = params.get("DRIVER_SPLIT_META_EXT", "-m.h5")
      raw_ext = params.get("DRIVER_SPLIT_RAW_EXT", "-r.h5")
      meta_name = meta_ext % name if "%s" in meta_ext else name + meta_ext
      raw_name = raw_ext % name if "%s" in raw_ext else name + raw_ext
      enc_meta_ext = encode_filename(meta_ext)
      enc_raw_ext = encode_filename(raw_ext)

    # Create a new file using default properties
    self.name = name

    # Encode the filename in case it is unicode
    encname = encode_filename(name)

    # These fields can be seen from Python.
    self._v_new = None  # this will be computed later
    # """Is this file going to be created from scratch?"""

    self._isPTFile = True  # assume a PyTables file by default
    # """Does this HDF5 file have a PyTables format?"""

    assert pymode in ('r', 'r+', 'a', 'w'), ("an invalid mode string ``%s`` "
           "passed the ``check_file_access()`` test; "
           "please report this to the authors" % pymode)

    image = params.get('DRIVER_CORE_IMAGE')
    if image:
      if driver != "H5FD_CORE":
        warnings.warn("The DRIVER_CORE_IMAGE parameter will be ignored by "
                      "the '%s' driver" % driver)
      elif not PyBytes_Check(image):
        raise TypeError("The DRIVER_CORE_IMAGE must be a string of bytes")

    # After the following check we can be quite sure
    # that the file or directory exists and permissions are right.
    if driver == "H5FD_SPLIT":
      for n in meta_name, raw_name:
        check_file_access(n, pymode)
    else:
      backing_store = params.get("DRIVER_CORE_BACKING_STORE", 1)
      if driver != "H5FD_CORE" or backing_store:
        check_file_access(name, pymode)

    # Should a new file be created?
    if image:
      exists = True
    elif driver == "H5FD_SPLIT":
      exists = os.path.exists(meta_name) and os.path.exists(raw_name)
    else:
      exists = os.path.exists(name)
    self._v_new = not (pymode in ('r', 'r+') or (pymode == 'a' and exists))

    user_block_size = params.get("USER_BLOCK_SIZE", 0)
    if user_block_size and not self._v_new:
        warnings.warn("The HDF5 file already esists: the USER_BLOCK_SIZE "
                      "will be ignored")
    elif user_block_size:
      user_block_size = int(user_block_size)
      is_pow_of_2 = ((user_block_size & (user_block_size - 1)) == 0)
      if user_block_size < 512 or not is_pow_of_2:
        raise ValueError("The USER_BLOCK_SIZE must be a power od 2 greather "
                         "than 512 or zero")

      # File creation property list
      create_plist = H5Pcreate(H5P_FILE_CREATE)
      err = H5Pset_userblock(create_plist, user_block_size)
      if err < 0:
        H5Pclose(create_plist)
        raise HDF5ExtError("Unable to set the user block size")

    # File access property list
    access_plist = H5Pcreate(H5P_FILE_ACCESS)

    # Set parameters for chunk cache
    H5Pset_cache(access_plist, 0,
                 params["CHUNK_CACHE_NELMTS"],
                 params["CHUNK_CACHE_SIZE"],
                 params["CHUNK_CACHE_PREEMPT"])

    # Set the I/O driver
    if driver == "H5FD_SEC2":
      err = H5Pset_fapl_sec2(access_plist)
    elif driver == "H5FD_DIRECT":
      if not H5_HAVE_DIRECT_DRIVER:
        H5Pclose(create_plist)
        H5Pclose(access_plist)
        raise RuntimeError("The H5FD_DIRECT driver is not available")
      err = pt_H5Pset_fapl_direct(access_plist,
                                  params["DRIVER_DIRECT_ALIGNMENT"],
                                  params["DRIVER_DIRECT_BLOCK_SIZE"],
                                  params["DRIVER_DIRECT_CBUF_SIZE"])
    #elif driver == "H5FD_LOG":
    #  if "DRIVER_LOG_FILE" not in params:
    #    H5Pclose(access_plist)
    #    raise ValueError("The DRIVER_LOG_FILE parameter is required for "
    #                     "the H5FD_LOG driver")
    #  logfile_name = encode_filename(params["DRIVER_LOG_FILE"])
    #  err = H5Pset_fapl_log(access_plist,
    #                        <char*>logfile_name,
    #                        params["DRIVER_LOG_FLAGS"],
    #                        params["DRIVER_LOG_BUF_SIZE"])
    elif driver == "H5FD_WINDOWS":
      if not H5_HAVE_WINDOWS_DRIVER:
        H5Pclose(access_plist)
        H5Pclose(create_plist)
        raise RuntimeError("The H5FD_WINDOWS driver is not available")
      err = pt_H5Pset_fapl_windows(access_plist)
    elif driver == "H5FD_STDIO":
      err = H5Pset_fapl_stdio(access_plist)
    elif driver == "H5FD_CORE":
      err = H5Pset_fapl_core(access_plist,
                             params["DRIVER_CORE_INCREMENT"],
                             backing_store)
      if image:
        img_buf_len = len(image)
        img_buf_p = <void *>PyBytes_AsString(image)
        err = H5Pset_file_image(access_plist, img_buf_p, img_buf_len)
        if err < 0:
          H5Pclose(create_plist)
          H5Pclose(access_plist)
          raise HDF5ExtError("Unable to set the file image")

    #elif driver == "H5FD_FAMILY":
    #  H5Pset_fapl_family(access_plist,
    #                     params["DRIVER_FAMILY_MEMB_SIZE"],
    #                     fapl_id)
    #elif driver == "H5FD_MULTI":
    #  err = H5Pset_fapl_multi(access_plist, memb_map, memb_fapl, memb_name,
    #                          memb_addr, relax)
    elif driver == "H5FD_SPLIT":
      err = H5Pset_fapl_split(access_plist, enc_meta_ext, meta_plist_id,
                              enc_raw_ext, raw_plist_id)
    if err < 0:
      e = HDF5ExtError("Unable to set the file access property list")
      H5Pclose(create_plist)
      H5Pclose(access_plist)
      raise e

    if pymode == 'r':
      self.file_id = H5Fopen(encname, H5F_ACC_RDONLY, access_plist)
    elif pymode == 'r+':
      self.file_id = H5Fopen(encname, H5F_ACC_RDWR, access_plist)
    elif pymode == 'a':
      if exists:
        # A test for logging.
        ## H5Pset_sieve_buf_size(access_plist, 0)
        ## H5Pset_fapl_log (access_plist, "test.log", H5FD_LOG_LOC_WRITE, 0)
        self.file_id = H5Fopen(encname, H5F_ACC_RDWR, access_plist)
      else:
        self.file_id = H5Fcreate(encname, H5F_ACC_TRUNC, create_plist,
                                 access_plist)
    elif pymode == 'w':
      self.file_id = H5Fcreate(encname, H5F_ACC_TRUNC, create_plist,
                               access_plist)

    if self.file_id < 0:
        e = HDF5ExtError("Unable to open/create file '%s'" % name)
        H5Pclose(create_plist)
        H5Pclose(access_plist)
        raise e

    H5Pclose(create_plist)
    H5Pclose(access_plist)

    # Set the cache size
    set_cache_size(self.file_id, params["METADATA_CACHE_SIZE"])

    # Set the maximum number of threads for Blosc
    set_blosc_max_threads(params["MAX_BLOSC_THREADS"])
    set_blosc2_max_threads(params["MAX_BLOSC_THREADS"])

  # XXX: add the possibility to pass a pre-allocated buffer
  def get_file_image(self):
    """Retrieves an in-memory image of an existing, open HDF5 file.

    .. versionadded:: 3.0

    """

    cdef ssize_t size = 0
    cdef size_t buf_len = 0
    cdef bytes image
    cdef char* cimage

    self.flush()

    # retrieve the size of the buffer for the file image
    size = H5Fget_file_image(self.file_id, NULL, buf_len)
    if size < 0:
      raise HDF5ExtError("Unable to retrieve the size of the buffer for the "
                         "file image.  Plese note that not all drivers "
                         "provide support for image files.")

    # allocate the memory buffer
    image = PyBytes_FromStringAndSize(NULL, size)
    if not image:
      raise RuntimeError("Unable to allecote meomory fir the file image")

    cimage = image
    buf_len = size
    size = H5Fget_file_image(self.file_id, <void*>cimage, buf_len)
    if size < 0:
      raise HDF5ExtError("Unable to retrieve the file image. "
                         "Plese note that not all drivers provide support "
                         "for image files.")

    return image

  def get_filesize(self):
    """Returns the size of an HDF5 file.

    The returned size is that of the entire file, as opposed to only
    the HDF5 portion of the file. I.e., size includes the user block,
    if any, the HDF5 portion of the file, and any data that may have
    been appended beyond the data written through the HDF5 Library.

    .. versionadded:: 3.0

    """

    cdef herr_t err = 0
    cdef hsize_t size = 0

    err = H5Fget_filesize(self.file_id, &size)
    if err < 0:
      raise HDF5ExtError("Unable to retrieve the HDF5 file size")

    return size

  def get_userblock_size(self):
    """Retrieves the size of a user block.

    .. versionadded:: 3.0

    """

    cdef herr_t err = 0
    cdef hsize_t size = 0
    cdef hid_t create_plist

    create_plist = H5Fget_create_plist(self.file_id)
    if create_plist < 0:
      raise HDF5ExtError("Unable to get the creation property list")

    err = H5Pget_userblock(create_plist, &size)
    if err < 0:
      H5Pclose(create_plist)
      raise HDF5ExtError("unable to retrieve the user block size")

    H5Pclose(create_plist)

    return size

  # Accessor definitions
  def _get_file_id(self):
    return self.file_id

  def fileno(self):
    """Return the underlying OS integer file descriptor.

    This is needed for lower-level file interfaces, such as the ``fcntl``
    module.

    """

    cdef void *file_handle
    cdef uintptr_t *descriptor
    cdef herr_t err
    err = H5Fget_vfd_handle(self.file_id, H5P_DEFAULT, &file_handle)
    if err < 0:
      raise HDF5ExtError(
        "Problems getting file descriptor for file ``%s``" % self.name)
    # Convert the 'void *file_handle' into an 'int *descriptor'
    descriptor = <uintptr_t *>file_handle
    return descriptor[0]


  def _flush_file(self, scope):
    # Close the file
    H5Fflush(self.file_id, scope)


  def _close_file(self):
    # Close the file
    H5Fclose( self.file_id )
    self.file_id = 0    # Means file closed


  # This method is moved out of scope, until we provide code to delete
  # the memory booked by this extension types
  def __dealloc__(self):
    cdef int ret
    if self.file_id > 0:
      # Close the HDF5 file because user didn't do that!
      ret = H5Fclose(self.file_id)
      if ret < 0:
        raise HDF5ExtError("Problems closing the file '%s'" % self.name)


cdef class AttributeSet:
  cdef object name

  def _g_new(self, node):
    self.name = node._v_name

  def _g_list_attr(self, node):
    """Return a tuple with the attribute list"""
    a = Aiterate(node._v_objectid)
    return a


  def _g_setattr(self, node, name, object value):
    """Save Python or NumPy objects as HDF5 attributes.

    Scalar Python objects, scalar NumPy & 0-dim NumPy objects will all be
    saved as H5T_SCALAR type.  N-dim NumPy objects will be saved as H5T_ARRAY
    type.

    """

    cdef int ret
    cdef hid_t dset_id, type_id
    cdef hsize_t *dims
    cdef ndarray ndv
    cdef object byteorder, rabyteorder, baseatom
    cdef char* cname = NULL
    cdef bytes encoded_name
    cdef int cset = H5T_CSET_DEFAULT

    encoded_name = name.encode('utf-8')
    # get the C pointer
    cname = encoded_name

    # The dataset id of the node
    dset_id = node._v_objectid

    # Convert a NumPy scalar into a NumPy 0-dim ndarray
    if isinstance(value, np.generic):
      value = np.array(value)

    # Check if value is a NumPy ndarray and of a supported type
    if (isinstance(value, np.ndarray) and
        value.dtype.kind in ('V', 'S', 'b', 'i', 'u', 'f', 'c')):
      # get a contiguous array: fixes #270 and gh-176
      #value = np.ascontiguousarray(value)
      value = value.copy()
      if value.dtype.kind == 'V':
        description, rabyteorder = descr_from_dtype(value.dtype, ptparams=node._v_file.params)
        byteorder = byteorders[rabyteorder]
        type_id = create_nested_type(description, byteorder)
        # Make sure the value is consistent with offsets of the description
        value = value.astype(description._v_dtype)
      else:
        # Get the associated native HDF5 type of the scalar type
        baseatom = Atom.from_dtype(value.dtype.base)
        byteorder = byteorders[value.dtype.byteorder]
        type_id = atom_to_hdf5_type(baseatom, byteorder)
      # Get dimensionality info
      ndv = <ndarray>value
      dims = npy_malloc_dims(PyArray_NDIM(ndv), PyArray_DIMS(ndv))
      # Actually write the attribute
      ret = H5ATTRset_attribute(dset_id, cname, type_id,
                                PyArray_NDIM(ndv), dims, PyArray_BYTES(ndv))
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in node:\n %s." %
                           (name, self._v_node))
      # Release resources
      free(<void *>dims)
      H5Tclose(type_id)
    else:
      # Object cannot be natively represented in HDF5.
      if (isinstance(value, np.ndarray) and
          value.dtype.kind == 'U' and
          value.shape == ()):
        value = value[()].encode('utf-8')
        cset = H5T_CSET_UTF8
      else:
        # Convert this object to a null-terminated string
        # (binary pickles are not supported at this moment)
        value = pickle.dumps(value, 0)

      ret = H5ATTRset_attribute_string(dset_id, cname, value, len(value), cset)
      if ret < 0:
        raise HDF5ExtError("Can't set attribute '%s' in node:\n %s." %
                           (name, self._v_node))


  # Get attributes
  def _g_getattr(self, node, attrname):
    """Get HDF5 attributes and retrieve them as NumPy objects.

    H5T_SCALAR types will be retrieved as scalar NumPy.
    H5T_ARRAY types will be retrieved as ndarray NumPy objects.

    """

    cdef hsize_t *dims
    cdef H5T_class_t class_id
    cdef size_t type_size
    cdef hid_t mem_type, dset_id, type_id, native_type
    cdef int rank, ret, enumtype
    cdef void *rbuf
    cdef char *str_value
    cdef char **str_values = NULL
    cdef ndarray ndvalue
    cdef object shape, stype_atom, shape_atom, retvalue
    cdef int i, nelements
    cdef char* cattrname = NULL
    cdef bytes encoded_attrname
    cdef int cset = H5T_CSET_DEFAULT

    encoded_attrname = attrname.encode('utf-8')
    # Get the C pointer
    cattrname = encoded_attrname

    # The dataset id of the node
    dset_id = node._v_objectid
    dims = NULL

    ret = H5ATTRget_type_ndims(dset_id, cattrname, &type_id, &class_id,
                               &type_size, &rank )
    if ret < 0:
      raise HDF5ExtError("Can't get type info on attribute %s in node %s." %
                         (attrname, self.name))

    # Call a fast function for scalar values and typical class types
    if (rank == 0 and class_id == H5T_STRING):
      type_size = H5ATTRget_attribute_string(dset_id, cattrname, &str_value,
                                             &cset)
      if type_size == 0:
        if cset == H5T_CSET_UTF8:
          retvalue = np.str_('')
        else:
          retvalue = np.bytes_(b'')

      elif cset == H5T_CSET_UTF8:
        retvalue = PyUnicode_DecodeUTF8(str_value, type_size, NULL)
        retvalue = np.str_(retvalue)
      else:
        retvalue = PyBytes_FromStringAndSize(str_value, type_size)
        # AV: oct 2012
        # since now we use the string size got form HDF5 we have to strip
        # trailing zeros used for padding.
        # The entire process is quite odd but due to a bug (??) in the way
        # numpy arrays are pickled in python 3 we can't assume that
        # strlen(attr_value) is the actual length of the attibute
        # and np.bytes_(attr_value) can give a truncated pickle sting
        retvalue = retvalue.rstrip(b'\x00')
        retvalue = np.bytes_(retvalue)     # bytes
      # Important to release attr_value, because it has been malloc'ed!
      if str_value:
        free(str_value)
      H5Tclose(type_id)
      return retvalue
    elif (rank == 0 and class_id in (H5T_BITFIELD, H5T_INTEGER, H5T_FLOAT)):
      dtype_ = get_dtype_scalar(type_id, class_id, type_size)
      if dtype_ is None:
        warnings.warn("Unsupported type for attribute '%s' in node '%s'. "
                      "Offending HDF5 class: %d" % (attrname, self.name,
                                                    class_id), DataTypeWarning)
        self._v_unimplemented.append(attrname)
        return None
      shape = ()
    else:
      # General case

      # Get the dimensional info
      dims = <hsize_t *>malloc(rank * sizeof(hsize_t))
      ret = H5ATTRget_dims(dset_id, cattrname, dims)
      if ret < 0:
        raise HDF5ExtError("Can't get dims info on attribute %s in node %s." %
                           (attrname, self.name))
      shape = getshape(rank, dims)
      # dims is not needed anymore
      free(<void *> dims)

      # Get the NumPy dtype from the type_id
      try:
        stype_, shape_ = hdf5_to_np_ext_type(type_id, pure_numpy_types=True, ptparams=node._v_file.params)
        dtype_ = np.dtype(stype_, shape_)
      except TypeError:
        if class_id == H5T_STRING and H5Tis_variable_str(type_id):
          nelements = H5ATTRget_attribute_vlen_string_array(dset_id, cattrname,
                                                            &str_values, &cset)
          if nelements < 0:
            raise HDF5ExtError("Can't read attribute %s in node %s." %
                               (attrname, self.name))

          # The following generator expressions do not work with Cython 0.15.1
          if cset == H5T_CSET_UTF8:
            #retvalue = np.fromiter(
            #  PyUnicode_DecodeUTF8(<char*>str_values[i],
            #                        strlen(<char*>str_values[i]),
            #                        NULL)
            #    for i in range(nelements), "O8")
            retvalue = np.array([
              PyUnicode_DecodeUTF8(<char*>str_values[i],
                                    strlen(<char*>str_values[i]),
                                    NULL)
                for i in range(nelements)], "O8")

          else:
            #retvalue = np.fromiter(
            #  <char*>str_values[i] for i in range(nelements), "O8")
            retvalue = np.array(
              [<char*>str_values[i] for i in range(nelements)], "O8")
          retvalue.shape = shape

          # Important to release attr_value, because it has been malloc'ed!
          for i in range(nelements):
            free(str_values[i])
          free(str_values)

          return retvalue

        # This class is not supported. Instead of raising a TypeError, issue a
        # warning explaining the problem. This will allow to continue browsing
        # native HDF5 files, while informing the user about the problem.
        warnings.warn("Unsupported type for attribute '%s' in node '%s'. "
                      "Offending HDF5 class: %d" % (attrname, self.name,
                                                    class_id), DataTypeWarning)
        self._v_unimplemented.append(attrname)
        return None

    # Get the container for data
    ndvalue = np.empty(dtype=dtype_, shape=shape)
    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(ndvalue)
    # Actually read the attribute from disk
    ret = H5ATTRget_attribute(dset_id, cattrname, type_id, rbuf)
    if ret < 0:
      raise HDF5ExtError("Attribute %s exists in node %s, but can't get it." %
                         (attrname, self.name))
    H5Tclose(type_id)

    if rank > 0:    # multidimensional case
      retvalue = ndvalue
    else:
      retvalue = ndvalue[()]   # 0-dim ndarray becomes a NumPy scalar

    return retvalue


  def _g_remove(self, node, attrname):
    cdef int ret
    cdef hid_t dset_id
    cdef char *cattrname = NULL
    cdef bytes encoded_attrname

    encoded_attrname = attrname.encode('utf-8')
    # Get the C pointer
    cattrname = encoded_attrname

    # The dataset id of the node
    dset_id = node._v_objectid

    ret = H5Adelete(dset_id, cattrname)
    if ret < 0:
      raise HDF5ExtError("Attribute '%s' exists in node '%s', but cannot be "
                         "deleted." % (attrname, self.name))


cdef class Node:
  # Instance variables declared in .pxd

  def _g_new(self, where, name, init):
    self.name = name
    # """The name of this node in its parent group."""
    self.parent_id = where._v_objectid
    # """The identifier of the parent group."""

  def _g_delete(self, parent):
    cdef int ret
    cdef bytes encoded_name

    encoded_name = self.name.encode('utf-8')

    # Delete this node
    ret = H5Ldelete(parent._v_objectid, encoded_name, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("problems deleting the node ``%s``" % self.name)
    return ret

  def __dealloc__(self):
    self.parent_id = 0

  def _get_obj_info(self):
    cdef herr_t ret = 0
    cdef H5O_info_t oinfo

    ret = H5Oget_info(self._v_objectid, &oinfo)
    if ret < 0:
      raise HDF5ExtError("Unable to get object info for '%s'" %
                         self. _v_pathname)

    return ObjInfo(oinfo.addr, oinfo.rc)

  def _get_obj_timestamps(self):
    cdef herr_t ret = 0
    cdef H5O_info_t oinfo

    ret = H5Oget_info(self._v_objectid, &oinfo)
    if ret < 0:
      raise HDF5ExtError("Unable to get object info for '%s'" %
                         self. _v_pathname)

    return ObjTimestamps(oinfo.atime, oinfo.mtime, oinfo.ctime,
                         oinfo.btime)


cdef class Group(Node):
  cdef hid_t   group_id

  def _g_create(self):
    cdef hid_t ret
    cdef bytes encoded_name

    encoded_name = self.name.encode('utf-8')

    # @TODO: set property list --> utf-8

    # Create a new group
    ret = H5Gcreate(self.parent_id, encoded_name, H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("Can't create the group %s." % self.name)
    self.group_id = ret
    return self.group_id

  def _g_open(self):
    cdef hid_t ret
    cdef bytes encoded_name

    encoded_name = self.name.encode('utf-8')

    ret = H5Gopen(self.parent_id, encoded_name, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("Can't open the group: '%s'." % self.name)
    self.group_id = ret
    return self.group_id

  def _g_get_objinfo(self, object h5name):
    """Check whether 'name' is a children of 'self' and return its type."""

    cdef int ret
    cdef object node_type
    cdef bytes encoded_name
    cdef char *cname

    encoded_name = h5name.encode('utf-8')
    # Get the C pointer
    cname = encoded_name

    ret = get_linkinfo(self.group_id, cname)
    if ret == -2 or ret == H5L_TYPE_ERROR:
      node_type = "NoSuchNode"
    elif ret == H5L_TYPE_SOFT:
      node_type = "SoftLink"
    elif ret == H5L_TYPE_EXTERNAL:
      node_type = "ExternalLink"
    elif ret == H5L_TYPE_HARD:
        ret = get_objinfo(self.group_id, cname)
        if ret == -2:
          node_type = "NoSuchNode"
        elif ret == H5O_TYPE_UNKNOWN:
          node_type = "Unknown"
        elif ret == H5O_TYPE_GROUP:
          node_type = "Group"
        elif ret == H5O_TYPE_DATASET:
          node_type = "Leaf"
        elif ret == H5O_TYPE_NAMED_DATATYPE:
          node_type = "NamedType"              # Not supported yet
        #else H5O_TYPE_LINK:
        #    # symbolic link
        #    raise RuntimeError('unexpected object type')
        else:
          node_type = "Unknown"
    return node_type

  def _g_list_group(self, parent):
    """Return a tuple with the groups and the leaves hanging from self."""

    cdef bytes encoded_name

    encoded_name = self.name.encode('utf-8')

    return Giterate(parent._v_objectid, self._v_objectid, encoded_name)


  def _g_get_gchild_attr(self, group_name, attr_name):
    """Return an attribute of a child `Group`.

    If the attribute does not exist, ``None`` is returned.

    """

    cdef hid_t gchild_id
    cdef object retvalue
    cdef bytes encoded_group_name
    cdef bytes encoded_attr_name

    encoded_group_name = group_name.encode('utf-8')
    encoded_attr_name = attr_name.encode('utf-8')

    # Open the group
    retvalue = None  # Default value
    gchild_id = H5Gopen(self.group_id, encoded_group_name, H5P_DEFAULT)
    if gchild_id < 0:
      raise HDF5ExtError("Non-existing node ``%s`` under ``%s``" %
                         (group_name, self._v_pathname))
    retvalue = get_attribute_string_or_none(gchild_id, encoded_attr_name)
    # Close child group
    H5Gclose(gchild_id)

    return retvalue


  def _g_get_lchild_attr(self, leaf_name, attr_name):
    """Return an attribute of a child `Leaf`.

    If the attribute does not exist, ``None`` is returned.

    """

    cdef hid_t leaf_id
    cdef object retvalue
    cdef bytes encoded_leaf_name
    cdef bytes encoded_attr_name

    encoded_leaf_name = leaf_name.encode('utf-8')
    encoded_attr_name = attr_name.encode('utf-8')

    # Open the dataset
    leaf_id = H5Dopen(self.group_id, encoded_leaf_name, H5P_DEFAULT)
    if leaf_id < 0:
      raise HDF5ExtError("Non-existing node ``%s`` under ``%s``" %
                         (leaf_name, self._v_pathname))
    retvalue = get_attribute_string_or_none(leaf_id, encoded_attr_name)
    # Close the dataset
    H5Dclose(leaf_id)
    return retvalue


  def _g_flush_group(self):
    # Close the group
    H5Fflush(self.group_id, H5F_SCOPE_GLOBAL)


  def _g_close_group(self):
    cdef int ret

    ret = H5Gclose(self.group_id)
    if ret < 0:
      raise HDF5ExtError("Problems closing the Group %s" % self.name)
    self.group_id = 0  # indicate that this group is closed


  def _g_move_node(self, hid_t oldparent, oldname, hid_t newparent, newname,
                   oldpathname, newpathname):
    cdef int ret
    cdef bytes encoded_oldname, encoded_newname

    encoded_oldname = oldname.encode('utf-8')
    encoded_newname = newname.encode('utf-8')

    ret = H5Lmove(oldparent, encoded_oldname, newparent, encoded_newname,
                  H5P_DEFAULT, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("Problems moving the node %s to %s" %
                         (oldpathname, newpathname) )
    return ret



cdef class Leaf(Node):
  # Instance variables declared in .pxd

  def _get_storage_size(self):
      return H5Dget_storage_size(self.dataset_id)

  def _get_obj_track_times(self):
    """Get track_times boolean for dataset

    Uses H5Pget_obj_track_times to determine if the dataset was
    created with the track_times property.  If the leaf is not a
    dataset, this will fail with HDF5ExtError.

    The track times dataset creation property does not seem to survive
    closing and reopening as of HDF5 1.8.17.  Currently, it may be
    more accurate to test whether the ctime for the dataset is 0:
    track_times = (leaf._get_obj_timestamps().ctime == 0)
    """
    cdef:
      hbool_t track_times = True

    if self.dataset_id < 0:
      raise ValueError('Invalid dataset id %s' % self.dataset_id)

    plist_id = H5Dget_create_plist(self.dataset_id)
    if plist_id < 0:
      raise HDF5ExtError("Could not get dataset creation property list "
                         "from dataset id %s" % self.dataset_id)

    try:
      # Get track_times boolean for dataset
      if H5Pget_obj_track_times(plist_id, &track_times) < 0:
        raise HDF5ExtError("Could not get dataset track_times property "
                           "from dataset id %s" % self.dataset_id)
    finally:
      H5Pclose(plist_id)

    return bool(track_times)

  def _g_new(self, where, name, init):
    if init:
      # Put this info to 0 just when the class is initialized
      self.dataset_id = -1
      self.type_id = -1
      self.base_type_id = -1
      self.disk_type_id = -1
    super()._g_new(where, name, init)

  def _g_chunk_info(self, ndarray coords):
    """Get storage information about chunk at `coords`.

    Return ``(filter_mask, offset, size)``, where items are ``None`` if the
    chunk is missing.

    """
    cdef herr_t ret
    cdef hsize_t *offset
    cdef unsigned filter_mask
    cdef haddr_t addr
    cdef hsize_t size

    # Get the pointer to the buffer data area of the coords array
    with nogil:
      offset = <hsize_t *>PyArray_DATA(coords)
      ret = H5Dget_chunk_info_by_coord(self.dataset_id, offset,
                                       &filter_mask, &addr, &size)
    if ret < 0:
      raise HDF5ExtError("Problems getting chunk info for ``%s``"
                         % self._v_pathname)
    return ((filter_mask, addr, size) if addr != HADDR_UNDEF
            else (None, None, None))

  def _g_read_chunk(self, ndarray coords, ndarray out):
    """Read the raw chunk at `coords` (into `out`).

    Return a new array of bytes if `out` is ``None``, `out` itself otherwise.
    Return ``None`` if the chunk is missing.

    """
    cdef ndarray rarr
    cdef herr_t ret
    cdef hsize_t *offset
    cdef uint32_t filters = 0
    cdef void *rbuf

    _, addr, size = self._g_chunk_info(coords)
    if addr is None:
      return None  # missing chunk
    if out is not None and len(out) < size:
      raise ValueError(f"Output buffer is too short: {len(out)} < {size}")

    rarr = np.empty((size,), dtype='u1') if out is None else out
    with nogil:
      rbuf = PyArray_DATA(rarr)
      offset = <hsize_t *>PyArray_DATA(coords)
      ret = H5Dread_chunk(self.dataset_id, H5P_DEFAULT, offset,
                          &filters, rbuf)
    if ret < 0:
      raise HDF5ExtError("Problems reading chunk from ``%s``"
                         % self._v_pathname)
    return rarr

  def _g_write_chunk(self, ndarray coords, ndarray data, uint32_t filters):
    """Write the raw `data` to the chunk in `coords`.

    The `filters` mask indicates which filters of the pipeline have not been
    used to create the `data`.

    """
    cdef herr_t ret
    cdef hsize_t *offset
    cdef size_t data_size
    cdef void *wbuf

    data_size = data.size
    with nogil:
      wbuf = PyArray_DATA(data)
      offset = <hsize_t *>PyArray_DATA(coords)
      ret = H5Dwrite_chunk(self.dataset_id, H5P_DEFAULT, filters,
                           offset, data_size, wbuf)
    if ret < 0:
      raise HDF5ExtError("Problems writing chunk to ``%s``"
                         % self._v_pathname)

  cdef _get_type_ids(self):
    """Get the disk and native HDF5 types associated with this leaf.

    It is guaranteed that both disk and native types are not the same
    descriptor (so that it is safe to close them separately).

    """

    cdef hid_t disk_type_id, native_type_id

    disk_type_id = H5Dget_type(self.dataset_id)
    native_type_id = get_native_type(disk_type_id)
    return disk_type_id, native_type_id

  cdef _convert_time64(self, ndarray nparr, int sense):
    """Converts a NumPy of Time64 elements between NumPy and HDF5 formats.

    NumPy to HDF5 conversion is performed when 'sense' is 0.  Otherwise, HDF5
    to NumPy conversion is performed.  The conversion is done in place,
    i.e. 'nparr' is modified.

    """

    cdef void *t64buf
    cdef long byteoffset, bytestride, nelements
    cdef hsize_t nrecords

    byteoffset = 0   # NumPy objects doesn't have an offset
    if (<object>nparr).shape == ():
      # 0-dim array does contain *one* element
      nrecords = 1
      bytestride = 8
    else:
      nrecords = len(nparr)
      bytestride = PyArray_STRIDE(nparr, 0)  # supports multi-dimensional recarray
    nelements = <size_t>nparr.size // nrecords
    t64buf = PyArray_DATA(nparr)

    conv_float64_timeval32(
      t64buf, byteoffset, bytestride, nrecords, nelements, sense)

  # can't do since cdef'd

  def _g_truncate(self, hsize_t size):
    """Truncate a Leaf to `size` nrows."""

    cdef hsize_t ret

    ret = truncate_dset(self.dataset_id, self.maindim, size)
    if ret < 0:
      raise HDF5ExtError("Problems truncating the leaf: %s" % self)

    classname = self.__class__.__name__
    if classname in ('EArray', 'CArray'):
      # Update the new dimensionality
      self.dims[self.maindim] = size
      # Update the shape
      shape = list(self.shape)
      shape[self.maindim] = SizeType(size)
      self.shape = tuple(shape)
    elif classname in ('Table', 'VLArray'):
      self.nrows = size
    else:
      raise ValueError("Unexpected classname: %s" % classname)

  def _g_flush(self):
    # Flush the dataset (in fact, the entire buffers in file!)
    if self.dataset_id >= 0:
        H5Fflush(self.dataset_id, H5F_SCOPE_GLOBAL)

  def _g_close(self):
    # Close dataset in HDF5 space
    # Release resources
    if self.type_id >= 0:
      H5Tclose(self.type_id)
    if self.disk_type_id >= 0:
      H5Tclose(self.disk_type_id)
    if self.base_type_id >= 0:
      H5Tclose(self.base_type_id)
    if self.dataset_id >= 0:
      H5Dclose(self.dataset_id)


cdef void* _array_data(ndarray arr):
    # When the object is not a 0-d ndarray and its strides == 0, that
    # means that the array does not contain actual data
    cdef npy_intp i, ndim

    ndim = PyArray_NDIM(arr)
    if ndim == 0:
        return PyArray_DATA(arr)
    for i in range(ndim):
        if PyArray_STRIDE(arr, i) > 0:
            return PyArray_DATA(arr)
    return NULL

def _supports_opt_blosc2_read_write(byteorder, filter_list, file_mode):
    if len(filter_list) == 1:  # Blosc2 must be the only filter
      opt_write = ((byteorder == sys.byteorder)
                   and ((filter_list[0] or "").startswith("blosc2")))
    else:
      opt_write = False
    # For reading, Windows does not support re-opening a file twice
    # in not read-only mode (for good reason), so we cannot use the
    # blosc2 opt
    opt_read = (opt_write
                and ((platform.system().lower() != 'windows') or
                     (file_mode == 'r')))
    return (opt_read, opt_write)

cdef class Array(Leaf):
  # Instance variables declared in .pxd

  def _create_array(self, ndarray nparr, object title, object atom):
    cdef int i
    cdef herr_t ret
    cdef void *rbuf
    cdef bytes complib, version, class_
    cdef object dtype_, atom_, shape
    cdef ndarray dims
    cdef bytes encoded_title, encoded_name
    cdef H5T_cset_t cset = H5T_CSET_ASCII

    encoded_title = title.encode('utf-8')
    encoded_name = self.name.encode('utf-8')

    # Get the HDF5 type associated with this numpy type
    shape = (<object>nparr).shape
    if atom is None or atom.shape == ():
      dtype_ = nparr.dtype.base
      atom_ = Atom.from_dtype(dtype_)
    else:
      atom_ = atom
      shape = shape[:-len(atom_.shape)]
    self.disk_type_id = atom_to_hdf5_type(atom_, self.byteorder)
    if self.disk_type_id < 0:
      raise HDF5ExtError(
        "Problems creating the %s: invalid disk type ID for atom %s" % (
            self.__class__.__name__, atom_))

    # Allocate space for the dimension axis info and fill it
    dims = np.array(shape, dtype=np.intp)
    self.rank = len(shape)
    self.dims = npy_malloc_dims(self.rank, <npy_intp *>PyArray_DATA(dims))
    rbuf = _array_data(nparr)

    # Blosc2 optimized operations cannot be used (no chunking nor filters).
    self.blosc2_support_read = False
    self.blosc2_support_wirte = False

    # Save the array
    complib = (self.filters.complib or '').encode('utf-8')
    version = self._v_version.encode('utf-8')
    class_ = self._c_classid.encode('utf-8')
    self.dataset_id = H5ARRAYmake(self.parent_id, encoded_name, version,
                                  self.rank, self.dims,
                                  self.extdim, self.disk_type_id, NULL, NULL,
                                  self.filters.complevel, complib,
                                  self.filters.shuffle_bitshuffle,
                                  self.filters.fletcher32,
                                  self._want_track_times,
                                  rbuf)
    if self.dataset_id < 0:
      raise HDF5ExtError("Problems creating the %s." % self.__class__.__name__)

    if self._v_file.params['PYTABLES_SYS_ATTRS']:
      cset = H5T_CSET_UTF8
      # Set the conforming array attributes
      H5ATTRset_attribute_string(self.dataset_id, "CLASS", class_,
                                 len(class_), cset)
      H5ATTRset_attribute_string(self.dataset_id, "VERSION", version,
                                 len(version), cset)
      H5ATTRset_attribute_string(self.dataset_id, "TITLE", encoded_title,
                                 len(encoded_title), cset)

    # Get the native type (so that it is HDF5 who is the responsible to deal
    # with non-native byteorders on-disk)
    self.type_id = get_native_type(self.disk_type_id)

    return self.dataset_id, shape, atom_


  def _create_carray(self, object title):
    cdef int i
    cdef herr_t ret
    cdef void *rbuf
    cdef bytes complib, version, class_
    cdef ndarray dflts
    cdef void *fill_data
    cdef ndarray extdim
    cdef object atom
    cdef bytes encoded_title, encoded_name

    encoded_title = title.encode('utf-8')
    encoded_name = self.name.encode('utf-8')

    atom = self.atom
    self.disk_type_id = atom_to_hdf5_type(atom, self.byteorder)

    self.rank = len(self.shape)
    self.dims = malloc_dims(self.shape)
    if self.chunkshape:
      self.dims_chunk = malloc_dims(self.chunkshape)

    # Decide whether Blosc2 optimized operations can be used.
    (self.blosc2_support_read, self.blosc2_support_write) = (
        _supports_opt_blosc2_read_write(self.byteorder, [self.filters.complib],
                                        self._v_file.mode))

    rbuf = NULL   # The data pointer. We don't have data to save initially
    # Encode strings
    complib = (self.filters.complib or '').encode('utf-8')
    version = self._v_version.encode('utf-8')
    class_ = self._c_classid.encode('utf-8')

    # Get the fill values
    if isinstance(atom.dflt, np.ndarray) or atom.dflt:
      dflts = np.array(atom.dflt, dtype=atom.dtype)
      fill_data = PyArray_DATA(dflts)
    else:
      dflts = np.zeros((), dtype=atom.dtype)
      fill_data = NULL
    if atom.shape == ():
      # The default is preferred as a scalar value instead of 0-dim array
      atom.dflt = dflts[()]
    else:
      atom.dflt = dflts

    cdef hsize_t blocksize = int(os.environ.get("PT_DEFAULT_B2_BLOCKSIZE", "0"))
    # Create the CArray/EArray
    self.dataset_id = H5ARRAYOmake(self.parent_id, encoded_name, version,
                                  self.rank, self.dims, self.extdim,
                                  self.disk_type_id, self.dims_chunk,
                                  blocksize, fill_data,
                                  self.filters.complevel, complib,
                                  self.filters.shuffle_bitshuffle,
                                  self.filters.fletcher32,
                                  self._want_track_times,
                                  rbuf)
    if self.dataset_id < 0:
      raise HDF5ExtError("Problems creating the %s." % self.__class__.__name__)

    if self._v_file.params['PYTABLES_SYS_ATTRS']:
      # Set the conforming array attributes
      H5ATTRset_attribute_string(self.dataset_id, "CLASS", class_,
                                 len(class_), H5T_CSET_ASCII)
      H5ATTRset_attribute_string(self.dataset_id, "VERSION", version,
                                 len(version), H5T_CSET_ASCII)
      H5ATTRset_attribute_string(self.dataset_id, "TITLE", encoded_title,
                                 len(encoded_title), H5T_CSET_ASCII)
      if self.extdim >= 0:
        extdim = <ndarray>np.array([self.extdim], dtype="int32")
        # Attach the EXTDIM attribute in case of enlargeable arrays
        H5ATTRset_attribute(self.dataset_id, "EXTDIM", H5T_NATIVE_INT,
                            0, NULL, PyArray_BYTES(extdim))

    # Get the native type (so that it is HDF5 who is the responsible to deal
    # with non-native byteorders on-disk)
    self.type_id = get_native_type(self.disk_type_id)

    return self.dataset_id


  def _open_array(self):
    cdef size_t type_size, type_precision
    cdef H5T_class_t class_id
    cdef char cbyteorder[11]  # "irrelevant" fits easily here
    cdef int i
    cdef int extdim
    cdef herr_t ret
    cdef object shape, chunkshapes, atom
    cdef int fill_status
    cdef ndarray dflts
    cdef void *fill_data
    cdef bytes encoded_name
    cdef str byteorder

    encoded_name = self.name.encode('utf-8')

    # Open the dataset
    self.dataset_id = H5Dopen(self.parent_id, encoded_name, H5P_DEFAULT)
    if self.dataset_id < 0:
      raise HDF5ExtError("Non-existing node ``%s`` under ``%s``" %
                         (self.name, self._v_parent._v_pathname))
    # Get the datatype handles
    self.disk_type_id, self.type_id = self._get_type_ids()
    # Get the atom for this type
    atom = atom_from_hdf5_type(self.type_id)

    # Get the rank for this array object
    if H5ARRAYget_ndims(self.dataset_id, &self.rank) < 0:
      raise HDF5ExtError("Problems getting ndims!")
    # Allocate space for the dimension axis info
    self.dims = <hsize_t *>malloc(self.rank * sizeof(hsize_t))
    self.maxdims = <hsize_t *>malloc(self.rank * sizeof(hsize_t))
    # Get info on dimensions, class and type (of base class)
    ret = H5ARRAYget_info(self.dataset_id, self.disk_type_id,
                          self.dims, self.maxdims,
                          &class_id, cbyteorder)
    if ret < 0:
      raise HDF5ExtError("Unable to get array info.")

    byteorder = cstr_to_pystr(cbyteorder)

    # Get the extendable dimension (if any)
    self.extdim = -1  # default is non-extensible Array
    for i from 0 <= i < self.rank:
      if self.maxdims[i] == <hsize_t>-1:
        self.extdim = i
        break

    # Get the shape as a python tuple
    shape = getshape(self.rank, self.dims)

    # Allocate space for the dimension chunking info
    self.dims_chunk = <hsize_t *>malloc(self.rank * sizeof(hsize_t))
    if H5ARRAYget_chunkshape(self.dataset_id, self.rank, self.dims_chunk) < 0:
      # The Array class is not chunked!
      chunkshapes = None
      # Blosc2 optimized operations cannot be used (no chunking nor filters).
      self.blosc2_support_read = False
      self.blosc2_support_write = False
    else:
      # Get the chunkshape as a python tuple
      chunkshapes = getshape(self.rank, self.dims_chunk)
      # Decide whether Blosc2 optimized operations can be used.
      filters = get_filters(self.parent_id, self.name) or {}
      (self.blosc2_support_read, self.blosc2_support_write) = (
          _supports_opt_blosc2_read_write(byteorder, list(filters),
                                          self._v_file.mode))

    # object arrays should not be read directly into memory
    if atom.dtype != object:
      # Get the fill value
      dflts = np.zeros((), dtype=atom.dtype)
      fill_data = PyArray_DATA(dflts)
      H5ARRAYget_fill_value(self.dataset_id, self.type_id,
                            &fill_status, fill_data);
      if fill_status == H5D_FILL_VALUE_UNDEFINED:
        # This can only happen with datasets created with other libraries
        # than PyTables.
        dflts = None
      if dflts is not None and atom.shape == ():
        # The default is preferred as a scalar value instead of 0-dim array
        atom.dflt = dflts[()]
      else:
        atom.dflt = dflts

    # Get the byteorder
    self.byteorder = correct_byteorder(atom.type, byteorder)

    return self.dataset_id, atom, shape, chunkshapes


  def _append(self, ndarray nparr):
    cdef int ret, extdim
    cdef hsize_t *dims_arr
    cdef void *rbuf
    cdef object shape

    if self.atom.kind == "reference":
      raise ValueError("Cannot append to the reference types")

    # Allocate space for the dimension axis info
    dims_arr = npy_malloc_dims(self.rank, PyArray_DIMS(nparr))
    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(nparr)
    # Convert some NumPy types to HDF5 before storing.
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 0)

    # Append the records
    extdim = self.extdim
    with nogil:
        ret = H5ARRAYappend_records(self.dataset_id, self.type_id, self.rank,
                                    self.dims, dims_arr, extdim, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems appending the elements")

    free(dims_arr)
    # Update the new dimensionality
    shape = list(self.shape)
    shape[self.extdim] = SizeType(self.dims[self.extdim])
    self.shape = tuple(shape)

  def _read_array(self, hsize_t start, hsize_t stop, hsize_t step,
                 ndarray nparr):
    cdef herr_t ret
    cdef void *rbuf
    cdef hsize_t nrows
    cdef int extdim
    cdef size_t item_size = H5Tget_size(self.type_id)
    cdef void * refbuf = NULL

    # Number of rows to read
    nrows = get_len_of_range(start, stop, step)

    # Get the pointer to the buffer data area
    if self.atom.kind == "reference":
      refbuf = malloc(nrows * item_size)
      rbuf = refbuf
    else:
      rbuf = PyArray_DATA(nparr)

    if hasattr(self, "extdim"):
      extdim = self.extdim
    else:
      extdim = -1

    # Do the physical read
    with nogil:
        ret = H5ARRAYread(self.dataset_id, self.type_id, start, nrows, step,
                          extdim, rbuf)

    try:
      if ret < 0:
        raise HDF5ExtError("Problems reading the array data.")

      # Get the pointer to the buffer data area
      if self.atom.kind == "reference":
        load_reference(self.dataset_id, <hobj_ref_t *>rbuf, item_size, nparr)
    finally:
      if refbuf:
        free(refbuf)
        refbuf = NULL

    if self.atom.kind == 'time':
      # Swap the byteorder by hand (this is not currently supported by HDF5)
      if H5Tget_order(self.type_id) != platform_byteorder:
        nparr.byteswap(True)

    # Convert some HDF5 types to NumPy after reading.
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 1)

    return


  def _g_read_slice(self, ndarray startl, ndarray stopl, ndarray stepl,
                   ndarray nparr):
    cdef herr_t ret
    cdef hsize_t *start
    cdef hsize_t *stop
    cdef hsize_t *step
    cdef void *rbuf
    cdef size_t item_size = H5Tget_size(self.type_id)
    cdef void * refbuf = NULL

    # Get the pointer to the buffer data area of startl, stopl and stepl arrays
    start = <hsize_t *>PyArray_DATA(startl)
    stop = <hsize_t *>PyArray_DATA(stopl)
    step = <hsize_t *>PyArray_DATA(stepl)

    # Get the pointer to the buffer data area
    if self.atom.kind == "reference":
      refbuf = malloc(nparr.size * item_size)
      rbuf = refbuf
    else:
      rbuf = PyArray_DATA(nparr)

    cdef bytes fname = self._v_file.filename.encode('utf8')
    cdef char *filename = fname
    # Do the physical read
    with nogil:
        ret = H5ARRAYOreadSlice(filename, self.blosc2_support_read, self.dataset_id, self.type_id,
                                start, stop, step, rbuf)
    try:
      if ret < 0:
        raise HDF5ExtError("Internal error reading the elements "
                           "(H5ARRAYOreadSlice returned errorcode %i)" % ret)

      # Get the pointer to the buffer data area
      if self.atom.kind == "reference":
        load_reference(self.dataset_id, <hobj_ref_t *>rbuf, item_size, nparr)
    finally:
      if refbuf:
        free(refbuf)
        refbuf = NULL

    if self.atom.kind == 'time':
      # Swap the byteorder by hand (this is not currently supported by HDF5)
      if H5Tget_order(self.type_id) != platform_byteorder:
        nparr.byteswap(True)

    # Convert some HDF5 types to NumPy after reading
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 1)

    return


  def _g_read_coords(self, ndarray coords, ndarray nparr):
    """Read coordinates in an already created NumPy array."""

    cdef herr_t ret
    cdef hid_t space_id
    cdef hid_t mem_space_id
    cdef hsize_t size
    cdef void *rbuf
    cdef object mode
    cdef size_t item_size = H5Tget_size(self.type_id)
    cdef void * refbuf = NULL

    # Get the dataspace handle
    space_id = H5Dget_space(self.dataset_id)
    # Create a memory dataspace handle
    size = nparr.size
    mem_space_id = H5Screate_simple(1, &size, NULL)

    # Select the dataspace to be read
    H5Sselect_elements(space_id, H5S_SELECT_SET,
                       <size_t>size, <hsize_t *>PyArray_DATA(coords))

    # Get the pointer to the buffer data area
    if self.atom.kind == "reference":
      refbuf = malloc(nparr.size * item_size)
      rbuf = refbuf
    else:
      rbuf = PyArray_DATA(nparr)

    # Do the actual read
    with nogil:
        ret = H5Dread(self.dataset_id, self.type_id, mem_space_id, space_id,
                      H5P_DEFAULT, rbuf)

    try:
      if ret < 0:
        raise HDF5ExtError("Problems reading the array data.")

      # Get the pointer to the buffer data area
      if self.atom.kind == "reference":
        load_reference(self.dataset_id, <hobj_ref_t *>rbuf, item_size, nparr)
    finally:
      if refbuf:
        free(refbuf)
        refbuf = NULL

    # Terminate access to the memory dataspace
    H5Sclose(mem_space_id)
    # Terminate access to the dataspace
    H5Sclose(space_id)

    if self.atom.kind == 'time':
      # Swap the byteorder by hand (this is not currently supported by HDF5)
      if H5Tget_order(self.type_id) != platform_byteorder:
        nparr.byteswap(True)

    # Convert some HDF5 types to NumPy after reading
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 1)

    return


  def perform_selection(self, space_id, start, count, step, idx, mode):
    """Performs a selection using start/count/step in the given axis.

    All other axes have their full range selected.  The selection is
    added to the current `space_id` selection using the given mode.

    Note: This is a backport from the h5py project.

    """

    cdef int select_mode
    cdef ndarray start_, count_, step_
    cdef hsize_t *startp
    cdef hsize_t *countp
    cdef hsize_t *stepp

    # Build arrays for the selection parameters
    startl, countl, stepl = [], [], []
    for i, x in enumerate(self.shape):
      if i != idx:
        startl.append(0)
        countl.append(x)
        stepl.append(1)
      else:
        startl.append(start)
        countl.append(count)
        stepl.append(step)
    start_ = np.array(startl, dtype="i8")
    count_ = np.array(countl, dtype="i8")
    step_ = np.array(stepl, dtype="i8")

    # Get the pointers to array data
    startp = <hsize_t *>PyArray_DATA(start_)
    countp = <hsize_t *>PyArray_DATA(count_)
    stepp = <hsize_t *>PyArray_DATA(step_)

    # Do the actual selection
    select_modes = {"AND": H5S_SELECT_AND, "NOTB": H5S_SELECT_NOTB}
    assert mode in select_modes
    select_mode = select_modes[mode]
    H5Sselect_hyperslab(space_id, <H5S_seloper_t>select_mode,
                        startp, stepp, countp, NULL)

  def _g_read_selection(self, object selection, ndarray nparr):
    """Read a selection in an already created NumPy array."""

    cdef herr_t ret
    cdef hid_t space_id
    cdef hid_t mem_space_id
    cdef hsize_t size
    cdef void *rbuf
    cdef object mode
    cdef size_t item_size = H5Tget_size(self.type_id)
    cdef void * refbuf = NULL

    # Get the dataspace handle
    space_id = H5Dget_space(self.dataset_id)
    # Create a memory dataspace handle
    size = nparr.size
    mem_space_id = H5Screate_simple(1, &size, NULL)

    # Select the dataspace to be read
    # Start by selecting everything
    H5Sselect_all(space_id)
    # Now refine with outstanding selections
    for args in selection:
      self.perform_selection(space_id, *args)

    # Get the pointer to the buffer data area
    if self.atom.kind == "reference":
      refbuf = malloc(nparr.size * item_size)
      rbuf = refbuf
    else:
      rbuf = PyArray_DATA(nparr)

    # Do the actual read
    with nogil:
        ret = H5Dread(self.dataset_id, self.type_id, mem_space_id, space_id,
                      H5P_DEFAULT, rbuf)

    try:
      if ret < 0:
        raise HDF5ExtError("Problems reading the array data.")

      # Get the pointer to the buffer data area
      if self.atom.kind == "reference":
        load_reference(self.dataset_id, <hobj_ref_t *>rbuf, item_size, nparr)
    finally:
      if refbuf:
        free(refbuf)
        refbuf = NULL

    # Terminate access to the memory dataspace
    H5Sclose(mem_space_id)
    # Terminate access to the dataspace
    H5Sclose(space_id)

    if self.atom.kind == 'time':
      # Swap the byteorder by hand (this is not currently supported by HDF5)
      if H5Tget_order(self.type_id) != platform_byteorder:
        nparr.byteswap(True)

    # Convert some HDF5 types to NumPy after reading
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 1)

    return


  def _g_write_slice(self, ndarray startl, ndarray stepl, ndarray countl,
                    ndarray nparr):
    """Write a slice in an already created NumPy array."""

    cdef int ret
    cdef void *rbuf
    cdef void *temp
    cdef hsize_t *start
    cdef hsize_t *step
    cdef hsize_t *count

    if self.atom.kind == "reference":
      raise ValueError("Cannot write reference types yet")
    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(nparr)
    # Get the start, step and count values
    start = <hsize_t *>PyArray_DATA(startl)
    step = <hsize_t *>PyArray_DATA(stepl)
    count = <hsize_t *>PyArray_DATA(countl)

    # Convert some NumPy types to HDF5 before storing.
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 0)

    # Modify the elements:
    with nogil:
        ret = H5ARRAYwrite_records(self.dataset_id, self.type_id, self.rank,
                                    start, step, count, rbuf)

    if ret < 0:
      raise HDF5ExtError("Internal error modifying the elements "
                         "(H5ARRAYwrite_records returned errorcode %i)" % ret)

    return


  def _g_write_coords(self, ndarray coords, ndarray nparr):
    """Write a selection in an already created NumPy array."""

    cdef herr_t ret
    cdef hid_t space_id
    cdef hid_t mem_space_id
    cdef hsize_t size
    cdef void *rbuf
    cdef object mode

    if self.atom.kind == "reference":
      raise ValueError("Cannot write reference types yet")
    # Get the dataspace handle
    space_id = H5Dget_space(self.dataset_id)
    # Create a memory dataspace handle
    size = nparr.size
    mem_space_id = H5Screate_simple(1, &size, NULL)

    # Select the dataspace to be written
    H5Sselect_elements(space_id, H5S_SELECT_SET,
                       <size_t>size, <hsize_t *>PyArray_DATA(coords))

    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(nparr)

    # Convert some NumPy types to HDF5 before storing.
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 0)

    # Do the actual write
    with nogil:
        ret = H5Dwrite(self.dataset_id, self.type_id, mem_space_id, space_id,
                       H5P_DEFAULT, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems writing the array data.")

    # Terminate access to the memory dataspace
    H5Sclose(mem_space_id)
    # Terminate access to the dataspace
    H5Sclose(space_id)

    return


  def _g_write_selection(self, object selection, ndarray nparr):
    """Write a selection in an already created NumPy array."""

    cdef herr_t ret
    cdef hid_t space_id
    cdef hid_t mem_space_id
    cdef hsize_t size
    cdef void *rbuf
    cdef object mode

    if self.atom.kind == "reference":
      raise ValueError("Cannot write reference types yet")
    # Get the dataspace handle
    space_id = H5Dget_space(self.dataset_id)
    # Create a memory dataspace handle
    size = nparr.size
    mem_space_id = H5Screate_simple(1, &size, NULL)

    # Select the dataspace to be written
    # Start by selecting everything
    H5Sselect_all(space_id)
    # Now refine with outstanding selections
    for args in selection:
      self.perform_selection(space_id, *args)

    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(nparr)

    # Convert some NumPy types to HDF5 before storing.
    if self.atom.type == 'time64':
      self._convert_time64(nparr, 0)

    # Do the actual write
    with nogil:
        ret = H5Dwrite(self.dataset_id, self.type_id, mem_space_id, space_id,
                       H5P_DEFAULT, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems writing the array data.")

    # Terminate access to the memory dataspace
    H5Sclose(mem_space_id)
    # Terminate access to the dataspace
    H5Sclose(space_id)

    return


  def __dealloc__(self):
    if self.dims:
      free(<void *>self.dims)
    if self.maxdims:
      free(<void *>self.maxdims)
    if self.dims_chunk:
      free(self.dims_chunk)


cdef class VLArray(Leaf):
  # Instance variables
  cdef hsize_t nrecords

  def _create_array(self, object title):
    cdef int rank
    cdef hsize_t *dims
    cdef herr_t ret
    cdef void *rbuf
    cdef bytes complib, version, class_
    cdef object type_, itemsize, atom, scatom
    cdef bytes encoded_title, encoded_name
    cdef H5T_cset_t cset = H5T_CSET_ASCII

    encoded_title = title.encode('utf-8')
    encoded_name = self.name.encode('utf-8')

    atom = self.atom
    if not hasattr(atom, 'size'):  # it is a pseudo-atom
      atom = atom.base

    # Get the HDF5 type of the *scalar* atom
    scatom = atom.copy(shape=())
    self.base_type_id = atom_to_hdf5_type(scatom, self.byteorder)
    if self.base_type_id < 0:
      raise HDF5ExtError(
        "Problems creating the %s: invalid base type ID for atom %s" % (
            self.__class__.__name__, scatom))

    # Allocate space for the dimension axis info
    rank = len(atom.shape)
    dims = malloc_dims(atom.shape)

    rbuf = NULL   # We don't have data to save initially

    # Encode strings
    complib = (self.filters.complib or '').encode('utf-8')
    version = self._v_version.encode('utf-8')
    class_ = self._c_classid.encode('utf-8')

    # Create the vlarray
    self.dataset_id = H5VLARRAYmake(self.parent_id, encoded_name, version,
                                    rank, dims, self.base_type_id,
                                    self.chunkshape[0], rbuf,
                                    self.filters.complevel, complib,
                                    self.filters.shuffle_bitshuffle,
                                    self.filters.fletcher32,
                                    self._want_track_times, rbuf)
    if dims:
      free(<void *>dims)
    if self.dataset_id < 0:
      raise HDF5ExtError("Problems creating the VLArray.")
    self.nrecords = 0  # Initialize the number of records saved

    if self._v_file.params['PYTABLES_SYS_ATTRS']:
      cset = H5T_CSET_UTF8
      # Set the conforming array attributes
      H5ATTRset_attribute_string(self.dataset_id, "CLASS", class_,
                                 len(class_), cset)
      H5ATTRset_attribute_string(self.dataset_id, "VERSION", version,
                                 len(version), cset)
      H5ATTRset_attribute_string(self.dataset_id, "TITLE", encoded_title,
                                 len(encoded_title), cset)

    # Get the datatype handles
    self.disk_type_id, self.type_id = self._get_type_ids()

    return self.dataset_id


  def _open_array(self):
    cdef char cbyteorder[11]  # "irrelevant" fits easily here
    cdef int i, enumtype
    cdef int rank
    cdef herr_t ret
    cdef hsize_t nrecords, chunksize
    cdef object shape, type_
    cdef bytes encoded_name
    cdef str byteorder

    encoded_name = self.name.encode('utf-8')

    # Open the dataset
    self.dataset_id = H5Dopen(self.parent_id, encoded_name, H5P_DEFAULT)
    if self.dataset_id < 0:
      raise HDF5ExtError("Non-existing node ``%s`` under ``%s``" %
                         (self.name, self._v_parent._v_pathname))
    # Get the datatype handles
    self.disk_type_id, self.type_id = self._get_type_ids()
    # Get the atom for this type
    atom = atom_from_hdf5_type(self.type_id)

    # Get info on dimensions & types (of base class)
    H5VLARRAYget_info(self.dataset_id, self.disk_type_id, &nrecords,
                      cbyteorder)

    byteorder = cstr_to_pystr(cbyteorder)

    # Get some properties of the atomic type
    self._atomicdtype = atom.dtype
    self._atomictype = atom.type
    self._atomicshape = atom.shape
    self._atomicsize = atom.size

    # Get the byteorder
    self.byteorder = correct_byteorder(atom.type, byteorder)

    # Get the chunkshape (VLArrays are unidimensional entities)
    H5ARRAYget_chunkshape(self.dataset_id, 1, &chunksize)

    self.nrecords = nrecords  # Initialize the number of records saved
    return self.dataset_id, SizeType(nrecords), (SizeType(chunksize),), atom


  def _append(self, ndarray nparr, int nobjects):
    cdef int ret
    cdef void *rbuf

    # Get the pointer to the buffer data area
    if nobjects:
      rbuf = PyArray_DATA(nparr)
      # Convert some NumPy types to HDF5 before storing.
      if self.atom.type == 'time64':
        self._convert_time64(nparr, 0)
    else:
      rbuf = NULL

    # Append the records:
    with nogil:
        ret = H5VLARRAYappend_records(self.dataset_id, self.type_id,
                                      nobjects, self.nrecords, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems appending the records.")

    self.nrecords = self.nrecords + 1

  def _modify(self, hsize_t nrow, ndarray nparr, int nobjects):
    cdef int ret
    cdef void *rbuf

    # Get the pointer to the buffer data area
    rbuf = PyArray_DATA(nparr)
    if nobjects:
      # Convert some NumPy types to HDF5 before storing.
      if self.atom.type == 'time64':
        self._convert_time64(nparr, 0)

    # Append the records:
    with nogil:
        ret = H5VLARRAYmodify_records(self.dataset_id, self.type_id,
                                      nrow, nobjects, rbuf)

    if ret < 0:
      raise HDF5ExtError("Problems modifying the record.")

    return nobjects

  # Because the size of each "row" is unknown, there is no easy way to
  # calculate this value
  def _get_memory_size(self):
    cdef hid_t space_id
    cdef hsize_t size
    cdef herr_t ret

    if self.nrows == 0:
      size = 0
    else:
      # Get the dataspace handle
      space_id = H5Dget_space(self.dataset_id)
      # Return the size of the entire dataset
      ret = H5Dvlen_get_buf_size(self.dataset_id, self.type_id, space_id,
                                 &size)
      if ret < 0:
        size = -1

      # Terminate access to the dataspace
      H5Sclose(space_id)

    return size

  def _read_array(self, hsize_t start, hsize_t stop, hsize_t step):
    cdef int i
    cdef size_t vllen
    cdef herr_t ret
    cdef hvl_t *rdata
    cdef hsize_t nrows
    cdef hid_t space_id
    cdef hid_t mem_space_id
    cdef object buf, nparr, shape, datalist

    # Compute the number of rows to read
    nrows = get_len_of_range(start, stop, step)
    if start + nrows > self.nrows:
      raise HDF5ExtError(
        "Asking for a range of rows exceeding the available ones!.",
        h5bt=False)

    # Now, read the chunk of rows
    with nogil:
        # Allocate the necessary memory for keeping the row handlers
        rdata = <hvl_t *>malloc(<size_t>nrows*sizeof(hvl_t))
        # Get the dataspace handle
        space_id = H5Dget_space(self.dataset_id)
        # Create a memory dataspace handle
        mem_space_id = H5Screate_simple(1, &nrows, NULL)
        # Select the data to be read
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, &start, &step, &nrows,
                            NULL)
        # Do the actual read
        ret = H5Dread(self.dataset_id, self.type_id, mem_space_id, space_id,
                      H5P_DEFAULT, rdata)

    if ret < 0:
      raise HDF5ExtError(
        "VLArray._read_array: Problems reading the array data.")

    datalist = []
    for i in range(<long long>nrows):
      # Number of atoms in row
      vllen = rdata[i].len
      # Get the pointer to the buffer data area
      if vllen > 0:
        # Create a buffer to keep this info. It is important to do a
        # copy, because we will dispose the buffer memory later on by
        # calling the H5Dvlen_reclaim. PyByteArray_FromStringAndSize does this.
        buf = PyByteArray_FromStringAndSize(<char *>rdata[i].p,
                                            vllen*self._atomicsize)
      else:
        # Case where there is info with zero lentgh
        buf = None
      # Compute the shape for the read array
      shape = list(self._atomicshape)
      shape.insert(0, vllen)  # put the length at the beginning of the shape
      nparr = np.ndarray(
        buffer=buf, dtype=self._atomicdtype.base, shape=shape)
      # Set the writeable flag for this ndarray object
      nparr.flags.writeable = True
      if self.atom.kind == 'time':
        # Swap the byteorder by hand (this is not currently supported by HDF5)
        if H5Tget_order(self.type_id) != platform_byteorder:
          nparr.byteswap(True)
      # Convert some HDF5 types to NumPy after reading.
      if self.atom.type == 'time64':
        self._convert_time64(nparr, 1)
      # Append this array to the output list
      datalist.append(nparr)

    # Release resources
    # Reclaim all the (nested) VL data
    ret = H5Dvlen_reclaim(self.type_id, mem_space_id, H5P_DEFAULT, rdata)
    if ret < 0:
      raise HDF5ExtError("VLArray._read_array: error freeing the data buffer.")
    # Terminate access to the memory dataspace
    H5Sclose(mem_space_id)
    # Terminate access to the dataspace
    H5Sclose(space_id)
    # Free the amount of row pointers to VL row data
    free(rdata)

    return datalist


  def get_row_size(self, row):
    """Return the total size in bytes of all the elements contained in a given row."""

    cdef hid_t space_id
    cdef hsize_t size
    cdef herr_t ret

    cdef hsize_t offset[1]
    cdef hsize_t count[1]

    if row >= self.nrows:
      raise HDF5ExtError(
        "Asking for a range of rows exceeding the available ones!.",
        h5bt=False)

    # Get the dataspace handle
    space_id = H5Dget_space(self.dataset_id)

    offset[0] = row
    count[0] = 1

    ret = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    if ret < 0:
      size = -1

    ret = H5Dvlen_get_buf_size(self.dataset_id, self.type_id, space_id, &size)
    if ret < 0:
      size = -1

    # Terminate access to the dataspace
    H5Sclose(space_id)

    return size


cdef class UnImplemented(Leaf):

  def _open_unimplemented(self):
    cdef object shape
    cdef char cbyteorder[11]  # "irrelevant" fits easily here
    cdef bytes encoded_name
    cdef str byteorder

    encoded_name = self.name.encode('utf-8')

    # Get info on dimensions
    shape = H5UIget_info(self.parent_id, encoded_name, cbyteorder)
    shape = tuple(map(SizeType, shape))
    self.dataset_id = H5Dopen(self.parent_id, encoded_name, H5P_DEFAULT)
    byteorder = cstr_to_pystr(cbyteorder)

    return (shape, byteorder, self.dataset_id)

  def _g_close(self):
    H5Dclose(self.dataset_id)


## Local Variables:
## mode: python
## py-indent-offset: 2
## tab-width: 2
## fill-column: 78
## End:
