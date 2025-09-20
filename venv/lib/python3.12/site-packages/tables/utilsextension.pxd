########################################################################
#
#       License: BSD
#       Created: March 03, 2008
#       Author:  Francesc Alted - faltet@pytables.com
#
#       $Id: definitions.pyd 1018 2005-06-20 09:43:34Z faltet $
#
########################################################################

"""
These are declarations for functions in utilsextension.pyx that have to
be shared with other extensions.
"""

from numpy cimport ndarray

from .definitions cimport hsize_t, hid_t, hobj_ref_t


cdef hsize_t *malloc_dims(object)
cdef hid_t get_native_type(hid_t) nogil
cdef str cstr_to_pystr(const char*)
cdef int load_reference(hid_t dataset_id, hobj_ref_t *refbuf, size_t item_size, ndarray nparr) except -1
