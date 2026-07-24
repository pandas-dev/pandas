########################################################################
#
# License: BSD
# Created: November 25, 2009
# Author:  Francesc Alted - faltet@pytables.com
#
# $Id$
#
########################################################################

"""Cython functions and classes for supporting links in HDF5."""

from .exceptions import HDF5ExtError

from libc.stdlib cimport malloc, free
from libc.string cimport strlen
from cpython.unicode cimport PyUnicode_DecodeUTF8

from .definitions cimport (
    H5P_DEFAULT,
    hid_t,
    herr_t,
    hbool_t,
    int64_t,
    H5T_cset_t,
    haddr_t,
)
from .hdf5extension cimport Node
from .utilsextension cimport cstr_to_pystr

#----------------------------------------------------------------------

# External declarations

cdef extern from "H5Lpublic.h" nogil:

  ctypedef enum H5L_type_t:
    H5L_TYPE_ERROR = (-1),       # Invalid link type id
    H5L_TYPE_HARD = 0,           # Hard link id
    H5L_TYPE_SOFT = 1,           # Soft link id
    H5L_TYPE_EXTERNAL = 64,      # External link id
    H5L_TYPE_MAX = 255           # Maximum link type id

  # Information struct for link (for H5Lget_info)
  cdef union _add_u:
    haddr_t address              # Address hard link points to
    size_t val_size              # Size of a soft link or UD link value

  ctypedef struct H5L_info_t:
    H5L_type_t     type          # Type of link
    hbool_t        corder_valid  # Indicate if creation order is valid
    int64_t        corder        # Creation order
    H5T_cset_t     cset          # Character set of link name
    _add_u         u             # Size of a soft link or UD link value

  # Operations with links
  herr_t H5Lcreate_hard(
    hid_t obj_loc_id, char *obj_name, hid_t link_loc_id, char *link_name,
    hid_t lcpl_id, hid_t lapl_id)

  herr_t H5Lcreate_soft(
    char *target_path, hid_t link_loc_id, char *link_name,
    hid_t lcpl_id, hid_t lapl_id)

  herr_t H5Lcreate_external(
    char *file_name, char *object_name, hid_t link_loc_id, char *link_name,
    hid_t lcpl_id, hid_t lapl_id)

  herr_t H5Lget_info(
    hid_t link_loc_id, char *link_name, H5L_info_t *link_buff,
    hid_t lapl_id)

  herr_t H5Lget_val(
    hid_t link_loc_id, char *link_name, void *linkval_buff, size_t size,
    hid_t lapl_id)

  herr_t H5Lunpack_elink_val(
    char *ext_linkval, size_t link_size, unsigned *flags,
    const char **filename, const char **obj_path)

  herr_t H5Lcopy(
    hid_t src_loc_id, char *src_name, hid_t dest_loc_id, char *dest_name,
    hid_t lcpl_id, hid_t lapl_id)


#----------------------------------------------------------------------

# Helper functions

def _get_link_class(parent_id, name):
    """Guess the link class."""

    cdef herr_t ret
    cdef H5L_info_t link_buff
    cdef H5L_type_t link_type

    ret = H5Lget_info(parent_id, name, &link_buff, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to get info about link")

    link_type = link_buff.type
    if link_type == H5L_TYPE_SOFT:
      return "SoftLink"
    elif link_type == H5L_TYPE_EXTERNAL:
      return "ExternalLink"
    # elif link_type == H5L_TYPE_HARD:
    #   return "HardLink"
    else:
      return "UnImplemented"




def _g_create_hard_link(parentnode, str name, targetnode):
  """Create a hard link in the file."""

  cdef herr_t ret
  cdef bytes encoded_name = name.encode('utf-8')
  cdef bytes encoded_v_name = targetnode._v_name.encode('utf-8')

  ret = H5Lcreate_hard(targetnode._v_parent._v_objectid, encoded_v_name,
                       parentnode._v_objectid, <char*>encoded_name,
                       H5P_DEFAULT, H5P_DEFAULT)
  if ret < 0:
    raise HDF5ExtError("failed to create HDF5 hard link")




#----------------------------------------------------------------------

# Public classes

cdef class Link(Node):
  """Extension class from which all link extensions inherits."""

  def _g_copy(self, newparent, newname, recursive, _log=True, **kwargs):
    """Private part for the _f_copy() method."""

    cdef herr_t ret
    cdef object stats
    cdef bytes encoded_name, encoded_newname

    encoded_name = self.name.encode('utf-8')
    encoded_newname = newname.encode('utf-8')

    # @TODO: set property list --> utf-8
    ret = H5Lcopy(self.parent_id, encoded_name, newparent._v_objectid,
                  encoded_newname, H5P_DEFAULT, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to copy HDF5 link")

    # Update statistics if needed.
    stats = kwargs.get('stats', None)
    if stats is not None:
      stats['links'] += 1

    return newparent._v_file.get_node(newparent, newname)


cdef class SoftLink(Link):
  """Extension class representing a soft link."""

  def _g_create(self):
    """Create the link in file."""

    cdef herr_t ret
    cdef bytes encoded_name = self.name.encode('utf-8')
    cdef bytes encoded_target = self.target.encode('utf-8')

    ret = H5Lcreate_soft(encoded_target, self.parent_id, encoded_name,
                         H5P_DEFAULT, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to create HDF5 soft link")

    return 0  # Object ID is zero'ed, as HDF5 does not assign one for links

  def _g_open(self):
    """Open the link in file."""

    cdef herr_t ret
    cdef H5L_info_t link_buff
    cdef size_t val_size
    cdef char *clinkval
    cdef bytes encoded_name

    encoded_name = self.name.encode('utf-8')

    ret = H5Lget_info(self.parent_id, encoded_name, &link_buff, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to get info about soft link")

    val_size = link_buff.u.val_size
    clinkval = <char *>malloc(val_size)

    ret = H5Lget_val(self.parent_id, encoded_name, clinkval, val_size,
                     H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to get target value")

    self.target = PyUnicode_DecodeUTF8(clinkval, strlen(clinkval), NULL)

    # Release resources
    free(clinkval)
    return 0  # Object ID is zero'ed, as HDF5 does not assign one for links


cdef class ExternalLink(Link):
  """Extension class representing an external link."""

  def _g_create(self):
    """Create the link in file."""

    cdef herr_t ret
    cdef bytes encoded_name, encoded_filename, encoded_target

    encoded_name = self.name.encode('utf-8')

    filename, target = self._get_filename_node()
    encoded_filename = filename.encode('utf-8')
    encoded_target = target.encode('utf-8')

    ret = H5Lcreate_external(<char*>encoded_filename, <char*>encoded_target,
                             self.parent_id, <char*>encoded_name,
                             H5P_DEFAULT, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to create HDF5 external link")

    return 0  # Object ID is zero'ed, as HDF5 does not assign one for links

  def _g_open(self):
    """Open the link in file."""

    cdef herr_t ret
    cdef H5L_info_t link_buff
    cdef size_t val_size
    cdef char *clinkval
    cdef char *cfilename
    cdef char *c_obj_path
    cdef unsigned flags
    cdef bytes encoded_name
    cdef str filename, obj_path

    encoded_name = self.name.encode('utf-8')

    ret = H5Lget_info(self.parent_id, encoded_name, &link_buff, H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to get info about external link")

    val_size = link_buff.u.val_size
    clinkval = <char *>malloc(val_size)

    ret = H5Lget_val(self.parent_id, encoded_name, clinkval, val_size,
                     H5P_DEFAULT)
    if ret < 0:
      raise HDF5ExtError("failed to get target value")

    ret = H5Lunpack_elink_val(clinkval, val_size, &flags,
                              <const char **>&cfilename,
                              <const char **>&c_obj_path)
    if ret < 0:
      raise HDF5ExtError("failed to unpack external link value")

    filename = cstr_to_pystr(cfilename)
    obj_path = cstr_to_pystr(c_obj_path)

    self.target = filename+':'+obj_path

    # Release resources
    free(clinkval)
    return 0  # Object ID is zero'ed, as HDF5 does not assign one for links

  def _get_obj_info(self):
    # ExternalLink do not have ObjectId. Hardcode addr and rc to 0, 1
    return 0, 1


## Local Variables:
## mode: python
## py-indent-offset: 2
## tab-width: 2
## fill-column: 78
## End:
