"""Here is defined the AttributeSet class."""

from __future__ import annotations

import re
import pickle
import warnings
from typing import Any, Literal, TYPE_CHECKING
from collections.abc import Callable

import numpy as np

from . import hdf5extension
from .path import check_attribute_name
from .utils import SizeType
from .filters import Filters
from .registry import class_name_dict
from .undoredo import attr_to_shadow
from .exceptions import ClosedNodeError, FiltersWarning, PerformanceWarning

if TYPE_CHECKING:
    from .node import Node


# System attributes
SYS_ATTRS = [
    "CLASS",
    "VERSION",
    "TITLE",
    "NROWS",
    "EXTDIM",
    "ENCODING",
    "PYTABLES_FORMAT_VERSION",
    "FLAVOR",
    "FILTERS",
    "AUTO_INDEX",
    "DIRTY",
    "NODE_TYPE",
    "NODE_TYPE_VERSION",
    "PSEUDOATOM",
]
# Prefixes of other system attributes
SYS_ATTRS_PREFIXES = ["FIELD_"]
# RO_ATTRS will be disabled and let the user modify them if they
# want to. The user is still not allowed to remove or rename
# system attributes. Francesc Alted 2004-12-19
# Read-only attributes:
# RO_ATTRS = ["CLASS", "FLAVOR", "VERSION", "NROWS", "EXTDIM",
#             "PYTABLES_FORMAT_VERSION", "FILTERS",
#             "NODE_TYPE", "NODE_TYPE_VERSION"]
# RO_ATTRS = []

# The next attributes are not meant to be copied during a Node copy process
SYS_ATTRS_NOTTOBECOPIED = [
    "CLASS",
    "VERSION",
    "TITLE",
    "NROWS",
    "EXTDIM",
    "PYTABLES_FORMAT_VERSION",
    "FILTERS",
    "ENCODING",
]
# Attributes forced to be copied during node copies
FORCE_COPY_CLASS = ["CLASS", "VERSION"]
# Regular expression for column default values.
_field_fill_re = re.compile("^FIELD_[0-9]+_FILL$")
# Regular expression for fixing old pickled filters.
_old_filters_re = re.compile(rb"\(([ic])tables\.Leaf\n")
# Fixed version of the previous string.
_new_filters_sub = rb"(\1tables.filters\n"


def issysattrname(name: str) -> bool:
    """Check if a name is a system attribute or not."""
    return bool(
        name in SYS_ATTRS
        or np.prod([name.startswith(prefix) for prefix in SYS_ATTRS_PREFIXES])
    )


class AttributeSet(hdf5extension.AttributeSet):
    r"""Container for the HDF5 attributes of a Node.

    This class provides methods to create new HDF5 node attributes,
    and to get, rename or delete existing ones.

    Like in Group instances (see :ref:`GroupClassDescr`), AttributeSet
    instances make use of the *natural naming* convention, i.e. you can
    access the attributes on disk as if they were normal Python
    attributes of the AttributeSet instance.

    This offers the user a very convenient way to access HDF5 node
    attributes. However, for this reason and in order not to pollute the
    object namespace, one can not assign *normal* attributes to
    AttributeSet instances, and their members use names which start by
    special prefixes as happens with Group objects.

    .. rubric:: Notes on native and pickled attributes

    The values of most basic types are saved as HDF5 native data in the
    HDF5 file.  This includes Python bool, int, float, complex and str
    (but not long nor unicode) values, as well as their NumPy scalar
    versions and homogeneous or *structured* NumPy arrays of them.  When
    read, these values are always loaded as NumPy scalar or array
    objects, as needed.

    For that reason, attributes in native HDF5 files will always be
    mapped into NumPy objects.  Specifically, a multidimensional
    attribute will be mapped into a multidimensional ndarray and a
    scalar will be mapped into a NumPy scalar object (for example, a
    scalar H5T_NATIVE_LLONG will be read and returned as a numpy.int64
    scalar).

    However, other kinds of values are serialized using pickle, so you
    only will be able to correctly retrieve them using a Python-aware
    HDF5 library.  Thus, if you want to save Python scalar values and
    make sure you are able to read them with generic HDF5 tools, you
    should make use of *scalar or homogeneous/structured array NumPy
    objects* (for example, numpy.int64(1) or numpy.array([1, 2, 3],
    dtype='int16')).

    One more advice: because of the various potential difficulties in
    restoring a Python object stored in an attribute, you may end up
    getting a pickle string where a Python object is expected. If this
    is the case, you may wish to run pickle.loads() on that string to
    get an idea of where things went wrong, as shown in this example::

        >>> import os, tempfile
        >>> import tables as tb
        >>>
        >>> class MyClass:
        ...     foo = 'bar'
        ...
        >>> myObject = MyClass()  # save object of custom class in HDF5 attr
        >>> h5fname = tempfile.mktemp(suffix='.h5')
        >>> h5f = tb.open_file(h5fname, 'w')
        >>> h5f.root._v_attrs.obj = myObject  # store the object
        >>> print(h5f.root._v_attrs.obj.foo)  # retrieve it
        bar
        >>> h5f.close()
        >>>
        >>> del MyClass, myObject  # delete class of object and reopen file
        >>> h5f = tb.open_file(h5fname, 'r')
        >>> print(repr(h5f.root._v_attrs.obj))
        b'ccopy_reg\\n_reconstructor...
        >>> import pickle  # let's unpickle that to see what went wrong
        >>> pickle.loads(h5f.root._v_attrs.obj)
        Traceback (most recent call last):
        ...
        AttributeError: Can't get attribute 'MyClass' ...
        >>> # So the problem was not in the stored object,
        ... # but in the *environment* where it was restored.
        ... h5f.close()
        >>> os.remove(h5fname)


    .. rubric:: Notes on AttributeSet methods

    Note that this class overrides the __getattr__(), __setattr__(),
    __delattr__() and __dir__() special methods.  This allows you to
    read, assign or delete attributes on disk by just using the next
    constructs::

        leaf.attrs.myattr = 'str attr'    # set a string (native support)
        leaf.attrs.myattr2 = 3            # set an integer (native support)
        leaf.attrs.myattr3 = [3, (1, 2)]  # a generic object (Pickled)
        attrib = leaf.attrs.myattr        # get the attribute ``myattr``
        del leaf.attrs.myattr             # delete the attribute ``myattr``

    In addition, the dictionary-like __getitem__(), __setitem__() and
    __delitem__() methods are available, so you may write things like
    this::

        for name in node._v_attrs._f_list():
            print("name: %s, value: %s" % (name, node._v_attrs[name]))

    Use whatever idiom you prefer to access the attributes.

    Finally, on interactive python sessions you may get autocompletions of
    attributes named as *valid python identifiers* by pressing the  `[Tab]`
    key, or to use the dir() global function.

    If an attribute is set on a target node that already has a large
    number of attributes, a PerformanceWarning will be issued.


    .. rubric:: AttributeSet attributes

    .. attribute:: _v_attrnames

        A list with all attribute names.

    .. attribute:: _v_attrnamessys

        A list with system attribute names.

    .. attribute:: _v_attrnamesuser

        A list with user attribute names.

    .. attribute:: _v_unimplemented

        A list of attribute names with unimplemented native HDF5 types.

    """

    def _g_getnode(self) -> Node:
        return self._v__nodefile._get_node(self._v__nodepath)

    @property
    def _v_node(self) -> Node:
        """:class:`Node` instance this attribute set is associated with."""
        return self._g_getnode()

    def __init__(self, node: Node) -> None:
        """Create the basic structures to keep the attribute information.

        Reads all the HDF5 attributes (if any) on disk for the node "node".

        Parameters
        ----------
        node
            The parent node

        """
        # Refuse to create an instance of an already closed node
        if not node._v_isopen:
            raise ClosedNodeError("the node for attribute set is closed")

        dict_ = self.__dict__

        self._g_new(node)
        dict_["_v__nodefile"] = node._v_file
        dict_["_v__nodepath"] = node._v_pathname
        dict_["_v_attrnames"] = self._g_list_attr(node)
        # The list of unimplemented attribute names
        dict_["_v_unimplemented"] = []

        # Get the file version format. This is an optimization
        # in order to avoid accessing it too much.
        try:
            format_version = node._v_file.format_version
        except AttributeError:
            parsed_version = None
        else:
            if format_version == "unknown":
                parsed_version = None
            else:
                parsed_version = tuple(map(int, format_version.split(".")))
        dict_["_v__format_version"] = parsed_version
        # Split the attribute list in system and user lists
        dict_["_v_attrnamessys"] = []
        dict_["_v_attrnamesuser"] = []
        for attr in self._v_attrnames:
            # put the attributes on the local dictionary to allow
            # tab-completion
            self.__getattr__(attr)
            if issysattrname(attr):
                self._v_attrnamessys.append(attr)
            else:
                self._v_attrnamesuser.append(attr)

        # Sort the attributes
        self._v_attrnames.sort()
        self._v_attrnamessys.sort()
        self._v_attrnamesuser.sort()

    def _g_update_node_location(self, node: Node) -> None:
        """Update the location information about the associated `node`."""
        dict_ = self.__dict__
        dict_["_v__nodefile"] = node._v_file
        dict_["_v__nodepath"] = node._v_pathname
        # hdf5extension operations:
        self._g_new(node)

    def _f_list(
        self, attrset: Literal["all", "sys", "user"] = "user"
    ) -> list[str]:
        """Get a list of attribute names.

        The attrset string selects the attribute set to be used.  A
        'user' value returns only user attributes (this is the default).
        A 'sys' value returns only system attributes.  Finally, 'all'
        returns both system and user attributes.

        """
        if attrset == "user":
            return self._v_attrnamesuser[:]
        elif attrset == "sys":
            return self._v_attrnamessys[:]
        elif attrset == "all":
            return self._v_attrnames[:]

    def __dir__(self) -> list[str]:
        """Autocomplete only children named as valid python identifiers.

        Only PY3 supports this special method.
        """
        return list(
            {
                c
                for c in super().__dir__() + self._v_attrnames
                if c.isidentifier()
            }
        )

    def __getattr__(self, name: str) -> Any:
        """Get the attribute named "name"."""
        # If attribute does not exist, raise AttributeError
        if name not in self._v_attrnames:
            raise AttributeError(
                f"Attribute {name!r} does not exist "
                f"in node: {self._v__nodepath!r}"
            )

        # Read the attribute from disk. This is an optimization to read
        # quickly system attributes that are _string_ values, but it
        # takes care of other types as well as for example NROWS for
        # Tables and EXTDIM for EArrays
        format_version = self._v__format_version
        value = self._g_getattr(self._v_node, name)

        # Check whether the value is pickled
        # Pickled values always seems to end with a "."
        maybe_pickled = (
            isinstance(value, np.generic)  # NumPy scalar?
            and value.dtype.type == np.bytes_  # string type?
            and value.itemsize > 0
            and value.endswith(b".")
        )

        if maybe_pickled and value in [b"0", b"0."]:
            # Workaround for a bug in many versions of Python (starting
            # somewhere after Python 2.6.1).  See ticket #253.
            retval = value
        elif (
            maybe_pickled
            and _field_fill_re.match(name)
            and format_version == (1, 5)
        ):
            # This format was used during the first 1.2 releases, just
            # for string defaults.
            try:
                retval = pickle.loads(value)
                retval = np.array(retval)
            except ImportError:
                retval = None  # signal error avoiding exception
        elif (
            maybe_pickled
            and name == "FILTERS"
            and format_version is not None
            and format_version < (2, 0)
        ):
            # This is a big hack, but we don't have other way to recognize
            # pickled filters of PyTables 1.x files.
            value = _old_filters_re.sub(_new_filters_sub, value, 1)
            retval = pickle.loads(value)  # pass unpickling errors through
        elif maybe_pickled:
            try:
                retval = pickle.loads(value)
            # except cPickle.UnpicklingError:
            # It seems that pickle may raise other errors than UnpicklingError
            # Perhaps it would be better just an "except:" clause?
            # except (cPickle.UnpicklingError, ImportError):
            # Definitely (see SF bug #1254636)
            except UnicodeDecodeError:
                # Object maybe pickled on python 2 and unpickled on python 3.
                # encoding='bytes' was added in python 3.4 to resolve this.
                # However 'bytes' mangles class attributes as they are
                # unplicked as bytestrings. Hence try 'latin1' first.
                # Ref: http://bugs.python.org/issue6784
                try:
                    retval = pickle.loads(value, encoding="latin1")
                except TypeError:
                    try:
                        retval = pickle.loads(value, encoding="bytes")
                    except Exception:
                        retval = value
                except Exception:
                    retval = value
            except Exception:
                # catch other unpickling errors:
                # ivb (2005-09-07): It is too hard to tell
                # whether the unpickling failed
                # because of the string not being a pickle one at all,
                # because of a malformed pickle string,
                # or because of some other problem in object reconstruction,
                # thus making inconvenient even the issuing of a warning here.
                # The documentation contains a note on this issue,
                # explaining how the user can tell where the problem was.
                retval = value
            # Additional check for allowing a workaround for #307
            if isinstance(retval, str) and retval == "":
                retval = np.array(retval)[()]
        elif (
            name == "FILTERS"
            and format_version is not None
            and format_version >= (2, 0)
        ):
            try:
                retval = Filters._unpack(value)
            except ValueError:
                warnings.warn(FiltersWarning("Failed parsing FILTERS key"))
                retval = None
        elif name == "TITLE" and not isinstance(value, str):
            retval = value.decode("utf-8")
        elif (
            issysattrname(name)
            and isinstance(value, (bytes, str))
            and not isinstance(value, str)
            and not _field_fill_re.match(name)
        ):
            # system attributes should always be str
            # python 3, bytes and not "FIELD_[0-9]+_FILL"
            retval = value.decode("utf-8")
        else:
            retval = value

        # Put this value in local directory
        self.__dict__[name] = retval
        return retval

    def _g__setattr(self, name: str, value: Any) -> None:
        """Set a PyTables attribute.

        Sets a (maybe new) PyTables attribute with the specified `name`
        and `value`.  If the attribute already exists, it is simply
        replaced.

        It does not log the change.

        """
        # Save this attribute to disk
        # (overwriting an existing one if needed)
        stvalue = value
        if issysattrname(name):
            if name in ["EXTDIM", "AUTO_INDEX", "DIRTY", "NODE_TYPE_VERSION"]:
                stvalue = np.array(value, dtype=np.int32)
                value = stvalue[()]
            elif name == "NROWS":
                stvalue = np.array(value, dtype=SizeType)
                value = stvalue[()]
            elif (
                name == "FILTERS"
                and self._v__format_version is not None
                and self._v__format_version >= (2, 0)
            ):
                stvalue = value._pack()
                # value will remain as a Filters instance here
        # Convert value from a Python scalar into a NumPy scalar
        # (only in case it has not been converted yet)
        # Fixes ticket #59
        if stvalue is value and type(value) in (
            bool,
            bytes,
            int,
            float,
            complex,
            str,
            np.str_,
        ):
            # Additional check for allowing a workaround for #307
            if isinstance(value, str) and len(value) == 0:
                stvalue = np.array("")
            else:
                stvalue = np.array(value)
            value = stvalue[()]

        self._g_setattr(self._v_node, name, stvalue)

        # New attribute or value. Introduce it into the local
        # directory
        self.__dict__[name] = value

        # Finally, add this attribute to the list if not present
        attrnames = self._v_attrnames
        if name not in attrnames:
            attrnames.append(name)
            attrnames.sort()
            if issysattrname(name):
                attrnamessys = self._v_attrnamessys
                attrnamessys.append(name)
                attrnamessys.sort()
            else:
                attrnamesuser = self._v_attrnamesuser
                attrnamesuser.append(name)
                attrnamesuser.sort()

    def __setattr__(self, name: str, value: Any) -> None:
        """Set a PyTables attribute.

        Sets a (maybe new) PyTables attribute with the specified `name`
        and `value`.  If the attribute already exists, it is simply
        replaced.

        A ``ValueError`` is raised when the name starts with a reserved
        prefix or contains a ``/``.  A `NaturalNameWarning` is issued if
        the name is not a valid Python identifier.  A
        `PerformanceWarning` is issued when the recommended maximum
        number of attributes in a node is going to be exceeded.

        """
        nodefile = self._v__nodefile
        attrnames = self._v_attrnames

        # Check for name validity
        check_attribute_name(name)

        nodefile._check_writable()

        # Check if there are too many attributes.
        max_node_attrs = nodefile.params["MAX_NODE_ATTRS"]
        if len(attrnames) >= max_node_attrs:
            warnings.warn(
                """\
node ``%s`` is exceeding the recommended maximum number of attributes (%d);\
be ready to see PyTables asking for *lots* of memory and possibly slow I/O"""
                % (self._v__nodepath, max_node_attrs),
                PerformanceWarning,
            )

        undo_enabled = nodefile.is_undo_enabled()
        # Log old attribute removal (if any).
        if undo_enabled and (name in attrnames):
            self._g_del_and_log(name)

        # Set the attribute.
        self._g__setattr(name, value)

        # Log new attribute addition.
        if undo_enabled:
            self._g_log_add(name)

    def _g_log_add(self, name: str) -> None:
        self._v__nodefile._log("ADDATTR", self._v__nodepath, name)

    def _g_del_and_log(self, name: str) -> None:
        nodefile = self._v__nodefile
        node_pathname = self._v__nodepath
        # Log *before* moving to use the right shadow name.
        nodefile._log("DELATTR", node_pathname, name)
        attr_to_shadow(nodefile, node_pathname, name)

    def _g__delattr(self, name: str) -> None:
        """Delete a PyTables attribute.

        Deletes the specified existing PyTables attribute.

        It does not log the change.

        """
        # Delete the attribute from disk
        self._g_remove(self._v_node, name)

        # Delete the attribute from local lists
        self._v_attrnames.remove(name)
        if name in self._v_attrnamessys:
            self._v_attrnamessys.remove(name)
        else:
            self._v_attrnamesuser.remove(name)

        # Delete the attribute from the local directory
        # closes (#1049285)
        del self.__dict__[name]

    def __delattr__(self, name: str) -> None:
        """Delete a PyTables attribute.

        Deletes the specified existing PyTables attribute from the
        attribute set.  If a nonexistent or system attribute is
        specified, an ``AttributeError`` is raised.

        """
        nodefile = self._v__nodefile

        # Check if attribute exists
        if name not in self._v_attrnames:
            raise AttributeError(
                "Attribute ('%s') does not exist in node '%s'"
                % (name, self._v__nodepath)
            )

        nodefile._check_writable()

        # Remove the PyTables attribute or move it to shadow.
        if nodefile.is_undo_enabled():
            self._g_del_and_log(name)
        else:
            self._g__delattr(name)

    def __getitem__(self, name: str) -> Any:
        """Implement a dictionary like interface for `__getattr__()`."""
        try:
            return self.__getattr__(name)
        except AttributeError:
            # Capture the AttributeError and re-raise a KeyError one
            raise KeyError(
                "Attribute ('%s') does not exist in node '%s'"
                % (name, self._v__nodepath)
            )

    def __setitem__(self, name: str, value: Any) -> None:
        """Implement a dictionary like interface for `__setattr__()`."""
        self.__setattr__(name, value)

    def __delitem__(self, name: str) -> None:
        """Implement a dictionary like interface for `__delattr__()`."""
        try:
            self.__delattr__(name)
        except AttributeError:
            # Capture the AttributeError and re-raise a KeyError one
            raise KeyError(
                "Attribute ('%s') does not exist in node '%s'"
                % (name, self._v__nodepath)
            )

    def __contains__(self, name: str) -> bool:
        """Return True if the set contains an attribute with the specified name.

        A true value is returned if the attribute set has an attribute
        with the given name, false otherwise.

        """
        return name in self._v_attrnames

    def _f_rename(self, oldattrname: str, newattrname: str) -> None:
        """Rename an attribute from oldattrname to newattrname."""
        if oldattrname == newattrname:
            # Do nothing
            return

        # First, fetch the value of the oldattrname
        attrvalue = getattr(self, oldattrname)

        # Now, create the new attribute
        setattr(self, newattrname, attrvalue)

        # Finally, remove the old attribute
        delattr(self, oldattrname)

    def _g_copy(
        self,
        newset: AttributeSet,
        set_attr: Callable[[str, Any], None] | None = None,
        copyclass: bool = False,
    ) -> None:
        """Copy set attributes.

        Copies all user and allowed system PyTables attributes to the
        given attribute set, replacing the existing ones.

        You can specify a *bound* method of the destination set that
        will be used to set its attributes.  Else, its `_g__setattr`
        method will be used.

        Changes are logged depending on the chosen setting method.  The
        default setting method does not log anything.

        .. versionchanged:: 3.0
           The *newSet* parameter has been renamed into *newset*.

        .. versionchanged:: 3.0
           The *copyClass* parameter has been renamed into *copyclass*.

        """
        copysysattrs = newset._v__nodefile.params["PYTABLES_SYS_ATTRS"]
        if set_attr is None:
            set_attr = newset._g__setattr

        for attrname in self._v_attrnamesuser:
            # Do not copy the unimplemented attributes.
            if attrname not in self._v_unimplemented:
                set_attr(attrname, getattr(self, attrname))
        # Copy the system attributes that we are allowed to.
        if copysysattrs:
            for attrname in self._v_attrnamessys:
                if (
                    (attrname not in SYS_ATTRS_NOTTOBECOPIED)
                    # Do not copy the FIELD_ attributes in tables as this can
                    # be really *slow* (don't know exactly the reason).
                    # See #304.
                    and not attrname.startswith("FIELD_")
                ):
                    set_attr(attrname, getattr(self, attrname))
            # Copy CLASS and VERSION attributes if requested
            if copyclass:
                for attrname in FORCE_COPY_CLASS:
                    if attrname in self._v_attrnamessys:
                        set_attr(attrname, getattr(self, attrname))

    def _f_copy(self, where: Node) -> None:
        """Copy attributes to the where node.

        Copies all user and certain system attributes to the given where
        node (a Node instance - see :ref:`NodeClassDescr`), replacing
        the existing ones.

        """
        # AttributeSet must be defined in order to define a Node.
        # However, we need to know Node here.
        # Using class_name_dict avoids a circular import.
        if not isinstance(where, class_name_dict["Node"]):
            raise TypeError(f"destination object is not a node: {where!r}")
        self._g_copy(where._v_attrs, where._v_attrs.__setattr__)

    def _g_close(self) -> None:
        # Nothing will be done here, as the existing instance is completely
        # operative now.
        pass

    def __str__(self) -> str:
        """Return the string representation for the object."""
        # The pathname
        pathname = self._v__nodepath
        # Get this class name
        classname = self.__class__.__name__
        # The attribute names
        attrnumber = sum(1 for _ in self._v_attrnames)
        return f"{pathname}._v_attrs ({classname}), {attrnumber} attributes"

    def __repr__(self) -> str:
        """Detailed string representation for this object."""
        # print additional info only if there are attributes to show
        attrnames = list(self._v_attrnames)
        if attrnames:
            rep = [f"{attr} := {getattr(self, attr)!r}" for attr in attrnames]
            return f"{self!s}:\n   [" + ",\n    ".join(rep) + "]"
        else:
            return str(self)


class NotLoggedAttributeSet(AttributeSet):
    """Attribut set without automatic logging."""

    def _g_log_add(self, name: str) -> None:
        pass

    def _g_del_and_log(self, name: str) -> None:
        self._g__delattr(name)
