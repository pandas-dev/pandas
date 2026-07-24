"""Create PyTables files and the object tree.

This module support importing generic HDF5 files, on top of which
PyTables files are created, read or extended. If a file exists, an
object tree mirroring their hierarchical structure is created in memory.
File class offer methods to traverse the tree, as well as to create new
nodes.

"""

from __future__ import annotations

import os
import atexit
import weakref
import datetime
import warnings
from typing import Any, Literal
from pathlib import Path
from collections import defaultdict
from collections.abc import Callable, Generator, Iterator

import numpy as np
import numexpr as ne
import numpy.typing as npt

from . import (
    hdf5extension,
    utilsextension,
    parameters,
    undoredo,
    linkextension,
    lrucacheextension,
)
from .atom import Atom
from .leaf import Leaf
from .link import SoftLink, ExternalLink
from .node import Node, NotLoggedMixin
from .path import join_path, split_path
from .array import Array
from .group import Group, RootGroup, TransactionGroupG, TransactionG, MarkG
from .table import Table
from .utils import detect_number_of_cores
from .carray import CArray
from .earray import EArray
from .flavor import flavor_of, array_as_internal
from .filters import Filters
from .vlarray import VLArray
from .registry import get_class_by_name
from .exceptions import (
    ClosedFileError,
    FileModeError,
    NodeError,
    NoSuchNodeError,
    UnclosedFileWarning,
    UndoRedoError,
    ClosedNodeError,
    PerformanceWarning,
)
from .description import (
    Description,
    IsDescription,
    UInt8Col,
    StringCol,
    descr_from_dtype,
    dtype_from_descr,
)

# format_version = "1.0"  # Initial format
# format_version = "1.1"  # Changes in ucl compression
# format_version = "1.2"  # Support for enlargeable arrays and VLA's
#                         # 1.2 was introduced in PyTables 0.8
# format_version = "1.3"  # Support for indexes in Tables
#                         # 1.3 was introduced in PyTables 0.9
# format_version = "1.4"  # Support for multidimensional attributes
#                         # 1.4 was introduced in PyTables 1.1
# format_version = "1.5"  # Support for persistent defaults in tables
#                         # 1.5 was introduced in PyTables 1.2
# format_version = "1.6"  # Support for NumPy objects and new flavors for
#                         # objects.
#                         # 1.6 was introduced in pytables 1.3
# format_version = "2.0"  # Pickles are not used anymore in system attrs
#                         # 2.0 was introduced in PyTables 2.0
format_version = "2.1"  # Numeric and numarray flavors are gone.

compatible_formats = []  # Old format versions we can read
#                        # Empty means that we support all the old formats


class _FileRegistry:
    def __init__(self) -> None:
        self._name_mapping: dict[str, set[File]] = defaultdict(set)
        self._handlers: set[File] = set()

    @property
    def filenames(self) -> list[str]:
        return list(self._name_mapping)

    @property
    def handlers(self) -> set[File]:
        # return set(self._handlers)  # return a copy
        return self._handlers

    def __len__(self) -> int:
        return len(self._handlers)

    def __contains__(self, filename: str) -> bool:
        return filename in self.filenames

    def add(self, handler: File) -> None:
        self._name_mapping[handler.filename].add(handler)
        self._handlers.add(handler)

    def remove(self, handler: File) -> None:
        filename = handler.filename
        self._name_mapping[filename].remove(handler)
        # remove enpty keys
        if not self._name_mapping[filename]:
            del self._name_mapping[filename]
        self._handlers.remove(handler)

    def get_handlers_by_name(self, filename: str) -> set[File]:
        # return set(self._name_mapping[filename])  # return a copy
        return self._name_mapping[filename]

    def close_all(self) -> None:
        handlers = list(self._handlers)  # make a copy
        for fileh in handlers:
            msg = f"Closing remaining open file: {fileh.filename}"
            warnings.warn(UnclosedFileWarning(msg))
            fileh.close()


# Dict of opened files (keys are filenames and values filehandlers)
_open_files = _FileRegistry()

# Opcodes for do-undo actions
_op_to_code = {
    "MARK": 0,
    "CREATE": 1,
    "REMOVE": 2,
    "MOVE": 3,
    "ADDATTR": 4,
    "DELATTR": 5,
}

_code_to_op = ["MARK", "CREATE", "REMOVE", "MOVE", "ADDATTR", "DELATTR"]


# Paths and names for hidden nodes related with transactions.
_trans_version = "1.0"

_trans_group_parent = "/"
_trans_group_name = "_p_transactions"
_trans_group_path = join_path(_trans_group_parent, _trans_group_name)

_action_log_parent = _trans_group_path
_action_log_name = "actionlog"
_action_log_path = join_path(_action_log_parent, _action_log_name)

_trans_parent = _trans_group_path
_trans_name = "t%d"  # %d -> transaction number
_trans_path = join_path(_trans_parent, _trans_name)

_mark_parent = _trans_path
_mark_name = "m%d"  # %d -> mark number
_mark_path = join_path(_mark_parent, _mark_name)

_shadow_parent = _mark_path
_shadow_name = "a%d"  # %d -> action number
_shadow_path = join_path(_shadow_parent, _shadow_name)


def _checkfilters(filters: Filters) -> None:
    if not (filters is None or isinstance(filters, Filters)):
        raise TypeError(
            f"filter parameter has to be None or a Filter "
            f"instance and the passed type is: '{type(filters)}'"
        )


def copy_file(
    srcfilename: str, dstfilename: str, overwrite: bool = False, **kwargs
) -> None:
    """Copy the content of a PyTables file to another.

    This function allows you to copy an existing PyTables file named
    srcfilename to another file called dstfilename. The source file
    must exist and be readable. The destination file can be
    overwritten in place if existing by asserting the overwrite
    argument.

    This function is a shorthand for the :meth:`File.copy_file` method,
    which acts on an already opened file. kwargs takes keyword
    arguments used to customize the copying process. See the
    documentation of :meth:`File.copy_file` for a description of those
    arguments.

    """
    # Open the source file.
    srcfileh = open_file(srcfilename, mode="r")

    try:
        # Copy it to the destination file.
        srcfileh.copy_file(dstfilename, overwrite=overwrite, **kwargs)
    finally:
        # Close the source file.
        srcfileh.close()


hdf5_version_str = utilsextension.get_hdf5_version()
hdf5_version_tup = tuple(map(int, hdf5_version_str.split("-")[0].split(".")))
_FILE_OPEN_POLICY = "strict" if hdf5_version_tup < (1, 8, 7) else "default"


def open_file(
    filename: str,
    mode: Literal["r", "w", "a", "r+"] = "r",
    title: str = "",
    root_uep: str = "/",
    filters: Filters | None = None,
    **kwargs,
) -> File:
    """Open a PyTables (or generic HDF5) file and return a File object.

    Parameters
    ----------
    filename : str
        The name of the file (supports environment variable expansion).
        It is suggested that file names have any of the .h5, .hdf or
        .hdf5 extensions, although this is not mandatory.
    mode : str
        The mode to open the file. It can be one of the
        following:

            * *'r'*: Read-only; no data can be modified.
            * *'w'*: Write; a new file is created (an existing file
              with the same name would be deleted).
            * *'a'*: Append; an existing file is opened for reading and
              writing, and if the file does not exist it is created.
            * *'r+'*: It is similar to 'a', but the file must already
              exist.

    title : str
        If the file is to be created, a TITLE string attribute will be
        set on the root group with the given value. Otherwise, the
        title will be read from disk, and this will not have any effect.
    root_uep : str
        The root User Entry Point. This is a group in the HDF5 hierarchy
        which will be taken as the starting point to create the object
        tree. It can be whatever existing group in the file, named by
        its HDF5 path. If it does not exist, an HDF5ExtError is issued.
        Use this if you do not want to build the *entire* object tree,
        but rather only a *subtree* of it.

        .. versionchanged:: 3.0
           The *rootUEP* parameter has been renamed into *root_uep*.

    filters : Filters
        An instance of the Filters (see :ref:`FiltersClassDescr`) class
        that provides information about the desired I/O filters
        applicable to the leaves that hang directly from the *root group*,
        unless other filter properties are specified for these leaves.
        Besides, if you do not specify filter properties for child groups,
        they will inherit these ones, which will in turn propagate to
        child nodes.

    Notes
    -----
    In addition, it recognizes the (lowercase) names of parameters
    present in :file:`tables/parameters.py` as additional keyword
    arguments.
    See :ref:`parameter_files` for a detailed info on the supported
    parameters.

    .. note::

        If you need to deal with a large number of nodes in an
        efficient way, please see :ref:`LRUOptim` for more info and
        advices about the integrated node cache engine.

    """
    filename = os.fspath(filename)
    # XXX filename normalization ??

    # Check already opened files
    if _FILE_OPEN_POLICY == "strict":
        # This policy does not allow to open the same file multiple times
        # even in read-only mode
        if filename in _open_files:
            raise ValueError(
                "The file '%s' is already opened.  "
                "Please close it before reopening.  "
                "HDF5 v.%s, FILE_OPEN_POLICY = '%s'"
                % (
                    filename,
                    utilsextension.get_hdf5_version(),
                    _FILE_OPEN_POLICY,
                )
            )
    else:
        for filehandle in _open_files.get_handlers_by_name(filename):
            omode = filehandle.mode
            # 'r' is incompatible with everything except 'r' itself
            if mode == "r" and omode != "r":
                raise ValueError(
                    "The file '%s' is already opened, but "
                    "not in read-only mode (as requested)." % filename
                )
            # 'a' and 'r+' are compatible with everything except 'r'
            elif mode in ("a", "r+") and omode == "r":
                raise ValueError(
                    "The file '%s' is already opened, but "
                    "in read-only mode.  Please close it before "
                    "reopening in append mode." % filename
                )
            # 'w' means that we want to destroy existing contents
            elif mode == "w":
                raise ValueError(
                    "The file '%s' is already opened.  Please "
                    "close it before reopening in write mode." % filename
                )

    # Finally, create the File instance, and return it
    return File(filename, mode, title, root_uep, filters, **kwargs)


# A dumb class that doesn't keep anything at all
class _NoCache:
    def __len__(self) -> int:
        return 0

    def __contains__(self, key: Any) -> bool:
        return False

    def __iter__(self) -> Iterator:
        return iter([])

    def __setitem__(self, key: Any, value: Any) -> None:
        pass

    __marker = object()

    def pop(self, key, d=__marker):
        if d is not self.__marker:
            return d
        raise KeyError(key)


class _DictCache(dict):
    def __init__(self, nslots: int) -> None:
        if nslots < 1:
            raise ValueError("Invalid number of slots: %d" % nslots)
        self.nslots = nslots
        super().__init__()

    def __setitem__(self, key: Any, value: Any) -> None:
        # Check if we are running out of space
        if len(self) > self.nslots:
            warnings.warn(
                "the dictionary of node cache is exceeding the recommended "
                "maximum number (%d); be ready to see PyTables asking for "
                "*lots* of memory and possibly slow I/O." % (self.nslots),
                PerformanceWarning,
            )
        super().__setitem__(key, value)


class NodeManager:
    """Node manager."""

    def __init__(self, nslots: int = 64, node_factory=None) -> None:
        super().__init__()

        self.registry = weakref.WeakValueDictionary()

        if nslots > 0:
            cache = lrucacheextension.NodeCache(nslots)
        elif nslots == 0:
            cache = _NoCache()
        else:
            # nslots < 0
            cache = _DictCache(-nslots)

        self.cache = cache

        # node_factory(node_path)
        self.node_factory = node_factory

    def register_node(self, node: Node, key: str | None) -> None:
        """Register a node."""
        if key is None:
            key = node._v_pathname

        if key in self.registry:
            if not self.registry[key]._v_isopen:
                del self.registry[key]
                self.registry[key] = node
            elif self.registry[key] is not node:
                raise RuntimeError(
                    "trying to register a node with an "
                    "existing key: ``%s``" % key
                )
        else:
            self.registry[key] = node

    def cache_node(self, node: Node, key: str | None = None) -> None:
        """Create a node."""
        if key is None:
            key = node._v_pathname

        self.register_node(node, key)
        if key in self.cache:
            oldnode = self.cache.pop(key)
            if oldnode is not node and oldnode._v_isopen:
                raise RuntimeError(
                    "trying to cache a node with an "
                    "existing key: ``%s``" % key
                )

        self.cache[key] = node

    def get_node(self, key: str) -> Node:
        """Return a node matching the input key."""
        node = self.cache.pop(key, None)
        if node is not None:
            if node._v_isopen:
                self.cache_node(node, key)
                return node
            else:
                # this should not happen
                warnings.warn("a closed node found in the cache: ``%s``" % key)

        if key in self.registry:
            node = self.registry[key]
            if node is None:
                # this should not happen since WeakValueDictionary drops all
                # dead weakrefs
                warnings.warn(
                    "None is stored in the registry for key: " "``%s``" % key
                )
            elif node._v_isopen:
                self.cache_node(node, key)
                return node
            else:
                # this should not happen
                warnings.warn(
                    "a closed node found in the registry: " "``%s``" % key
                )
                del self.registry[key]
                node = None

        if self.node_factory:
            node = self.node_factory(key)
            self.cache_node(node, key)

        return node

    def rename_node(self, oldkey: str, newkey: str) -> None:
        """Rename a node."""
        for cache in (self.cache, self.registry):
            if oldkey in cache:
                node = cache.pop(oldkey)
                cache[newkey] = node

    def drop_from_cache(self, nodepath: str) -> None:
        """Remove the node from cache."""
        # Remove the node from the cache.
        self.cache.pop(nodepath, None)

    def drop_node(self, node: Node, check_unregistered: bool = True) -> None:
        """Drop the `node`.

        Remove the node from the cache and, if it has no more references,
        close it.

        """
        # Remove all references to the node.
        nodepath = node._v_pathname

        self.drop_from_cache(nodepath)

        if nodepath in self.registry:
            if not node._v_isopen:
                del self.registry[nodepath]
        elif check_unregistered:
            # If the node is not in the registry (this should never happen)
            # we close it forcibly since it is not ensured that the __del__
            # method is called for object that are still alive when the
            # interpreter is shut down
            if node._v_isopen:
                warnings.warn(
                    "dropping a node that is not in the registry: "
                    "``%s``" % nodepath
                )

                node._g_pre_kill_hook()
                node._f_close()

    def flush_nodes(self) -> None:
        """Flush all nodes."""
        # Only iter on the nodes in the registry since nodes in the cache
        # should always have an entry in the registry
        closed_keys = []
        for path, node in list(self.registry.items()):
            if not node._v_isopen:
                closed_keys.append(path)
            elif "/_i_" not in path:  # Indexes are not necessary to be flushed
                if isinstance(node, Leaf):
                    node.flush()

        for path in closed_keys:
            # self.cache.pop(path, None)
            if path in self.cache:
                warnings.warn("closed node the cache: ``%s``" % path)
                self.cache.pop(path, None)
            self.registry.pop(path)

    @staticmethod
    def _close_nodes(
        nodepaths: list[str], get_node: Callable[[str], Group | Node]
    ) -> None:
        for nodepath in nodepaths:
            try:
                node = get_node(nodepath)
            except KeyError:
                pass
            else:
                if not node._v_isopen or node._v__deleting:
                    continue

                try:
                    # Avoid descendent nodes to also iterate over
                    # their descendents, which are already to be
                    # closed by this loop.
                    if hasattr(node, "_f_get_child"):
                        node._g_close()
                    else:
                        node._f_close()
                    del node
                except ClosedNodeError:
                    # import traceback
                    # type_, value, tb = sys.exc_info()
                    # exception_dump = ''.join(
                    #     traceback.format_exception(type_, value, tb))
                    # warnings.warn(
                    #     "A '%s' exception occurred trying to close a node "
                    #     "that was supposed to be open.\n"
                    #     "%s" % (type_.__name__, exception_dump))
                    pass

    def close_subtree(self, prefix: str = "/") -> None:
        """Close a sub-tree of nodes."""
        if not prefix.endswith("/"):
            prefix = prefix + "/"

        cache = self.cache
        registry = self.registry

        # Ensure tables are closed before their indices
        paths = [
            path
            for path in cache
            if path.startswith(prefix) and "/_i_" not in path
        ]
        self._close_nodes(paths, cache.pop)

        # Close everything else (i.e. indices)
        paths = [path for path in cache if path.startswith(prefix)]
        self._close_nodes(paths, cache.pop)

        # Ensure tables are closed before their indices
        paths = [
            path
            for path in registry
            if path.startswith(prefix) and "/_i_" not in path
        ]
        self._close_nodes(paths, registry.pop)

        # Close everything else (i.e. indices)
        paths = [path for path in registry if path.startswith(prefix)]
        self._close_nodes(paths, registry.pop)

    def shutdown(self) -> None:
        """Shutdown the node manager."""
        registry = self.registry
        cache = self.cache

        # self.close_subtree('/')

        keys = list(cache)  # copy
        for key in keys:
            node = cache.pop(key)
            if node._v_isopen:
                registry.pop(node._v_pathname, None)
                node._f_close()

        while registry:
            key, node = registry.popitem()
            if node._v_isopen:
                node._f_close()


class File(hdf5extension.File):
    """The in-memory representation of a PyTables file.

    An instance of this class is returned when a PyTables file is
    opened with the :func:`tables.open_file` function. It offers methods
    to manipulate (create, rename, delete...) nodes and handle their
    attributes, as well as methods to traverse the object tree.
    The *user entry point* to the object tree attached to the HDF5 file
    is represented in the root_uep attribute.
    Other attributes are available.

    File objects support an *Undo/Redo mechanism* which can be enabled
    with the :meth:`File.enable_undo` method. Once the Undo/Redo
    mechanism is enabled, explicit *marks* (with an optional unique
    name) can be set on the state of the database using the
    :meth:`File.mark`
    method. There are two implicit marks which are always available:
    the initial mark (0) and the final mark (-1).  Both the identifier
    of a mark and its name can be used in *undo* and *redo* operations.

    Hierarchy manipulation operations (node creation, movement and
    removal) and attribute handling operations (setting and deleting)
    made after a mark can be undone by using the :meth:`File.undo`
    method, which returns the database to the state of a past mark.
    If undo() is not followed by operations that modify the hierarchy
    or attributes, the :meth:`File.redo` method can be used to return
    the database to the state of a future mark. Else, future states of
    the database are forgotten.

    Note that data handling operations can not be undone nor redone by
    now. Also, hierarchy manipulation operations on nodes that do not
    support the Undo/Redo mechanism issue an UndoRedoWarning *before*
    changing the database.

    The Undo/Redo mechanism is persistent between sessions and can
    only be disabled by calling the :meth:`File.disable_undo` method.

    File objects can also act as context managers when using the with
    statement introduced in Python 2.5.  When exiting a context, the
    file is automatically closed.

    Parameters
    ----------
    filename : str
        The name of the file (supports environment variable expansion).
        It is suggested that file names have any of the .h5, .hdf or
        .hdf5 extensions, although this is not mandatory.

    mode : str
        The mode to open the file. It can be one of the
        following:

            * *'r'*: Read-only; no data can be modified.
            * *'w'*: Write; a new file is created (an existing file
              with the same name would be deleted).
            * *'a'*: Append; an existing file is opened for reading
              and writing, and if the file does not exist it is created.
            * *'r+'*: It is similar to 'a', but the file must already
              exist.

    title : str
        If the file is to be created, a TITLE string attribute will be
        set on the root group with the given value. Otherwise, the
        title will be read from disk, and this will not have any effect.

    root_uep : str
        The root User Entry Point. This is a group in the HDF5 hierarchy
        which will be taken as the starting point to create the object
        tree. It can be whatever existing group in the file, named by
        its HDF5 path. If it does not exist, an HDF5ExtError is issued.
        Use this if you do not want to build the *entire* object tree,
        but rather only a *subtree* of it.

        .. versionchanged:: 3.0
           The *rootUEP* parameter has been renamed into *root_uep*.

    filters : Filters
        An instance of the Filters (see :ref:`FiltersClassDescr`) class that
        provides information about the desired I/O filters applicable to the
        leaves that hang directly from the *root group*, unless other filter
        properties are specified for these leaves. Besides, if you do not
        specify filter properties for child groups, they will inherit these
        ones, which will in turn propagate to child nodes.

    Notes
    -----
    In addition, it recognizes the (lowercase) names of parameters
    present in :file:`tables/parameters.py` as additional keyword
    arguments.
    See :ref:`parameter_files` for a detailed info on the supported
    parameters.


    .. rubric:: File attributes

    .. attribute:: filename

        The name of the opened file.

    .. attribute:: format_version

        The PyTables version number of this file.

    .. attribute:: isopen

        True if the underlying file is open, false otherwise.

    .. attribute:: mode

        The mode in which the file was opened.

    .. attribute:: root

        The *root* of the object tree hierarchy (a Group instance).

    .. attribute:: root_uep

        The UEP (user entry point) group name in the file (see
        the :func:`open_file` function).

        .. versionchanged:: 3.0
           The *rootUEP* attribute has been renamed into *root_uep*.

    """

    # The top level kinds. Group must go first!
    _node_kinds = ("Group", "Leaf", "Link", "Unknown")

    @property
    def title(self) -> str:
        """Title of the root group in the file."""
        return self.root._v_title

    @title.setter
    def title(self, title: str) -> None:
        self.root._v_title = title

    @title.deleter
    def title(self) -> None:
        del self.root._v_title

    @property
    def filters(self) -> Filters:
        """Filter properties for the root group.

        See :ref:`FiltersClassDescr`).
        """
        return self.root._v_filters

    @filters.setter
    def filters(self, filters: Filters) -> None:
        self.root._v_filters = filters

    @filters.deleter
    def filters(self) -> None:
        del self.root._v_filters

    def __init__(
        self,
        filename: str,
        mode: Literal["r", "w", "a", "r+"] = "r",
        title: str = "",
        root_uep: str = "/",
        filters: Filters | None = None,
        **kwargs,
    ) -> None:

        self.filename = os.fspath(filename)
        """The name of the opened file."""

        self.mode = mode
        """The mode in which the file was opened."""

        if mode not in ("r", "r+", "a", "w"):
            raise ValueError(
                "invalid mode string ``%s``. Allowed modes are: "
                "'r', 'r+', 'a' and 'w'" % mode
            )

        # Get all the parameters in parameter file(s)
        params = {
            k: v
            for k, v in parameters.__dict__.items()
            if k.isupper() and not k.startswith("_")
        }
        # Update them with possible keyword arguments
        if [k for k in kwargs if k.isupper()]:
            warnings.warn(
                "The use of uppercase keyword parameters is " "deprecated",
                DeprecationWarning,
            )

        kwargs = {k.upper(): v for k, v in kwargs.items()}
        params.update(kwargs)

        # If MAX_ * _THREADS is not set yet, set it to the number of cores
        # on this machine.

        if params["MAX_NUMEXPR_THREADS"] is None:
            params["MAX_NUMEXPR_THREADS"] = detect_number_of_cores()

        if params["MAX_BLOSC_THREADS"] is None:
            params["MAX_BLOSC_THREADS"] = detect_number_of_cores()

        self.params = params

        # Now, it is time to initialize the File extension
        self._g_new(filename, mode, **params)

        # Check filters and set PyTables format version for new files.
        new = self._v_new
        if new:
            _checkfilters(filters)
            self.format_version = format_version
            """The PyTables version number of this file."""

        # The node manager must be initialized before the root group
        # initialization but the node_factory attribute is set onl later
        # because it is a bound method of the root grop itself.
        node_cache_slots = params["NODE_CACHE_SLOTS"]
        self._node_manager = NodeManager(nslots=node_cache_slots)

        # For the moment Undo/Redo is not enabled.
        self._undoEnabled = False

        # Set the flag to indicate that the file has been opened.
        # It must be set before opening the root group
        # to allow some basic access to its attributes.
        self.isopen = 1
        """True if the underlying file os open, False otherwise."""

        # Append the name of the file to the global dict of files opened.
        _open_files.add(self)

        # Set the number of times this file has been opened to 1
        self._open_count = 1

        # Get the root group from this file
        self.root = root = self.__get_root_group(root_uep, title, filters)
        """The *root* of the object tree hierarchy (a Group instance)."""
        # Complete the creation of the root node
        # (see the explanation in ``RootGroup.__init__()``).
        root._g_post_init_hook()
        self._node_manager.node_factory = self.root._g_load_child

        # Save the PyTables format version for this file.
        if new:
            if params["PYTABLES_SYS_ATTRS"]:
                root._v_attrs._g__setattr(
                    "PYTABLES_FORMAT_VERSION", format_version
                )

        # If the file is old, and not opened in "read-only" mode,
        # check if it has a transaction log
        if not new and self.mode != "r" and _trans_group_path in self:
            # It does. Enable the undo.
            self.enable_undo()

        # Set the maximum number of threads for Numexpr
        ne.set_vml_num_threads(params["MAX_NUMEXPR_THREADS"])

    def __get_root_group(
        self, root_uep: str | None, title: str, filters: Filters
    ) -> RootGroup:
        """Return a Group instance which can act as the root group in a tree.

        If file is opened in "r", "r+" or "a" mode, and the file already
        exists, this method dynamically builds a python object tree
        emulating the structure present on file.

        """
        self._v_objectid = self._get_file_id()

        if root_uep in [None, ""]:
            root_uep = "/"
        # Save the User Entry Point in a variable class
        self.root_uep = root_uep

        new = self._v_new

        # Get format version *before* getting the object tree
        if not new:
            # Firstly, get the PyTables format version for this file
            self.format_version = utilsextension.read_f_attr(
                self._v_objectid, "PYTABLES_FORMAT_VERSION"
            )
            if not self.format_version:
                # PYTABLES_FORMAT_VERSION attribute is not present
                self.format_version = "unknown"
                self._isPTFile = False
            elif not isinstance(self.format_version, str):
                # system attributes should always be str
                self.format_version = self.format_version.decode("utf-8")

        # Create new attributes for the root Group instance and
        # create the object tree
        return RootGroup(self, root_uep, title=title, new=new, filters=filters)

    def _get_or_create_path(
        self, path: Node | str, create: bool
    ) -> Group | Node | RootGroup:
        """Get the given `path` or create it if `create` is true.

        If `create` is true, `path` *must* be a string path and not a
        node, otherwise a `TypeError`will be raised.

        """
        if create:
            path = path._v_pathname if hasattr(path, "_v_pathname") else path
            return self._create_path(path)
        else:
            return self.get_node(path)

    def _create_path(self, path: str) -> Group:
        """Create the groups needed for the `path` to exist.

        The group associated with the given `path` is returned.

        """
        if not hasattr(path, "split"):
            raise TypeError("when creating parents, parent must be a path")

        if path == "/":
            return self.root

        parent, create_group = self.root, self.create_group
        for pcomp in path.split("/")[1:]:
            try:
                child = parent._f_get_child(pcomp)
            except NoSuchNodeError:
                child = create_group(parent, pcomp)
            parent = child
        return parent

    def create_group(
        self,
        where: Group | str,
        name: str,
        title: str = "",
        filters: Filters | None = None,
        createparents: bool = False,
    ) -> Group:
        """Create a new group.

        Parameters
        ----------
        where : str or Group
            The parent group from which the new group will hang. It can be a
            path string (for example '/level1/leaf5'), or a Group instance
            (see :ref:`GroupClassDescr`).
        name : str
            The name of the new group.
        title : str, optional
            A description for this node (it sets the TITLE HDF5 attribute on
            disk).
        filters : Filters
            An instance of the Filters class (see :ref:`FiltersClassDescr`)
            that provides information about the desired I/O filters applicable
            to the leaves that hang directly from this new group (unless other
            filter properties are specified for these leaves). Besides, if you
            do not specify filter properties for its child groups, they will
            inherit these ones.
        createparents : bool
            Whether to create the needed groups for the parent
            path to exist (not done by default).

        See Also
        --------
        Group : for more information on groups

        """
        parentnode = self._get_or_create_path(where, createparents)
        _checkfilters(filters)
        return Group(parentnode, name, title=title, new=True, filters=filters)

    def create_table(
        self,
        where: Group | str,
        name: str,
        description: (
            dict | type[IsDescription] | Description | npt.DTypeLike | None
        ) = None,
        title: str = "",
        filters: Filters | None = None,
        expectedrows: int = 10_000,
        chunkshape: int | tuple[int] | None = None,
        byteorder: str | None = None,
        createparents: bool = False,
        obj: np.ndarray | None = None,
        track_times: bool = True,
    ) -> Table:
        """Create a new table with the given name in where location.

        Parameters
        ----------
        where : str or Group
            The parent group from which the new table will hang. It can be a
            path string (for example '/level1/leaf5'), or a Group instance
            (see :ref:`GroupClassDescr`).
        name : str
            The name of the new table.
        description : Description
            This is an object that describes the table, i.e. how
            many columns it has, their names, types, shapes, etc.  It
            can be any of the following:

                * *A user-defined class*: This should inherit from the
                  IsDescription class (see :ref:`IsDescriptionClassDescr`)
                  where table fields are specified.
                * *A dictionary*: For example, when you do not know
                  beforehand which structure your table will have).
                * *A Description instance*: You can use the description
                  attribute of another table to create a new one with the
                  same structure.
                * *A NumPy dtype*: A completely general structured NumPy
                  dtype.
                * *A NumPy (structured) array instance*: The dtype of
                  this structured array will be used as the description.
                  Also, in case the array has actual data, it will be
                  injected into the newly created table.

            .. versionchanged:: 3.0
               The *description* parameter can be None (default) if *obj* is
               provided.  In that case the structure of the table is deduced
               by *obj*.

        title : str
            A description for this node (it sets the TITLE HDF5 attribute
            on disk).
        filters : Filters
            An instance of the Filters class (see :ref:`FiltersClassDescr`)
            that provides information about the desired I/O filters to be
            applied during the life of this object.
        expectedrows : int
            A user estimate of the number of records that will be in the table.
            If not provided, the default value is EXPECTED_ROWS_TABLE (see
            :file:`tables/parameters.py`). If you plan to create a bigger
            table try providing a guess; this will optimize the HDF5 B-Tree
            creation and management process time and memory used.
        chunkshape
            The shape of the data chunk to be read or written in a
            single HDF5 I/O operation. Filters are applied to those
            chunks of data. The rank of the chunkshape for tables must
            be 1. If None, a sensible value is calculated based on the
            expectedrows parameter (which is recommended).
        byteorder : str
            The byteorder of data *on disk*, specified as 'little' or 'big'.
            If this is not specified, the byteorder is that of the platform,
            unless you passed an array as the description, in which case
            its byteorder will be used.
        createparents : bool
            Whether to create the needed groups for the parent path to exist
            (not done by default).
        obj : python object
            The recarray to be saved.  Accepted types are NumPy record
            arrays.

            The *obj* parameter is optional and it can be provided in
            alternative to the *description* parameter.
            If both *obj* and *description* are provided they must
            be consistent with each other.

            .. versionadded:: 3.0

        track_times
            Whether time data associated with the leaf are recorded (object
            access time, raw data modification time, metadata change time,
            object birth time); default True.  Semantics of these times
            depend on their implementation in the HDF5 library: refer to
            documentation of the H5O_info_t data structure.  As of HDF5
            1.8.15, only ctime (metadata change time) is implemented.

            .. versionadded:: 3.4.3

        See Also
        --------
        Table : for more information on tables

        """
        if obj is not None:
            if not isinstance(obj, np.ndarray):
                raise TypeError("invalid obj parameter %r" % obj)

            descr, _ = descr_from_dtype(obj.dtype, ptparams=self.params)
            if (
                description is not None
                and dtype_from_descr(description, ptparams=self.params)
                != obj.dtype
            ):
                raise TypeError(
                    "the desctiption parameter is not consistent "
                    "with the data type of the obj parameter"
                )
            elif description is None:
                description = descr

        parentnode = self._get_or_create_path(where, createparents)
        if description is None:
            raise ValueError("invalid table description: None")
        _checkfilters(filters)

        ptobj = Table(
            parentnode,
            name,
            description=description,
            title=title,
            filters=filters,
            expectedrows=expectedrows,
            chunkshape=chunkshape,
            byteorder=byteorder,
            track_times=track_times,
        )

        if obj is not None:
            ptobj.append(obj)

        return ptobj

    def create_array(
        self,
        where: Group | str,
        name: str,
        obj: np.ndarray | None = None,
        title: str = "",
        byteorder: str | None = None,
        createparents: bool = False,
        atom: Atom | None = None,
        shape: tuple[int, ...] | None = None,
        track_times: bool = True,
    ) -> Array:
        """Create a new array.

        Parameters
        ----------
        where : str or Group
            The parent group from which the new array will hang. It can be a
            path string (for example '/level1/leaf5'), or a Group instance
            (see :ref:`GroupClassDescr`).
        name : str
            The name of the new array
        obj : python object
            The array or scalar to be saved.  Accepted types are NumPy
            arrays and scalars, as well as native Python sequences and
            scalars, provided that values are regular (i.e. they are
            not like ``[[1,2],2]``) and homogeneous (i.e. all the
            elements are of the same type).

            Also, objects that have some of their dimensions equal to 0
            are not supported (use an EArray node (see
            :ref:`EArrayClassDescr`) if you want to store an array with
            one of its dimensions equal to 0).

            .. versionchanged:: 3.0
               The *Object parameter has been renamed into *obj*.*

        title : str
            A description for this node (it sets the TITLE HDF5 attribute on
            disk).
        byteorder : str
            The byteorder of the data *on disk*, specified as 'little' or
            'big'.  If this is not specified, the byteorder is that of the
            given object.
        createparents : bool, optional
            Whether to create the needed groups for the parent path to exist
            (not done by default).
        atom : Atom
            An Atom (see :ref:`AtomClassDescr`) instance representing
            the *type* and *shape* of the atomic objects to be saved.

            .. versionadded:: 3.0

        shape : tuple of ints
            The shape of the stored array.

            .. versionadded:: 3.0

        track_times
            Whether time data associated with the leaf are recorded (object
            access time, raw data modification time, metadata change time,
            object birth time); default True.  Semantics of these times
            depend on their implementation in the HDF5 library: refer to
            documentation of the H5O_info_t data structure.  As of HDF5
            1.8.15, only ctime (metadata change time) is implemented.

            .. versionadded:: 3.4.3

        See Also
        --------
        Array : for more information on arrays
        create_table : for more information on the rest of parameters

        """
        if obj is None:
            if atom is None or shape is None:
                raise TypeError(
                    "if the obj parameter is not specified "
                    "(or None) then both the atom and shape "
                    "parametes should be provided."
                )
            else:
                # Making strides=(0,...) below is a trick to create the
                # array fast and without memory consumption
                dflt = np.zeros((), dtype=atom.dtype)
                obj = np.ndarray(
                    shape,
                    dtype=atom.dtype,
                    buffer=dflt,
                    strides=(0,) * len(shape),
                )
        else:
            flavor = flavor_of(obj)
            # use a temporary object because converting obj at this stage
            # breaks some test. This is solution performs a double,
            # potentially expensive, conversion of the obj parameter.
            _obj = array_as_internal(obj, flavor)

            if shape is not None and shape != _obj.shape:
                raise TypeError("the shape parameter do not match obj.shape")

            if atom is not None and atom.dtype != _obj.dtype:
                raise TypeError(
                    "the atom parameter is not consistent with "
                    "the data type of the obj parameter"
                )

        parentnode = self._get_or_create_path(where, createparents)
        return Array(
            parentnode,
            name,
            obj=obj,
            title=title,
            byteorder=byteorder,
            track_times=track_times,
        )

    def create_carray(
        self,
        where: Group | str,
        name: str,
        atom: Atom | None = None,
        shape: tuple[int, ...] | None = None,
        title: str = "",
        filters: Filters | None = None,
        chunkshape: int | tuple[int, ...] | None = None,
        byteorder: str | None = None,
        createparents: bool = False,
        obj: np.ndarray | None = None,
        track_times: bool = True,
    ) -> CArray:
        """Create a new chunked array.

        Parameters
        ----------
        where : str or Group
            The parent group from which the new array will hang. It can
            be a path string (for example '/level1/leaf5'), or a Group
            instance (see :ref:`GroupClassDescr`).
        name : str
            The name of the new array
        atom : Atom
            An Atom (see :ref:`AtomClassDescr`) instance representing
            the *type* and *shape* of the atomic objects to be saved.

            .. versionchanged:: 3.0
               The *atom* parameter can be None (default) if *obj* is
               provided.

        shape : tuple
            The shape of the new array.

            .. versionchanged:: 3.0
               The *shape* parameter can be None (default) if *obj* is
               provided.

        title : str, optional
            A description for this node (it sets the TITLE HDF5 attribute
            on disk).
        filters : Filters, optional
            An instance of the Filters class (see :ref:`FiltersClassDescr`)
            that provides information about the desired I/O filters to
            be applied during the life of this object.
        chunkshape : tuple or number or None, optional
            The shape of the data chunk to be read or written in a
            single HDF5 I/O operation.  Filters are applied to those
            chunks of data.  The dimensionality of chunkshape must be
            the same as that of shape.  If None, a sensible value is
            calculated (which is recommended).
        byteorder : str, optional
            The byteorder of the data *on disk*, specified as 'little'
            or 'big'.  If this is not specified, the byteorder is that
            of the given object.
        createparents : bool, optional
            Whether to create the needed groups for the parent path to
            exist (not done by default).
        obj : python object
            The array or scalar to be saved.  Accepted types are NumPy
            arrays and scalars, as well as native Python sequences and
            scalars, provided that values are regular (i.e. they are
            not like ``[[1,2],2]``) and homogeneous (i.e. all the
            elements are of the same type).

            Also, objects that have some of their dimensions equal to 0
            are not supported. Please use an EArray node (see
            :ref:`EArrayClassDescr`) if you want to store an array with
            one of its dimensions equal to 0.

            The *obj* parameter is optional and it can be provided in
            alternative to the *atom* and *shape* parameters.
            If both *obj* and *atom* and/or *shape* are provided they must
            be consistent with each other.

            .. versionadded:: 3.0

        track_times
            Whether time data associated with the leaf are recorded (object
            access time, raw data modification time, metadata change time,
            object birth time); default True.  Semantics of these times
            depend on their implementation in the HDF5 library: refer to
            documentation of the H5O_info_t data structure.  As of HDF5
            1.8.15, only ctime (metadata change time) is implemented.

            .. versionadded:: 3.4.3

        See Also
        --------
        CArray : for more information on chunked arrays

        """
        if obj is not None:
            flavor = flavor_of(obj)
            obj = array_as_internal(obj, flavor)

            if shape is not None and shape != obj.shape:
                raise TypeError("the shape parameter do not match obj.shape")
            else:
                shape = obj.shape

            if atom is not None and atom.dtype != obj.dtype:
                raise TypeError(
                    "the 'atom' parameter is not consistent with "
                    "the data type of the 'obj' parameter"
                )
            elif atom is None:
                atom = Atom.from_dtype(obj.dtype)
        else:
            if atom is None and shape is None:
                raise TypeError(
                    "the 'atom' and 'shape' parameters or the 'obj' parameter "
                    "must be provided"
                )

        parentnode = self._get_or_create_path(where, createparents)
        _checkfilters(filters)
        ptobj = CArray(
            parentnode,
            name,
            atom=atom,
            shape=shape,
            title=title,
            filters=filters,
            chunkshape=chunkshape,
            byteorder=byteorder,
            track_times=track_times,
        )

        if obj is not None:
            ptobj[...] = obj

        return ptobj

    def create_earray(
        self,
        where: Group | str,
        name: str,
        atom: Atom | None = None,
        shape: tuple[int, ...] | None = None,
        title: str = "",
        filters: Filters | None = None,
        expectedrows: int = 1000,
        chunkshape: int | tuple[int, ...] | None = None,
        byteorder: str | None = None,
        createparents: bool = False,
        obj: np.ndarray | None = None,
        track_times: bool = True,
    ) -> EArray:
        """Create a new enlargeable array.

        Parameters
        ----------
        where : str or Group
            The parent group from which the new array will hang. It can be a
            path string (for example '/level1/leaf5'), or a Group instance
            (see :ref:`GroupClassDescr`).
        name : str
            The name of the new array
        atom : Atom
            An Atom (see :ref:`AtomClassDescr`) instance representing the
            *type* and *shape* of the atomic objects to be saved.

            .. versionchanged:: 3.0
               The *atom* parameter can be None (default) if *obj* is
               provided.

        shape : tuple
            The shape of the new array.  One (and only one) of the shape
            dimensions *must* be 0.  The dimension being 0 means that the
            resulting EArray object can be extended along it.  Multiple
            enlargeable dimensions are not supported right now.

            .. versionchanged:: 3.0
               The *shape* parameter can be None (default) if *obj* is
               provided.

        title : str, optional
            A description for this node (it sets the TITLE HDF5 attribute on
            disk).
        filters : Filters, optional
            An instance of the Filters class (see :ref:`FiltersClassDescr`)
            that provides information about the desired I/O filters to
            be applied during the life of this object.
        expectedrows : int, optional
            A user estimate about the number of row elements that will be added
            to the growable dimension in the EArray node.  If not provided, the
            default value is EXPECTED_ROWS_EARRAY (see tables/parameters.py).
            If you plan to create either a much smaller or a much bigger array
            try providing a guess; this will optimize the HDF5 B-Tree creation
            and management process time and the amount of memory used.
        chunkshape : tuple, numeric, or None, optional
            The shape of the data chunk to be read or written in a single HDF5
            I/O operation.  Filters are applied to those chunks of data.  The
            dimensionality of chunkshape must be the same as that of shape
            (beware: no dimension should be 0 this time!).  If None, a sensible
            value is calculated based on the expectedrows parameter (which is
            recommended).
        byteorder : str, optional
            The byteorder of the data *on disk*, specified as 'little' or
            'big'. If this is not specified, the byteorder is that of the
            platform.
        createparents : bool, optional
            Whether to create the needed groups for the parent path to exist
            (not done by default).
        obj : python object
            The array or scalar to be saved.  Accepted types are NumPy
            arrays and scalars, as well as native Python sequences and
            scalars, provided that values are regular (i.e. they are
            not like ``[[1,2],2]``) and homogeneous (i.e. all the
            elements are of the same type).

            The *obj* parameter is optional and it can be provided in
            alternative to the *atom* and *shape* parameters.
            If both *obj* and *atom* and/or *shape* are provided they must
            be consistent with each other.

            .. versionadded:: 3.0

        track_times
            Whether time data associated with the leaf are recorded (object
            access time, raw data modification time, metadata change time,
            object birth time); default True.  Semantics of these times
            depend on their implementation in the HDF5 library: refer to
            documentation of the H5O_info_t data structure.  As of HDF5
            1.8.15, only ctime (metadata change time) is implemented.

            .. versionadded:: 3.4.3

        See Also
        --------
        EArray : for more information on enlargeable arrays

        """
        if obj is not None:
            flavor = flavor_of(obj)
            obj = array_as_internal(obj, flavor)

            earray_shape = (0,) + obj.shape[1:]

            if shape is not None and shape != earray_shape:
                raise TypeError(
                    "the shape parameter is not compatible " "with obj.shape."
                )
            else:
                shape = earray_shape

            if atom is not None and atom.dtype != obj.dtype:
                raise TypeError(
                    "the atom parameter is not consistent with "
                    "the data type of the obj parameter"
                )
            elif atom is None:
                atom = Atom.from_dtype(obj.dtype)

        parentnode = self._get_or_create_path(where, createparents)
        _checkfilters(filters)
        ptobj = EArray(
            parentnode,
            name,
            atom=atom,
            shape=shape,
            title=title,
            filters=filters,
            expectedrows=expectedrows,
            chunkshape=chunkshape,
            byteorder=byteorder,
            track_times=track_times,
        )

        if obj is not None:
            ptobj.append(obj)

        return ptobj

    def create_vlarray(
        self,
        where: Group | str,
        name: str,
        atom: Atom | None = None,
        title: str = "",
        filters: Filters | None = None,
        expectedrows: int | None = None,
        chunkshape: int | tuple[int, ...] | None = None,
        byteorder: str | None = None,
        createparents: bool = False,
        obj: np.ndarray | None = None,
        track_times: bool = True,
    ) -> VLArray:
        """Create a new variable-length array.

        Parameters
        ----------
        where : str or Group
            The parent group from which the new array will hang. It can
            be a path string (for example '/level1/leaf5'), or a Group
            instance (see :ref:`GroupClassDescr`).
        name : str
            The name of the new array
        atom : Atom
            An Atom (see :ref:`AtomClassDescr`) instance representing
            the *type* and *shape* of the atomic objects to be saved.

            .. versionchanged:: 3.0
               The *atom* parameter can be None (default) if *obj* is
               provided.

        title : str, optional
            A description for this node (it sets the TITLE HDF5 attribute
            on disk).
        filters : Filters
            An instance of the Filters class (see :ref:`FiltersClassDescr`)
            that provides information about the desired I/O filters to
            be applied during the life of this object.
        expectedrows : int, optional
            A user estimate about the number of row elements that will
            be added to the growable dimension in the `VLArray` node.
            If not provided, the default value is ``EXPECTED_ROWS_VLARRAY``
            (see ``tables/parameters.py``).  If you plan to create either
            a much smaller or a much bigger `VLArray` try providing a guess;
            this will optimize the HDF5 B-Tree creation and management
            process time and the amount of memory used.

            .. versionadded:: 3.0

        chunkshape : int or tuple of int, optional
            The shape of the data chunk to be read or written in a
            single HDF5 I/O operation. Filters are applied to those
            chunks of data. The dimensionality of chunkshape must be 1.
            If None, a sensible value is calculated (which is recommended).
        byteorder : str, optional
            The byteorder of the data *on disk*, specified as 'little' or
            'big'. If this is not specified, the byteorder is that of the
            platform.
        createparents : bool, optional
            Whether to create the needed groups for the parent path to
            exist (not done by default).
        obj : python object
            The array or scalar to be saved.  Accepted types are NumPy
            arrays and scalars, as well as native Python sequences and
            scalars, provided that values are regular (i.e. they are
            not like ``[[1,2],2]``) and homogeneous (i.e. all the
            elements are of the same type).

            The *obj* parameter is optional and it can be provided in
            alternative to the *atom* parameter.
            If both *obj* and *atom* and are provided they must
            be consistent with each other.

            .. versionadded:: 3.0

        track_times
            Whether time data associated with the leaf are recorded (object
            access time, raw data modification time, metadata change time,
            object birth time); default True.  Semantics of these times
            depend on their implementation in the HDF5 library: refer to
            documentation of the H5O_info_t data structure.  As of HDF5
            1.8.15, only ctime (metadata change time) is implemented.

            .. versionadded:: 3.4.3

        See Also
        --------
        VLArray : for more informationon variable-length arrays

        .. versionchanged:: 3.0
           The *expectedsizeinMB* parameter has been replaced by
           *expectedrows*.

        """
        if obj is not None:
            flavor = flavor_of(obj)
            obj = array_as_internal(obj, flavor)

            if atom is not None and atom.dtype != obj.dtype:
                raise TypeError(
                    "the atom parameter is not consistent with "
                    "the data type of the obj parameter"
                )
            if atom is None:
                atom = Atom.from_dtype(obj.dtype)
        elif atom is None:
            raise ValueError("atom parameter cannot be None")

        parentnode = self._get_or_create_path(where, createparents)
        _checkfilters(filters)
        ptobj = VLArray(
            parentnode,
            name,
            atom=atom,
            title=title,
            filters=filters,
            expectedrows=expectedrows,
            chunkshape=chunkshape,
            byteorder=byteorder,
            track_times=track_times,
        )

        if obj is not None:
            ptobj.append(obj)

        return ptobj

    def create_hard_link(
        self,
        where: Node | str,
        name: str,
        target: Node | str,
        createparents: bool = False,
    ) -> Group | Leaf:
        """Create a hard link.

        Create a hard link to a `target` node with the given `name` in
        `where` location.  `target` can be a node object or a path
        string.  If `createparents` is true, the intermediate groups
        required for reaching `where` are created (the default is not
        doing so).

        The returned node is a regular `Group` or `Leaf` instance.

        """
        targetnode = self.get_node(target)
        parentnode = self._get_or_create_path(where, createparents)
        linkextension._g_create_hard_link(parentnode, name, targetnode)
        # Refresh children names in link's parent node
        parentnode._g_add_children_names()
        # Return the target node
        return self.get_node(parentnode, name)

    def create_soft_link(
        self,
        where: Node | str,
        name: str,
        target: Node | str,
        createparents: bool = False,
    ) -> SoftLink:
        """Create a soft link (aka symbolic link) to a `target` node.

        Create a soft link (aka symbolic link) to a `target` nodewith
        the given `name` in `where` location.  `target` can be a node
        object or a path string.  If `createparents` is true, the
        intermediate groups required for reaching `where` are created.

        (the default is not doing so).

        The returned node is a SoftLink instance.  See the SoftLink
        class (in :ref:`SoftLinkClassDescr`) for more information on
        soft links.

        """
        if not isinstance(target, str):
            if hasattr(target, "_v_pathname"):  # quacks like a Node
                target = target._v_pathname
            else:
                raise ValueError(
                    "`target` has to be a string or a node object"
                )
        parentnode = self._get_or_create_path(where, createparents)
        slink = SoftLink(parentnode, name, target)
        # Refresh children names in link's parent node
        parentnode._g_add_children_names()
        return slink

    def create_external_link(
        self,
        where: Node | str,
        name: str,
        target: Node | str,
        createparents: bool = False,
    ) -> ExternalLink:
        """Create an external link.

        Create an external link to a *target* node with the given *name*
        in *where* location.  *target* can be a node object in another
        file or a path string in the form 'file:/path/to/node'.  If
        *createparents* is true, the intermediate groups required for
        reaching *where* are created (the default is not doing so).

        The returned node is an :class:`ExternalLink` instance.

        """
        if not isinstance(target, str):
            if hasattr(target, "_v_pathname"):  # quacks like a Node
                target = target._v_file.filename + ":" + target._v_pathname
            else:
                raise ValueError(
                    "`target` has to be a string or a node object"
                )
        elif target.find(":/") == -1:
            raise ValueError("`target` must expressed as 'file:/path/to/node'")
        parentnode = self._get_or_create_path(where, createparents)
        elink = ExternalLink(parentnode, name, target)
        # Refresh children names in link's parent node
        parentnode._g_add_children_names()
        return elink

    def _get_node(self, nodepath: str) -> Node | RootGroup:
        # The root node is always at hand.
        if nodepath == "/":
            return self.root

        node = self._node_manager.get_node(nodepath)
        assert node is not None, "unable to instantiate node ``%s``" % nodepath

        return node

    def get_node(
        self,
        where: Node | str,
        name: str | None = None,
        classname: str | None = None,
    ) -> Node:
        """Get the node under where with the given name.

        Parameters
        ----------
        where : str or Node
            This can be a path string leading to a node or a Node instance (see
            :ref:`NodeClassDescr`). If no name is specified, that node is
            returned.

            .. note::

                If where is a Node instance from a different file than the one
                on which this function is called, the returned node will also
                be from that other file.

        name : str, optional
            If a name is specified, this must be a string with the name of
            a node under where.  In this case the where argument can only
            lead to a Group (see :ref:`GroupClassDescr`) instance (else a
            TypeError is raised). The node called name under the group
            where is returned.
        classname : str, optional
            If the classname argument is specified, it must be the name of
            a class derived from Node (e.g. Table). If the node is found but it
            is not an instance of that class, a NoSuchNodeError is also raised.

        Notes
        -----
        If the node to be returned does not exist, a NoSuchNodeError is
        raised. Please note that hidden nodes are also considered.

        """
        self._check_open()

        if isinstance(where, Node):
            where._g_check_open()

            basepath = where._v_pathname
            nodepath = join_path(basepath, name or "") or "/"
            node = where._v_file._get_node(nodepath)
        elif isinstance(where, (str, np.str_)):
            if not where.startswith("/"):
                raise NameError("``where`` must start with a slash ('/')")

            basepath = where
            nodepath = join_path(basepath, name or "") or "/"
            node = self._get_node(nodepath)
        else:
            raise TypeError(f"``where`` must be a string or a node: {where!r}")

        # Finally, check whether the desired node is an instance
        # of the expected class.
        if classname:
            class_ = get_class_by_name(classname)
            if not isinstance(node, class_):
                npathname = node._v_pathname
                nclassname = node.__class__.__name__
                # This error message is right since it can never be shown
                # for ``classname in [None, 'Node']``.
                raise NoSuchNodeError(
                    "could not find a ``%s`` node at ``%s``; "
                    "instead, a ``%s`` node has been found there"
                    % (classname, npathname, nclassname)
                )

        return node

    def is_visible_node(self, path: str) -> bool:
        """Return True if the node under `path` is visible.

        If the node does not exist, a NoSuchNodeError is raised.

        """
        # ``util.isvisiblepath()`` is still recommended for internal use.
        return self.get_node(path)._f_isvisible()

    def rename_node(
        self,
        where: Node | str,
        newname: str,
        name: str | None = None,
        overwrite: bool = False,
    ) -> None:
        """Change the name of the node specified by where and name to newname.

        Parameters
        ----------
        where, name
            These arguments work as in
            :meth:`File.get_node`, referencing the node to be acted upon.
        newname : str
            The new name to be assigned to the node (a string).
        overwrite : bool
            Whether to recursively remove a node with the same
            newname if it already exists (not done by default).

        """
        obj = self.get_node(where, name=name)
        obj._f_rename(newname, overwrite)

    def move_node(
        self,
        where: Node | str,
        newparent: Group | str | None = None,
        newname: str | None = None,
        name: str | None = None,
        overwrite: bool = False,
        createparents: bool = False,
    ) -> None:
        """Move the node specified by where and name to newparent/newname.

        Parameters
        ----------
        where, name : path
            These arguments work as in
            :meth:`File.get_node`, referencing the node to be acted upon.
        newparent
            The destination group the node will be moved into (a
            path name or a Group instance). If it is
            not specified or None, the current parent
            group is chosen as the new parent.
        newname
            The new name to be assigned to the node in its
            destination (a string). If it is not specified or
            None, the current name is chosen as the
            new name.
        overwrite : bool, optional
            Whether to allow moving over another node, in which case that
            node is recursively removed before moving (not done by default).
        createparents : bool, optional
            If True, any necessary parents of newparent will be created.
            Defaults to False.

        Notes
        -----
        The other arguments work as in :meth:`Node._f_move`.

        """
        obj = self.get_node(where, name=name)
        obj._f_move(newparent, newname, overwrite, createparents)

    def copy_node(
        self,
        where: Node | str,
        newparent: Group | str | None = None,
        newname: str | None = None,
        name: str | None = None,
        overwrite: bool = False,
        recursive: bool = False,
        createparents: bool = False,
        **kwargs,
    ) -> Node:
        """Copy the node specified by where and name to newparent/newname.

        Parameters
        ----------
        where : str
            These arguments work as in
            :meth:`File.get_node`, referencing the node to be acted
            upon.
        newparent : str or Group
            The destination group that the node will be copied
            into (a path name or a Group
            instance). If not specified or None, the
            current parent group is chosen as the new parent.
        newname : str
            The name to be assigned to the new copy in its
            destination (a string).  If it is not specified or
            None, the current name is chosen as the
            new name.
        name : str
            These arguments work as in
            :meth:`File.get_node`, referencing the node to be acted
            upon.
        overwrite : bool, optional
            If True, the destination group will be overwritten if it already
            exists.  Defaults to False.
        recursive : bool, optional
            If True, all descendant nodes of the node to be acted upon are
            recursively copied.  Defaults to False.
        createparents : bool, optional
            If True, any necessary parents of newparent will be created.
            Defaults to False.
        kwargs
           Additional keyword arguments can be used to customize the copying
           process.  See the documentation of :meth:`Group._f_copy`
           for a description of those arguments.

        Returns
        -------
        node : Node
            The newly created copy of the source node (i.e. the destination
            node).  See :meth:`.Node._f_copy` for further details on the
            semantics of copying nodes.

        """
        obj = self.get_node(where, name=name)
        if obj._v_depth == 0 and newparent and not newname:
            npobj = self.get_node(newparent)
            if obj._v_file is not npobj._v_file:
                # Special case for copying file1:/ --> file2:/path
                self.root._f_copy_children(
                    npobj, overwrite=overwrite, recursive=recursive, **kwargs
                )
                return npobj
            else:
                raise OSError(
                    "You cannot copy a root group over the same file"
                )
        return obj._f_copy(
            newparent, newname, overwrite, recursive, createparents, **kwargs
        )

    def remove_node(
        self,
        where: Node | str,
        name: str | None = None,
        recursive: bool = False,
    ) -> None:
        """Remove the object node *name* under *where* location.

        Parameters
        ----------
        where, name
            These arguments work as in
            :meth:`File.get_node`, referencing the node to be acted upon.
        recursive : bool
            If not supplied or false, the node will be removed
            only if it has no children; if it does, a
            NodeError will be raised. If supplied
            with a true value, the node and all its descendants will be
            completely removed.

        """
        obj = self.get_node(where, name=name)
        obj._f_remove(recursive)

    def get_node_attr(
        self, where: Node | str, attrname: str, name: str | None = None
    ) -> Any:
        """Get a PyTables attribute from the given node.

        Parameters
        ----------
        where, name
            These arguments work as in :meth:`File.get_node`, referencing the
            node to be acted upon.
        attrname
            The name of the attribute to retrieve.  If the named
            attribute does not exist, an AttributeError is raised.

        """
        obj = self.get_node(where, name=name)
        return obj._f_getattr(attrname)

    def set_node_attr(
        self,
        where: Node | str,
        attrname: str,
        attrvalue: Any,
        name: str | None = None,
    ) -> None:
        """Set a PyTables attribute for the given node.

        Parameters
        ----------
        where, name
            These arguments work as in
            :meth:`File.get_node`, referencing the node to be acted upon.
        attrname
            The name of the attribute to set.
        attrvalue
            The value of the attribute to set. Any kind of Python
            object (like strings, ints, floats, lists, tuples, dicts,
            small NumPy objects ...) can be stored as an attribute.
            However, if necessary, pickle is automatically used so as
            to serialize objects that you might want to save.
            See the :class:`AttributeSet` class for details.

        Notes
        -----
        If the node already has a large number of attributes, a
        PerformanceWarning is issued.

        """
        obj = self.get_node(where, name=name)
        obj._f_setattr(attrname, attrvalue)

    def del_node_attr(
        self, where: Node | str, attrname: str, name: str | None = None
    ) -> None:
        """Delete a PyTables attribute from the given node.

        Parameters
        ----------
        where, name
            These arguments work as in :meth:`File.get_node`, referencing the
            node to be acted upon.
        attrname
            The name of the attribute to delete.  If the named
            attribute does not exist, an AttributeError is raised.

        """
        obj = self.get_node(where, name=name)
        obj._f_delattr(attrname)

    def copy_node_attrs(
        self, where: Node | str, dstnode: Node | str, name: str | None = None
    ) -> None:
        """Copy PyTables attributes from one node to another.

        Parameters
        ----------
        where, name
            These arguments work as in :meth:`File.get_node`, referencing the
            node to be acted upon.
        dstnode
            The destination node where the attributes will be copied to. It can
            be a path string or a Node instance (see :ref:`NodeClassDescr`).

        """
        srcobject = self.get_node(where, name=name)
        dstobject = self.get_node(dstnode)
        srcobject._v_attrs._f_copy(dstobject)

    def copy_children(
        self,
        srcgroup: str,
        dstgroup: str,
        overwrite: bool = False,
        recursive: bool = False,
        createparents: bool = False,
        **kwargs,
    ) -> None:
        """Copy the children of a group into another group.

        Parameters
        ----------
        srcgroup : str
            The group to copy from.
        dstgroup : str
            The destination group.
        overwrite : bool, optional
            If True, the destination group will be overwritten if it already
            exists.  Defaults to False.
        recursive : bool, optional
            If True, all descendant nodes of srcgroup are recursively copied.
            Defaults to False.
        createparents : bool, optional
            If True, any necessary parents of dstgroup will be created.
            Defaults to False.
        kwargs : dict
           Additional keyword arguments can be used to customize the copying
           process.  See the documentation of :meth:`Group._f_copy_children`
           for a description of those arguments.

        """
        srcgroup = self.get_node(srcgroup)  # Does the source node exist?
        self._check_group(srcgroup)  # Is it a group?

        srcgroup._f_copy_children(
            dstgroup, overwrite, recursive, createparents, **kwargs
        )

    def copy_file(
        self, dstfilename: str, overwrite: bool = False, **kwargs
    ) -> None:
        """Copy the contents of this file to dstfilename.

        Parameters
        ----------
        dstfilename : str
            A path string indicating the name of the destination file. If
            it already exists, the copy will fail with an IOError, unless
            the overwrite argument is true.
        overwrite : bool, optional
            If true, the destination file will be overwritten if it already
            exists.  In this case, the destination file must be closed, or
            errors will occur.  Defaults to False.
        kwargs
            Additional keyword arguments discussed below.

        Notes
        -----
        Additional keyword arguments may be passed to customize the
        copying process. For instance, title and filters may be changed,
        user attributes may be or may not be copied, data may be
        sub-sampled, stats may be collected, etc. Arguments unknown to
        nodes are simply ignored. Check the documentation for copying
        operations of nodes to see which options they support.

        In addition, it recognizes the names of parameters present in
        :file:`tables/parameters.py` as additional keyword arguments.
        See :ref:`parameter_files` for a detailed info on the supported
        parameters.

        Copying a file usually has the beneficial side effect of
        creating a more compact and cleaner version of the original
        file.

        """
        self._check_open()

        # Check that we are not treading our own shoes
        if Path(self.filename).resolve() == Path(dstfilename).resolve():
            raise OSError("You cannot copy a file over itself")

        # Compute default arguments.
        # These are *not* passed on.
        filters = kwargs.pop("filters", None)
        if filters is None:
            # By checking the HDF5 attribute, we avoid setting filters
            # in the destination file if not explicitly set in the
            # source file.  Just by assigning ``self.filters`` we would
            # not be able to tell.
            filters = getattr(self.root._v_attrs, "FILTERS", None)
        copyuserattrs = kwargs.get("copyuserattrs", True)
        title = kwargs.pop("title", self.title)

        if Path(dstfilename).is_file() and not overwrite:
            raise OSError(
                f"file ``{dstfilename}`` already exists; you may want to "
                f"use the ``overwrite`` argument"
            )

        # Create destination file, overwriting it.
        dstfileh = open_file(
            dstfilename, mode="w", title=title, filters=filters, **kwargs
        )

        try:
            # Maybe copy the user attributes of the root group.
            if copyuserattrs:
                self.root._v_attrs._f_copy(dstfileh.root)

            # Copy the rest of the hierarchy.
            self.root._f_copy_children(dstfileh.root, recursive=True, **kwargs)
        finally:
            dstfileh.close()

    def list_nodes(
        self, where: Node | str, classname: str | None = None
    ) -> list[Node]:
        """Return a *list* with children nodes hanging from where.

        This is a list-returning version of :meth:`File.iter_nodes`.

        """
        group = self.get_node(where)  # Does the parent exist?
        self._check_group(group)  # Is it a group?

        return group._f_list_nodes(classname)

    def iter_nodes(
        self, where: Node | str, classname: str | None = None
    ) -> Generator[Node]:
        """Iterate over children nodes hanging from where.

        Parameters
        ----------
        where
            This argument works as in :meth:`File.get_node`, referencing the
            node to be acted upon.
        classname
            If the name of a class derived from
            Node (see :ref:`NodeClassDescr`) is supplied, only instances of
            that class (or subclasses of it) will be returned.

        Notes
        -----
        The returned nodes are alphanumerically sorted by their name.
        This is an iterator version of :meth:`File.list_nodes`.

        """
        group = self.get_node(where)  # Does the parent exist?
        self._check_group(group)  # Is it a group?

        return group._f_iter_nodes(classname)

    def __contains__(self, path: str) -> bool:
        """Return True if there is a node with the specified path.

        Returns True if the file has a node with the given path (a
        string), False otherwise.

        """
        try:
            self.get_node(path)
        except NoSuchNodeError:
            return False
        else:
            return True

    def __iter__(self) -> Generator[Node]:
        """Recursively iterate over the nodes in the tree.

        This is equivalent to calling :meth:`File.walk_nodes` with no
        arguments.

        Examples
        --------
        ::

            # Recursively list all the nodes in the object tree.
            h5file = tables.open_file('vlarray1.h5')
            print("All nodes in the object tree:")
            for node in h5file:
                print(node)

        """
        return self.walk_nodes("/")

    def walk_nodes(
        self, where: Group | str = "/", classname: str | None = None
    ) -> Generator[Node]:
        """Recursively iterate over nodes hanging from where.

        Parameters
        ----------
        where : str or Group, optional
            If supplied, the iteration starts from (and includes)
            this group. It can be a path string or a
            Group instance (see :ref:`GroupClassDescr`).
        classname
            If the name of a class derived from
            Node (see :ref:`GroupClassDescr`) is supplied, only instances of
            that class (or subclasses of it) will be returned.

        Notes
        -----
        This version iterates over the leaves in the same group in order
        to avoid having a list referencing to them and thus, preventing
        the LRU cache to remove them after their use.

        Examples
        --------
        ::

            # Recursively print all the nodes hanging from '/detector'.
            print("Nodes hanging from group '/detector':")
            for node in h5file.walk_nodes('/detector', classname='EArray'):
                print(node)

        """
        class_ = get_class_by_name(classname)

        if class_ is Group:  # only groups
            yield from self.walk_groups(where)
        elif class_ is Node:  # all nodes
            yield self.get_node(where)
            for group in self.walk_groups(where):
                yield from self.iter_nodes(group)
        else:  # only nodes of the named type
            for group in self.walk_groups(where):
                yield from self.iter_nodes(group, classname)

    def walk_groups(self, where: Group | str = "/") -> Generator[Group]:
        """Recursively iterate over groups (not leaves) hanging from where.

        The where group itself is listed first (preorder), then each of its
        child groups (following an alphanumerical order) is also traversed,
        following the same procedure.  If where is not supplied, the root
        group is used.

        The where argument can be a path string
        or a Group instance (see :ref:`GroupClassDescr`).

        """
        group = self.get_node(where)  # Does the parent exist?
        self._check_group(group)  # Is it a group?
        return group._f_walk_groups()

    def _check_open(self) -> None:
        """Check the state of the file.

        If the file is closed, a `ClosedFileError` is raised.

        """
        if not self.isopen:
            raise ClosedFileError("the file object is closed")

    def _iswritable(self) -> bool:
        """Return True if the file is writable."""
        return self.mode in ("w", "a", "r+")

    def _check_writable(self) -> None:
        """Check whether the file is writable.

        If the file is not writable, a `FileModeError` is raised.

        """
        if not self._iswritable():
            raise FileModeError("the file is not writable")

    def _check_group(self, node: Group) -> None:
        # `node` must already be a node.
        if not isinstance(node, Group):
            raise TypeError(f"node ``{node._v_pathname}`` is not a group")

    def is_undo_enabled(self) -> bool:
        """Return True uf the Undo/Redo mechanism is enabled.

        Returns True if the Undo/Redo mechanism has been enabled for
        this file, False otherwise. Please note that this mechanism is
        persistent, so a newly opened PyTables file may already have
        Undo/Redo support enabled.

        """
        self._check_open()
        return self._undoEnabled

    def _check_undo_enabled(self) -> None:
        if not self._undoEnabled:
            raise UndoRedoError("Undo/Redo feature is currently disabled!")

    def _create_transaction_group(self) -> TransactionGroupG:
        tgroup = TransactionGroupG(
            self.root,
            _trans_group_name,
            "Transaction information container",
            new=True,
        )
        # The format of the transaction container.
        tgroup._v_attrs._g__setattr("FORMATVERSION", _trans_version)
        return tgroup

    def _create_transaction(
        self, troot: TransactionGroupG, tid: int
    ) -> TransactionG:
        return TransactionG(
            troot, _trans_name % tid, f"Transaction number {tid}", new=True
        )

    def _create_mark(self, trans: TransactionG, mid: int) -> MarkG:
        return MarkG(trans, _mark_name % mid, "Mark number %d" % mid, new=True)

    def enable_undo(self, filters: Filters = Filters(complevel=1)) -> None:
        """Enable the Undo/Redo mechanism.

        This operation prepares the database for undoing and redoing
        modifications in the node hierarchy. This
        allows :meth:`File.mark`, :meth:`File.undo`, :meth:`File.redo` and
        other methods to be called.

        The filters argument, when specified,
        must be an instance of class Filters (see :ref:`FiltersClassDescr`) and
        is meant for setting the compression values for the action log. The
        default is having compression enabled, as the gains in terms of
        space can be considerable. You may want to disable compression if
        you want maximum speed for Undo/Redo operations.

        Calling this method when the Undo/Redo mechanism is already
        enabled raises an UndoRedoError.

        """
        maxundo = self.params["MAX_UNDO_PATH_LENGTH"]

        class ActionLog(NotLoggedMixin, Table):
            pass

        class ActionLogDesc(IsDescription):
            opcode = UInt8Col(pos=0)
            arg1 = StringCol(maxundo, pos=1, dflt=b"")
            arg2 = StringCol(maxundo, pos=2, dflt=b"")

        self._check_open()

        # Enabling several times is not allowed to avoid the user having
        # the illusion that a new implicit mark has been created
        # when calling enable_undo for the second time.

        if self.is_undo_enabled():
            raise UndoRedoError("Undo/Redo feature is already enabled!")

        self._markers: dict[str, int] = {}
        self._seqmarkers: list[int] = []
        self._nmarks = 0
        self._curtransaction = 0
        self._curmark = -1  # No marks yet

        # Get the Group for keeping user actions
        try:
            tgroup = self.get_node(_trans_group_path)
        except NodeError:
            # The file is going to be changed.
            self._check_writable()

            # A transaction log group does not exist. Create it
            tgroup = self._create_transaction_group()

            # Create a transaction.
            self._trans = self._create_transaction(
                tgroup, self._curtransaction
            )

            # Create an action log
            self._actionlog = ActionLog(
                tgroup,
                _action_log_name,
                ActionLogDesc,
                "Action log",
                filters=filters,
            )

            # Create an implicit mark
            self._actionlog.append([(_op_to_code["MARK"], str(0), "")])
            self._nmarks += 1
            self._seqmarkers.append(0)  # current action is 0

            # Create a group for mark 0
            self._create_mark(self._trans, 0)
            # Initialize the marker pointer
            self._curmark = int(self._nmarks - 1)
            # Initialize the action pointer
            self._curaction = self._actionlog.nrows - 1
        else:
            # The group seems to exist already
            # Get the default transaction
            self._trans = tgroup._f_get_child(
                _trans_name % self._curtransaction
            )
            # Open the action log and go to the end of it
            self._actionlog = tgroup.actionlog
            for row in self._actionlog:
                if row["opcode"] == _op_to_code["MARK"]:
                    name = row["arg2"].decode("utf-8")
                    self._markers[name] = self._nmarks
                    self._seqmarkers.append(row.nrow)
                    self._nmarks += 1
            # Get the current mark and current action
            self._curmark = int(self._actionlog.attrs.CURMARK)
            self._curaction = self._actionlog.attrs.CURACTION

        # The Undo/Redo mechanism has been enabled.
        self._undoEnabled = True

    def disable_undo(self) -> None:
        """Disable the Undo/Redo mechanism.

        Disabling the Undo/Redo mechanism leaves the database in the
        current state and forgets past and future database states. This
        makes :meth:`File.mark`, :meth:`File.undo`, :meth:`File.redo` and other
        methods fail with an UndoRedoError.

        Calling this method when the Undo/Redo mechanism is already
        disabled raises an UndoRedoError.

        """
        self._check_open()

        if not self.is_undo_enabled():
            raise UndoRedoError("Undo/Redo feature is already disabled!")

        # The file is going to be changed.
        self._check_writable()

        del self._markers
        del self._seqmarkers
        del self._curmark
        del self._curaction
        del self._curtransaction
        del self._nmarks
        del self._actionlog
        # Recursively delete the transaction group
        tnode = self.get_node(_trans_group_path)
        tnode._g_remove(recursive=1)

        # The Undo/Redo mechanism has been disabled.
        self._undoEnabled = False

    def mark(self, name: str | None = None) -> int:
        """Mark the state of the database.

        Creates a mark for the current state of the database. A unique (and
        immutable) identifier for the mark is returned. An optional name (a
        string) can be assigned to the mark. Both the identifier of a mark and
        its name can be used in :meth:`File.undo` and :meth:`File.redo`
        operations. When the name has already been used for another mark,
        an UndoRedoError is raised.

        This method can only be called when the Undo/Redo mechanism has been
        enabled. Otherwise, an UndoRedoError is raised.

        """
        self._check_open()
        self._check_undo_enabled()

        if name is None:
            name = ""
        else:
            if not isinstance(name, str):
                raise TypeError(
                    "Only strings are allowed as mark names. "
                    "You passed object: '%s'" % name
                )
            if name in self._markers:
                raise UndoRedoError(
                    "Name '%s' is already used as a marker "
                    "name. Try another one." % name
                )

            # The file is going to be changed.
            self._check_writable()

            self._markers[name] = self._curmark + 1

        # Create an explicit mark
        # Insert the mark in the action log
        self._log("MARK", str(self._curmark + 1), name)
        self._curmark += 1
        self._nmarks = self._curmark + 1
        self._seqmarkers.append(self._curaction)
        # Create a group for the current mark
        self._create_mark(self._trans, self._curmark)
        return self._curmark

    def _log(self, action: str, *args) -> None:
        """Log an action.

        The `action` must be an all-uppercase string identifying it.
        Arguments must also be strings.

        This method should be called once the action has been completed.

        This method can only be called when the Undo/Redo mechanism has
        been enabled.  Otherwise, an `UndoRedoError` is raised.

        """
        assert self.is_undo_enabled()

        maxundo = self.params["MAX_UNDO_PATH_LENGTH"]
        # Check whether we are at the end of the action log or not
        if self._curaction != self._actionlog.nrows - 1:
            # We are not, so delete the trailing actions
            self._actionlog.remove_rows(
                self._curaction + 1, self._actionlog.nrows
            )
            # Reset the current marker group
            mnode = self.get_node(
                _mark_path % (self._curtransaction, self._curmark)
            )
            mnode._g_reset()
            # Delete the marker groups with backup objects
            for mark in range(self._curmark + 1, self._nmarks):
                mnode = self.get_node(
                    _mark_path % (self._curtransaction, mark)
                )
                mnode._g_remove(recursive=1)
            # Update the new number of marks
            self._nmarks = self._curmark + 1
            self._seqmarkers = self._seqmarkers[: self._nmarks]

        if action not in _op_to_code:  # INTERNAL
            raise UndoRedoError(
                "Action ``%s`` not in ``_op_to_code`` "
                "dictionary: %r" % (action, _op_to_code)
            )

        arg1 = ""
        arg2 = ""
        if len(args) <= 1:
            arg1 = args[0]
        elif len(args) <= 2:
            arg1 = args[0]
            arg2 = args[1]
        else:  # INTERNAL
            raise UndoRedoError(
                "Too many parameters for action log: " "%r"
            ).with_traceback(args)
        if len(arg1) > maxundo or len(arg2) > maxundo:  # INTERNAL
            raise UndoRedoError(
                "Parameter arg1 or arg2 is too long: "
                "(%r, %r)" % (arg1, arg2)
            )
        # print("Logging-->", (action, arg1, arg2))
        self._actionlog.append(
            [(_op_to_code[action], arg1.encode("utf-8"), arg2.encode("utf-8"))]
        )
        self._curaction += 1

    def _get_mark_id(self, mark: int | str) -> int:
        """Get an integer markid from a mark sequence number or name."""
        if isinstance(mark, int):
            markid = mark
        elif isinstance(mark, str):
            if mark not in self._markers:
                lmarkers = sorted(self._markers)
                raise UndoRedoError(
                    "The mark that you have specified has not "
                    "been found in the internal marker list: "
                    "%r" % lmarkers
                )
            markid = self._markers[mark]
        else:
            raise TypeError(
                "Parameter mark can only be an integer or a "
                "string, and you passed a type <%s>" % type(mark)
            )
        # print("markid, self._nmarks:", markid, self._nmarks)
        return markid

    def _get_final_action(self, markid: int) -> int:
        """Get the action to go.

        It does not touch the self private attributes

        """
        if markid > self._nmarks - 1:
            # The required mark is beyond the end of the action log
            # The final action is the last row
            return self._actionlog.nrows
        elif markid <= 0:
            # The required mark is the first one
            # return the first row
            return 0

        return self._seqmarkers[markid]

    def _doundo(self, finalaction: int, direction: int) -> None:
        """Undo/Redo actions up to final action in the specified direction."""
        if direction < 0:
            actionlog = self._actionlog[finalaction + 1 : self._curaction + 1][
                ::-1
            ]
        else:
            actionlog = self._actionlog[self._curaction : finalaction]

        # Uncomment this for debugging
        #         print("curaction, finalaction, direction", \
        #               self._curaction, finalaction, direction)
        for i in range(len(actionlog)):
            if actionlog["opcode"][i] != _op_to_code["MARK"]:
                # undo/redo the action
                if direction > 0:
                    # Uncomment this for debugging
                    # print("redo-->", \
                    #       _code_to_op[actionlog['opcode'][i]],\
                    #       actionlog['arg1'][i],\
                    #       actionlog['arg2'][i])
                    undoredo.redo(
                        self,
                        # _code_to_op[actionlog['opcode'][i]],
                        # The next is a workaround for python < 2.5
                        _code_to_op[int(actionlog["opcode"][i])],
                        actionlog["arg1"][i].decode("utf8"),
                        actionlog["arg2"][i].decode("utf8"),
                    )
                else:
                    # Uncomment this for debugging
                    # print("undo-->", \
                    #       _code_to_op[actionlog['opcode'][i]],\
                    #       actionlog['arg1'][i].decode('utf8'),\
                    #       actionlog['arg2'][i].decode('utf8'))
                    undoredo.undo(
                        self,
                        # _code_to_op[actionlog['opcode'][i]],
                        # The next is a workaround for python < 2.5
                        _code_to_op[int(actionlog["opcode"][i])],
                        actionlog["arg1"][i].decode("utf8"),
                        actionlog["arg2"][i].decode("utf8"),
                    )
            else:
                if direction > 0:
                    self._curmark = int(actionlog["arg1"][i])
                else:
                    self._curmark = int(actionlog["arg1"][i]) - 1
                    # Protection against negative marks
                    if self._curmark < 0:
                        self._curmark = 0
            self._curaction += direction

    def undo(self, mark: int | str | None = None) -> None:
        """Go to a past state of the database.

        Returns the database to the state associated with the specified mark.
        Both the identifier of a mark and its name can be used. If the mark is
        omitted, the last created mark is used. If there are no past
        marks, or the specified mark is not older than the current one, an
        UndoRedoError is raised.

        This method can only be called when the Undo/Redo mechanism
        has been enabled. Otherwise, an UndoRedoError
        is raised.

        """
        self._check_open()
        self._check_undo_enabled()

        #         print("(pre)UNDO: (curaction, curmark) = (%s,%s)" % \
        #               (self._curaction, self._curmark))
        if mark is None:
            markid = self._curmark
            # Correction if we are settled on top of a mark
            opcode = self._actionlog.cols.opcode
            if opcode[self._curaction] == _op_to_code["MARK"]:
                markid -= 1
        else:
            # Get the mark ID number
            markid = self._get_mark_id(mark)
        # Get the final action ID to go
        finalaction = self._get_final_action(markid)
        if finalaction > self._curaction:
            raise UndoRedoError(
                "Mark ``%s`` is newer than the current mark. "
                "Use `redo()` or `goto()` instead." % (mark,)
            )

        # The file is going to be changed.
        self._check_writable()

        # Try to reach this mark by unwinding actions in the log
        self._doundo(finalaction - 1, -1)
        if self._curaction < self._actionlog.nrows - 1:
            self._curaction += 1
        self._curmark = int(self._actionlog.cols.arg1[self._curaction])

    #         print("(post)UNDO: (curaction, curmark) = (%s,%s)" % \
    #               (self._curaction, self._curmark))

    def redo(self, mark: int | str | None = None) -> None:
        """Go to a future state of the database.

        Returns the database to the state associated with the specified
        mark.  Both the identifier of a mark and its name can be used.
        If the `mark` is omitted, the next created mark is used.  If
        there are no future marks, or the specified mark is not newer
        than the current one, an UndoRedoError is raised.

        This method can only be called when the Undo/Redo mechanism has
        been enabled.  Otherwise, an UndoRedoError is raised.

        """
        self._check_open()
        self._check_undo_enabled()

        #         print("(pre)REDO: (curaction, curmark) = (%s, %s)" % \
        #               (self._curaction, self._curmark))
        if self._curaction >= self._actionlog.nrows - 1:
            # We are at the end of log, so no action
            return

        if mark is None:
            mark = self._curmark + 1
        elif mark == -1:
            mark = int(self._nmarks)  # Go beyond the mark bounds up to the end
        # Get the mark ID number
        markid = self._get_mark_id(mark)
        finalaction = self._get_final_action(markid)
        if finalaction < self._curaction + 1:
            raise UndoRedoError(
                "Mark ``%s`` is older than the current mark. "
                "Use `redo()` or `goto()` instead." % (mark,)
            )

        # The file is going to be changed.
        self._check_writable()

        # Get the final action ID to go
        self._curaction += 1

        # Try to reach this mark by redoing the actions in the log
        self._doundo(finalaction, 1)
        # Increment the current mark only if we are not at the end of marks
        if self._curmark < self._nmarks - 1:
            self._curmark += 1
        if self._curaction > self._actionlog.nrows - 1:
            self._curaction = self._actionlog.nrows - 1

        # print("(post)REDO: (curaction, curmark) = (%s,%s)" % \
        #       (self._curaction, self._curmark))

    def goto(self, mark: int | str) -> None:
        """Go to a specific mark of the database.

        Returns the database to the state associated with the specified mark.
        Both the identifier of a mark and its name can be used.

        This method can only be called when the Undo/Redo mechanism has been
        enabled. Otherwise, an UndoRedoError is raised.

        """
        self._check_open()
        self._check_undo_enabled()

        if mark == -1:  # Special case
            mark = self._nmarks  # Go beyond the mark bounds up to the end
        # Get the mark ID number
        markid = self._get_mark_id(mark)
        finalaction = self._get_final_action(markid)
        if finalaction < self._curaction:
            self.undo(mark)
        else:
            self.redo(mark)

    def get_current_mark(self) -> int:
        """Get the identifier of the current mark.

        Returns the identifier of the current mark. This can be used
        to know the state of a database after an application crash, or to
        get the identifier of the initial implicit mark after a call
        to :meth:`File.enable_undo`.

        This method can only be called when the Undo/Redo mechanism
        has been enabled. Otherwise, an UndoRedoError
        is raised.

        """
        self._check_open()
        self._check_undo_enabled()
        return self._curmark

    def _shadow_name(self) -> tuple[Node, str]:
        """Compute and return a shadow name.

        Computes the current shadow name according to the current
        transaction, mark and action.  It returns a tuple with the
        shadow parent node and the name of the shadow in it.

        """
        parent = self.get_node(
            _shadow_parent % (self._curtransaction, self._curmark)
        )
        name = _shadow_name % (self._curaction,)

        return (parent, name)

    def flush(self) -> None:
        """Flush all the alive leaves in the object tree."""
        self._check_open()

        # Flush the cache to disk
        self._node_manager.flush_nodes()
        self._flush_file(0)  # 0 means local scope, 1 global (virtual) scope

    def close(self) -> None:
        """Flush all the alive leaves in object tree and close the file."""
        # If the file is already closed, return immediately
        if not self.isopen:
            return

        # If this file has been opened more than once, decrease the
        # counter and return
        if self._open_count > 1:
            self._open_count -= 1
            return

        filename = self.filename

        if self._undoEnabled and self._iswritable():
            # Save the current mark and current action
            self._actionlog.attrs._g__setattr("CURMARK", self._curmark)
            self._actionlog.attrs._g__setattr("CURACTION", self._curaction)

        # Close all loaded nodes.
        self.root._f_close()

        self._node_manager.shutdown()

        # Post-conditions
        assert (
            len(self._node_manager.cache) == 0
        ), "cached nodes remain after closing: %s" % list(
            self._node_manager.cache
        )

        # No other nodes should have been revived.
        assert (
            len(self._node_manager.registry) == 0
        ), "alive nodes remain after closing: %s" % list(
            self._node_manager.registry
        )

        # Close the file
        self._close_file()

        # After the objects are disconnected, destroy the
        # object dictionary using the brute force ;-)
        # This should help to the garbage collector
        self.__dict__.clear()

        # Set the flag to indicate that the file is closed
        self.isopen = 0

        # Restore the filename attribute that is used by _FileRegistry
        self.filename = filename

        # Delete the entry from the registry of opened files
        _open_files.remove(self)

    def __enter__(self) -> File:
        """Enter a context and return the same file."""
        return self

    def __exit__(self, *exc_info) -> bool:
        """Exit a context and close the file."""
        self.close()
        return False  # do not hide exceptions

    def __str__(self) -> str:
        """Return a short string representation of the object tree.

        Examples
        --------
        ::

            >>> import tables
            >>> f = tables.open_file('tables/tests/Tables_lzo2.h5')
            >>> print(f)
            tables/tests/Tables_lzo2.h5 (File) 'Table Benchmark'
            Last modif.: '...'
            Object Tree:
            / (RootGroup) 'Table Benchmark'
            /tuple0 (Table(100,)lzo(1)) 'This is the table title'
            /group0 (Group) ''
            /group0/tuple1 (Table(100,)lzo(1)) 'This is the table title'
            /group0/group1 (Group) ''
            /group0/group1/tuple2 (Table(100,)lzo(1)) 'This is the table title'
            /group0/group1/group2 (Group) ''
            >>> f.close()

        """
        if not self.isopen:
            return "<closed File>"

        # Print all the nodes (Group and Leaf objects) on object tree
        try:
            date = datetime.datetime.fromtimestamp(
                Path(self.filename).stat().st_mtime, datetime.UTC
            ).isoformat(timespec="seconds")
        except OSError:
            # in-memory file
            date = "<in-memory file>"
        lines = [
            f"{self.filename} (File) {self.title!r}",
            f"Last modif.: {date!r}",
            "Object Tree: ",
        ]

        for group in self.walk_groups("/"):
            lines.append(f"{group}")
            for kind in self._node_kinds[1:]:
                for node in self.list_nodes(group, kind):
                    lines.append(f"{node}")
        return "\n".join(lines) + "\n"

    def __repr__(self) -> str:
        """Return a detailed string representation of the object tree."""
        if not self.isopen:
            return "<closed File>"

        # Print all the nodes (Group and Leaf objects) on object tree
        lines = [
            f"File(filename={self.filename!s}, title={self.title!r}, "
            f"mode={self.mode!r}, root_uep={self.root_uep!r}, "
            f"filters={self.filters!r})"
        ]
        for group in self.walk_groups("/"):
            lines.append(f"{group}")
            for kind in self._node_kinds[1:]:
                for node in self.list_nodes(group, kind):
                    lines.append(f"{node!r}")
        return "\n".join(lines) + "\n"

    def _update_node_locations(self, oldpath: str, newpath: str) -> None:
        """Update location information of nodes under `oldpath`.

        This only affects *already loaded* nodes.

        """
        oldprefix = oldpath + "/"  # root node can not be renamed, anyway
        oldprefix_len = len(oldprefix)

        # Update alive and dead descendents.
        for cache in [self._node_manager.cache, self._node_manager.registry]:
            for nodepath in list(cache):
                if nodepath.startswith(oldprefix) and nodepath != oldprefix:
                    nodesuffix = nodepath[oldprefix_len:]
                    newnodepath = join_path(newpath, nodesuffix)
                    newnodeppath = split_path(newnodepath)[0]
                    descendent_node = self._get_node(nodepath)
                    descendent_node._g_update_location(newnodeppath)


# If a user hits ^C during a run, it is wise to gracefully close the
# opened files.
atexit.register(_open_files.close_all)
