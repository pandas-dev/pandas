"""Here is defined the Group class."""

import os
import weakref
import warnings

from .misc.proxydict import ProxyDict
from . import hdf5extension
from . import utilsextension
from .registry import class_id_dict
from .exceptions import (NodeError, NoSuchNodeError, NaturalNameWarning,
                         PerformanceWarning)
from .filters import Filters
from .registry import get_class_by_name
from .path import check_name_validity, join_path, isvisiblename
from .node import Node, NotLoggedMixin
from .leaf import Leaf
from .unimplemented import UnImplemented, Unknown

from .link import Link, SoftLink, ExternalLink


obversion = "1.0"


class _ChildrenDict(ProxyDict):
    def _get_value_from_container(self, container, key):
        return container._f_get_child(key)


class Group(hdf5extension.Group, Node):
    """Basic PyTables grouping structure.

    Instances of this class are grouping structures containing *child*
    instances of zero or more groups or leaves, together with
    supporting metadata. Each group has exactly one *parent* group.

    Working with groups and leaves is similar in many ways to working
    with directories and files, respectively, in a Unix filesystem.
    As with Unix directories and files, objects in the object tree are
    often described by giving their full (or absolute) path names.
    This full path can be specified either as a string (like in
    '/group1/group2') or as a complete object path written in *natural
    naming* schema (like in file.root.group1.group2).

    A collateral effect of the *natural naming* schema is that the
    names of members in the Group class and its instances must be
    carefully chosen to avoid colliding with existing children node
    names.  For this reason and to avoid polluting the children
    namespace all members in a Group start with some reserved prefix,
    like _f_ (for public methods), _g_ (for private ones), _v_ (for
    instance variables) or _c_ (for class variables). Any attempt to
    create a new child node whose name starts with one of these
    prefixes will raise a ValueError exception.

    Another effect of natural naming is that children named after
    Python keywords or having names not valid as Python identifiers
    (e.g.  class, $a or 44) can not be accessed using the node.child
    syntax. You will be forced to use node._f_get_child(child) to
    access them (which is recommended for programmatic accesses).

    You will also need to use _f_get_child() to access an existing
    child node if you set a Python attribute in the Group with the
    same name as that node (you will get a NaturalNameWarning when
    doing this).

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.
    name : str
        The name of this node in its parent group.
    title
        The title for this group
    new
        If this group is new or has to be read from disk
    filters : Filters
        A Filters instance


    .. versionchanged:: 3.0
       *parentNode* renamed into *parentnode*

    Notes
    -----
    The following documentation includes methods that are automatically
    called when a Group instance is accessed in a special way.

    For instance, this class defines the __setattr__, __getattr__,
    __delattr__ and __dir__ methods, and they set, get and delete
    *ordinary Python attributes* as normally intended. In addition to that,
    __getattr__ allows getting *child nodes* by their name for the sake of
    easy interaction on the command line, as long as there is no Python
    attribute with the same name. Groups also allow the interactive
    completion (when using readline) of the names of child nodes.
    For instance::

        # get a Python attribute
        nchild = group._v_nchildren

        # Add a Table child called 'table' under 'group'.
        h5file.create_table(group, 'table', myDescription)
        table = group.table          # get the table child instance
        group.table = 'foo'          # set a Python attribute

        # (PyTables warns you here about using the name of a child node.)
        foo = group.table            # get a Python attribute
        del group.table              # delete a Python attribute
        table = group.table          # get the table child instance again

    Additionally, on interactive python sessions you may get autocompletions
    of children named as *valid python identifiers* by pressing the  `[Tab]`
    key, or to use the dir() global function.

    .. rubric:: Group attributes

    The following instance variables are provided in addition to those
    in Node (see :ref:`NodeClassDescr`):

    .. attribute:: _v_children

        Dictionary with all nodes hanging from this group.

    .. attribute:: _v_groups

        Dictionary with all groups hanging from this group.

    .. attribute:: _v_hidden

        Dictionary with all hidden nodes hanging from this group.

    .. attribute:: _v_leaves

        Dictionary with all leaves hanging from this group.

    .. attribute:: _v_links

        Dictionary with all links hanging from this group.

    .. attribute:: _v_unknown

        Dictionary with all unknown nodes hanging from this group.

    """

    # Class identifier.
    _c_classid = 'GROUP'

    # Children containers that should be loaded only in a lazy way.
    # These are documented in the ``Group._g_add_children_names`` method.
    _c_lazy_children_attrs = (
        '__members__', '_v_children', '_v_groups', '_v_leaves',
        '_v_links', '_v_unknown', '_v_hidden')

    # `_v_nchildren` is a direct read-only shorthand
    # for the number of *visible* children in a group.
    def _g_getnchildren(self):
        """The number of children hanging from this group."""
        return len(self._v_children)

    _v_nchildren = property(_g_getnchildren)

    # `_v_filters` is a direct read-write shorthand for the ``FILTERS``
    # attribute with the default `Filters` instance as a default value.
    def _g_getfilters(self):
        filters = getattr(self._v_attrs, 'FILTERS', None)
        if filters is None:
            filters = Filters()
        return filters

    def _g_setfilters(self, value):
        if not isinstance(value, Filters):
            raise TypeError(
                f"value is not an instance of `Filters`: {value!r}")
        self._v_attrs.FILTERS = value

    def _g_delfilters(self):
        del self._v_attrs.FILTERS

    _v_filters = property(
        _g_getfilters, _g_setfilters, _g_delfilters,
        """Default filter properties for child nodes.

        You can (and are encouraged to) use this property to get, set and
        delete the FILTERS HDF5 attribute of the group, which stores a Filters
        instance (see :ref:`FiltersClassDescr`). When the group has no such
        attribute, a default Filters instance is used.
        """)

    def __init__(self, parentnode, name,
                 title="", new=False, filters=None,
                 _log=True):

        # Remember to assign these values in the root group constructor
        # if it does not use this one!

        # First, set attributes belonging to group objects.

        self._v_version = obversion
        """The object version of this group."""

        self._v_new = new
        """Is this the first time the node has been created?"""

        self._v_new_title = title
        """New title for this node."""

        self._v_new_filters = filters
        """New default filter properties for child nodes."""

        self._v_max_group_width = parentnode._v_file.params['MAX_GROUP_WIDTH']
        """Maximum number of children on each group before warning the user.

        .. versionchanged:: 3.0
           The *_v_maxGroupWidth* attribute has been renamed into
           *_v_max_group_width*.

        """

        # Finally, set up this object as a node.
        super().__init__(parentnode, name, _log)

    def _g_post_init_hook(self):
        if self._v_new:
            if self._v_file.params['PYTABLES_SYS_ATTRS']:
                # Save some attributes for the new group on disk.
                set_attr = self._v_attrs._g__setattr
                # Set the title, class and version attributes.
                set_attr('TITLE', self._v_new_title)
                set_attr('CLASS', self._c_classid)
                set_attr('VERSION', self._v_version)

                # Set the default filter properties.
                newfilters = self._v_new_filters
                if newfilters is None:
                    # If no filters have been passed in the constructor,
                    # inherit them from the parent group, but only if they
                    # have been inherited or explicitly set.
                    newfilters = getattr(
                        self._v_parent._v_attrs, 'FILTERS', None)
                if newfilters is not None:
                    set_attr('FILTERS', newfilters)
        else:
            # If the file has PyTables format, get the VERSION attr
            if 'VERSION' in self._v_attrs._v_attrnamessys:
                self._v_version = self._v_attrs.VERSION
            else:
                self._v_version = "0.0 (unknown)"
            # We don't need to get more attributes from disk,
            # since the most important ones are defined as properties.

    def __del__(self):
        if (self._v_isopen and
            self._v_pathname in self._v_file._node_manager.registry and
                '_v_children' in self.__dict__):
            # The group is going to be killed.  Rebuild weak references
            # (that Python cancelled just before calling this method) so
            # that they are still usable if the object is revived later.
            selfref = weakref.ref(self)
            self._v_children.containerref = selfref
            self._v_groups.containerref = selfref
            self._v_leaves.containerref = selfref
            self._v_links.containerref = selfref
            self._v_unknown.containerref = selfref
            self._v_hidden.containerref = selfref

        super().__del__()

    def _g_get_child_group_class(self, childname):
        """Get the class of a not-yet-loaded group child.

        `childname` must be the name of a *group* child.

        """

        childCID = self._g_get_gchild_attr(childname, 'CLASS')
        if childCID is not None and not isinstance(childCID, str):
            childCID = childCID.decode('utf-8')

        if childCID in class_id_dict:
            return class_id_dict[childCID]  # look up group class
        else:
            return Group  # default group class

    def _g_get_child_leaf_class(self, childname, warn=True):
        """Get the class of a not-yet-loaded leaf child.

        `childname` must be the name of a *leaf* child.  If the child
        belongs to an unknown kind of leaf, or if its kind can not be
        guessed, `UnImplemented` will be returned and a warning will be
        issued if `warn` is true.

        """

        if self._v_file.params['PYTABLES_SYS_ATTRS']:
            childCID = self._g_get_lchild_attr(childname, 'CLASS')
            if childCID is not None and not isinstance(childCID, str):
                childCID = childCID.decode('utf-8')
        else:
            childCID = None

        if childCID in class_id_dict:
            return class_id_dict[childCID]  # look up leaf class
        else:
            # Unknown or no ``CLASS`` attribute, try a guess.
            childCID2 = utilsextension.which_class(self._v_objectid, childname)
            if childCID2 == 'UNSUPPORTED':
                if warn:
                    if childCID is None:
                        warnings.warn(
                            "leaf ``%s`` is of an unsupported type; "
                            "it will become an ``UnImplemented`` node"
                            % self._g_join(childname))
                    else:
                        warnings.warn(
                            ("leaf ``%s`` has an unknown class ID ``%s``; "
                             "it will become an ``UnImplemented`` node")
                            % (self._g_join(childname), childCID))
                return UnImplemented
            assert childCID2 in class_id_dict
            return class_id_dict[childCID2]  # look up leaf class

    def _g_add_children_names(self):
        """Add children names to this group taking into account their
        visibility and kind."""

        mydict = self.__dict__

        # The names of the lazy attributes
        mydict['__members__'] = members = []
        """The names of visible children nodes for readline-style completion.
        """
        mydict['_v_children'] = children = _ChildrenDict(self)
        """The number of children hanging from this group."""
        mydict['_v_groups'] = groups = _ChildrenDict(self)
        """Dictionary with all groups hanging from this group."""
        mydict['_v_leaves'] = leaves = _ChildrenDict(self)
        """Dictionary with all leaves hanging from this group."""
        mydict['_v_links'] = links = _ChildrenDict(self)
        """Dictionary with all links hanging from this group."""
        mydict['_v_unknown'] = unknown = _ChildrenDict(self)
        """Dictionary with all unknown nodes hanging from this group."""
        mydict['_v_hidden'] = hidden = _ChildrenDict(self)
        """Dictionary with all hidden nodes hanging from this group."""

        # Get the names of *all* child groups and leaves.
        (group_names, leaf_names, link_names, unknown_names) = \
            self._g_list_group(self._v_parent)

        # Separate groups into visible groups and hidden nodes,
        # and leaves into visible leaves and hidden nodes.
        for (childnames, childdict) in ((group_names, groups),
                                        (leaf_names, leaves),
                                        (link_names, links),
                                        (unknown_names, unknown)):

            for childname in childnames:
                # See whether the name implies that the node is hidden.
                # (Assigned values are entirely irrelevant.)
                if isvisiblename(childname):
                    # Visible node.
                    members.insert(0, childname)
                    children[childname] = None
                    childdict[childname] = None
                else:
                    # Hidden node.
                    hidden[childname] = None

    def _g_check_has_child(self, name):
        """Check whether 'name' is a children of 'self' and return its type."""

        # Get the HDF5 name matching the PyTables name.
        node_type = self._g_get_objinfo(name)
        if node_type == "NoSuchNode":
            raise NoSuchNodeError(
                "group ``%s`` does not have a child named ``%s``"
                % (self._v_pathname, name))
        return node_type

    def __iter__(self):
        """Iterate over the child nodes hanging directly from the group.

        This iterator is *not* recursive.

        Examples
        --------

        ::

            # Non-recursively list all the nodes hanging from '/detector'
            print("Nodes in '/detector' group:")
            for node in h5file.root.detector:
                print(node)

        """

        return self._f_iter_nodes()

    def __contains__(self, name):
        """Is there a child with that `name`?

        Returns a true value if the group has a child node (visible or
        hidden) with the given `name` (a string), false otherwise.

        """

        self._g_check_open()
        try:
            self._g_check_has_child(name)
        except NoSuchNodeError:
            return False
        return True

    def __getitem__(self, childname):
        """Return the (visible or hidden) child with that `name` ( a string).

        Raise IndexError if child not exist.
        """
        try:
            return self._f_get_child(childname)
        except NoSuchNodeError:
            raise IndexError(childname)

    def _f_walknodes(self, classname=None):
        """Iterate over descendant nodes.

        This method recursively walks *self* top to bottom (preorder),
        iterating over child groups in alphanumerical order, and yielding
        nodes.  If classname is supplied, only instances of the named class are
        yielded.

        If *classname* is Group, it behaves like :meth:`Group._f_walk_groups`,
        yielding only groups.  If you don't want a recursive behavior,
        use :meth:`Group._f_iter_nodes` instead.

        Examples
        --------

        ::

            # Recursively print all the arrays hanging from '/'
            print("Arrays in the object tree '/':")
            for array in h5file.root._f_walknodes('Array', recursive=True):
                print(array)

        """

        self._g_check_open()

        # For compatibility with old default arguments.
        if classname == '':
            classname = None

        if classname == "Group":
            # Recursive algorithm
            yield from self._f_walk_groups()
        else:
            for group in self._f_walk_groups():
                yield from group._f_iter_nodes(classname)

    def _g_join(self, name):
        """Helper method to correctly concatenate a name child object with the
        pathname of this group."""

        if name == "/":
            # This case can happen when doing copies
            return self._v_pathname
        return join_path(self._v_pathname, name)

    def _g_width_warning(self):
        """Issue a :exc:`PerformanceWarning` on too many children."""

        warnings.warn("""\
group ``%s`` is exceeding the recommended maximum number of children (%d); \
be ready to see PyTables asking for *lots* of memory and possibly slow I/O."""
                      % (self._v_pathname, self._v_max_group_width),
                      PerformanceWarning)

    def _g_refnode(self, childnode, childname, validate=True):
        """Insert references to a `childnode` via a `childname`.

        Checks that the `childname` is valid and does not exist, then
        creates references to the given `childnode` by that `childname`.
        The validation of the name can be omitted by setting `validate`
        to a false value (this may be useful for adding already existing
        nodes to the tree).

        """

        # Check for name validity.
        if validate:
            check_name_validity(childname)
            childnode._g_check_name(childname)

        # Check if there is already a child with the same name.
        # This can be triggered because of the user
        # (via node construction or renaming/movement).
        # Links are not checked here because they are copied and referenced
        # using ``File.get_node`` so they already exist in `self`.
        if (not isinstance(childnode, Link)) and childname in self:
            raise NodeError(
                "group ``%s`` already has a child node named ``%s``"
                % (self._v_pathname, childname))

        # Show a warning if there is an object attribute with that name.
        if childname in self.__dict__:
            warnings.warn(
                "group ``%s`` already has an attribute named ``%s``; "
                "you will not be able to use natural naming "
                "to access the child node"
                % (self._v_pathname, childname), NaturalNameWarning)

        # Check group width limits.
        if (len(self._v_children) + len(self._v_hidden) >=
                self._v_max_group_width):
            self._g_width_warning()

        # Update members information.
        # Insert references to the new child.
        # (Assigned values are entirely irrelevant.)
        if isvisiblename(childname):
            # Visible node.
            self.__members__.insert(0, childname)  # enable completion
            self._v_children[childname] = None  # insert node
            if isinstance(childnode, Unknown):
                self._v_unknown[childname] = None
            elif isinstance(childnode, Link):
                self._v_links[childname] = None
            elif isinstance(childnode, Leaf):
                self._v_leaves[childname] = None
            elif isinstance(childnode, Group):
                self._v_groups[childname] = None
        else:
            # Hidden node.
            self._v_hidden[childname] = None  # insert node

    def _g_unrefnode(self, childname):
        """Remove references to a node.

        Removes all references to the named node.

        """

        # This can *not* be triggered because of the user.
        assert childname in self, \
            ("group ``%s`` does not have a child node named ``%s``"
                % (self._v_pathname, childname))

        # Update members information, if needed
        if '_v_children' in self.__dict__:
            if childname in self._v_children:
                # Visible node.
                members = self.__members__
                member_index = members.index(childname)
                del members[member_index]  # disables completion

                del self._v_children[childname]  # remove node
                self._v_unknown.pop(childname, None)
                self._v_links.pop(childname, None)
                self._v_leaves.pop(childname, None)
                self._v_groups.pop(childname, None)
            else:
                # Hidden node.
                del self._v_hidden[childname]  # remove node

    def _g_move(self, newparent, newname):
        # Move the node to the new location.
        oldpath = self._v_pathname
        super()._g_move(newparent, newname)
        newpath = self._v_pathname

        # Update location information in children.  This node shouldn't
        # be affected since it has already been relocated.
        self._v_file._update_node_locations(oldpath, newpath)

    def _g_copy(self, newparent, newname, recursive, _log=True, **kwargs):
        # Compute default arguments.
        title = kwargs.get('title', self._v_title)
        filters = kwargs.get('filters', None)
        stats = kwargs.get('stats', None)

        # Fix arguments with explicit None values for backwards compatibility.
        if title is None:
            title = self._v_title
        # If no filters have been passed to the call, copy them from the
        # source group, but only if inherited or explicitly set.
        if filters is None:
            filters = getattr(self._v_attrs, 'FILTERS', None)

        # Create a copy of the object.
        new_node = Group(newparent, newname,
                         title, new=True, filters=filters, _log=_log)

        # Copy user attributes if needed.
        if kwargs.get('copyuserattrs', True):
            self._v_attrs._g_copy(new_node._v_attrs, copyclass=True)

        # Update statistics if needed.
        if stats is not None:
            stats['groups'] += 1

        if recursive:
            # Copy child nodes if a recursive copy was requested.
            # Some arguments should *not* be passed to children copy ops.
            kwargs = kwargs.copy()
            kwargs.pop('title', None)
            self._g_copy_children(new_node, **kwargs)

        return new_node

    def _g_copy_children(self, newparent, **kwargs):
        """Copy child nodes.

        Copies all nodes descending from this one into the specified
        `newparent`.  If the new parent has a child node with the same
        name as one of the nodes in this group, the copy fails with a
        `NodeError`, maybe resulting in a partial copy.  Nothing is
        logged.

        """

        # Recursive version of children copy.
        # for srcchild in self._v_children.itervalues():
        #     srcchild._g_copy_as_child(newparent, **kwargs)

        # Non-recursive version of children copy.
        use_hardlinks = kwargs.get('use_hardlinks', False)
        if use_hardlinks:
            address_map = kwargs.setdefault('address_map', {})

        parentstack = [(self, newparent)]  # [(source, destination), ...]
        while parentstack:
            (srcparent, dstparent) = parentstack.pop()

            if use_hardlinks:
                for srcchild in srcparent._v_children.values():
                    addr, rc = srcchild._get_obj_info()
                    if rc > 1 and addr in address_map:
                        where, name = address_map[addr][0]
                        localsrc = os.path.join(where, name)
                        dstparent._v_file.create_hard_link(dstparent,
                                                           srcchild.name,
                                                           localsrc)
                        address_map[addr].append(
                            (dstparent._v_pathname, srcchild.name)
                        )

                        # Update statistics if needed.
                        stats = kwargs.pop('stats', None)
                        if stats is not None:
                            stats['hardlinks'] += 1
                    else:
                        dstchild = srcchild._g_copy_as_child(dstparent,
                                                             **kwargs)
                        if isinstance(srcchild, Group):
                            parentstack.append((srcchild, dstchild))

                        if rc > 1:
                            address_map[addr] = [
                                (dstparent._v_pathname, srcchild.name)
                            ]
            else:
                for srcchild in srcparent._v_children.values():
                    dstchild = srcchild._g_copy_as_child(dstparent, **kwargs)
                    if isinstance(srcchild, Group):
                        parentstack.append((srcchild, dstchild))

    def _f_get_child(self, childname):
        """Get the child called childname of this group.

        If the child exists (be it visible or not), it is returned.  Else, a
        NoSuchNodeError is raised.

        Using this method is recommended over getattr() when doing programmatic
        accesses to children if childname is unknown beforehand or when its
        name is not a valid Python identifier.

        """

        self._g_check_open()

        self._g_check_has_child(childname)

        childpath = join_path(self._v_pathname, childname)
        return self._v_file._get_node(childpath)

    def _f_list_nodes(self, classname=None):
        """Return a *list* with children nodes.

        This is a list-returning version of :meth:`Group._f_iter_nodes()`.

        """

        return list(self._f_iter_nodes(classname))

    def _f_iter_nodes(self, classname=None):
        """Iterate over children nodes.

        Child nodes are yielded alphanumerically sorted by node name.  If the
        name of a class derived from Node (see :ref:`NodeClassDescr`) is
        supplied in the classname parameter, only instances of that class (or
        subclasses of it) will be returned.

        This is an iterator version of :meth:`Group._f_list_nodes`.

        """

        self._g_check_open()

        if not classname:
            # Returns all the children alphanumerically sorted
            for name in sorted(self._v_children):
                yield self._v_children[name]
        elif classname == 'Group':
            # Returns all the groups alphanumerically sorted
            for name in sorted(self._v_groups):
                yield self._v_groups[name]
        elif classname == 'Leaf':
            # Returns all the leaves alphanumerically sorted
            for name in sorted(self._v_leaves):
                yield self._v_leaves[name]
        elif classname == 'Link':
            # Returns all the links alphanumerically sorted
            for name in sorted(self._v_links):
                yield self._v_links[name]
        elif classname == 'IndexArray':
            raise TypeError(
                "listing ``IndexArray`` nodes is not allowed")
        else:
            class_ = get_class_by_name(classname)
            for childname, childnode in sorted(self._v_children.items()):
                if isinstance(childnode, class_):
                    yield childnode

    def _f_walk_groups(self):
        """Recursively iterate over descendent groups (not leaves).

        This method starts by yielding *self*, and then it goes on to
        recursively iterate over all child groups in alphanumerical order, top
        to bottom (preorder), following the same procedure.

        """

        self._g_check_open()

        stack = [self]
        yield self
        # Iterate over the descendants
        while stack:
            objgroup = stack.pop()
            groupnames = sorted(objgroup._v_groups)
            # Sort the groups before delivering. This uses the groups names
            # for groups in tree (in order to sort() can classify them).
            for groupname in groupnames:
                # TODO: check recursion
                stack.append(objgroup._v_groups[groupname])
                yield objgroup._v_groups[groupname]

    def __delattr__(self, name):
        """Delete a Python attribute called name.

        This method only provides a extra warning in case the user
        tries to delete a children node using __delattr__.

        To remove a children node from this group use
        :meth:`File.remove_node` or :meth:`Node._f_remove`. To delete
        a PyTables node attribute use :meth:`File.del_node_attr`,
        :meth:`Node._f_delattr` or :attr:`Node._v_attrs``.

        If there is an attribute and a child node with the same name,
        the child node will be made accessible again via natural naming.

        """

        try:
            super().__delattr__(name)  # nothing particular
        except AttributeError as ae:
            hint = " (use ``node._f_remove()`` if you want to remove a node)"
            raise ae.__class__(str(ae) + hint)

    def __dir__(self):
        """Autocomplete only children named as valid python identifiers.

        Only PY3 supports this special method.
        """
        subnods = [c for c in self._v_children if c.isidentifier()]
        return super().__dir__() + subnods

    def __getattr__(self, name):
        """Get a Python attribute or child node called name.
        If the node has a child node called name it is returned,
        else an AttributeError is raised.
        """

        if name in self._c_lazy_children_attrs:
            self._g_add_children_names()
            return self.__dict__[name]
        return self._f_get_child(name)

    def __setattr__(self, name, value):
        """Set a Python attribute called name with the given value.

        This method stores an *ordinary Python attribute* in the object. It
        does *not* store new children nodes under this group; for that, use the
        File.create*() methods (see the File class
        in :ref:`FileClassDescr`). It does *neither* store a PyTables node
        attribute; for that,
        use :meth:`File.set_node_attr`, :meth`:Node._f_setattr`
        or :attr:`Node._v_attrs`.

        If there is already a child node with the same name, a
        NaturalNameWarning will be issued and the child node will not be
        accessible via natural naming nor getattr(). It will still be available
        via :meth:`File.get_node`, :meth:`Group._f_get_child` and children
        dictionaries in the group (if visible).

        """

        # Show a warning if there is an child node with that name.
        #
        # ..note::
        #
        #   Using ``if name in self:`` is not right since that would
        #   require ``_v_children`` and ``_v_hidden`` to be already set
        #   when the very first attribute assignments are made.
        #   Moreover, this warning is only concerned about clashes with
        #   names used in natural naming, i.e. those in ``__members__``.
        #
        # ..note::
        #
        #   The check ``'__members__' in myDict`` allows attribute
        #   assignment to happen before calling `Group.__init__()`, by
        #   avoiding to look into the still not assigned ``__members__``
        #   attribute.  This allows subclasses to set up some attributes
        #   and then call the constructor of the superclass.  If the
        #   check above is disabled, that results in Python entering an
        #   endless loop on exit!

        mydict = self.__dict__
        if '__members__' in mydict and name in self.__members__:
            warnings.warn(
                "group ``%s`` already has a child node named ``%s``; "
                "you will not be able to use natural naming "
                "to access the child node"
                % (self._v_pathname, name), NaturalNameWarning)

        super().__setattr__(name, value)

    def _f_flush(self):
        """Flush this Group."""

        self._g_check_open()
        self._g_flush_group()

    def _g_close_descendents(self):
        """Close all the *loaded* descendent nodes of this group."""

        node_manager = self._v_file._node_manager
        node_manager.close_subtree(self._v_pathname)

    def _g_close(self):
        """Close this (open) group."""

        if self._v_isopen:
            # hdf5extension operations:
            #   Close HDF5 group.
            self._g_close_group()

        # Close myself as a node.
        super()._f_close()

    def _f_close(self):
        """Close this group and all its descendents.

        This method has the behavior described in :meth:`Node._f_close`.
        It should be noted that this operation closes all the nodes
        descending from this group.

        You should not need to close nodes manually because they are
        automatically opened/closed when they are loaded/evicted from
        the integrated LRU cache.

        """

        # If the group is already closed, return immediately
        if not self._v_isopen:
            return

        # First, close all the descendents of this group, unless a) the
        # group is being deleted (evicted from LRU cache) or b) the node
        # is being closed during an aborted creation, in which cases
        # this is not an explicit close issued by the user.
        if not (self._v__deleting or self._v_objectid is None):
            self._g_close_descendents()

        # When all the descendents have been closed, close this group.
        # This is done at the end because some nodes may still need to
        # be loaded during the closing process; thus this node must be
        # open until the very end.
        self._g_close()

    def _g_remove(self, recursive=False, force=False):
        """Remove (recursively if needed) the Group.

        This version correctly handles both visible and hidden nodes.

        """

        if self._v_nchildren > 0:
            if not (recursive or force):
                raise NodeError("group ``%s`` has child nodes; "
                                "please set `recursive` or `force` to true "
                                "to remove it"
                                % (self._v_pathname,))

            # First close all the descendents hanging from this group,
            # so that it is not possible to use a node that no longer exists.
            self._g_close_descendents()

        # Remove the node itself from the hierarchy.
        super()._g_remove(recursive, force)

    def _f_copy(self, newparent=None, newname=None,
                overwrite=False, recursive=False, createparents=False,
                **kwargs):
        """Copy this node and return the new one.

        This method has the behavior described in :meth:`Node._f_copy`.
        In addition, it recognizes the following keyword arguments:

        Parameters
        ----------
        title
            The new title for the destination. If omitted or None, the
            original title is used. This only applies to the topmost
            node in recursive copies.
        filters : Filters
            Specifying this parameter overrides the original filter
            properties in the source node. If specified, it must be an
            instance of the Filters class (see :ref:`FiltersClassDescr`).
            The default is to copy the filter properties from the source
            node.
        copyuserattrs
            You can prevent the user attributes from being copied by setting
            thisparameter to False. The default is to copy them.
        stats
            This argument may be used to collect statistics on the copy
            process. When used, it should be a dictionary with keys 'groups',
            'leaves', 'links' and 'bytes' having a numeric value. Their values
            willbe incremented to reflect the number of groups, leaves and
            bytes, respectively, that have been copied during the operation.

        """

        return super()._f_copy(
            newparent, newname,
            overwrite, recursive, createparents, **kwargs)

    def _f_copy_children(self, dstgroup, overwrite=False, recursive=False,
                         createparents=False, **kwargs):
        """Copy the children of this group into another group.

        Children hanging directly from this group are copied into dstgroup,
        which can be a Group (see :ref:`GroupClassDescr`) object or its
        pathname in string form. If createparents is true, the needed groups
        for the given destination group path to exist will be created.

        The operation will fail with a NodeError if there is a child node
        in the destination group with the same name as one of the copied
        children from this one, unless overwrite is true; in this case,
        the former child node is recursively removed before copying the
        later.

        By default, nodes descending from children groups of this node
        are not copied. If the recursive argument is true, all descendant
        nodes of this node are recursively copied.

        Additional keyword arguments may be passed to customize the
        copying process. For instance, title and filters may be changed,
        user attributes may be or may not be copied, data may be sub-sampled,
        stats may be collected, etc. Arguments unknown to nodes are simply
        ignored. Check the documentation for copying operations of nodes to
        see which options they support.

        """

        self._g_check_open()

        # `dstgroup` is used instead of its path to avoid accepting
        # `Node` objects when `createparents` is true.  Also, note that
        # there is no risk of creating parent nodes and failing later
        # because of destination nodes already existing.
        dstparent = self._v_file._get_or_create_path(dstgroup, createparents)
        self._g_check_group(dstparent)  # Is it a group?

        if not overwrite:
            # Abort as early as possible when destination nodes exist
            # and overwriting is not enabled.
            for childname in self._v_children:
                if childname in dstparent:
                    raise NodeError(
                        "destination group ``%s`` already has "
                        "a node named ``%s``; "
                        "you may want to use the ``overwrite`` argument"
                        % (dstparent._v_pathname, childname))

        use_hardlinks = kwargs.get('use_hardlinks', False)
        if use_hardlinks:
            address_map = kwargs.setdefault('address_map', {})

            for child in self._v_children.values():
                addr, rc = child._get_obj_info()
                if rc > 1 and addr in address_map:
                    where, name = address_map[addr][0]
                    localsrc = os.path.join(where, name)
                    dstparent._v_file.create_hard_link(dstparent, child.name,
                                                       localsrc)
                    address_map[addr].append(
                        (dstparent._v_pathname, child.name)
                    )

                    # Update statistics if needed.
                    stats = kwargs.pop('stats', None)
                    if stats is not None:
                        stats['hardlinks'] += 1
                else:
                    child._f_copy(dstparent, None, overwrite, recursive,
                                  **kwargs)
                    if rc > 1:
                        address_map[addr] = [
                            (dstparent._v_pathname, child.name)
                        ]
        else:
            for child in self._v_children.values():
                child._f_copy(dstparent, None, overwrite, recursive, **kwargs)

    def __str__(self):
        """Return a short string representation of the group.

        Examples
        --------

        ::

            >>> import tables
            >>> f = tables.open_file('tables/tests/Tables_lzo2.h5')
            >>> print(f.root.group0)
            /group0 (Group) ''
            >>> f.close()

        """

        return (f"{self._v_pathname} ({self.__class__.__name__}) "
                f"{self._v_title!r}")

    def __repr__(self):
        """Return a detailed string representation of the group.

        Examples
        --------

        ::

            >>> import tables
            >>> f = tables.open_file('tables/tests/Tables_lzo2.h5')
            >>> f.root.group0
            /group0 (Group) ''
              children := ['group1' (Group), 'tuple1' (Table)]
            >>> f.close()

        """

        rep = [
            f'{childname!r} ({child.__class__.__name__})'
            for (childname, child) in self._v_children.items()
        ]
        return f'{self!s}\n  children := [{", ".join(rep)}]'


# Special definition for group root
class RootGroup(Group):

    def __init__(self, ptfile, name, title, new, filters):
        mydict = self.__dict__

        # Set group attributes.
        self._v_version = obversion
        self._v_new = new
        if new:
            self._v_new_title = title
            self._v_new_filters = filters
        else:
            self._v_new_title = None
            self._v_new_filters = None

        # Set node attributes.
        self._v_file = ptfile
        self._v_isopen = True  # root is always open
        self._v_pathname = '/'
        self._v_name = '/'
        self._v_depth = 0
        self._v_max_group_width = ptfile.params['MAX_GROUP_WIDTH']
        self._v__deleting = False
        self._v_objectid = None  # later

        # Only the root node has the file as a parent.
        # Bypass __setattr__ to avoid the ``Node._v_parent`` property.
        mydict['_v_parent'] = ptfile
        ptfile._node_manager.register_node(self, '/')

        # hdf5extension operations (do before setting an AttributeSet):
        #   Update node attributes.
        self._g_new(ptfile, name, init=True)
        #   Open the node and get its object ID.
        self._v_objectid = self._g_open()

        # Set disk attributes and read children names.
        #
        # This *must* be postponed because this method needs the root node
        # to be created and bound to ``File.root``.
        # This is an exception to the rule, handled by ``File.__init()__``.
        #
        # self._g_post_init_hook()

    def _g_load_child(self, childname):
        """Load a child node from disk.

        The child node `childname` is loaded from disk and an adequate
        `Node` object is created and returned.  If there is no such
        child, a `NoSuchNodeError` is raised.

        """

        if self._v_file.root_uep != "/":
            childname = join_path(self._v_file.root_uep, childname)
        # Is the node a group or a leaf?
        node_type = self._g_check_has_child(childname)

        # Nodes that HDF5 report as H5G_UNKNOWN
        if node_type == 'Unknown':
            return Unknown(self, childname)

        # Guess the PyTables class suited to the node,
        # build a PyTables node and return it.
        if node_type == "Group":
            if self._v_file.params['PYTABLES_SYS_ATTRS']:
                ChildClass = self._g_get_child_group_class(childname)
            else:
                # Default is a Group class
                ChildClass = Group
            return ChildClass(self, childname, new=False)
        elif node_type == "Leaf":
            ChildClass = self._g_get_child_leaf_class(childname, warn=True)
            # Building a leaf may still fail because of unsupported types
            # and other causes.
            # return ChildClass(self, childname)  # uncomment for debugging
            try:
                return ChildClass(self, childname)
            except Exception as exc:  # XXX
                warnings.warn(
                    "problems loading leaf ``%s``::\n\n"
                    "  %s\n\n"
                    "The leaf will become an ``UnImplemented`` node."
                    % (self._g_join(childname), exc))
                # If not, associate an UnImplemented object to it
                return UnImplemented(self, childname)
        elif node_type == "SoftLink":
            return SoftLink(self, childname)
        elif node_type == "ExternalLink":
            return ExternalLink(self, childname)
        else:
            return UnImplemented(self, childname)

    def _f_rename(self, newname):
        raise NodeError("the root node can not be renamed")

    def _f_move(self, newparent=None, newname=None, createparents=False):
        raise NodeError("the root node can not be moved")

    def _f_remove(self, recursive=False):
        raise NodeError("the root node can not be removed")


class TransactionGroupG(NotLoggedMixin, Group):
    _c_classid = 'TRANSGROUP'

    def _g_width_warning(self):
        warnings.warn("""\
the number of transactions is exceeding the recommended maximum (%d);\
be ready to see PyTables asking for *lots* of memory and possibly slow I/O"""
                      % (self._v_max_group_width,), PerformanceWarning)


class TransactionG(NotLoggedMixin, Group):
    _c_classid = 'TRANSG'

    def _g_width_warning(self):
        warnings.warn("""\
transaction ``%s`` is exceeding the recommended maximum number of marks (%d);\
be ready to see PyTables asking for *lots* of memory and possibly slow I/O"""
                      % (self._v_pathname, self._v_max_group_width),
                      PerformanceWarning)


class MarkG(NotLoggedMixin, Group):
    # Class identifier.
    _c_classid = 'MARKG'

    import re
    _c_shadow_name_re = re.compile(r'^a[0-9]+$')

    def _g_width_warning(self):
        warnings.warn("""\
mark ``%s`` is exceeding the recommended maximum action storage (%d nodes);\
be ready to see PyTables asking for *lots* of memory and possibly slow I/O"""
                      % (self._v_pathname, self._v_max_group_width),
                      PerformanceWarning)

    def _g_reset(self):
        """Empty action storage (nodes and attributes).

        This method empties all action storage kept in this node: nodes
        and attributes.

        """

        # Remove action storage nodes.
        for child in list(self._v_children.values()):
            child._g_remove(True, True)

        # Remove action storage attributes.
        attrs = self._v_attrs
        shname = self._c_shadow_name_re
        for attrname in attrs._v_attrnamesuser[:]:
            if shname.match(attrname):
                attrs._g__delattr(attrname)
