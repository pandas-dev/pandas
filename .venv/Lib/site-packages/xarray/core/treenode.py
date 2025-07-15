from __future__ import annotations

import collections
import sys
from collections.abc import Iterator, Mapping
from pathlib import PurePosixPath
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    TypeVar,
)

from xarray.core.types import Self
from xarray.core.utils import Frozen, is_dict_like

if TYPE_CHECKING:
    from xarray.core.types import T_DataArray


class InvalidTreeError(Exception):
    """Raised when user attempts to create an invalid tree in some way."""


class NotFoundInTreeError(ValueError):
    """Raised when operation can't be completed because one node is not part of the expected tree."""


class NodePath(PurePosixPath):
    """Represents a path from one node to another within a tree."""

    def __init__(self, *pathsegments):
        if sys.version_info >= (3, 12):
            super().__init__(*pathsegments)
        else:
            super().__new__(PurePosixPath, *pathsegments)
        if self.drive:
            raise ValueError("NodePaths cannot have drives")

        if self.root not in ["/", ""]:
            raise ValueError(
                'Root of NodePath can only be either "/" or "", with "" meaning the path is relative.'
            )
        # TODO should we also forbid suffixes to avoid node names with dots in them?


Tree = TypeVar("Tree", bound="TreeNode")


class TreeNode(Generic[Tree]):
    """
    Base class representing a node of a tree, with methods for traversing and altering the tree.

    This class stores no data, it has only parents and children attributes, and various methods.

    Stores child nodes in an dict, ensuring that equality checks between trees
    and order of child nodes is preserved (since python 3.7).

    Nodes themselves are intrinsically unnamed (do not possess a ._name attribute), but if the node has a parent you can
    find the key it is stored under via the .name property.

    The .parent attribute is read-only: to replace the parent using public API you must set this node as the child of a
    new parent using `new_parent.children[name] = child_node`, or to instead detach from the current parent use
    `child_node.orphan()`.

    This class is intended to be subclassed by DataTree, which will overwrite some of the inherited behaviour,
    in particular to make names an inherent attribute, and allow setting parents directly. The intention is to mirror
    the class structure of xarray.Variable & xarray.DataArray, where Variable is unnamed but DataArray is (optionally)
    named.

    Also allows access to any other node in the tree via unix-like paths, including upwards referencing via '../'.

    (This class is heavily inspired by the anytree library's NodeMixin class.)

    """

    _parent: Tree | None
    _children: dict[str, Tree]

    def __init__(self, children: Mapping[str, Tree] | None = None):
        """Create a parentless node."""
        self._parent = None
        self._children = {}

        if children:
            # shallow copy to avoid modifying arguments in-place (see GH issue #9196)
            self.children = {name: child.copy() for name, child in children.items()}

    @property
    def parent(self) -> Tree | None:
        """Parent of this node."""
        return self._parent

    @parent.setter
    def parent(self: Tree, new_parent: Tree) -> None:
        raise AttributeError(
            "Cannot set parent attribute directly, you must modify the children of the other node instead using dict-like syntax"
        )

    def _set_parent(
        self, new_parent: Tree | None, child_name: str | None = None
    ) -> None:
        # TODO is it possible to refactor in a way that removes this private method?

        if new_parent is not None and not isinstance(new_parent, TreeNode):
            raise TypeError(
                "Parent nodes must be of type DataTree or None, "
                f"not type {type(new_parent)}"
            )

        old_parent = self._parent
        if new_parent is not old_parent:
            self._check_loop(new_parent)
            self._detach(old_parent)
            self._attach(new_parent, child_name)

    def _check_loop(self, new_parent: Tree | None) -> None:
        """Checks that assignment of this new parent will not create a cycle."""
        if new_parent is not None:
            if new_parent is self:
                raise InvalidTreeError(
                    f"Cannot set parent, as node {self} cannot be a parent of itself."
                )

            if self._is_descendant_of(new_parent):
                raise InvalidTreeError(
                    "Cannot set parent, as intended parent is already a descendant of this node."
                )

    def _is_descendant_of(self, node: Tree) -> bool:
        return any(n is self for n in node.parents)

    def _detach(self, parent: Tree | None) -> None:
        if parent is not None:
            self._pre_detach(parent)
            parents_children = parent.children
            parent._children = {
                name: child
                for name, child in parents_children.items()
                if child is not self
            }
            self._parent = None
            self._post_detach(parent)

    def _attach(self, parent: Tree | None, child_name: str | None = None) -> None:
        if parent is not None:
            if child_name is None:
                raise ValueError(
                    "To directly set parent, child needs a name, but child is unnamed"
                )

            self._pre_attach(parent, child_name)
            parentchildren = parent._children
            assert not any(child is self for child in parentchildren), (
                "Tree is corrupt."
            )
            parentchildren[child_name] = self
            self._parent = parent
            self._post_attach(parent, child_name)
        else:
            self._parent = None

    def orphan(self) -> None:
        """Detach this node from its parent."""
        self._set_parent(new_parent=None)

    @property
    def children(self: Tree) -> Mapping[str, Tree]:
        """Child nodes of this node, stored under a mapping via their names."""
        return Frozen(self._children)

    @children.setter
    def children(self: Tree, children: Mapping[str, Tree]) -> None:
        self._check_children(children)
        children = {**children}

        old_children = self.children
        del self.children
        try:
            self._pre_attach_children(children)
            for name, child in children.items():
                child._set_parent(new_parent=self, child_name=name)
            self._post_attach_children(children)
            assert len(self.children) == len(children)
        except Exception:
            # if something goes wrong then revert to previous children
            self.children = old_children
            raise

    @children.deleter
    def children(self) -> None:
        # TODO this just detaches all the children, it doesn't actually delete them...
        children = self.children
        self._pre_detach_children(children)
        for child in self.children.values():
            child.orphan()
        assert len(self.children) == 0
        self._post_detach_children(children)

    @staticmethod
    def _check_children(children: Mapping[str, Tree]) -> None:
        """Check children for correct types and for any duplicates."""
        if not is_dict_like(children):
            raise TypeError(
                "children must be a dict-like mapping from names to node objects"
            )

        seen = set()
        for name, child in children.items():
            if not isinstance(child, TreeNode):
                raise TypeError(
                    f"Cannot add object {name}. It is of type {type(child)}, "
                    "but can only add children of type DataTree"
                )

            childid = id(child)
            if childid not in seen:
                seen.add(childid)
            else:
                raise InvalidTreeError(
                    f"Cannot add same node {name} multiple times as different children."
                )

    def __repr__(self) -> str:
        return f"TreeNode(children={dict(self._children)})"

    def _pre_detach_children(self: Tree, children: Mapping[str, Tree]) -> None:
        """Method call before detaching `children`."""
        pass

    def _post_detach_children(self: Tree, children: Mapping[str, Tree]) -> None:
        """Method call after detaching `children`."""
        pass

    def _pre_attach_children(self: Tree, children: Mapping[str, Tree]) -> None:
        """Method call before attaching `children`."""
        pass

    def _post_attach_children(self: Tree, children: Mapping[str, Tree]) -> None:
        """Method call after attaching `children`."""
        pass

    def copy(self, *, inherit: bool = True, deep: bool = False) -> Self:
        """
        Returns a copy of this subtree.

        Copies this node and all child nodes.

        If `deep=True`, a deep copy is made of each of the component variables.
        Otherwise, a shallow copy of each of the component variable is made, so
        that the underlying memory region of the new datatree is the same as in
        the original datatree.

        Parameters
        ----------
        inherit : bool
            Whether inherited coordinates defined on parents of this node should
            also be copied onto the new tree. Only relevant if the `parent` of
            this node is not yet, and "Inherited coordinates" appear in its
            repr.
        deep : bool
            Whether each component variable is loaded into memory and copied onto
            the new object. Default is False.

        Returns
        -------
        object : DataTree
            New object with dimensions, attributes, coordinates, name, encoding,
            and data of this node and all child nodes copied from original.

        See Also
        --------
        xarray.Dataset.copy
        pandas.DataFrame.copy
        """
        return self._copy_subtree(inherit=inherit, deep=deep)

    def _copy_subtree(
        self, inherit: bool, deep: bool = False, memo: dict[int, Any] | None = None
    ) -> Self:
        """Copy entire subtree recursively."""
        new_tree = self._copy_node(inherit=inherit, deep=deep, memo=memo)
        for name, child in self.children.items():
            # TODO use `.children[name] = ...` once #9477 is implemented
            new_tree._set(
                name, child._copy_subtree(inherit=False, deep=deep, memo=memo)
            )
        return new_tree

    def _copy_node(
        self, inherit: bool, deep: bool = False, memo: dict[int, Any] | None = None
    ) -> Self:
        """Copy just one node of a tree"""
        new_empty_node = type(self)()
        return new_empty_node

    def __copy__(self) -> Self:
        return self._copy_subtree(inherit=True, deep=False)

    def __deepcopy__(self, memo: dict[int, Any] | None = None) -> Self:
        return self._copy_subtree(inherit=True, deep=True, memo=memo)

    def _iter_parents(self: Tree) -> Iterator[Tree]:
        """Iterate up the tree, starting from the current node's parent."""
        node: Tree | None = self.parent
        while node is not None:
            yield node
            node = node.parent

    def iter_lineage(self: Tree) -> tuple[Tree, ...]:
        """Iterate up the tree, starting from the current node."""
        from warnings import warn

        warn(
            "`iter_lineage` has been deprecated, and in the future will raise an error."
            "Please use `parents` from now on.",
            DeprecationWarning,
            stacklevel=2,
        )
        return (self, *self.parents)

    @property
    def lineage(self: Tree) -> tuple[Tree, ...]:
        """All parent nodes and their parent nodes, starting with the closest."""
        from warnings import warn

        warn(
            "`lineage` has been deprecated, and in the future will raise an error."
            "Please use `parents` from now on.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.iter_lineage()

    @property
    def parents(self: Tree) -> tuple[Tree, ...]:
        """All parent nodes and their parent nodes, starting with the closest."""
        return tuple(self._iter_parents())

    @property
    def ancestors(self: Tree) -> tuple[Tree, ...]:
        """All parent nodes and their parent nodes, starting with the most distant."""

        from warnings import warn

        warn(
            "`ancestors` has been deprecated, and in the future will raise an error."
            "Please use `parents`. Example: `tuple(reversed(node.parents))`",
            DeprecationWarning,
            stacklevel=2,
        )
        return (*reversed(self.parents), self)

    @property
    def root(self: Tree) -> Tree:
        """Root node of the tree"""
        node = self
        while node.parent is not None:
            node = node.parent
        return node

    @property
    def is_root(self) -> bool:
        """Whether this node is the tree root."""
        return self.parent is None

    @property
    def is_leaf(self) -> bool:
        """
        Whether this node is a leaf node.

        Leaf nodes are defined as nodes which have no children.
        """
        return self.children == {}

    @property
    def leaves(self: Tree) -> tuple[Tree, ...]:
        """
        All leaf nodes.

        Leaf nodes are defined as nodes which have no children.
        """
        return tuple(node for node in self.subtree if node.is_leaf)

    @property
    def siblings(self: Tree) -> dict[str, Tree]:
        """
        Nodes with the same parent as this node.
        """
        if self.parent:
            return {
                name: child
                for name, child in self.parent.children.items()
                if child is not self
            }
        else:
            return {}

    @property
    def subtree(self: Tree) -> Iterator[Tree]:
        """
        Iterate over all nodes in this tree, including both self and all descendants.

        Iterates breadth-first.

        See Also
        --------
        DataTree.subtree_with_keys
        DataTree.descendants
        group_subtrees
        """
        # https://en.wikipedia.org/wiki/Breadth-first_search#Pseudocode
        queue = collections.deque([self])
        while queue:
            node = queue.popleft()
            yield node
            queue.extend(node.children.values())

    @property
    def subtree_with_keys(self: Tree) -> Iterator[tuple[str, Tree]]:
        """
        Iterate over relative paths and node pairs for all nodes in this tree.

        Iterates breadth-first.

        See Also
        --------
        DataTree.subtree
        DataTree.descendants
        group_subtrees
        """
        queue = collections.deque([(NodePath(), self)])
        while queue:
            path, node = queue.popleft()
            yield str(path), node
            queue.extend((path / name, child) for name, child in node.children.items())

    @property
    def descendants(self: Tree) -> tuple[Tree, ...]:
        """
        Child nodes and all their child nodes.

        Returned in depth-first order.

        See Also
        --------
        DataTree.subtree
        """
        all_nodes = tuple(self.subtree)
        this_node, *descendants = all_nodes
        return tuple(descendants)

    @property
    def level(self: Tree) -> int:
        """
        Level of this node.

        Level means number of parent nodes above this node before reaching the root.
        The root node is at level 0.

        Returns
        -------
        level : int

        See Also
        --------
        depth
        width
        """
        return len(self.parents)

    @property
    def depth(self: Tree) -> int:
        """
        Maximum level of this tree.

        Measured from the root, which has a depth of 0.

        Returns
        -------
        depth : int

        See Also
        --------
        level
        width
        """
        return max(node.level for node in self.root.subtree)

    @property
    def width(self: Tree) -> int:
        """
        Number of nodes at this level in the tree.

        Includes number of immediate siblings, but also "cousins" in other branches and so-on.

        Returns
        -------
        depth : int

        See Also
        --------
        level
        depth
        """
        return len([node for node in self.root.subtree if node.level == self.level])

    def _pre_detach(self: Tree, parent: Tree) -> None:
        """Method call before detaching from `parent`."""
        pass

    def _post_detach(self: Tree, parent: Tree) -> None:
        """Method call after detaching from `parent`."""
        pass

    def _pre_attach(self: Tree, parent: Tree, name: str) -> None:
        """Method call before attaching to `parent`."""
        pass

    def _post_attach(self: Tree, parent: Tree, name: str) -> None:
        """Method call after attaching to `parent`."""
        pass

    def get(self: Tree, key: str, default: Tree | None = None) -> Tree | None:
        """
        Return the child node with the specified key.

        Only looks for the node within the immediate children of this node,
        not in other nodes of the tree.
        """
        if key in self.children:
            return self.children[key]
        else:
            return default

    # TODO `._walk` method to be called by both `_get_item` and `_set_item`

    def _get_item(self: Tree, path: str | NodePath) -> Tree | T_DataArray:
        """
        Returns the object lying at the given path.

        Raises a KeyError if there is no object at the given path.
        """
        if isinstance(path, str):
            path = NodePath(path)

        if path.root:
            current_node = self.root
            root, *parts = list(path.parts)
        else:
            current_node = self
            parts = list(path.parts)

        for part in parts:
            if part == "..":
                if current_node.parent is None:
                    raise KeyError(f"Could not find node at {path}")
                else:
                    current_node = current_node.parent
            elif part in ("", "."):
                pass
            elif current_node.get(part) is None:
                raise KeyError(f"Could not find node at {path}")
            else:
                current_node = current_node.get(part)
        return current_node

    def _set(self: Tree, key: str, val: Tree) -> None:
        """
        Set the child node with the specified key to value.

        Counterpart to the public .get method, and also only works on the immediate node, not other nodes in the tree.
        """
        new_children = {**self.children, key: val}
        self.children = new_children

    def _set_item(
        self: Tree,
        path: str | NodePath,
        item: Tree | T_DataArray,
        new_nodes_along_path: bool = False,
        allow_overwrite: bool = True,
    ) -> None:
        """
        Set a new item in the tree, overwriting anything already present at that path.

        The given value either forms a new node of the tree or overwrites an
        existing item at that location.

        Parameters
        ----------
        path
        item
        new_nodes_along_path : bool
            If true, then if necessary new nodes will be created along the
            given path, until the tree can reach the specified location.
        allow_overwrite : bool
            Whether or not to overwrite any existing node at the location given
            by path.

        Raises
        ------
        KeyError
            If node cannot be reached, and new_nodes_along_path=False.
            Or if a node already exists at the specified path, and allow_overwrite=False.
        """
        if isinstance(path, str):
            path = NodePath(path)

        if not path.name:
            raise ValueError("Can't set an item under a path which has no name")

        if path.root:
            # absolute path
            current_node = self.root
            root, *parts, name = path.parts
        else:
            # relative path
            current_node = self
            *parts, name = path.parts

        if parts:
            # Walk to location of new node, creating intermediate node objects as we go if necessary
            for part in parts:
                if part == "..":
                    if current_node.parent is None:
                        # We can't create a parent if `new_nodes_along_path=True` as we wouldn't know what to name it
                        raise KeyError(f"Could not reach node at path {path}")
                    else:
                        current_node = current_node.parent
                elif part in ("", "."):
                    pass
                elif part in current_node.children:
                    current_node = current_node.children[part]
                elif new_nodes_along_path:
                    # Want child classes (i.e. DataTree) to populate tree with their own types
                    new_node = type(self)()
                    current_node._set(part, new_node)
                    current_node = current_node.children[part]
                else:
                    raise KeyError(f"Could not reach node at path {path}")

        if name in current_node.children:
            # Deal with anything already existing at this location
            if allow_overwrite:
                current_node._set(name, item)
            else:
                raise KeyError(f"Already a node object at path {path}")
        else:
            current_node._set(name, item)

    def __delitem__(self: Tree, key: str) -> None:
        """Remove a child node from this tree object."""
        if key in self.children:
            child = self._children[key]
            del self._children[key]
            child.orphan()
        else:
            raise KeyError(key)

    def same_tree(self, other: Tree) -> bool:
        """True if other node is in the same tree as this node."""
        return self.root is other.root


AnyNamedNode = TypeVar("AnyNamedNode", bound="NamedNode")


def _validate_name(name: str | None) -> None:
    if name is not None:
        if not isinstance(name, str):
            raise TypeError("node name must be a string or None")
        if "/" in name:
            raise ValueError("node names cannot contain forward slashes")


class NamedNode(TreeNode, Generic[Tree]):
    """
    A TreeNode which knows its own name.

    Implements path-like relationships to other nodes in its tree.
    """

    _name: str | None
    _parent: Tree | None
    _children: dict[str, Tree]

    def __init__(self, name=None, children=None):
        super().__init__(children=children)
        _validate_name(name)
        self._name = name

    @property
    def name(self) -> str | None:
        """The name of this node."""
        return self._name

    @name.setter
    def name(self, name: str | None) -> None:
        if self.parent is not None:
            raise ValueError(
                "cannot set the name of a node which already has a parent. "
                "Consider creating a detached copy of this node via .copy() "
                "on the parent node."
            )
        _validate_name(name)
        self._name = name

    def __repr__(self, level=0):
        repr_value = "\t" * level + self.__str__() + "\n"
        for child in self.children:
            repr_value += self.get(child).__repr__(level + 1)
        return repr_value

    def __str__(self) -> str:
        name_repr = repr(self.name) if self.name is not None else ""
        return f"NamedNode({name_repr})"

    def _post_attach(self, parent: Self, name: str) -> None:
        """Ensures child has name attribute corresponding to key under which it has been stored."""
        _validate_name(name)  # is this check redundant?
        self._name = name

    def _copy_node(
        self, inherit: bool, deep: bool = False, memo: dict[int, Any] | None = None
    ) -> Self:
        """Copy just one node of a tree"""
        new_node = super()._copy_node(inherit=inherit, deep=deep, memo=memo)
        new_node._name = self.name
        return new_node

    @property
    def path(self) -> str:
        """Return the file-like path from the root to this node."""
        if self.is_root:
            return "/"
        else:
            root, *ancestors = tuple(reversed(self.parents))
            # don't include name of root because (a) root might not have a name & (b) we want path relative to root.
            names = [*(node.name for node in ancestors), self.name]
            return "/" + "/".join(names)

    def relative_to(self: NamedNode, other: NamedNode) -> str:
        """
        Compute the relative path from this node to node `other`.

        If other is not in this tree, or it's otherwise impossible, raise a ValueError.
        """
        if not self.same_tree(other):
            raise NotFoundInTreeError(
                "Cannot find relative path because nodes do not lie within the same tree"
            )

        this_path = NodePath(self.path)
        if any(other.path == parent.path for parent in (self, *self.parents)):
            return str(this_path.relative_to(other.path))
        else:
            common_ancestor = self.find_common_ancestor(other)
            path_to_common_ancestor = other._path_to_ancestor(common_ancestor)
            return str(
                path_to_common_ancestor / this_path.relative_to(common_ancestor.path)
            )

    def find_common_ancestor(self, other: NamedNode) -> NamedNode:
        """
        Find the first common ancestor of two nodes in the same tree.

        Raise ValueError if they are not in the same tree.
        """
        if self is other:
            return self

        other_paths = [op.path for op in other.parents]
        for parent in (self, *self.parents):
            if parent.path in other_paths:
                return parent

        raise NotFoundInTreeError(
            "Cannot find common ancestor because nodes do not lie within the same tree"
        )

    def _path_to_ancestor(self, ancestor: NamedNode) -> NodePath:
        """Return the relative path from this node to the given ancestor node"""

        if not self.same_tree(ancestor):
            raise NotFoundInTreeError(
                "Cannot find relative path to ancestor because nodes do not lie within the same tree"
            )
        if ancestor.path not in [a.path for a in (self, *self.parents)]:
            raise NotFoundInTreeError(
                "Cannot find relative path to ancestor because given node is not an ancestor of this node"
            )

        parents_paths = [parent.path for parent in (self, *self.parents)]
        generation_gap = list(parents_paths).index(ancestor.path)
        path_upwards = "../" * generation_gap if generation_gap > 0 else "."
        return NodePath(path_upwards)


class TreeIsomorphismError(ValueError):
    """Error raised if two tree objects do not share the same node structure."""


def group_subtrees(
    *trees: AnyNamedNode,
) -> Iterator[tuple[str, tuple[AnyNamedNode, ...]]]:
    """Iterate over subtrees grouped by relative paths in breadth-first order.

    `group_subtrees` allows for applying operations over all nodes of a
    collection of DataTree objects with nodes matched by their relative paths.

    Example usage::

        outputs = {}
        for path, (node_a, node_b) in group_subtrees(tree_a, tree_b):
            outputs[path] = f(node_a, node_b)
        tree_out = DataTree.from_dict(outputs)

    Parameters
    ----------
    *trees : Tree
        Trees to iterate over.

    Yields
    ------
    A tuple of the relative path and corresponding nodes for each subtree in the
    inputs.

    Raises
    ------
    TreeIsomorphismError
        If trees are not isomorphic, i.e., they have different structures.

    See also
    --------
    DataTree.subtree
    DataTree.subtree_with_keys
    """
    if not trees:
        raise TypeError("must pass at least one tree object")

    # https://en.wikipedia.org/wiki/Breadth-first_search#Pseudocode
    queue = collections.deque([(NodePath(), trees)])

    while queue:
        path, active_nodes = queue.popleft()

        # yield before raising an error, in case the caller chooses to exit
        # iteration early
        yield str(path), active_nodes

        first_node = active_nodes[0]
        if any(
            sibling.children.keys() != first_node.children.keys()
            for sibling in active_nodes[1:]
        ):
            path_str = "root node" if not path.parts else f"node {str(path)!r}"
            child_summary = " vs ".join(
                str(list(node.children)) for node in active_nodes
            )
            raise TreeIsomorphismError(
                f"children at {path_str} do not match: {child_summary}"
            )

        for name in first_node.children:
            child_nodes = tuple(node.children[name] for node in active_nodes)
            queue.append((path / name, child_nodes))


def zip_subtrees(
    *trees: AnyNamedNode,
) -> Iterator[tuple[AnyNamedNode, ...]]:
    """Zip together subtrees aligned by relative path."""
    for _, nodes in group_subtrees(*trees):
        yield nodes
