#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

import contextlib
import os
from collections.abc import Iterator, MutableMapping
from typing import TYPE_CHECKING

import numpy as np

import blosc2
from blosc2.dict_store import DictStore
from blosc2.schunk import SChunk

if TYPE_CHECKING:
    from blosc2.c2array import C2Array
    from blosc2.ndarray import NDArray


class vlmetaProxy(MutableMapping):
    """Proxy for SChunk.vlmeta to control access and slicing.

    - Ensures `vlmeta[:]` returns a dict of {name: value} using decoded values.
    - Enforces TreeStore read-only mode for set/del operations.
    - Delegates iteration and length to the underlying vlmeta object.
    """

    def __init__(self, tstore: "TreeStore", inner_vlmeta):
        self._tstore = tstore
        self._inner = inner_vlmeta

    def __setitem__(self, key, value):
        if self._tstore.mode == "r":
            raise ValueError("TreeStore is in read-only mode")

        # Ensure the vlmeta SChunk is persisted before any write operation.
        # This handles the case where vlmeta is being created lazily.
        # Use DictStore's methods directly to bypass TreeStore's vlmeta filtering
        if not DictStore.__contains__(self._tstore, self._tstore._vlmeta_key):
            DictStore.__setitem__(self._tstore, self._tstore._vlmeta_key, self._tstore._vlmeta)

        # Support bulk set via [:]
        if isinstance(key, slice):
            if key.start is None and key.stop is None:
                # Merge/update existing values instead of replacing
                for k, v in value.items():
                    self._inner[k] = v
                # Persist once after bulk update
                self._tstore._persist_vlmeta()
                return
            raise NotImplementedError("Slicing is not supported, unless [:]")

        self._inner[key] = value
        # Persist changes in the embed store snapshot
        self._tstore._persist_vlmeta()

    def __getitem__(self, key):
        # Support bulk get via [:]
        if isinstance(key, slice):
            if key.start is None and key.stop is None:
                # Build a Python dict to ensure keys are str and values decoded
                return {name: self._inner[name] for name in self._inner}
            raise NotImplementedError("Slicing is not supported, unless [:]")
        return self._inner[key]

    def __delitem__(self, key):
        if self._tstore.mode == "r":
            raise ValueError("TreeStore is in read-only mode")
        self._inner.__delitem__(key)
        # Persist changes in the embed store snapshot
        self._tstore._persist_vlmeta()

    def __iter__(self):
        return iter(self._inner)

    def __len__(self):
        return len(self._inner)


class TreeStore(DictStore):
    """
    A hierarchical tree-based storage container for Blosc2 data.

    Extends :class:`blosc2.DictStore` with strict hierarchical key validation
    and tree traversal capabilities. Keys must follow a hierarchical structure
    using '/' as separator and always start with '/'. If user passes a key
    that doesn't start with '/', it will be automatically added.

    It supports the same arguments as :class:`blosc2.DictStore`.

    Parameters
    ----------
    localpath : str
        Local path for the directory (`.b2d`) or file (`.b2z`); other extensions
        are not supported. If a directory is specified, it will be treated as
        a Blosc2 directory format (B2DIR). If a file is specified, it
        will be treated as a Blosc2 zip format (B2ZIP).
    mode : str, optional
        File mode ('r', 'w', 'a'). Default is 'a'.
    tmpdir : str or None, optional
        Temporary directory to use when working with `.b2z` files. If None,
        a system temporary directory will be managed. Default is None.
    cparams : dict or None, optional
        Compression parameters for the internal embed store.
        If None, the default Blosc2 parameters are used.
    dparams : dict or None, optional
        Decompression parameters for the internal embed store.
        If None, the default Blosc2 parameters are used.
    storage : blosc2.Storage or None, optional
        Storage properties for the internal embed store.
        If None, the default Blosc2 storage properties are used.
    threshold : int, optional
        Threshold for the array size (bytes) to be kept in the embed store.
        If the *compressed* array size is below this threshold, it will be
        stored in the embed store instead of as a separate file. If None,
        in-memory arrays are stored in the embed store and on-disk arrays
        are stored as separate files.
        C2Array objects will always be stored in the embed store,
        regardless of their size.

    Examples
    --------
    >>> tstore = TreeStore(localpath="my_tstore.b2z", mode="w")
    >>> # Create a hierarchy. Data is stored in leaf nodes.
    >>> # Structural nodes like /child0 and /child0/child1 are created automatically.
    >>> tstore["/child0/leaf1"] = np.array([1, 2, 3])
    >>> tstore["/child0/child1/leaf2"] = np.array([4, 5, 6])
    >>> tstore["/child0/child2"] = np.array([7, 8, 9])
    >>>
    >>> # Walk the tree structure
    >>> for path, children, nodes in tstore.walk("/child0"):
    ...     print(f"Path: {path}, Children: {sorted(children)}, Nodes: {sorted(nodes)}")
    Path: /child0, Children: ['/child0/child1'], Nodes: ['/child0/child2', '/child0/leaf1']
    Path: /child0/child1, Children: [], Nodes: ['/child0/child1/leaf2']
    >>>
    >>> # Get a subtree view
    >>> subtree = tstore.get_subtree("/child0")
    >>> sorted(list(subtree.keys()))
    ['/child1/leaf2', '/child2', '/leaf1']

    Notes
    -----
    The TreeStore is still experimental and subject to change.
    Please report any issues you may find.
    """

    # For some reason, we had to revert the explicit parametrisation of the
    # constructor to make benchmarks wrok fine again.
    def __init__(self, *args, _from_parent_store=None, **kwargs):
        """Initialize TreeStore with subtree support.

        It supports the same arguments as :class:`blosc2.DictStore`.
        """
        if _from_parent_store is not None:
            # This is a subtree view, copy state from parent
            self.__dict__.update(_from_parent_store.__dict__)
        else:
            super().__init__(*args, **kwargs)
            self.subtree_path = ""  # Empty string means full tree

    def _is_vlmeta_key(self, key: str) -> bool:
        """Check if a key is a vlmeta key that should be hidden from regular access."""
        return key.endswith("/__vlmeta__")

    def _translate_key_to_full(self, key: str) -> str:
        """Translate subtree-relative key to full tree key."""
        if not self.subtree_path:
            return key
        if key == "/":
            return self.subtree_path
        else:
            return self.subtree_path + key

    def _translate_key_from_full(self, full_key: str) -> str | None:
        """Translate full tree key to subtree-relative key."""
        if not self.subtree_path:
            return full_key
        if full_key == self.subtree_path:
            return "/"
        elif full_key.startswith(self.subtree_path + "/"):
            return full_key[len(self.subtree_path) :]
        else:
            # Key is not within this subtree
            return None

    def _validate_key(self, key: str) -> str:
        """Validate and normalize hierarchical key structure.

        Parameters
        ----------
        key : str
            The key to validate and normalize.

        Returns
        -------
        normalized_key : str
            The normalized key with leading '/' added if missing.

        Raises
        ------
        ValueError
            If key doesn't follow hierarchical rules.
        """
        if not isinstance(key, str):
            raise ValueError(f"Key must be a string, got {type(key)}")

        # Auto-add leading '/' if missing
        if not key.startswith("/"):
            key = "/" + key

        if key != "/" and key.endswith("/"):
            raise ValueError(f"Key cannot end with '/' (except for root), got: {key}")

        if "//" in key:
            raise ValueError(f"Key cannot contain empty path segments '//', got: {key}")

        # Additional validation for special characters that might cause issues
        invalid_chars = ["\0", "\n", "\r", "\t"]
        for char in invalid_chars:
            if char in key:
                raise ValueError(f"Key cannot contain invalid character {char!r}, got: {key}")

        return key

    def __setitem__(self, key: str, value: blosc2.Array | SChunk) -> None:
        """Add a node with hierarchical key validation.

        Parameters
        ----------
        key : str
            Hierarchical node key.
        value : np.ndarray or blosc2.NDArray or blosc2.C2Array or blosc2.SChunk
            to store.

        Raises
        ------
        ValueError
            If key doesn't follow hierarchical structure rules, if trying to
            assign to a structural path that already has children, or if trying
            to add a child to a path that already contains data.
        """
        key = self._validate_key(key)

        # Check if this key already has children (is a structural subtree)
        children = self.get_children(key)
        if children:
            raise ValueError(
                f"Cannot assign array to structural path '{key}' that already has children: {children}"
            )

        # Check if we're trying to add a child to a path that already has data
        # Extract parent path from the key
        if key != "/":
            parent_path = "/".join(key.split("/")[:-1])
            if not parent_path:  # Handle case where parent is root
                parent_path = "/"

            full_parent_key = self._translate_key_to_full(parent_path)
            if super().__contains__(full_parent_key):
                raise ValueError(
                    f"Cannot add child '{key}' to path '{parent_path}' that already contains data"
                )

        full_key = self._translate_key_to_full(key)
        super().__setitem__(full_key, value)

    def __getitem__(self, key: str) -> "NDArray | C2Array | SChunk | TreeStore":
        """Retrieve a node or subtree view.

        If the key points to a subtree (intermediate path with children),
        returns a TreeStore view of that subtree. If the key points to
        a final node (leaf), returns the stored array or schunk.

        Parameters
        ----------
        key : str
            Hierarchical node key.

        Returns
        -------
        out : blosc2.NDArray or blosc2.C2Array or blosc2.SChunk or TreeStore
            The stored array/chunk if key is a leaf node, or a TreeStore subtree view
            if key is an intermediate path with children.

        Raises
        ------
        KeyError
            If key is not found.
        ValueError
            If key doesn't follow hierarchical structure rules.
        """
        key = self._validate_key(key)
        if self._is_vlmeta_key(key):
            raise KeyError(f"Key '{key}' not found; vlmeta keys are not directly accessible.")

        full_key = self._translate_key_to_full(key)

        # Check if this key has children (is a subtree)
        children = self.get_children(key)

        # Check if the key exists as an actual data node
        key_exists_as_data = super().__contains__(full_key)

        if children:
            # If it has children, return a subtree view
            return self.get_subtree(key)
        elif key_exists_as_data:
            # If no children but exists as data, it's a leaf node - get the actual data
            return super().__getitem__(full_key)
        else:
            # Key doesn't exist at all
            raise KeyError(f"Key '{key}' not found")

    def __delitem__(self, key: str) -> None:
        """Remove a node or subtree.

        If the key points to a subtree (intermediate path with children),
        removes all nodes in that subtree recursively. If the key points to a final
        node (leaf), removes only that node.

        Parameters
        ----------
        key : str
            Hierarchical node key.

        Raises
        ------
        KeyError
            If key is not found.
        ValueError
            If key doesn't follow hierarchical structure rules.
        """
        key = self._validate_key(key)

        if self._is_vlmeta_key(key):
            raise KeyError(f"Key '{key}' not found; vlmeta keys are not directly accessible.")

        # Check if the key exists (either as data or as a structural node with descendants)
        full_key = self._translate_key_to_full(key)
        key_exists_as_data = super().__contains__(full_key)
        descendants = self.get_descendants(key)

        if not key_exists_as_data and not descendants:
            raise KeyError(f"Key '{key}' not found")

        # Collect all keys to delete (leaf nodes only, since structural nodes don't exist as data)
        keys_to_delete = []

        # If the key itself has data, include it
        if key_exists_as_data:
            keys_to_delete.append(key)

        # Add all descendant leaf nodes (only those that actually exist as data)
        for descendant in descendants:
            full_descendant_key = self._translate_key_to_full(descendant)
            if super().__contains__(full_descendant_key):
                keys_to_delete.append(descendant)

        # Delete all data keys in the subtree
        for k in keys_to_delete:
            full_key_to_delete = self._translate_key_to_full(k)
            super().__delitem__(full_key_to_delete)

    def __contains__(self, key: str) -> bool:
        """Check if a key exists.

        Parameters
        ----------
        key : str
            Hierarchical node key.

        Returns
        -------
        exists : bool
            True if key exists, False otherwise.
        """
        try:
            key = self._validate_key(key)
            if self._is_vlmeta_key(key):
                return False
            full_key = self._translate_key_to_full(key)
            return super().__contains__(full_key)
        except ValueError:
            return False

    def keys(self):
        """Return all keys in the current subtree view."""
        if not self.subtree_path:
            all_keys = set(super().keys())
        else:
            all_keys = set()
            for full_key in super().keys():  # noqa: SIM118
                relative_key = self._translate_key_from_full(full_key)
                if relative_key is not None:
                    all_keys.add(relative_key)

        # Filter out vlmeta keys
        all_keys = {key for key in all_keys if not self._is_vlmeta_key(key)}

        # Also include structural paths (intermediate nodes that have children but no data)
        structural_keys = set()
        for key in all_keys:
            # For each leaf key, add all its parent paths
            parts = key.split("/")[1:]  # Remove empty first element from split
            current_path = ""
            for part in parts[:-1]:  # Exclude the leaf itself
                current_path = current_path + "/" + part if current_path else "/" + part
                if current_path and current_path != "/" and current_path not in all_keys:
                    structural_keys.add(current_path)

        return all_keys | structural_keys

    def __iter__(self) -> Iterator[str]:
        """Iterate over keys, excluding vlmeta keys."""
        return iter(self.keys())

    def items(self) -> Iterator[tuple[str, "NDArray | C2Array | SChunk | TreeStore"]]:
        """Return key-value pairs in the current subtree view."""
        for key in self.keys():
            yield key, self[key]

    def get_children(self, path: str) -> list[str]:
        """Get direct children of a given path.

        Parameters
        ----------
        path : str
            The parent path to get children for.

        Returns
        -------
        children : list[str]
            List of direct child paths.
        """
        path = self._validate_key(path)

        if path == "/":
            prefix = "/"
        else:
            prefix = path + "/"

        prefix_len = len(prefix)
        children_names = set()

        for key in self.keys():
            if self._is_vlmeta_key(key):
                continue  # Should be already filtered by self.keys(), but for safety
            if key.startswith(prefix):
                # e.g. key = /hierarchy/level1/data, prefix = /hierarchy/
                # rest = level1/data
                rest = key[prefix_len:]
                # child_name = level1
                child_name = rest.split("/")[0]
                children_names.add(child_name)

        if path == "/":
            return sorted(["/" + name for name in children_names])
        else:
            return sorted([path + "/" + name for name in children_names])

    def get_descendants(self, path: str) -> list[str]:
        """Get all descendants of a given path.

        Parameters
        ----------
        path : str
            The parent path to get descendants for.

        Returns
        -------
        descendants : list[str]
            List of all descendant paths.
        """
        path = self._validate_key(path)

        if path == "/":
            prefix = "/"
        else:
            prefix = path + "/"

        descendants = set()

        # Get all leaf nodes under this path
        for key in self.keys():
            if self._is_vlmeta_key(key):
                continue  # Should be already filtered by self.keys(), but for safety
            if key.startswith(prefix) and key != path:
                descendants.add(key)

        return sorted(descendants)

    def walk(self, path: str = "/", topdown: bool = True) -> Iterator[tuple[str, list[str], list[str]]]:
        """Walk the tree structure.

        Similar to os.walk(), this visits all structural nodes in the hierarchy,
        yielding information about each level. Returns relative names, not full paths.

        Parameters
        ----------
        path : str, optional
            The root path to start walking from. Default is "/".
        topdown : bool, optional
            If True (default), traverse top-down (yield parent before children).
            If False, traverse bottom-up (yield children before parent), mimicking os.walk(topdown=False).

        Yields
        ------
        path : str
            Current path being walked.
        children : list[str]
            List of child directory names (structural nodes that have descendants).
            These are just the names, not full paths.
        nodes : list[str]
            List of leaf node names (nodes that contain data).
            These are just the names, not full paths.

        Examples
        --------
        >>> for path, children, nodes in tstore.walk("/child0", topdown=True):
        ...     print(f"Path: {path}, Children: {children}, Nodes: {nodes}")
        """
        path = self._validate_key(path)

        # Get all direct children of this path
        direct_children = self.get_children(path)

        # Separate children into directories (have descendants) and leaf nodes
        children_dirs = []
        leaf_nodes = []

        for child in direct_children:
            child_descendants = self.get_descendants(child)
            if child_descendants:
                # Extract just the name from the full path
                child_name = child.split("/")[-1]
                children_dirs.append(child_name)
            else:
                # Extract just the name from the full path
                child_name = child.split("/")[-1]
                leaf_nodes.append(child_name)

        # Validate and normalize names to ensure robustness
        # 1) Enforce that returned names are simple (no '/')
        children_dirs = [
            name for name in children_dirs if isinstance(name, str) and "/" not in name and name != ""
        ]
        leaf_nodes = [
            name for name in leaf_nodes if isinstance(name, str) and "/" not in name and name != ""
        ]

        # 2) Ensure leaf nodes correspond to actual data nodes in the underlying store
        valid_leaf_nodes: list[str] = []
        for name in leaf_nodes:
            # Compose subtree-relative child path
            child_rel_path = path + "/" + name if path != "/" else "/" + name
            # Translate to full key in the backing store and verify it's a data node
            full_key = self._translate_key_to_full(child_rel_path)
            if super().__contains__(full_key):
                valid_leaf_nodes.append(name)
        leaf_nodes = valid_leaf_nodes

        if topdown:
            # Yield current level first (pre-order)
            yield path, children_dirs, leaf_nodes

        # Recursively walk child directories (structural nodes)
        for child in direct_children:
            child_descendants = self.get_descendants(child)
            if child_descendants:
                yield from self.walk(child, topdown=topdown)

        if not topdown:
            # Yield current level after children (post-order)
            yield path, children_dirs, leaf_nodes

    def get_subtree(self, path: str) -> "TreeStore":
        """Create a subtree view with the specified path as root.

        Parameters
        ----------
        path : str
            The path that will become the root of the subtree view (relative to current subtree,
            will be normalized to start with '/' if missing).

        Returns
        -------
        subtree : TreeStore
            A new TreeStore instance that presents the subtree as if `path` were the root.

        Examples
        --------
        >>> tstore["/child0/child1/data"] = np.array([1, 2, 3])
        >>> tstore["/child0/child1/grandchild"] = np.array([4, 5, 6])
        >>> subtree = tstore.get_subtree("/child0/child1")
        >>> list(subtree.keys())
        ['/data', '/grandchild']
        >>> subtree["/grandchild"][:]
        array([4, 5, 6])

        Notes
        -----
        This is equivalent to `tstore[path]` when path is a structural path.
        """
        path = self._validate_key(path)
        full_path = self._translate_key_to_full(path)

        # Create a new TreeStore instance that shares the same underlying storage
        # but with a different subtree_path
        subtree = TreeStore(_from_parent_store=self)
        subtree.subtree_path = full_path

        return subtree

    @property
    def vlmeta(self) -> MutableMapping:
        """Access variable-length metadata for the TreeStore or current subtree.

        Returns a proxy to the vlmeta attribute of an internal SChunk stored at
        '/__vlmeta__' for the root tree, or '<subtree_path>/__vlmeta__' for subtrees.
        The SChunk is created on-demand if it doesn't exist.

        Notes
        -----
        The metadata is stored as vlmeta of an internal SChunk, ensuring robust
        serialization and persistence. This mirrors SChunk.vlmeta behavior, with
        additional guarantees:
        - Bulk get via `[:]` always returns a dict with string keys and decoded values.
        - Read-only protection is enforced at the TreeStore level.
        - Each subtree has its own independent vlmeta storage.
        """
        # Create vlmeta key based on subtree_path
        if not self.subtree_path:
            # Root tree uses global vlmeta
            vlmeta_key = "/__vlmeta__"
        else:
            # Subtree uses path-specific vlmeta: <subtree_path>/__vlmeta__
            vlmeta_key = f"{self.subtree_path}/__vlmeta__"

        # Use super().__contains__ to bypass our own filtering logic
        if super().__contains__(vlmeta_key):
            # Load the current snapshot from the store to ensure freshness
            self._vlmeta = super().__getitem__(vlmeta_key)
        else:
            # Create a new, empty SChunk in memory. It will be persisted on first write.
            self._vlmeta = blosc2.SChunk()

        # Store the key for _persist_vlmeta method
        self._vlmeta_key = vlmeta_key

        # Return a fresh proxy that wraps the latest inner vlmeta
        return vlmetaProxy(self, self._vlmeta.vlmeta)

    def _persist_vlmeta(self) -> None:
        """Persist current vlmeta SChunk into the store.

        This is needed because the EmbedStore keeps a serialized snapshot of
        stored objects; mutating the in-memory SChunk does not automatically
        update the snapshot. We emulate an update by deleting and re-adding
        the object in the embed store.
        """
        if hasattr(self, "_vlmeta_key"):
            vlmeta_key = self._vlmeta_key
            # Only embedded case is expected; handle it safely.
            if hasattr(self, "_estore") and vlmeta_key in self._estore:
                # Replace the stored snapshot
                with contextlib.suppress(KeyError):
                    del self._estore[vlmeta_key]
                self._estore[vlmeta_key] = self._vlmeta


if __name__ == "__main__":
    # Example usage
    localpath = "example_tstore.b2z"

    with TreeStore(localpath, mode="w") as tstore:
        # Create a hierarchical structure.
        # Note: data is stored in leaf nodes, not structural nodes.
        tstore["/child0/data_node"] = np.array([1, 2, 3])
        tstore["/child0/child1/data_node"] = np.array([4, 5, 6])
        tstore["/child0/child2"] = np.array([7, 8, 9])
        tstore["/child0/child1/grandchild"] = np.array([10, 11, 12])
        tstore["/other"] = np.array([13, 14, 15])

        print("TreeStore keys:", sorted(tstore.keys()))

        # Test subtree view
        root_subtree = tstore["/child0"]
        root_subtree.vlmeta["foo"] = "bar"
        print("Subtree keys:", sorted(root_subtree.keys()))
        print("Subtree vlmeta:", root_subtree.vlmeta)

        # Walk the tree
        for path, children, nodes in root_subtree.walk("/"):
            print(f"Path: {path}, Children: {children}, Nodes: {nodes}")

    # Clean up
    if os.path.exists(localpath):
        os.remove(localpath)
