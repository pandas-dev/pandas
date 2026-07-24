#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import contextlib
import os
from collections.abc import Iterator, MutableMapping
from typing import TYPE_CHECKING

import numpy as np

import blosc2
from blosc2.dict_store import DictStore

if TYPE_CHECKING:
    from blosc2.c2array import C2Array
    from blosc2.ndarray import NDArray
    from blosc2.schunk import SChunk


class vlmetaProxy(MutableMapping):
    """Proxy for SChunk.vlmeta to control access and slicing.

    - Ensures `vlmeta[:]` returns a dict of {name: value} using decoded values.
    - Enforces TreeStore read-only mode for set/del operations.
    - Delegates iteration and length to the underlying vlmeta object.
    """

    def __init__(self, tstore: TreeStore, inner_vlmeta):
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
        Local path for the directory-backed store or compact zip-backed file.
        A ``.b2z`` suffix selects the zip-backed format. Existing directories,
        and new paths not ending in ``.b2z``, use Blosc2 directory format
        (B2DIR); a ``.b2d`` suffix is recommended for these directory-backed
        stores. Existing files are treated as Blosc2 zip format (B2ZIP).
    mode : str, optional
        File mode ('r', 'w', 'a'). Default is 'a'.
    tmpdir : str or None, optional
        Temporary directory to use when working with `.b2z` files. If None,
        a temporary directory is created in the same directory as the `.b2z`
        file, so that unpacked data stays on the same filesystem. Default is None.
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
        Default is 0, meaning values are persisted as external files by
        default.  C2Array objects are always stored in the embed store
        regardless of this setting.

    Examples
    --------
    Store plain arrays in a hierarchy:

    >>> tstore = TreeStore(localpath="my_tstore.b2z", mode="w")
    >>> # Data lives in leaf nodes; structural nodes are created automatically.
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

    Mix NDArrays and CTables in the same bundle:

    >>> import dataclasses
    >>> @dataclasses.dataclass
    ... class Row:
    ...     x: int = 0
    ...     y: float = 0.0
    >>> table = blosc2.CTable(Row)
    >>> _ = table.append(Row(x=1, y=1.5))
    >>> _ = table.append(Row(x=2, y=3.0))
    >>> with blosc2.TreeStore("bundle.b2z", mode="w") as ts:
    ...     ts["/data/array"] = blosc2.arange(5)
    ...     ts["/data/table"] = table
    >>> with blosc2.open("bundle.b2z", mode="r") as ts:
    ...     print(sorted(ts.keys()))
    ...     arr = ts["/data/array"]
    ...     tbl = ts["/data/table"]
    ...     print(type(tbl).__name__, len(tbl))
    ['/data', '/data/array', '/data/table']
    CTable 2

    """

    # For some reason, we had to revert the explicit parametrisation of the
    # constructor to make benchmarks working again.
    def __init__(self, *args, _from_parent_store=None, **kwargs):
        """Initialize TreeStore with subtree support.

        It supports the same arguments as :class:`blosc2.DictStore`.
        """
        if _from_parent_store is not None:
            # This is a subtree view, copy state from parent.
            # Mark it as closed so DictStore.__del__ does not attempt to pack
            # or clean up the shared backing store when this ephemeral view
            # is garbage-collected.
            self.__dict__.update(_from_parent_store.__dict__)
            self._closed = True
        else:
            # Call initialization and mark this storage as a b2tree object
            super().__init__(*args, **kwargs, _storage_meta={"b2tree": {"version": 1}})

            self.subtree_path = ""  # Empty string means full tree
            self._inline_handles: list = []  # inline object handles opened from this store
            self._known_object_roots_cache: set[str] | None = None
            self._effective_object_roots_cache: tuple[str, set[str]] | None = None

    # ------------------------------------------------------------------
    # Object registry helpers
    # ------------------------------------------------------------------

    def _objects_registry(self) -> dict:
        """Return the object registry dict stored in the embed-store SChunk vlmeta."""
        try:
            reg = self._estore._store.vlmeta.get("_object_registry")
            return dict(reg) if reg else {}
        except Exception:
            return {}

    def _invalidate_object_roots_cache(self) -> None:
        """Invalidate cached object-root views."""
        self._known_object_roots_cache = None
        self._effective_object_roots_cache = None

    def _register_object(self, full_key: str, *, kind: str, version: int, layout: str) -> None:
        """Register *full_key* as an object root in the persistent registry."""
        try:
            reg = self._objects_registry()
            reg[full_key] = {"kind": kind, "version": version, "layout": layout}
            self._estore._store.vlmeta["_object_registry"] = reg
            self._invalidate_object_roots_cache()
        except Exception:
            pass  # best-effort

    def _unregister_object(self, full_key: str) -> None:
        """Remove *full_key* from the object registry."""
        try:
            reg = self._objects_registry()
            reg.pop(full_key, None)
            self._estore._store.vlmeta["_object_registry"] = reg
            self._invalidate_object_roots_cache()
        except Exception:
            pass

    def _object_info(self, full_key: str) -> dict | None:
        """Return registry metadata for *full_key*, or ``None`` if not registered."""
        return self._objects_registry().get(full_key)

    def _object_roots(self) -> set:
        """Return all registered object-root full keys."""
        return set(self._objects_registry().keys())

    def _probed_object_roots(self) -> set:
        """Return object roots discovered from physical CTable manifests."""
        roots = set()
        candidates = set(self.map_tree.keys()) | set(self._estore.keys())
        for key in candidates:
            if not key.endswith("/_meta"):
                continue
            root = key[: -len("/_meta")]
            if not root:
                # A manifest at /_meta marks the TreeStore itself as a CTable
                # backing store; it is not an inline object root to collapse.
                continue
            if self._probe_object_info(root) is not None:
                roots.add(root)
        return roots

    def _known_object_roots(self) -> set:
        """Return registered plus physically probed object-root full keys."""
        if self._known_object_roots_cache is not None:
            return set(self._known_object_roots_cache)

        registered = self._object_roots()
        # Fast path: when registry is non-empty, avoid costly full-store probing.
        if registered:
            roots = registered
        else:
            roots = self._probed_object_roots()
        self._known_object_roots_cache = set(roots)
        return set(roots)

    def _effective_object_roots(self) -> set:
        """Object root keys relative to the current view (subtree or root)."""
        current_subtree_path = self.subtree_path or ""
        if (
            self._effective_object_roots_cache is not None
            and self._effective_object_roots_cache[0] == current_subtree_path
        ):
            return set(self._effective_object_roots_cache[1])

        all_roots = self._known_object_roots()
        if not self.subtree_path:
            self._effective_object_roots_cache = (current_subtree_path, set(all_roots))
            return set(all_roots)
        result = set()
        for full_key in all_roots:
            relative = self._translate_key_from_full(full_key)
            if relative is not None:
                result.add(relative)
        self._effective_object_roots_cache = (current_subtree_path, result)
        return set(result)

    def _is_object_internal_key(self, key: str, object_roots: set[str] | None = None) -> bool:
        """Return ``True`` when *key* (subtree-relative) is inside an object root."""
        roots = self._effective_object_roots() if object_roots is None else object_roots
        return any(key != root and key.startswith(root + "/") for root in roots)

    def _probe_object_info(self, full_key: str) -> dict | None:
        """Probe the physical store for a CTable manifest at *full_key*/_meta.

        Used as a fallback for stores written before the registry was introduced.
        """
        meta_full_key = full_key + "/_meta"
        if meta_full_key not in self.map_tree and meta_full_key not in self._estore:
            return None
        try:
            from blosc2.dict_store import DictStore

            meta_obj = DictStore.__getitem__(self, meta_full_key)
            if not isinstance(meta_obj, blosc2.SChunk):
                return None
            kind = meta_obj.vlmeta.get("kind")
            if isinstance(kind, bytes):
                kind = kind.decode()
            if kind == "ctable":
                return {"kind": "ctable", "version": 1, "layout": "inline-tree-subtree"}
        except Exception:
            pass
        return None

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

    def __setitem__(
        self, key: str, value: blosc2.Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray | blosc2.CTable
    ) -> None:
        """Add a node with hierarchical key validation.

        Parameters
        ----------
        key : str
            Hierarchical node key.
        value : np.ndarray or blosc2.NDArray or blosc2.C2Array or blosc2.SChunk or blosc2.CTable
            to store.

        Raises
        ------
        ValueError
            If key doesn't follow hierarchical structure rules, if trying to
            assign to a structural path that already has children, or if trying
            to add a child to a path that already contains data.

        Examples
        --------
        Store an NDArray and a CTable together:

        >>> import dataclasses
        >>> @dataclasses.dataclass
        ... class Row:
        ...     x: int = 0
        >>> t = blosc2.CTable(Row)
        >>> _ = t.append(Row(x=10))
        >>> with blosc2.TreeStore("store.b2z", mode="w") as ts:
        ...     ts["/arr"] = blosc2.zeros(5, dtype="i4")
        ...     ts["/table"] = t   # CTable stored inline

        Replacing an existing object root requires an explicit delete first::

            del ts["/table"]
            ts["/table"] = new_table
        """
        key = self._validate_key(key)

        # --- CTable: store as inline subtree object ---
        if isinstance(value, blosc2.CTable):
            self._set_ctable_object(key, value)
            return

        # Block writes to object internals
        if self._is_object_internal_key(key):
            raise ValueError(
                f"Cannot write to '{key}': it is an internal component of an object root. "
                f"Use the object's own API to modify it."
            )

        # Block overwriting an existing object root with a plain value
        full_key = self._translate_key_to_full(key)
        if (self._object_info(full_key) or self._probe_object_info(full_key)) is not None:
            raise ValueError(
                f"'{key}' is an object root (e.g. CTable). "
                f"Delete it first with `del ts['{key}']` before assigning a new value."
            )

        # Check if this key already has children (is a structural subtree)
        children = self.get_children(key)
        if children:
            raise ValueError(
                f"Cannot assign array to structural path '{key}' that already has children: {children}"
            )

        # Check if we're trying to add a child to a path that already has data
        if key != "/":
            parent_path = "/".join(key.split("/")[:-1])
            if not parent_path:  # Handle case where parent is root
                parent_path = "/"

            full_parent_key = self._translate_key_to_full(parent_path)
            if super().__contains__(full_parent_key):
                raise ValueError(
                    f"Cannot add child '{key}' to path '{parent_path}' that already contains data"
                )

        super().__setitem__(full_key, value)

    def _set_ctable_object(self, key: str, value: blosc2.CTable) -> None:
        """Materialise a CTable inline into this store at *key*."""
        if self.mode == "r":
            raise ValueError("TreeStore is in read-only mode")

        full_key = self._translate_key_to_full(key)

        # Raise if already exists as object root (no silent replace)
        if (self._object_info(full_key) or self._probe_object_info(full_key)) is not None:
            raise ValueError(
                f"'{key}' already exists as an object root. Delete it first with `del ts['{key}']`."
            )

        # Raise if already exists as data leaf
        if super().__contains__(full_key):
            raise ValueError(f"'{key}' already exists as a data leaf. Delete it first.")

        # Raise if key is inside an existing object root
        if self._is_object_internal_key(key):
            raise ValueError(f"Cannot assign to '{key}': it is inside an existing object root.")

        # Raise if key already has structural children
        children = self.get_children(key)
        if children:
            raise ValueError(
                f"Cannot assign CTable to '{key}': structural children already exist: {children}."
            )

        value._save_to_treestore(self, full_key)
        self._register_object(full_key, kind="ctable", version=1, layout="inline-tree-subtree")
        self._modified = True

    def __getitem__(
        self, key: str
    ) -> NDArray | C2Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray | blosc2.CTable | TreeStore:
        """Retrieve a node, object, or subtree view.

        If the key is a registered object root (e.g. CTable) returns that object.
        If the key is a structural intermediate path returns a subtree view.
        If the key is a leaf returns the stored array/schunk.

        Examples
        --------
        >>> import dataclasses
        >>> @dataclasses.dataclass
        ... class Row:
        ...     x: int = 0
        >>> t = blosc2.CTable(Row)
        >>> _ = t.append(Row(x=42))
        >>> with blosc2.TreeStore("store.b2z", mode="w") as ts:
        ...     ts["/arr"] = blosc2.zeros(3, dtype="i4")
        ...     ts["/group/val"] = blosc2.ones(2, dtype="f4")
        ...     ts["/table"] = t
        >>> with blosc2.open("store.b2z", mode="r") as ts:
        ...     arr = ts["/arr"]            # NDArray leaf
        ...     sub = ts["/group"]           # TreeStore subtree view
        ...     tbl = ts["/table"]           # CTable object
        ...     print(type(arr).__name__, type(sub).__name__, type(tbl).__name__)
        NDArray TreeStore CTable
        """
        key = self._validate_key(key)
        if self._is_vlmeta_key(key):
            raise KeyError(f"Key '{key}' not found; vlmeta keys are not directly accessible.")

        full_key = self._translate_key_to_full(key)

        # --- Object root dispatch (registry first, then probe fallback) ---
        info = self._object_info(full_key)
        if info is None:
            info = self._probe_object_info(full_key)
        if info is not None and info["kind"] == "ctable":
            ctable = blosc2.CTable._open_from_treestore(self, full_key)
            self._inline_handles.append(ctable)
            return ctable

        # Check if the key exists as an actual data node
        key_exists_as_data = super().__contains__(full_key)

        # Check if this key has children (is a structural subtree)
        children = self.get_children(key)

        if children:
            return self.get_subtree(key)
        elif key_exists_as_data:
            return super().__getitem__(full_key)
        else:
            raise KeyError(f"Key '{key}' not found")

    def __delitem__(self, key: str) -> None:
        """Remove a node, object root, or subtree.

        If *key* is a registered object root, all its physical leaves and the
        registry entry are removed.  If *key* has children, all descendants are
        removed recursively.  Object internals cannot be deleted directly.

        Examples
        --------
        >>> import dataclasses
        >>> @dataclasses.dataclass
        ... class Row:
        ...     x: int = 0
        >>> t = blosc2.CTable(Row)
        >>> _ = t.append(Row(x=1))
        >>> with blosc2.TreeStore("store.b2z", mode="w") as ts:
        ...     ts["/arr"] = blosc2.zeros(3, dtype="i4")
        ...     ts["/table"] = t
        ...     del ts["/table"]          # removes all CTable leaves + registry entry
        ...     print("/table" in ts)
        False
        """
        key = self._validate_key(key)

        if self._is_vlmeta_key(key):
            raise KeyError(f"Key '{key}' not found; vlmeta keys are not directly accessible.")

        full_key = self._translate_key_to_full(key)

        # --- Object root deletion ---
        if (self._object_info(full_key) or self._probe_object_info(full_key)) is not None:
            self._delete_object_subtree(full_key)
            return

        # Block direct deletion of object internals
        if self._is_object_internal_key(key):
            raise ValueError(
                f"Cannot delete '{key}': it is an internal component of an object root. "
                f"Delete the object root itself."
            )

        # Regular node / subtree deletion
        key_exists_as_data = super().__contains__(full_key)
        descendants = self.get_descendants(key)
        prefix = full_key + "/" if full_key != "/" else "/"
        object_roots_to_delete = sorted(
            [root for root in self._known_object_roots() if root.startswith(prefix)],
            key=len,
            reverse=True,
        )

        if not key_exists_as_data and not descendants and not object_roots_to_delete:
            raise KeyError(f"Key '{key}' not found")

        keys_to_delete = []
        if key_exists_as_data:
            keys_to_delete.append(key)
        for descendant in descendants:
            full_desc = self._translate_key_to_full(descendant)
            if super().__contains__(full_desc):
                keys_to_delete.append(descendant)

        for object_root in object_roots_to_delete:
            self._delete_object_subtree(object_root)

        for k in keys_to_delete:
            full_desc = self._translate_key_to_full(k)
            if super().__contains__(full_desc):
                super().__delitem__(full_desc)

        # Remove stale registry entries for any nested objects that were deleted as plain descendants.
        for root in list(self._object_roots()):
            if root.startswith(prefix):
                self._unregister_object(root)

    def _delete_object_subtree(self, full_key: str) -> None:
        """Delete all physical leaves under *full_key* and unregister it."""
        prefix = full_key + "/"
        # Remove from map_tree
        for k in list(self.map_tree.keys()):
            if k == full_key or k.startswith(prefix):
                filepath = self.map_tree.pop(k)
                full_path = os.path.join(self.working_dir, filepath)
                if os.path.exists(full_path) and not os.path.isdir(full_path):
                    os.remove(full_path)
        # Remove any embedded entries
        for k in list(self._estore.keys()):
            if k == full_key or k.startswith(prefix):
                import contextlib

                with contextlib.suppress(KeyError):
                    del self._estore[k]
        # Remove leftover directory (e.g. _indexes)
        table_dir = os.path.join(self.working_dir, full_key.lstrip("/"))
        if os.path.isdir(table_dir):
            import shutil

            shutil.rmtree(table_dir, ignore_errors=True)
        self._unregister_object(full_key)
        self._modified = True

    def __contains__(self, key: str) -> bool:
        """Check if a key exists (includes object roots, excludes object internals).

        Examples
        --------
        >>> import dataclasses
        >>> @dataclasses.dataclass
        ... class Row:
        ...     x: int = 0
        >>> t = blosc2.CTable(Row)
        >>> _ = t.append(Row(x=7))
        >>> with blosc2.TreeStore("store.b2z", mode="w") as ts:
        ...     ts["/arr"] = blosc2.zeros(2, dtype="i4")
        ...     ts["/table"] = t
        ...     print("/table" in ts)       # object root: True
        ...     print("/table/_meta" in ts) # internal key: False
        ...     print("/arr" in ts)         # normal leaf: True
        True
        False
        True
        """
        try:
            key = self._validate_key(key)
            if self._is_vlmeta_key(key):
                return False
            object_roots = self._effective_object_roots()
            if self._is_object_internal_key(key, object_roots):
                return False
            full_key = self._translate_key_to_full(key)
            return (
                super().__contains__(full_key)
                or self._object_info(full_key) is not None
                or self._probe_object_info(full_key) is not None
            )
        except ValueError:
            return False

    def keys(self):
        """Return all keys in the current subtree view.

        Object root keys (e.g. CTable) are included as single entries.
        Object-internal keys are hidden from normal traversal.

        Examples
        --------
        >>> import dataclasses
        >>> @dataclasses.dataclass
        ... class Row:
        ...     x: int = 0
        >>> t = blosc2.CTable(Row)
        >>> _ = t.append(Row(x=1))
        >>> with blosc2.TreeStore("store.b2z", mode="w") as ts:
        ...     ts["/arr"] = blosc2.zeros(3, dtype="i4")
        ...     ts["/group/val"] = blosc2.ones(2, dtype="f4")
        ...     ts["/table"] = t
        ...     print(sorted(ts.keys()))
        ['/arr', '/group', '/group/val', '/table']
        """
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

        # Filter out object-internal keys
        object_roots = self._effective_object_roots()
        all_keys = {key for key in all_keys if not self._is_object_internal_key(key, object_roots)}

        # Add object roots (they are not stored as DictStore keys themselves)
        # Build structural paths from both data leaves and object root keys
        all_with_roots = all_keys | object_roots
        structural_keys = set()
        for key in all_with_roots:
            parts = key.split("/")[1:]  # Remove empty first element from split
            current_path = ""
            for part in parts[:-1]:  # Exclude the leaf itself
                current_path = current_path + "/" + part if current_path else "/" + part
                if current_path and current_path != "/" and current_path not in all_with_roots:
                    structural_keys.add(current_path)

        return all_keys | structural_keys | object_roots

    def __iter__(self) -> Iterator[str]:
        """Iterate over keys, excluding vlmeta keys."""
        return iter(self.keys())

    def items(self) -> Iterator[tuple[str, NDArray | C2Array | SChunk | TreeStore]]:
        """Return key-value pairs in the current subtree view."""
        for key in self.keys():
            yield key, self[key]

    def values(
        self,
    ) -> Iterator[
        NDArray | C2Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray | blosc2.CTable | TreeStore
    ]:
        """Return values in the current subtree view, with object roots collapsed."""
        for key in self.keys():
            yield self[key]

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

        # 2) Ensure leaf nodes correspond to actual data nodes or object roots
        valid_leaf_nodes: list[str] = []
        for name in leaf_nodes:
            # Compose subtree-relative child path
            child_rel_path = path + "/" + name if path != "/" else "/" + name
            # Translate to full key in the backing store and verify it's a data node or object root
            full_key = self._translate_key_to_full(child_rel_path)
            if (
                super().__contains__(full_key)
                or self._object_info(full_key) is not None
                or self._probe_object_info(full_key) is not None
            ):
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

    def get_subtree(self, path: str) -> TreeStore:
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

        # Object roots cannot be navigated as subtrees
        if (self._object_info(full_path) or self._probe_object_info(full_path)) is not None:
            raise ValueError(
                f"'{path}' is an object root (e.g. CTable), not a TreeStore subtree. "
                f"Use ts['{path}'] to access the object."
            )

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
            if vlmeta_key in self.map_tree:
                filepath = self.map_tree[vlmeta_key]
                dest_path = os.path.join(self.working_dir, filepath)
                parent_dir = os.path.dirname(dest_path)
                if parent_dir and not os.path.exists(parent_dir):
                    os.makedirs(parent_dir, exist_ok=True)
                with open(dest_path, "wb") as f:
                    f.write(self._vlmeta.to_cframe())
            elif hasattr(self, "_estore") and vlmeta_key in self._estore:
                # Replace the stored snapshot
                with contextlib.suppress(KeyError):
                    del self._estore[vlmeta_key]
                self._estore[vlmeta_key] = self._vlmeta

    # ------------------------------------------------------------------
    # Lifecycle overrides (inline handle management)
    # ------------------------------------------------------------------

    def close(self) -> None:
        """Flush inline object handles then delegate to DictStore.close()."""
        if self._closed:
            return
        # Close any inline object handles (CTable etc.) before packing.
        for handle in list(getattr(self, "_inline_handles", [])):
            try:
                storage = getattr(handle, "_storage", None)
                if storage is not None:
                    handle.close()
            except Exception:
                pass
        if hasattr(self, "_inline_handles"):
            self._inline_handles.clear()
        super().close()

    def discard(self) -> None:
        """Discard without repacking; also discard inline handle storage."""
        if self._closed:
            return
        for handle in list(getattr(self, "_inline_handles", [])):
            try:
                storage = getattr(handle, "_storage", None)
                if storage is not None and hasattr(storage, "discard"):
                    storage.discard()
            except Exception:
                pass
        if hasattr(self, "_inline_handles"):
            self._inline_handles.clear()
        super().discard()


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
