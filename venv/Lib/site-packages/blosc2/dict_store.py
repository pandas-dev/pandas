#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import contextlib
import os
import shutil
import tempfile
import warnings
import zipfile
from typing import TYPE_CHECKING, Any

import numpy as np

import blosc2
from blosc2.c2array import C2Array
from blosc2.embed_store import EmbedStore
from blosc2.schunk import SChunk, process_opened_object

if TYPE_CHECKING:
    from collections.abc import Iterator, Set


class DictStore:
    """
    Dictionary-like storage for compressed Blosc2 objects.

    Manages arrays in a directory or zip-file backed format.

    Supports the following types:

    - blosc2.NDArray: n-dimensional arrays. When persisted externally they
      are stored as .b2nd files.
    - blosc2.SChunk: super-chunks. When persisted externally they are stored
      as .b2f files.
    - blosc2.BatchArray: batched variable-length containers. When persisted
      externally they are stored as .b2b files.
    - blosc2.C2Array: columnar containers. These are always kept inside the
      embedded store (never externalized).
    - numpy.ndarray: converted to blosc2.NDArray on assignment.

    Parameters
    ----------
    localpath : str
        Local path for the directory or zip file. A ``.b2z`` suffix selects the
        compact zip-backed format. Existing directories, and new paths not
        ending in ``.b2z``, use Blosc2 directory format (B2DIR); a ``.b2d``
        suffix is recommended for these directory-backed stores. Existing files
        are treated as Blosc2 zip format (B2ZIP).
    mode : str, optional
        File mode ('r', 'w', 'a'). Default is 'a'.
    mmap_mode : str or None, optional
        Memory mapping mode for read access. For now, only ``"r"`` is supported,
        and only when ``mode="r"``. Default is None.
    tmpdir : str or None, optional
        Temporary directory to use when working with ".b2z" files. If None,
        a temporary directory is created in the same directory as the ".b2z"
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
    threshold : int or None, optional
        Threshold (in bytes of uncompressed data) under which values are kept
        in the embedded store. Default is 0, meaning all values are persisted
        as external files by default.  C2Array objects are always stored in
        the embedded store regardless of this setting.
    locking : bool, optional
        Serialize accesses to a directory-backed (``.b2d``) store against
        other handles and other processes.  See the Notes below.  Not
        supported for zip stores nor together with `mmap_mode`.  Default is
        False.

    Notes
    -----
    Without ``locking``, DictStore is single-process, single-writer: the key
    maps are cached in Python, so mutations made through another handle are
    not seen until the store is reopened, and concurrent writers can corrupt
    each other's entries.

    With ``locking=True`` on *every* handle (or the ``BLOSC_LOCKING``
    environment variable), a directory-backed store can be shared across
    processes: whole mutations (external files plus key maps) run under an
    exclusive lock, and every access re-syncs the key maps, so readers follow
    keys added or removed by other processes.  Two caveats: a reader holding
    a value whose key another process deletes may get errors from that value
    afterwards, and a crash mid-mutation can leave a partial external file
    behind.  Not supported on network filesystems (NFS).

    A ``.b2z`` file needs no locking: it is safe to share read-only across
    processes, and :meth:`to_b2z` replaces it atomically, so readers see
    either the old or the new archive, never a torn one.  All of this also
    applies to :class:`blosc2.TreeStore`, which builds on DictStore.

    External persistence uses the following file extensions: .b2nd for
    NDArray, .b2f for SChunk, and .b2b for BatchArray. These suffixes are a
    naming convention for newly written leaves; when reopening an existing
    store, leaf typing is resolved from object metadata instead of trusting
    the suffix alone.

    Examples
    --------
    >>> dstore = DictStore(localpath="my_dstore.b2z", mode="w")
    >>> dstore["/node1"] = np.array([1, 2, 3])
    >>> dstore["/node2"] = blosc2.ones(2)
    >>> arr_external = blosc2.arange(3, urlpath="ext_node3.b2nd", mode="w")
    >>> dstore["/dir1/node3"] = arr_external  # external file in dir1 (.b2nd)
    >>> schunk = blosc2.SChunk(chunksize=32)
    >>> schunk.append_data(b"abcd")
    1
    >>> dstore["/dir1/schunk1"] = schunk  # externalized as .b2f if above threshold
    >>> _ = dstore.to_b2z(filename="my_dstore.b2z")  # persist to the zip file; external files are copied in
    >>> print(sorted(dstore.keys()))
    ['/dir1/node3', '/dir1/schunk1', '/node1', '/node2']
    >>> print(dstore["/node1"][:])
    [1 2 3]
    """

    def __init__(
        self,
        localpath: os.PathLike[Any] | str | bytes,
        mode: str = "a",
        tmpdir: str | None = None,
        cparams: blosc2.CParams | None = None,
        dparams: blosc2.DParams | None = None,
        storage: blosc2.Storage | None = None,
        threshold: int | None = 0,
        *,
        mmap_mode: str | None = None,
        locking: bool = False,
        _storage_meta: dict | None = None,
    ):
        """
        See :class:`DictStore` for full documentation of parameters.
        """
        self.localpath = localpath if isinstance(localpath, str | bytes) else str(localpath)
        if mode not in ("r", "w", "a"):
            raise ValueError("For DictStore containers, mode must be 'r', 'w', or 'a'")
        if mmap_mode not in (None, "r"):
            raise ValueError("For DictStore containers, mmap_mode must be None or 'r'")
        if mmap_mode == "r" and mode != "r":
            raise ValueError("For DictStore containers, mmap_mode='r' requires mode='r'")
        if locking and mmap_mode is not None:
            raise ValueError("locking is not supported together with mmap_mode")

        self.mode = mode
        self.mmap_mode = mmap_mode
        self.threshold = threshold
        self.cparams = cparams or blosc2.CParams()
        self.dparams = dparams or blosc2.DParams()
        self.storage = storage or blosc2.Storage()
        self._locking = bool(locking)

        if _storage_meta:
            self.storage.meta = _storage_meta
        else:
            # Mark this storage as a b2dict object
            self.storage.meta = {"b2dict": {"version": 1}}

        self.offsets = {}
        self.map_tree = {}
        self._temp_dir_obj = None
        self._closed = False
        self._modified = False

        self._setup_paths_and_dirs(tmpdir)
        if locking and self.is_zip_store:
            raise ValueError("locking is not supported for zip (.b2z) stores; share them read-only")
        env = os.environ.get("BLOSC_LOCKING", "")
        env_locking = env not in ("", "0") and mmap_mode is None
        self._shared = (self._locking or env_locking) and not self.is_zip_store
        self._store_tick = 0

        if self.mode == "r":
            self._init_read_mode(self.dparams)
        else:
            self._init_write_append_mode(self.cparams, self.dparams, storage)
        self._store_tick = self._estore._backing_schunk.change_tick

    def _setup_paths_and_dirs(self, tmpdir: str | None):
        """Set up working directories and paths."""
        localpath_exists = os.path.exists(self.localpath)
        if localpath_exists:
            self.is_zip_store = os.path.isfile(self.localpath)
        elif self.localpath.endswith(".b2z"):
            self.is_zip_store = True
        elif self.localpath.endswith(".b2d"):
            self.is_zip_store = False
        else:
            # Default extensionless new stores to directory-backed layout.
            self.is_zip_store = False
        if self.is_zip_store:
            if self.mode == "r":
                # Read mode only needs working_dir as a path namespace for relative↔absolute
                # key arithmetic.  No files are ever written, so no real directory is needed.
                self._temp_dir_obj = None
                self.working_dir = os.path.splitext(os.path.abspath(self.localpath))[0]
            elif tmpdir is None:
                b2z_parent = os.path.dirname(os.path.abspath(self.localpath))
                self._temp_dir_obj = tempfile.TemporaryDirectory(dir=b2z_parent)
                self.working_dir = self._temp_dir_obj.name
            else:
                self._temp_dir_obj = None
                self.working_dir = tmpdir
                os.makedirs(tmpdir, exist_ok=True)
            self.b2z_path = self.localpath
        else:
            self.working_dir = self.localpath
            if self.mode in ("w", "a"):
                os.makedirs(self.working_dir, exist_ok=True)
            if self.localpath.endswith(".b2d"):
                self.b2z_path = self.localpath[:-4] + ".b2z"
            else:
                self.b2z_path = self.localpath + ".b2z"

        self.estore_path = os.path.join(self.working_dir, "embed.b2e")

    def _init_read_mode(self, dparams: blosc2.DParams | None = None):
        """Initialize store in read mode."""
        if not os.path.exists(self.localpath):
            raise FileNotFoundError(f"dir/zip file {self.localpath} does not exist.")

        if self.is_zip_store:
            self.offsets = self._get_zip_offsets()
            if "embed.b2e" not in self.offsets:
                raise FileNotFoundError("Embed file embed.b2e not found in store.")
            estore_offset = self.offsets["embed.b2e"]["offset"]
            schunk = blosc2.blosc2_ext.open(
                self.b2z_path,
                mode="r",
                offset=estore_offset,
                mmap_mode=self.mmap_mode,
                dparams=dparams,
            )
            self._update_map_tree_from_offsets()
        else:  # directory-backed store
            if not os.path.isdir(self.localpath):
                raise FileNotFoundError(f"Directory {self.localpath} does not exist for reading.")
            schunk = blosc2.blosc2_ext.open(
                self.estore_path,
                mode="r",
                offset=0,
                mmap_mode=self.mmap_mode,
                dparams=dparams,
                locking=self._locking,
            )
            self._update_map_tree()

        self._estore = EmbedStore(_from_schunk=schunk)
        self.storage.meta = self._estore.storage.meta

    @staticmethod
    def _logical_key_from_relpath(rel_path: str) -> str:
        """Map an external leaf path to its logical tree key."""
        rel_path = rel_path.replace(os.sep, "/")
        key = os.path.splitext(rel_path)[0]
        if not key.startswith("/"):
            key = "/" + key
        return key

    @staticmethod
    def _expected_ext_from_kind(kind: str) -> str:
        """Return the canonical write-time suffix for a supported external leaf kind."""
        if kind == "ndarray":
            return ".b2nd"
        if kind in ("batcharray", "listarray"):
            return ".b2b"
        return ".b2f"

    @classmethod
    def _opened_external_kind(
        cls,
        opened: blosc2.NDArray | SChunk | blosc2.ObjectArray | blosc2.BatchArray | C2Array | Any,
        rel_path: str,
    ) -> str | None:
        """Return the supported external leaf kind for an already opened object."""
        meta = getattr(opened, "schunk", opened).meta
        if "b2o" in meta and isinstance(opened, blosc2.NDArray):
            # Keep b2o carriers as NDArray external leaves during discovery.
            # Rehydrating them here can recurse when a lazy recipe points back
            # into the same DictStore via dictstore_key refs.
            kind = "ndarray"
            processed_name = type(opened).__name__
        else:
            processed = process_opened_object(opened)
            processed_name = type(processed).__name__
            if isinstance(processed, blosc2.BatchArray):
                kind = "batcharray"
            elif isinstance(processed, blosc2.ObjectArray):
                kind = "vlarray"
            elif isinstance(processed, blosc2.NDArray):
                kind = "ndarray"
            elif isinstance(processed, SChunk):
                kind = "schunk"
            elif processed_name == "ListArray":
                kind = "listarray"
            else:
                warnings.warn(
                    f"Ignoring unsupported Blosc2 object at '{rel_path}' during DictStore discovery: "
                    f"{processed_name}",
                    UserWarning,
                    stacklevel=2,
                )
                return None

        expected_ext = cls._expected_ext_from_kind(kind)
        found_ext = os.path.splitext(rel_path)[1]
        if found_ext != expected_ext:
            warnings.warn(
                f"External leaf '{rel_path}' uses extension '{found_ext}' but metadata resolves to "
                f"{processed_name}; expected '{expected_ext}'.",
                UserWarning,
                stacklevel=2,
            )
        return kind

    def _probe_external_leaf_path(self, rel_path: str) -> bool:
        """Return whether a working-dir file is a supported external leaf."""
        urlpath = os.path.join(self.working_dir, rel_path)
        try:
            opened = blosc2.blosc2_ext.open(
                urlpath,
                mode="r",
                offset=0,
                mmap_mode=self.mmap_mode,
                dparams=self.dparams,
            )
        except Exception:
            return False
        return self._opened_external_kind(opened, rel_path) is not None

    def _probe_external_leaf_offset(self, filepath: str) -> bool:
        """Return whether a zip member is a supported external leaf."""
        offset = self.offsets[filepath]["offset"]
        try:
            opened = blosc2.blosc2_ext.open(
                self.b2z_path,
                mode="r",
                offset=offset,
                mmap_mode=self.mmap_mode,
                dparams=self.dparams,
            )
        except Exception:
            return False
        return self._opened_external_kind(opened, filepath) is not None

    def _init_write_append_mode(
        self,
        cparams: blosc2.CParams | None,
        dparams: blosc2.DParams | None,
        storage: blosc2.Storage | None,
    ):
        """Initialize store in write/append mode."""
        if self.mode == "a" and os.path.exists(self.localpath):
            if self.is_zip_store:
                # When using an explicit tmpdir the directory may already contain
                # stale files from a previous open that was never closed.  Clear
                # it before extracting so we always start from a clean slate.
                if self._temp_dir_obj is None:
                    shutil.rmtree(self.working_dir, ignore_errors=True)
                    os.makedirs(self.working_dir, exist_ok=True)
                with zipfile.ZipFile(self.localpath, "r") as zf:
                    zf.extractall(self.working_dir)
            elif not os.path.isdir(self.working_dir):
                raise FileNotFoundError(f"Directory {self.working_dir} does not exist for reading.")

        if self._locking:
            # The embed store's frame lock doubles as the store-wide lock, so
            # its handle must participate; supply a Storage carrying the flag
            # (plus the urlpath/mode that EmbedStore takes from it)
            storage = storage or blosc2.Storage(contiguous=True)
            storage.locking = True
            if storage.urlpath is None:
                storage.urlpath = self.estore_path
                storage.mode = self.mode
        self._estore = EmbedStore(
            urlpath=self.estore_path,
            mode=self.mode,
            cparams=cparams,
            dparams=dparams,
            storage=storage,
            meta=self.storage.meta,
        )
        self._update_map_tree()

    def _update_map_tree(self):
        """Build map_tree from supported external leaves in working dir.

        Trust canonical external leaf suffixes on the fast path.  Fall back to
        metadata probing for legacy or manually renamed leaves with unusual
        suffixes, preserving discovery warnings and compatibility.
        """
        external_exts = {".b2nd", ".b2f", ".b2b"}
        for root, _, files in os.walk(self.working_dir):
            for file in files:
                filepath = os.path.join(root, file)
                if os.path.abspath(filepath) == os.path.abspath(self.estore_path):
                    continue
                rel_path = os.path.relpath(filepath, self.working_dir).replace(os.sep, "/")
                if os.path.splitext(rel_path)[1] in external_exts or self._probe_external_leaf_path(
                    rel_path
                ):
                    self.map_tree[self._logical_key_from_relpath(rel_path)] = rel_path

    def _update_map_tree_from_offsets(self):
        """Build map_tree from supported external leaves in a zip store.

        Zip-backed stores written by DictStore/TreeStore use canonical external
        leaf suffixes.  Trusting those suffixes avoids opening every member just
        to classify it, which is especially important for compact CTable stores
        with many columns.
        """
        external_exts = {".b2nd", ".b2f", ".b2b"}
        for filepath in self.offsets:
            if filepath == "embed.b2e":
                continue
            if os.path.splitext(filepath)[1] in external_exts or self._probe_external_leaf_offset(filepath):
                self.map_tree[self._logical_key_from_relpath(filepath)] = filepath

    def _annotate_external_value(
        self,
        key: str,
        value: blosc2.NDArray | SChunk | blosc2.ObjectArray | blosc2.BatchArray | C2Array,
    ):
        """Attach DictStore origin metadata so structured msgpack can preserve member identity."""
        value._blosc2_ref = blosc2.Ref.dictstore_key(self.localpath, key)
        return value

    @property
    def estore(self) -> EmbedStore:
        """Access the underlying EmbedStore."""
        return self._estore

    @contextlib.contextmanager
    def _mutation_bracket(self):
        """Make a whole store mutation atomic against other processes.

        The embed store's exclusive frame lock covers the ensemble (external
        leaves + key maps): every locked handle on the same store takes it for
        its own mutations, and its re-entrancy makes the nested embed-store
        operations free.  A no-op for non-shared stores.
        """
        with self._estore._write_bracket():
            self._sync_store()
            yield
            self._store_tick = self._estore._backing_schunk.change_tick

    def _bump_store_tick(self) -> None:
        """Signal an external-leaf mutation through the embed store.

        Mutations that only touch external files would otherwise be invisible
        to other handles, which watch the embed store's ``change_tick``.
        Must be called inside the mutation bracket.
        """
        if not self._shared:
            return
        try:
            tick = self._estore._store.vlmeta["dstore_tick"]
        except KeyError:
            tick = 0
        self._estore._store.vlmeta["dstore_tick"] = tick + 1

    def _sync_store(self) -> None:
        """Re-scan the external leaves if another handle changed the store.

        A no-op for non-shared stores; when nothing changed, it costs one
        staleness poll of the embed store.  The re-scan itself runs under the
        exclusive store lock so that it cannot observe the half-written files
        of an in-flight transaction.
        """
        if not self._shared:
            return
        self._estore._sync_metadata()
        if self._estore._backing_schunk.change_tick == self._store_tick:
            return
        with self._estore._backing_schunk.holding_lock():
            self._estore._sync_metadata()
            self.map_tree = {}
            self._update_map_tree()
            self._store_tick = self._estore._backing_schunk.change_tick

    @staticmethod
    def _value_nbytes(value: blosc2.Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray) -> int:
        if isinstance(value, blosc2.ObjectArray | blosc2.BatchArray):
            return value.schunk.nbytes
        return value.nbytes

    @staticmethod
    def _is_external_value(value: blosc2.Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray) -> bool:
        return isinstance(value, blosc2.NDArray | SChunk | blosc2.ObjectArray | blosc2.BatchArray) and bool(
            getattr(value, "urlpath", None)
        )

    @staticmethod
    def _external_ext(value: blosc2.Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray) -> str:
        if isinstance(value, blosc2.NDArray):
            return ".b2nd"
        if isinstance(value, blosc2.BatchArray):
            return ".b2b"
        return ".b2f"

    def __setitem__(
        self, key: str, value: blosc2.Array | SChunk | blosc2.ObjectArray | blosc2.BatchArray
    ) -> None:
        """Add a node to the DictStore."""
        self._modified = True
        if isinstance(value, np.ndarray):
            value = blosc2.asarray(value, cparams=self.cparams, dparams=self.dparams)
        with self._mutation_bracket():
            # Dict-like overwrite: drop any previous value under this key, in
            # either tier.  Otherwise behavior depends on the *size* of old and
            # new values: embedded keys refused overwrite while external ones
            # accepted it, and an embed->external overwrite left a stale
            # embedded value that resurrected after a delete.
            if key in self.map_tree:
                old_filepath = self.map_tree.pop(key)
                old_full_path = os.path.join(self.working_dir, old_filepath)
                if os.path.exists(old_full_path):
                    os.remove(old_full_path)
            if key in self._estore:
                del self._estore[key]
            # C2Array should always go to embed store; let estore handle it directly
            if isinstance(value, C2Array):
                self._estore[key] = value
                return
            exceeds_threshold = self.threshold is not None and self._value_nbytes(value) >= self.threshold
            external_file = self._is_external_value(value)
            if exceeds_threshold or (external_file and self.threshold is None):
                ext = self._external_ext(value)
                # Convert key to a proper file path within the tree directory
                rel_key = key.lstrip("/")
                dest_path = os.path.join(self.working_dir, rel_key + ext)

                # Ensure the parent directory exists
                parent_dir = os.path.dirname(dest_path)
                if parent_dir and not os.path.exists(parent_dir):
                    os.makedirs(parent_dir, exist_ok=True)

                # Save the value to the destination path
                if not external_file:
                    if isinstance(value, blosc2.NDArray) and "b2o" in value.schunk.meta:
                        carrier = blosc2.empty(
                            value.shape,
                            value.dtype,
                            chunks=value.chunks,
                            blocks=value.blocks,
                            cparams=value.cparams,
                            urlpath=dest_path,
                            mode="w",
                            meta={"b2o": value.schunk.meta["b2o"]},
                        )
                        for meta_key, meta_value in value.schunk.vlmeta[:].items():
                            carrier.schunk.vlmeta[meta_key] = meta_value
                    elif hasattr(value, "save"):
                        value.save(urlpath=dest_path)
                    else:
                        # SChunk, ObjectArray and BatchArray can all be persisted via their cframe.
                        with open(dest_path, "wb") as f:
                            f.write(value.to_cframe())
                else:
                    # This should be faster than using value.save() ?
                    shutil.copy2(value.urlpath, dest_path)

                # Store relative path from tree directory
                rel_path = os.path.relpath(dest_path, self.working_dir)
                # Normalize to forward slashes
                rel_path = rel_path.replace(os.sep, "/")
                self.map_tree[key] = rel_path
                self._bump_store_tick()
            else:
                if external_file:
                    # Embed a copy by using cframe
                    value = blosc2.from_cframe(value.to_cframe())
                self._estore[key] = value

    def __getitem__(
        self, key: str
    ) -> blosc2.NDArray | SChunk | blosc2.ObjectArray | blosc2.BatchArray | C2Array:
        """Retrieve a node from the DictStore."""
        self._sync_store()
        # Check map_tree first
        if key in self.map_tree:
            filepath = self.map_tree[key]
            if filepath in self.offsets:
                offset = self.offsets[filepath]["offset"]
                opened = blosc2.blosc2_ext.open(
                    self.b2z_path,
                    mode="r",
                    offset=offset,
                    mmap_mode=self.mmap_mode,
                    dparams=self.dparams,
                )
                return self._annotate_external_value(key, process_opened_object(opened))
            else:
                urlpath = os.path.join(self.working_dir, filepath)
                if os.path.exists(urlpath):
                    return self._annotate_external_value(
                        key,
                        blosc2.open(
                            urlpath,
                            mode="r" if self.mode == "r" else "a",
                            mmap_mode=self.mmap_mode if self.mode == "r" else None,
                            dparams=self.dparams,
                        ),
                    )
                else:
                    raise KeyError(f"File for key '{key}' not found in offsets or temporary directory.")

        # Fall back to EmbedStore
        return self._estore[key]

    def get(
        self, key: str, default: Any = None
    ) -> blosc2.NDArray | SChunk | blosc2.ObjectArray | blosc2.BatchArray | C2Array | Any:
        """Retrieve a node, or default if not found."""
        try:
            return self[key]
        except KeyError:
            return default

    def __delitem__(self, key: str) -> None:
        """Remove a node from the DictStore."""
        self._modified = True
        with self._mutation_bracket():
            if key in self.map_tree:
                # Remove from map_tree and delete the external file
                filepath = self.map_tree[key]
                del self.map_tree[key]

                # Delete the physical file if it exists
                full_path = os.path.join(self.working_dir, filepath)
                if os.path.exists(full_path):
                    os.remove(full_path)
                self._bump_store_tick()
            elif key in self._estore:
                del self._estore[key]
            else:
                raise KeyError(f"Key '{key}' not found")

    def __contains__(self, key: str) -> bool:
        """Check if a key exists."""
        self._sync_store()
        return key in self.map_tree or key in self._estore

    def __len__(self) -> int:
        """Return number of nodes."""
        self._sync_store()
        return len(self.map_tree) + len(self._estore)

    def __iter__(self) -> Iterator[str]:
        """Iterate over keys."""
        self._sync_store()
        yield from self.map_tree.keys()
        for key in self._estore:
            if key not in self.map_tree:
                yield key

    def keys(self) -> Set[str]:
        """Return all keys."""
        self._sync_store()
        return self.map_tree.keys() | self._estore.keys()

    def values(self) -> Iterator[blosc2.NDArray | SChunk | C2Array]:
        """Iterate over all values."""
        self._sync_store()
        # Get all unique keys from both map_tree and _estore, with map_tree taking precedence
        all_keys = set(self.map_tree.keys()) | set(self._estore.keys())

        for key in all_keys:
            if key in self.map_tree:
                filepath = self.map_tree[key]
                if self.is_zip_store:
                    if filepath in self.offsets:
                        offset = self.offsets[filepath]["offset"]
                        yield self._annotate_external_value(
                            key,
                            process_opened_object(
                                blosc2.blosc2_ext.open(
                                    self.b2z_path,
                                    mode="r",
                                    offset=offset,
                                    mmap_mode=self.mmap_mode,
                                    dparams=self.dparams,
                                )
                            ),
                        )
                else:
                    urlpath = os.path.join(self.working_dir, filepath)
                    yield self._annotate_external_value(
                        key,
                        blosc2.open(
                            urlpath,
                            mode="r" if self.mode == "r" else "a",
                            mmap_mode=self.mmap_mode if self.mode == "r" else None,
                            dparams=self.dparams,
                        ),
                    )
            elif key in self._estore:
                yield self._estore[key]

    def items(self) -> Iterator[tuple[str, blosc2.NDArray | SChunk | C2Array]]:
        """Iterate over (key, value) pairs."""
        self._sync_store()
        # Get all unique keys from both map_tree and _estore, with map_tree taking precedence
        all_keys = set(self.map_tree.keys()) | set(self._estore.keys())

        for key in all_keys:
            # Check map_tree first, then fall back to _estore
            if key in self.map_tree:
                filepath = self.map_tree[key]
                if self.is_zip_store:
                    if filepath in self.offsets:
                        offset = self.offsets[filepath]["offset"]
                        yield (
                            key,
                            self._annotate_external_value(
                                key,
                                process_opened_object(
                                    blosc2.blosc2_ext.open(
                                        self.b2z_path,
                                        mode="r",
                                        offset=offset,
                                        mmap_mode=self.mmap_mode,
                                        dparams=self.dparams,
                                    )
                                ),
                            ),
                        )
                else:
                    urlpath = os.path.join(self.working_dir, filepath)
                    yield (
                        key,
                        self._annotate_external_value(
                            key,
                            blosc2.open(
                                urlpath,
                                mode="r" if self.mode == "r" else "a",
                                mmap_mode=self.mmap_mode if self.mode == "r" else None,
                                dparams=self.dparams,
                            ),
                        ),
                    )
            elif key in self._estore:
                yield key, self._estore[key]

    def to_b2z(self, overwrite=False, filename=None) -> os.PathLike[Any] | str:
        """
        Serialize store contents to a compact ``.b2z`` file.

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite the existing b2z file if it exists. Default is False.
        filename : str, optional
            If provided, use this filename instead of the default b2z file path.
            Keyword use is recommended for clarity.

        Returns
        -------
        filename : str
            The absolute path to the created b2z file.

        Notes
        -----
        The file is written to a temporary sibling and moved onto ``filename``
        atomically: concurrent readers of an existing target see either the
        old archive or the new one, never a partial write. On Windows, the
        final replace fails if another process holds the target open.

        Examples
        --------
        Pack a directory-backed store into a zip store.  A ``.b2d`` suffix is
        recommended for directory-backed stores, but not required::

            with blosc2.DictStore("data.b2d", mode="w") as dstore:
                dstore["/values"] = np.arange(10)

            with blosc2.DictStore("data.b2d", mode="r") as dstore:
                dstore.to_b2z(filename="data.b2z", overwrite=True)

        ``filename`` can also be passed positionally::

            with blosc2.DictStore("data.b2d", mode="r") as dstore:
                dstore.to_b2z("copy.b2z", overwrite=True)
        """
        if isinstance(overwrite, str | os.PathLike) and filename is None:
            filename = overwrite
            overwrite = False

        if self.mode == "r" and self.is_zip_store:
            raise ValueError("Cannot call to_b2z() on a .b2z DictStore opened in read mode.")

        b2z_path = self.b2z_path if filename is None else filename
        b2z_path = os.fspath(b2z_path)
        if not b2z_path.endswith(".b2z"):
            raise ValueError("b2z_path must have a .b2z extension")

        if os.path.exists(b2z_path) and not overwrite:
            raise FileExistsError(f"'{b2z_path}' already exists. Use overwrite=True to overwrite.")

        # Gather all files except estore_path
        filepaths = []
        for root, _, files in os.walk(self.working_dir):
            for file in files:
                filepath = os.path.join(root, file)
                if os.path.abspath(filepath) != os.path.abspath(self.estore_path):
                    filepaths.append(filepath)

        # Sort filepaths by file size from largest to smallest
        filepaths.sort(key=os.path.getsize, reverse=True)

        # Build the zip in a temporary file on the same filesystem, then move it
        # atomically onto the target: concurrent readers of the target see either
        # the old zip or the new one, never a torn state.  (On Windows the final
        # replace fails if another process holds the target open.)
        fd, tmp_path = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(b2z_path)), suffix=".b2z.tmp")
        os.close(fd)
        try:
            with zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_STORED) as zf:
                # Write all files (except estore_path) first (sorted by size)
                for filepath in filepaths:
                    arcname = os.path.relpath(filepath, self.working_dir)
                    zf.write(filepath, arcname)
                # Write estore last
                if os.path.exists(self.estore_path):
                    arcname = os.path.relpath(self.estore_path, self.working_dir)
                    zf.write(self.estore_path, arcname)
            os.replace(tmp_path, b2z_path)
        except BaseException:
            with contextlib.suppress(OSError):
                os.remove(tmp_path)
            raise
        return os.path.abspath(b2z_path)

    def to_b2d(self, dirname=None, *, overwrite: bool = False) -> os.PathLike[Any] | str:
        """
        Serialize store contents to a directory-backed store.

        Parameters
        ----------
        dirname : str, optional
            If provided, use this directory instead of the default directory
            path.  A ``.b2d`` suffix is recommended for clarity, but not
            required.
        overwrite : bool, optional
            If True, overwrite the existing b2d directory if it exists.
            Default is False.

        Returns
        -------
        dirname : str
            The absolute path to the created directory-backed store.

        Examples
        --------
        Unpack a zip-backed store into a directory-backed store::

            with blosc2.DictStore("data.b2z", mode="r") as dstore:
                dstore.to_b2d("data.b2d", overwrite=True)

            with blosc2.DictStore("data.b2d", mode="r") as dstore:
                values = dstore["/values"][:]

        Copy an existing directory-backed store to another directory.  A
        ``.b2d`` suffix is recommended for directory-backed stores::

            with blosc2.DictStore("data.b2d", mode="r") as dstore:
                dstore.to_b2d("backup.b2d", overwrite=True)
        """
        b2d_path = self.localpath if dirname is None and not self.is_zip_store else dirname
        if b2d_path is None:
            b2d_path = (
                self.b2z_path[:-4] + ".b2d" if self.b2z_path.endswith(".b2z") else self.b2z_path + ".b2d"
            )
        b2d_path = os.fspath(b2d_path)

        target_path = os.path.abspath(b2d_path)
        source_path = os.path.abspath(self.working_dir)
        if not self.is_zip_store and target_path == source_path:
            return target_path

        if os.path.exists(target_path):
            if not overwrite:
                raise FileExistsError(f"'{target_path}' already exists. Use overwrite=True to overwrite.")
            if os.path.isdir(target_path):
                shutil.rmtree(target_path)
            else:
                os.remove(target_path)

        if self.is_zip_store and self.mode == "r":
            os.makedirs(target_path, exist_ok=True)
            with zipfile.ZipFile(self.b2z_path, "r") as zf:
                zf.extractall(target_path)
        else:
            shutil.copytree(self.working_dir, target_path)
        return target_path

    def _get_zip_offsets(self) -> dict[str, dict[str, int]]:
        """Get offset and length of all files in the zip archive."""
        self.offsets = {}  # Reset offsets
        with open(self.b2z_path, "rb") as f, zipfile.ZipFile(f) as zf:
            for info in zf.infolist():
                # info.header_offset points to the local file header
                # The actual file data starts after the header
                f.seek(info.header_offset)
                local_header = f.read(30)
                filename_len = int.from_bytes(local_header[26:28], "little")
                extra_len = int.from_bytes(local_header[28:30], "little")
                data_offset = info.header_offset + 30 + filename_len + extra_len
                self.offsets[info.filename] = {"offset": data_offset, "length": info.file_size}
        return self.offsets

    def close(self) -> None:
        """Persist changes and cleanup."""
        if self._closed:
            return
        self._closed = True

        # Repack estore
        # TODO: for some reason this is not working
        # if self.mode != "r":
        #     cframe = self._estore.to_cframe()
        #     with open(self._estore.urlpath, "wb") as f:
        #         f.write(cframe)

        if self.is_zip_store and self.mode in ("w", "a"):
            # Serialize to b2z file.
            self.to_b2z(overwrite=True)

        # Clean up temporary directory if we created it
        if self._temp_dir_obj is not None:
            self._temp_dir_obj.cleanup()

    def discard(self) -> None:
        """Clean up resources *without* repacking the .b2z file.

        Use this instead of :meth:`close` when the store was opened only for
        inspection and should be thrown away without persisting any changes
        back to the archive.
        """
        if self._closed:
            return
        self._closed = True
        if self._temp_dir_obj is not None:
            self._temp_dir_obj.cleanup()

    def __del__(self):
        """Ensure the temporary directory is removed and, if writes were made
        through this store's own API, the store is flushed back to the .b2z
        file.

        When no Python-level writes went through ``__setitem__`` / ``__delitem__``
        (``_modified`` is False), we skip ``to_b2z()`` to avoid repacking a
        potentially partial temp dir during garbage collection.  Explicit
        ``close()`` / ``__exit__`` always repacks regardless.
        """
        try:
            if not self._closed and self.is_zip_store and self.mode in ("w", "a") and not self._modified:
                # Skip repacking — discard is safe and avoids corrupting the
                # archive when the temp dir is torn down during GC.
                self.discard()
            else:
                self.close()
        except Exception:
            pass

    def __enter__(self):
        """Context manager enter."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
        # No need to handle exceptions, just close the DictStore
        return False


if __name__ == "__main__":
    # Example usage
    localpath = "example_dstore.b2z"
    if True:
        with DictStore(localpath, mode="w") as dstore:
            dstore["/node1"] = np.array([1, 2, 3])
            dstore["/node2"] = blosc2.ones(2)

            # Make /node3 an external file
            arr_external = blosc2.arange(3, urlpath="ext_node3.b2nd", mode="w")
            dstore["/dir1/node3"] = arr_external

            print("DictStore keys:", list(dstore.keys()))
            print("Node1 data:", dstore["/node1"][:])
            print("Node2 data:", dstore["/node2"][:])
            print("Node3 data (external):", dstore["/dir1/node3"][:])

            del dstore["/node1"]
            print("After deletion, keys:", list(dstore.keys()))

    # Open the stored zip file
    with DictStore(localpath, mode="r") as dstore_opened:
        print("Opened dstore keys:", list(dstore_opened.keys()))
        for key, value in dstore_opened.items():
            if isinstance(value, blosc2.NDArray):
                print(
                    f"Key: {key}, Shape: {value.shape}, Values: {value[:10] if len(value) > 3 else value[:]}"
                )
