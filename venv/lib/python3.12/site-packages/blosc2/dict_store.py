#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################
import os
import shutil
import tempfile
import zipfile
from collections.abc import Iterator, Set
from typing import Any

import numpy as np

import blosc2
from blosc2.c2array import C2Array
from blosc2.embed_store import EmbedStore
from blosc2.schunk import SChunk


class DictStore:
    """
    Directory-based storage for compressed data using Blosc2.
    Manages arrays in a directory (.b2d) or zip (.b2z) format.

    Supports the following types:

    - blosc2.NDArray: n-dimensional arrays. When persisted externally they
      are stored as .b2nd files.
    - blosc2.SChunk: super-chunks. When persisted externally they are stored
      as .b2f files.
    - blosc2.C2Array: columnar containers. These are always kept inside the
      embedded store (never externalized).
    - numpy.ndarray: converted to blosc2.NDArray on assignment.

    Parameters
    ----------
    localpath : str
        Local path for the directory (".b2d") or file (".b2z"); other extensions
        are not supported. If a directory is specified, it will be treated as
        a Blosc2 directory format (B2DIR). If a file is specified, it
        will be treated as a Blosc2 zip format (B2ZIP).
    mode : str, optional
        File mode ('r', 'w', 'a'). Default is 'a'.
    tmpdir : str or None, optional
        Temporary directory to use when working with ".b2z" files. If None,
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
    threshold : int or None, optional
        Threshold (in bytes of uncompressed data) under which values are kept
        in the embedded store. If None, in-memory arrays are stored in the
        embedded store and on-disk arrays are stored as separate files.
        C2Array objects will always be stored in the embedded store,
        regardless of their size.

    Examples
    --------
    >>> dstore = DictStore(localpath="my_dstore.b2z", mode="w")
    >>> dstore["/node1"] = np.array([1, 2, 3])  # goes to embed store
    >>> dstore["/node2"] = blosc2.ones(2)  # goes to embed store
    >>> arr_external = blosc2.arange(3, urlpath="ext_node3.b2nd", mode="w")
    >>> dstore["/dir1/node3"] = arr_external  # external file in dir1 (.b2nd)
    >>> schunk = blosc2.SChunk(chunksize=32)
    >>> schunk.append_data(b"abcd")
    4
    >>> dstore["/dir1/schunk1"] = schunk  # externalized as .b2f if above threshold
    >>> dstore.to_b2z()  # persist to the zip file; external files are copied in
    >>> print(sorted(dstore.keys()))
    ['/dir1/node3', '/dir1/schunk1', '/node1', '/node2']
    >>> print(dstore["/node1"][:]))
    array([1, 2, 3])

    Notes
    -----
    - The DictStore is still experimental and subject to change.
      Please report any issues you may find.
    - External persistence uses the following file extensions:
      .b2nd for NDArray and .b2f for SChunk.
    """

    def __init__(
        self,
        localpath: os.PathLike[Any] | str | bytes,
        mode: str = "a",
        tmpdir: str | None = None,
        cparams: blosc2.CParams | None = None,
        dparams: blosc2.DParams | None = None,
        storage: blosc2.Storage | None = None,
        threshold: int | None = 2**13,
    ):
        """
        See :class:`DictStore` for full documentation of parameters.
        """
        self.localpath = localpath if isinstance(localpath, (str, bytes)) else str(localpath)
        if not self.localpath.endswith((".b2z", ".b2d")):
            raise ValueError(f"localpath must have a .b2z or .b2d extension; you passed: {self.localpath}")
        if mode not in ("r", "w", "a"):
            raise ValueError("For DictStore containers, mode must be 'r', 'w', or 'a'")

        self.mode = mode
        self.threshold = threshold
        self.cparams = cparams or blosc2.CParams()
        self.dparams = dparams or blosc2.DParams()
        self.storage = storage or blosc2.Storage()

        self.offsets = {}
        self.map_tree = {}
        self._temp_dir_obj = None

        self._setup_paths_and_dirs(tmpdir)

        if self.mode == "r":
            self._init_read_mode(self.dparams)
        else:
            self._init_write_append_mode(self.cparams, self.dparams, storage)

    def _setup_paths_and_dirs(self, tmpdir: str | None):
        """Set up working directories and paths."""
        self.is_zip_store = self.localpath.endswith(".b2z")
        if self.is_zip_store:
            if tmpdir is None:
                self._temp_dir_obj = tempfile.TemporaryDirectory()
                self.working_dir = self._temp_dir_obj.name
            else:
                self.working_dir = tmpdir
                os.makedirs(tmpdir, exist_ok=True)
            self.b2z_path = self.localpath
        else:  # .b2d
            self.working_dir = self.localpath
            if self.mode in ("w", "a"):
                os.makedirs(self.working_dir, exist_ok=True)
            self.b2z_path = self.localpath[:-4] + ".b2z"

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
            schunk = blosc2.open(self.b2z_path, mode="r", offset=estore_offset, dparams=dparams)
            for filepath in self.offsets:
                if filepath.endswith((".b2nd", ".b2f")):
                    key = "/" + filepath[: -5 if filepath.endswith(".b2nd") else -4]
                    self.map_tree[key] = filepath
        else:  # .b2d
            if not os.path.isdir(self.localpath):
                raise FileNotFoundError(f"Directory {self.localpath} does not exist for reading.")
            schunk = blosc2.open(self.estore_path, mode="r", dparams=dparams)
            self._update_map_tree()

        self._estore = EmbedStore(_from_schunk=schunk)

    def _init_write_append_mode(
        self,
        cparams: blosc2.CParams | None,
        dparams: blosc2.DParams | None,
        storage: blosc2.Storage | None,
    ):
        """Initialize store in write/append mode."""
        if self.mode == "a" and os.path.exists(self.localpath):
            if self.is_zip_store:
                with zipfile.ZipFile(self.localpath, "r") as zf:
                    zf.extractall(self.working_dir)
            elif not os.path.isdir(self.working_dir):
                raise FileNotFoundError(f"Directory {self.working_dir} does not exist for reading.")

        self._estore = EmbedStore(
            urlpath=self.estore_path,
            mode=self.mode,
            cparams=cparams,
            dparams=dparams,
            storage=storage,
        )
        self._update_map_tree()

    def _update_map_tree(self):
        # Build map_tree from .b2nd and .b2f files in working dir
        for root, _, files in os.walk(self.working_dir):
            for file in files:
                filepath = os.path.join(root, file)
                if filepath.endswith((".b2nd", ".b2f")):
                    # Convert filename to key: remove extension and ensure starts with /
                    rel_path = os.path.relpath(filepath, self.working_dir)
                    # Normalize path separators to forward slashes for cross-platform consistency
                    rel_path = rel_path.replace(os.sep, "/")
                    if rel_path.endswith(".b2nd"):
                        key = rel_path[:-5]
                    elif rel_path.endswith(".b2f"):
                        key = rel_path[:-4]
                    else:
                        continue
                    if not key.startswith("/"):
                        key = "/" + key
                    self.map_tree[key] = rel_path

    @property
    def estore(self) -> EmbedStore:
        """Access the underlying EmbedStore."""
        return self._estore

    def __setitem__(self, key: str, value: np.ndarray | blosc2.NDArray | SChunk | C2Array) -> None:
        """Add a node to the DictStore."""
        if isinstance(value, np.ndarray):
            value = blosc2.asarray(value, cparams=self.cparams, dparams=self.dparams)
        # C2Array should always go to embed store; let estore handle it directly
        if isinstance(value, C2Array):
            self._estore[key] = value
            return
        exceeds_threshold = self.threshold is not None and value.nbytes >= self.threshold
        # Consider both NDArray and SChunk external files (have urlpath)
        external_file = isinstance(value, (blosc2.NDArray, SChunk)) and getattr(value, "urlpath", None)
        if exceeds_threshold or (external_file and self.threshold is None):
            # Choose extension based on type
            ext = ".b2f" if isinstance(value, SChunk) else ".b2nd"
            # Convert key to a proper file path within the tree directory
            rel_key = key.lstrip("/")
            dest_path = os.path.join(self.working_dir, rel_key + ext)

            # Ensure the parent directory exists
            parent_dir = os.path.dirname(dest_path)
            if parent_dir and not os.path.exists(parent_dir):
                os.makedirs(parent_dir, exist_ok=True)

            # Save the value to the destination path
            if not external_file:
                if hasattr(value, "save"):
                    value.save(urlpath=dest_path)
                else:
                    # An SChunk does not have a save() method
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
        else:
            if external_file:
                # Embed a copy by using cframe
                value = blosc2.from_cframe(value.to_cframe())
            self._estore[key] = value

    def __getitem__(self, key: str) -> blosc2.NDArray | SChunk | C2Array:
        """Retrieve a node from the DictStore."""
        # Check map_tree first
        if key in self.map_tree:
            filepath = self.map_tree[key]
            if filepath in self.offsets:
                offset = self.offsets[filepath]["offset"]
                return blosc2.open(self.b2z_path, mode="r", offset=offset, dparams=self.dparams)
            else:
                urlpath = os.path.join(self.working_dir, filepath)
                if os.path.exists(urlpath):
                    return blosc2.open(urlpath, mode="r" if self.mode == "r" else "a", dparams=self.dparams)
                else:
                    raise KeyError(f"File for key '{key}' not found in offsets or temporary directory.")

        # Fall back to EmbedStore
        return self._estore[key]

    def get(self, key: str, default: Any = None) -> blosc2.NDArray | SChunk | C2Array | Any:
        """Retrieve a node, or default if not found."""
        try:
            return self[key]
        except KeyError:
            return default

    def __delitem__(self, key: str) -> None:
        """Remove a node from the DictStore."""
        if key in self.map_tree:
            # Remove from map_tree and delete the external file
            filepath = self.map_tree[key]
            del self.map_tree[key]

            # Delete the physical file if it exists
            full_path = os.path.join(self.working_dir, filepath)
            if os.path.exists(full_path):
                os.remove(full_path)
        elif key in self._estore:
            del self._estore[key]
        else:
            raise KeyError(f"Key '{key}' not found")

    def __contains__(self, key: str) -> bool:
        """Check if a key exists."""
        return key in self.map_tree or key in self._estore

    def __len__(self) -> int:
        """Return number of nodes."""
        return len(self.map_tree) + len(self._estore)

    def __iter__(self) -> Iterator[str]:
        """Iterate over keys."""
        yield from self.map_tree.keys()
        for key in self._estore:
            if key not in self.map_tree:
                yield key

    def keys(self) -> Set[str]:
        """Return all keys."""
        return self.map_tree.keys() | self._estore.keys()

    def values(self) -> Iterator[blosc2.NDArray | SChunk | C2Array]:
        """Iterate over all values."""
        # Get all unique keys from both map_tree and _estore, with map_tree taking precedence
        all_keys = set(self.map_tree.keys()) | set(self._estore.keys())

        for key in all_keys:
            if key in self.map_tree:
                filepath = self.map_tree[key]
                if self.is_zip_store:
                    if filepath in self.offsets:
                        offset = self.offsets[filepath]["offset"]
                        yield blosc2.open(self.b2z_path, mode="r", offset=offset, dparams=self.dparams)
                else:
                    urlpath = os.path.join(self.working_dir, filepath)
                    yield blosc2.open(urlpath, mode="r" if self.mode == "r" else "a", dparams=self.dparams)
            elif key in self._estore:
                yield self._estore[key]

    def items(self) -> Iterator[tuple[str, blosc2.NDArray | SChunk | C2Array]]:
        """Iterate over (key, value) pairs."""
        # Get all unique keys from both map_tree and _estore, with map_tree taking precedence
        all_keys = set(self.map_tree.keys()) | set(self._estore.keys())

        for key in all_keys:
            # Check map_tree first, then fall back to _estore
            if key in self.map_tree:
                filepath = self.map_tree[key]
                if self.is_zip_store:
                    if filepath in self.offsets:
                        offset = self.offsets[filepath]["offset"]
                        yield key, blosc2.open(self.b2z_path, mode="r", offset=offset)
                else:
                    urlpath = os.path.join(self.working_dir, filepath)
                    yield key, blosc2.open(urlpath, mode="r" if self.mode == "r" else "a")
            elif key in self._estore:
                yield key, self._estore[key]

    def to_b2z(self, overwrite=False, filename=None) -> os.PathLike[Any] | str:
        """
        Serialize zip store contents to the b2z file.

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite the existing b2z file if it exists. Default is False.
        filename : str, optional
            If provided, use this filename instead of the default b2z file path.

        Returns
        -------
        filename : str
            The absolute path to the created b2z file.
        """
        if self.mode == "r":
            raise ValueError("Cannot call to_b2z() on a DictStore opened in read mode.")

        b2z_path = self.b2z_path if filename is None else filename
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
        filepaths.sort(key=lambda f: os.path.getsize(f), reverse=True)

        with zipfile.ZipFile(self.b2z_path, "w", zipfile.ZIP_STORED) as zf:
            # Write all files (except estore_path) first (sorted by size)
            for filepath in filepaths:
                arcname = os.path.relpath(filepath, self.working_dir)
                zf.write(filepath, arcname)
            # Write estore last
            if os.path.exists(self.estore_path):
                arcname = os.path.relpath(self.estore_path, self.working_dir)
                zf.write(self.estore_path, arcname)
        return os.path.abspath(self.b2z_path)

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
        # Repack estore
        # TODO: for some reason this is not working
        # if self.mode != "r":
        #     cframe = self._estore.to_cframe()
        #     with open(self._estore.urlpath, "wb") as f:
        #         f.write(cframe)

        if self.is_zip_store and self.mode in ("w", "a"):
            # Serialize to b2z file
            self.to_b2z(overwrite=True)

        # Clean up temporary directory if we created it
        if self._temp_dir_obj is not None:
            self._temp_dir_obj.cleanup()

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
