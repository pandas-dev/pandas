#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

import copy
from collections.abc import Iterator, KeysView
from typing import Any

import numpy as np

import blosc2
from blosc2.c2array import C2Array
from blosc2.schunk import SChunk

PROFILE = False  # Set to True to enable PROFILE prints in EmbedStore


class EmbedStore:
    """
    A dictionary-like container for storing NumPy/Blosc2 arrays (NDArray or SChunk) as nodes.

    For NumPy arrays, Blosc2 NDArrays (even if they live in external ``.b2nd`` files),
    and Blosc2 SChunk objects, the data is read and embedded into the store. For remote
    arrays (``C2Array``), only lightweight references (URL base and path) are stored.
    If you need a richer hierarchical container with optional external references, consider using
    `blosc2.TreeStore` or `blosc2.DictStore`.

    Parameters
    ----------
    urlpath : str or None, optional
        Path for persistent storage. Using a '.b2e' extension is recommended.
        If None, the embed store will be in memory only, which can be
        deserialized later using the :func:`blosc2.from_cframe` function.
    mode : str, optional
        File mode ('r', 'w', 'a'). Default is 'w'.
    cparams : dict or None, optional
        Compression parameters for nodes and the embed store itself.
        Default is None, which uses the default Blosc2 parameters.
    dparams : dict or None, optional
        Decompression parameters for nodes and the embed store itself.
        Default is None, which uses the default Blosc2 parameters.
    storage : blosc2.Storage or None, optional
        Storage properties for the embed store.  If passed, it will override
        the `urlpath` and `mode` parameters.
    chunksize : int, optional
        Size of chunks for the backing storage. Default is 1 MiB.

    Examples
    --------
    >>> estore = EmbedStore(urlpath="example_estore.b2e", mode="w")
    >>> estore["/node1"] = np.array([1, 2, 3])
    >>> estore["/node2"] = blosc2.ones(2)
    >>> estore["/node3"] = blosc2.arange(3, dtype="i4", urlpath="external_node3.b2nd", mode="w")
    >>> urlpath = blosc2.URLPath("@public/examples/ds-1d.b2nd", "https://cat2.cloud/demo")
    >>> estore["/node4"] = blosc2.open(urlpath, mode="r")
    >>> print(list(estore.keys()))
    ['/node1', '/node2', '/node3', '/node4']
    >>> print(estore["/node1"][:])
    [1 2 3]

    Notes
    -----
    The EmbedStore is still experimental and subject to change.
    Please report any issues you may find.
    """

    def __init__(
        self,
        urlpath: str | None = None,
        mode: str = "a",
        cparams: blosc2.CParams | None = None,
        dparams: blosc2.CParams | None = None,
        storage: blosc2.Storage | None = None,
        chunksize: int | None = 2**13,
        _from_schunk: SChunk | None = None,
    ):
        """Initialize EmbedStore."""

        # For some reason, the SChunk store cannot achieve the same compression ratio as the NDArray store,
        # although it is more efficient in terms of CPU usage.
        # Let's use the SChunk store by default and continue experimenting.
        self._schunk_store = True  # put this to False to use an NDArray instead of a SChunk
        self.urlpath = urlpath

        if _from_schunk is not None:
            self.cparams = _from_schunk.cparams
            self.dparams = _from_schunk.dparams
            self.mode = mode
            self._store = _from_schunk
            self._load_metadata()
            return

        self.mode = mode
        self.cparams = cparams or blosc2.CParams()
        # self.cparams.nthreads = 1  # for debugging purposes, use only one thread
        self.dparams = dparams or blosc2.DParams()
        # self.dparams.nthreads = 1  # for debugging purposes, use only one thread
        if storage is None:
            self.storage = blosc2.Storage(
                contiguous=True,
                urlpath=urlpath,
                mode=mode,
            )
        else:
            self.storage = storage

        if mode in ("r", "a") and urlpath:
            self._store = blosc2.blosc2_ext.open(urlpath, mode=mode, offset=0)
            self._load_metadata()
            return

        _cparams = copy.deepcopy(self.cparams)
        _cparams.typesize = 1  # ensure typesize is set to 1 for byte storage
        _storage = self.storage
        # Mark this storage as a b2embed object
        _storage.meta = {"b2embed": {"version": 1}}
        if self._schunk_store:
            self._store = blosc2.SChunk(
                chunksize=chunksize,
                data=None,
                cparams=_cparams,
                dparams=self.dparams,
                storage=_storage,
            )
        else:
            self._store = blosc2.zeros(
                chunksize,
                dtype=np.uint8,
                cparams=_cparams,
                dparams=self.dparams,
                storage=_storage,
            )
        self._embed_map: dict = {}
        self._current_offset = 0

    def _validate_key(self, key: str) -> None:
        """Validate node key."""
        if not isinstance(key, str):
            raise TypeError("Key must be a string.")
        if not key.startswith("/"):
            raise ValueError("Key must start with '/'.")
        if len(key) > 1 and key.endswith("/"):
            raise ValueError("Key cannot end with '/' unless it is the root key '/'.")
        if "//" in key:
            raise ValueError("Key cannot contain consecutive slashes '//'.")
        for char in (":", "\0", "\n", "\r", "\t"):
            if char in key:
                raise ValueError(f"Key cannot contain character: {char!r}")
        if key in self._embed_map:
            raise ValueError(f"Key '{key}' already exists in store.")

    def _ensure_capacity(self, needed_bytes: int) -> None:
        """Ensure backing storage has enough capacity."""
        required_size = self._current_offset + needed_bytes
        if required_size > self._store.shape[0]:
            new_size = max(required_size, int(self._store.shape[0] * 1.5))
            self._store.resize((new_size,))

    def __setitem__(self, key: str, value: blosc2.Array | SChunk) -> None:
        """Add a node to the embed store."""
        if self.mode == "r":
            raise ValueError("Cannot set items in read-only mode.")
        self._validate_key(key)
        if isinstance(value, C2Array):
            self._embed_map[key] = {"urlbase": value.urlbase, "path": value.path}
        else:
            if isinstance(value, np.ndarray):
                value = blosc2.asarray(value, cparams=self.cparams, dparams=self.dparams)
            serialized_data = value.to_cframe()
            data_len = len(serialized_data)
            if not self._schunk_store:
                self._ensure_capacity(data_len)
            offset = self._current_offset
            if self._schunk_store:
                self._store[offset : offset + data_len] = serialized_data
            else:
                self._store[offset : offset + data_len] = np.frombuffer(serialized_data, dtype=np.uint8)
            self._current_offset += data_len
            self._embed_map[key] = {"offset": offset, "length": data_len}
        self._save_metadata()

    def __getitem__(self, key: str) -> blosc2.NDArray | SChunk:
        """Retrieve a node from the embed store."""
        if key not in self._embed_map:
            raise KeyError(f"Key '{key}' not found in the embed store.")
        node_info = self._embed_map[key]
        urlbase = node_info.get("urlbase", None)
        if urlbase:
            urlpath = blosc2.URLPath(node_info["path"], urlbase=urlbase)
            return blosc2.open(urlpath, mode="r")
        offset = node_info["offset"]
        length = node_info["length"]
        serialized_data = bytes(self._store[offset : offset + length])
        # It is safer to copy data here, as the reference to the SChunk may disappear
        # Use from_cframe so we can deserialize either an NDArray or an SChunk
        return blosc2.from_cframe(serialized_data, copy=True)

    def get(self, key: str, default: Any = None) -> blosc2.NDArray | SChunk | Any:
        """Retrieve a node, or default if not found."""
        return self[key] if key in self._embed_map else default

    def __delitem__(self, key: str) -> None:
        """Remove a node from the embed store."""
        if key not in self._embed_map:
            raise KeyError(f"Key '{key}' not found in the embed store.")
        del self._embed_map[key]
        self._save_metadata()

    def __contains__(self, key: str) -> bool:
        """Check if a key exists."""
        return key in self._embed_map

    def __len__(self) -> int:
        """Return number of nodes."""
        return len(self._embed_map)

    def __iter__(self) -> Iterator[str]:
        """Iterate over keys."""
        return iter(self._embed_map)

    def keys(self) -> KeysView[str]:
        """Return all keys."""
        return self._embed_map.keys()

    def values(self) -> Iterator[blosc2.NDArray | SChunk]:
        """Iterate over all values."""
        for key in self._embed_map:
            yield self[key]

    def items(self) -> Iterator[tuple[str, blosc2.NDArray | SChunk]]:
        """Iterate over (key, value) pairs."""
        for key in self._embed_map:
            yield key, self[key]

    def _save_metadata(self) -> None:
        """Save embed store map to vlmeta."""
        metadata = {"embed_map": self._embed_map, "current_offset": self._current_offset}
        self._store.vlmeta["estore_metadata"] = metadata

    def _load_metadata(self) -> None:
        """Load embed store map from vlmeta."""
        if "estore_metadata" in self._store.vlmeta:
            metadata = self._store.vlmeta["estore_metadata"]
            self._embed_map = metadata["embed_map"]
            self._current_offset = metadata["current_offset"]
        else:
            self._embed_map = {}
            self._current_offset = 0

    def to_cframe(self) -> bytes:
        """Serialize embed store to CFrame format."""
        return self._store.to_cframe()

    def __enter__(self):
        """Context manager enter."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        # No need to close anything as SChunk/NDArray handles persistence automatically
        return False


def estore_from_cframe(cframe: bytes, copy: bool = False) -> EmbedStore:
    """
    Deserialize a CFrame to an EmbedStore object.

    Parameters
    ----------
    cframe : bytes
        CFrame data to deserialize.
    copy : bool, optional
        If True, copy the data. Default is False.

    Returns
    -------
    estore : EmbedStore
        The deserialized EmbedStore object.
    """
    schunk = blosc2.schunk_from_cframe(cframe, copy=copy)
    return EmbedStore(_from_schunk=schunk)


if __name__ == "__main__":
    # Example usage
    persistent = False
    if persistent:
        estore = EmbedStore(urlpath="example_estore.b2e", mode="w")  # , cparams=blosc2.CParams(clevel=0))
    else:
        estore = EmbedStore()  # , cparams=blosc2.CParams(clevel=0))
    # import pdb;  pdb.set_trace()
    estore["/node1"] = np.array([1, 2, 3])
    estore["/node2"] = blosc2.ones(2)
    urlpath = blosc2.URLPath("@public/examples/ds-1d.b2nd", "https://cat2.cloud/demo")
    arr_remote = blosc2.open(urlpath, mode="r")
    estore["/dir1/node3"] = arr_remote

    print("EmbedStore keys:", list(estore.keys()))
    print("Node1 data:", estore["/node1"][:])
    print("Node2 data:", estore["/node2"][:])
    print("Node3 data (remote):", estore["/dir1/node3"][:3])

    del estore["/node1"]
    print("After deletion, keys:", list(estore.keys()))

    # Reading back the estore
    if persistent:
        estore_read = EmbedStore(urlpath="example_estore.b2e", mode="r")
    else:
        estore_read = blosc2.from_cframe(estore.to_cframe())

    print("Read keys:", list(estore_read.keys()))
    for key, value in estore_read.items():
        print(
            f"shape of {key}: {value.shape}, dtype: {value.dtype}, map: {estore_read._embed_map[key]}, "
            f"values: {value[:10] if len(value) > 3 else value[:]}"
        )
